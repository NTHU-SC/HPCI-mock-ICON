// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

// libfyaml is not very clean...so we have to suppress some warnings

#if defined(__NVCOMPILER)
#  pragma diag_suppress unsigned_compare_with_zero
#elif defined(__GNUC__)
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Wpedantic"
#  pragma GCC diagnostic ignored "-Wall"
#  pragma GCC diagnostic ignored "-Wextra"
#endif
#include <libfyaml.h>
#if defined(__NVCOMPILER)
#  pragma diag_default unsigned_compare_with_zero
#elif defined(__GNUC__)
#  pragma GCC diagnostic pop
#endif

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "yac.h"
#include "utils_mci.h"
#include "config_yaml.h"
#include "mtime_calendar.h"
#include "geometry.h"
#include "io_utils.h"
#include "interp_method_avg.h"
#include "interp_method_callback.h"
#include "interp_method_check.h"
#include "interp_method_conserv.h"
#include "interp_method_creep.h"
#include "interp_method_file.h"
#include "interp_method_fixed.h"
#include "interp_method_hcsbb.h"
#include "interp_method_ncc.h"
#include "interp_method_nnn.h"
#include "interp_method_spmap.h"
#include "interp_stack_config.h"
#include "instance.h"
#include "fields.h"

typedef struct fy_document * fy_document_t;
typedef struct fy_node * fy_node_t;
typedef struct fy_node_pair * fy_node_pair_t;

enum {
  EMITTER_DEFAULT = FYECF_DEFAULT,
  EMITTER_JSON = FYECF_MODE_JSON,
  PARSER_DEFAULT = 0,
  PARSER_JSON_AUTO = FYPCF_JSON_AUTO,
  PARSER_JSON_FORCE = FYPCF_JSON_FORCE,
};

int const  YAC_YAML_EMITTER_DEFAULT = EMITTER_DEFAULT;
int const  YAC_YAML_EMITTER_JSON = EMITTER_JSON;
int const  YAC_YAML_PARSER_DEFAULT = PARSER_DEFAULT;
int const  YAC_YAML_PARSER_JSON_AUTO = PARSER_JSON_AUTO;
int const  YAC_YAML_PARSER_JSON_FORCE = PARSER_JSON_FORCE;

char const * yac_time_to_ISO(
  char const * time, enum yac_time_unit_type time_unit);

struct field_couple_buffer {
  struct {
    char const * comp_name;
    struct {
      char const ** name;
      size_t count;
    } grid;
    int lag;
  } src, tgt;
  char const * coupling_period;
  enum yac_reduction_type time_reduction;
  struct yac_interp_stack_config * interp_stack;
  struct  {
    char const * name;
    enum yac_weight_file_on_existing on_existing;
  } weight_file;
  int mapping_on_source;
  double scale_factor, scale_summand;

  struct field_couple_field_names {
    char const * src, * tgt;
  } * field_names;
  size_t num_field_names;

  char const ** src_mask_names;
  size_t num_src_mask_names;
  char const * tgt_mask_name;

  char const * yaxt_exchanger_name;

  int use_raw_exchange;
};

enum yaml_base_key_types {
  START_DATE,
  END_DATE,
  CALENDAR,
  TIMESTEP_UNIT,
  COUPLING,
  DEBUG,
};

enum yaml_couple_key_types {
  SOURCE_NAMES,
  SOURCE_COMPONENT,
  SOURCE_GRID,
  TARGET_NAMES,
  TARGET_COMPONENT,
  TARGET_GRID,
  FIELD,
  COUPLING_PERIOD,
  TIME_REDUCTION,
  SOURCE_LAG,
  TARGET_LAG,
  WEIGHT_FILE_DATA,
  WEIGHT_FILE_NAME,
  WEIGHT_FILE_ON_EXISTING,
  MAPPING_SIDE,
  SCALE_FACTOR,
  SCALE_SUMMAND,
  INTERPOLATION,
  SOURCE_MASK_NAME,
  SOURCE_MASK_NAMES,
  TARGET_MASK_NAME,
  YAXT_EXCHANGER_NAME,
  USE_RAW_EXCHANGE,
};

enum yaml_debug_key_types {
  GLOBAL_CONFIG,
  GLOBAL_DEFS,
  OUTPUT_GRIDS,
  MISSING_DEF
};

enum yaml_debug_sync_loc_key_types {
  SYNC_LOC_DEF_COMP,
  SYNC_LOC_SYNC_DEF,
  SYNC_LOC_ENDDEF,
  SYNC_LOC_COUNT, // has to be the last entry
};

enum yaml_comp_grid_names_key_types {
  COMP_NAME,
  GRID_NAMES,
};

enum yaml_weight_file_data_key_types {
  WEIGHT_FILE_DATA_NAME,
  WEIGHT_FILE_DATA_ON_EXISTING,
};

struct debug_config_file_buffer {
  struct debug_config_file {
    char const * name;
    enum yac_text_filetype type;
    int include_definitions;
  } config_file[SYNC_LOC_COUNT];
  char const * sync_loc_ref[SYNC_LOC_COUNT];
};

enum yaml_debug_output_grid_key_types{
  OUTPUT_GRID_GRID_NAME,
  OUTPUT_GRID_FILE_NAME
};

#define CAST_NAME_TYPE_PAIRS(...) (struct yac_name_type_pair[]) {__VA_ARGS__}
#define COUNT_NAME_TYPE_PAIRS(...) \
  sizeof(CAST_NAME_TYPE_PAIRS(__VA_ARGS__)) / sizeof(struct yac_name_type_pair)
#define DEF_NAME_TYPE_PAIR(NAME, TYPE) {.name = #NAME, .type = (int)(TYPE)}
#define DEF_NAME_TYPE_PAIRS(NAME, ...) \
  static const struct yac_name_type_pair NAME [] = {__VA_ARGS__}; \
  static const size_t num_ ## NAME = COUNT_NAME_TYPE_PAIRS(__VA_ARGS__);

DEF_NAME_TYPE_PAIRS(
  yaml_base_keys,
  DEF_NAME_TYPE_PAIR(start_date,     START_DATE),
  DEF_NAME_TYPE_PAIR(start_datetime, START_DATE),
  DEF_NAME_TYPE_PAIR(end_date,       END_DATE),
  DEF_NAME_TYPE_PAIR(end_datetime,   END_DATE),
  DEF_NAME_TYPE_PAIR(calendar,       CALENDAR),
  DEF_NAME_TYPE_PAIR(timestep_unit,  TIMESTEP_UNIT),
  DEF_NAME_TYPE_PAIR(coupling,       COUPLING),
  DEF_NAME_TYPE_PAIR(debug,          DEBUG))

DEF_NAME_TYPE_PAIRS(
  yaml_couple_keys,
  DEF_NAME_TYPE_PAIR(source,                  SOURCE_NAMES),
  DEF_NAME_TYPE_PAIR(src_component,           SOURCE_COMPONENT),
  DEF_NAME_TYPE_PAIR(src_grid,                SOURCE_GRID),
  DEF_NAME_TYPE_PAIR(target,                  TARGET_NAMES),
  DEF_NAME_TYPE_PAIR(tgt_component,           TARGET_COMPONENT),
  DEF_NAME_TYPE_PAIR(tgt_grid,                TARGET_GRID),
  DEF_NAME_TYPE_PAIR(field,                   FIELD),
  DEF_NAME_TYPE_PAIR(coupling_period,         COUPLING_PERIOD),
  DEF_NAME_TYPE_PAIR(time_reduction,          TIME_REDUCTION),
  DEF_NAME_TYPE_PAIR(src_lag,                 SOURCE_LAG),
  DEF_NAME_TYPE_PAIR(tgt_lag,                 TARGET_LAG),
  DEF_NAME_TYPE_PAIR(weight_file,             WEIGHT_FILE_DATA),
  DEF_NAME_TYPE_PAIR(weight_file_name,        WEIGHT_FILE_NAME),
  DEF_NAME_TYPE_PAIR(weight_file_on_existing, WEIGHT_FILE_ON_EXISTING),
  DEF_NAME_TYPE_PAIR(mapping_side,            MAPPING_SIDE),
  DEF_NAME_TYPE_PAIR(scale_factor,            SCALE_FACTOR),
  DEF_NAME_TYPE_PAIR(scale_summand,           SCALE_SUMMAND),
  DEF_NAME_TYPE_PAIR(interpolation,           INTERPOLATION),
  DEF_NAME_TYPE_PAIR(src_mask_name,           SOURCE_MASK_NAME),
  DEF_NAME_TYPE_PAIR(src_mask_names,          SOURCE_MASK_NAMES),
  DEF_NAME_TYPE_PAIR(tgt_mask_name,           TARGET_MASK_NAME),
  DEF_NAME_TYPE_PAIR(yaxt_exchanger,          YAXT_EXCHANGER_NAME),
  DEF_NAME_TYPE_PAIR(yaxt_exchanger_name,     YAXT_EXCHANGER_NAME),
  DEF_NAME_TYPE_PAIR(use_raw_exchange,        USE_RAW_EXCHANGE))

DEF_NAME_TYPE_PAIRS(
  yaml_debug_sync_loc_keys,
  DEF_NAME_TYPE_PAIR(def_comp, SYNC_LOC_DEF_COMP),
  DEF_NAME_TYPE_PAIR(DEF_COMP, SYNC_LOC_DEF_COMP),
  DEF_NAME_TYPE_PAIR(def_comps, SYNC_LOC_DEF_COMP),
  DEF_NAME_TYPE_PAIR(DEF_COMPS, SYNC_LOC_DEF_COMP),
  DEF_NAME_TYPE_PAIR(sync_def, SYNC_LOC_SYNC_DEF),
  DEF_NAME_TYPE_PAIR(SYNC_DEF, SYNC_LOC_SYNC_DEF),
  DEF_NAME_TYPE_PAIR(enddef, SYNC_LOC_ENDDEF),
  DEF_NAME_TYPE_PAIR(ENDDEF, SYNC_LOC_ENDDEF))

DEF_NAME_TYPE_PAIRS(
  yaml_debug_output_grid_keys,
  DEF_NAME_TYPE_PAIR(grid_name, OUTPUT_GRID_GRID_NAME),
  DEF_NAME_TYPE_PAIR(file_name, OUTPUT_GRID_FILE_NAME))

DEF_NAME_TYPE_PAIRS(
  bool_names,
  DEF_NAME_TYPE_PAIR(true,  true),
  DEF_NAME_TYPE_PAIR(TRUE,  true),
  DEF_NAME_TYPE_PAIR(yes,   true),
  DEF_NAME_TYPE_PAIR(YES,   true),
  DEF_NAME_TYPE_PAIR(false, false),
  DEF_NAME_TYPE_PAIR(FALSE, false),
  DEF_NAME_TYPE_PAIR(no,    false),
  DEF_NAME_TYPE_PAIR(NO,    false))

DEF_NAME_TYPE_PAIRS(
  timestep_units,
  DEF_NAME_TYPE_PAIR(millisecond, C_MILLISECOND),
  DEF_NAME_TYPE_PAIR(second,      C_SECOND),
  DEF_NAME_TYPE_PAIR(minute,      C_MINUTE),
  DEF_NAME_TYPE_PAIR(hour,        C_HOUR),
  DEF_NAME_TYPE_PAIR(day,         C_DAY),
  DEF_NAME_TYPE_PAIR(month,       C_MONTH),
  DEF_NAME_TYPE_PAIR(year,        C_YEAR),
  DEF_NAME_TYPE_PAIR(ISO_format,  C_ISO_FORMAT))

DEF_NAME_TYPE_PAIRS(
  time_operations,
  DEF_NAME_TYPE_PAIR(accumulate, TIME_ACCUMULATE),
  DEF_NAME_TYPE_PAIR(average,    TIME_AVERAGE),
  DEF_NAME_TYPE_PAIR(minimum,    TIME_MINIMUM),
  DEF_NAME_TYPE_PAIR(maximum,    TIME_MAXIMUM),
  DEF_NAME_TYPE_PAIR(none,       TIME_NONE))

DEF_NAME_TYPE_PAIRS(
  calendar_types,
  DEF_NAME_TYPE_PAIR(proleptic-gregorian, PROLEPTIC_GREGORIAN),
  DEF_NAME_TYPE_PAIR(360d, YEAR_OF_360_DAYS),
  DEF_NAME_TYPE_PAIR(365d, YEAR_OF_365_DAYS))

DEF_NAME_TYPE_PAIRS(
  mapping_sides,
  DEF_NAME_TYPE_PAIR(source, 1),
  DEF_NAME_TYPE_PAIR(target, 0))

DEF_NAME_TYPE_PAIRS(
  yaml_debug_keys,
  DEF_NAME_TYPE_PAIR(global_config, GLOBAL_CONFIG),
  DEF_NAME_TYPE_PAIR(output_grids, OUTPUT_GRIDS),
  DEF_NAME_TYPE_PAIR(missing_definition_is_fatal, MISSING_DEF))

DEF_NAME_TYPE_PAIRS(
  config_filetypes,
  DEF_NAME_TYPE_PAIR(yaml, YAC_TEXT_FILETYPE_YAML),
  DEF_NAME_TYPE_PAIR(YAML, YAC_TEXT_FILETYPE_YAML),
  DEF_NAME_TYPE_PAIR(json, YAC_TEXT_FILETYPE_JSON),
  DEF_NAME_TYPE_PAIR(JSON, YAC_TEXT_FILETYPE_JSON))

DEF_NAME_TYPE_PAIRS(
  role_types,
  DEF_NAME_TYPE_PAIR(target, TARGET),
  DEF_NAME_TYPE_PAIR(source, SOURCE),
  DEF_NAME_TYPE_PAIR(nothing, NOTHING))

DEF_NAME_TYPE_PAIRS(
  yaml_comp_grid_names_keys,
  DEF_NAME_TYPE_PAIR(component, COMP_NAME),
  DEF_NAME_TYPE_PAIR(grid,      GRID_NAMES))

DEF_NAME_TYPE_PAIRS(
  yaml_weight_file_data_keys,
  DEF_NAME_TYPE_PAIR(name,        WEIGHT_FILE_DATA_NAME),
  DEF_NAME_TYPE_PAIR(on_existing, WEIGHT_FILE_DATA_ON_EXISTING))

DEF_NAME_TYPE_PAIRS(
  weight_file_on_existing_types,
  DEF_NAME_TYPE_PAIR(error,     YAC_WEIGHT_FILE_ERROR),
  DEF_NAME_TYPE_PAIR(keep,      YAC_WEIGHT_FILE_KEEP),
  DEF_NAME_TYPE_PAIR(overwrite, YAC_WEIGHT_FILE_OVERWRITE))

enum interp_method_parameter_value_type{
  ENUM_PARAM,
  INT_PARAM,
  DBLE_PARAM,
  BOOL_PARAM,
  STR_PARAM,
  DEG_PARAM,
  MAP_PARAM,
  SEQ_PARAM,
};

struct interp_method_parameter_value;
typedef struct interp_method_parameter_value {

#if defined __NVCOMPILER && (__NVCOMPILER_MAJOR__ <= 23 || __NVCOMPILER_MAJOR__ == 24 && __NVCOMPILER_MINOR__ <= 3)
// Older versions of NVHPC have serious problems with unions that contain
// pointers that are not the first members: versions 23.7 and older fail with
// 'Internal compiler error. unhandled type', the newer versions produce code
// that fails at the runtime with
// 'Segmentation fault: address not mapped to object'.
  struct
#else
  union
#endif
  {
    int enum_value;
    int int_value;
    double dble_value;
    int bool_value;
    char const * str_value;
    struct interp_method_parameter_value * map_values;
    struct {
      struct interp_method_parameter_value * values;
      size_t count;
    } seq;
  } data;
  int is_set;
} interp_method_parameter_value;

struct interp_method_parameter {

  char const * name;
  enum interp_method_parameter_value_type type;

#if defined __NVCOMPILER && (__NVCOMPILER_MAJOR__ <= 23 || __NVCOMPILER_MAJOR__ == 24 && __NVCOMPILER_MINOR__ <= 3)
  // Older versions of NVHPC generate invalid LLVM IR for unions that contain
  // structures with doubles: an attempt to initialize a member of type double
  // that is not the first member of the respective structure will lead to a
  // type mismatch.
  struct
#else
  union
#endif
  {
    struct { // enum
      struct yac_name_type_pair const * valid_values;
      size_t num_valid_values;
    } enum_param;
    struct { // integer
      int valid_min, valid_max;
    } int_param;
    struct { // double
      double valid_min, valid_max;
    } dble_param;
    struct { // string
      size_t max_str_len;
    } str_param;
    struct { // bool
      int dummy;
    } bool_param;
    struct {
      struct interp_method_parameter * sub_params;
      size_t num_sub_params;
    } map_param;
    struct {
      struct interp_method_parameter * sub_param;
    } seq_param;
  } data;
  interp_method_parameter_value const default_value;
};

//! generates routine, which adds an interpolation based on a list of
//! interpolation parameters to an interpolation stack
#define DEF_INTERP_METHOD_ADD_FUNC(NAME, FUNC) \
  static void add_interp_method_ ## NAME ( \
    struct yac_interp_stack_config * interp_stack, \
    interp_method_parameter_value parameter_value, \
    char const * yaml_filename) { \
    char const * routine_name = "add_interp_method_" #NAME; \
    (void)parameter_value; \
    (void)yaml_filename; \
    (void)routine_name; \
    {FUNC} }
//! generates a routine, which gets interpolation parameters from
//! a interpolation stack entry
#define DEF_INTERP_METHOD_GET_FUNC(NAME, FUNC) \
  static void get_interp_method_ ## NAME ( \
    union yac_interp_stack_config_entry const * interp_stack_entry, \
    interp_method_parameter_value * parameter_value) { \
    char const * routine_name = "get_interp_method_" #NAME; \
    (void)interp_stack_entry; \
    (void)parameter_value; \
    (void)routine_name; \
    {FUNC} }

// macros for generating various types of interpolation parameters
#define DEF_ENUM_PARAM(NAME, DEFAULT, ...) \
  {.name = #NAME, \
   .type = ENUM_PARAM, \
   .data.enum_param = \
     {.valid_values = CAST_NAME_TYPE_PAIRS(__VA_ARGS__), \
      .num_valid_values = COUNT_NAME_TYPE_PAIRS(__VA_ARGS__)}, \
   .default_value.data.enum_value = (int)(DEFAULT)}
#define DEF_INT_PARAM(NAME, DEFAULT, VALID_MIN, VALID_MAX) \
  {.name = #NAME, \
   .type = INT_PARAM, \
   .data.int_param = \
     {.valid_min = (int)(VALID_MIN), \
      .valid_max = (int)(VALID_MAX)}, \
   .default_value.data.int_value = (int)(DEFAULT)}
#define DEF_DBLE_PARAM(NAME, DEFAULT, VALID_MIN, VALID_MAX) \
  {.name = #NAME, \
   .type = DBLE_PARAM, \
   .data.dble_param = \
     {.valid_min = (double)(VALID_MIN), \
      .valid_max = (double)(VALID_MAX)}, \
   .default_value.data.dble_value = (double)(DEFAULT)}
#define DEF_DEG_PARAM(NAME, DEFAULT, VALID_MIN, VALID_MAX) \
  {.name = #NAME, \
   .type = DEG_PARAM, \
   .data.dble_param = \
     {.valid_min = (double)(VALID_MIN), \
      .valid_max = (double)(VALID_MAX)}, \
   .default_value.data.dble_value = (double)(DEFAULT)}
#define DEF_BOOL_PARAM(NAME, DEFAULT) \
  {.name = #NAME, \
   .type = BOOL_PARAM, \
   .default_value.data.bool_value = (int)(DEFAULT)}
#define DEF_STR_PARAM(NAME, DEFAULT, MAX_STR_LEN) \
  {.name = #NAME, \
   .type = STR_PARAM, \
   .data.str_param.max_str_len = (MAX_STR_LEN), \
   .default_value.data.str_value = (DEFAULT)}
#define DEF_MAP_PARAM(NAME, SUBPARAMS) \
  {.name = #NAME, \
   .type = MAP_PARAM, \
   .data.map_param.sub_params = (SUBPARAMS ## _subparams), \
   .data.map_param.num_sub_params = (SUBPARAMS ## _subparams_count), \
   .default_value.data.map_values = NULL}
#define DEF_MAP_SUBPARAMS(NAME, ...) \
  static struct interp_method_parameter \
    NAME ## _subparams [] = {__VA_ARGS__}; \
  enum { \
    NAME ## _subparams_count = \
      sizeof(NAME ## _subparams) / sizeof(NAME ## _subparams[0])};
#define DEF_SEQ_PARAM(NAME, SUBPARAM) \
  {.name = #NAME, \
   .type = SEQ_PARAM, \
   .data.seq_param.sub_param = (SUBPARAM), \
   .default_value.data.seq.values = NULL, \
   .default_value.data.seq.count = 0}

// definition of some assert macros
#define YAML_ASSERT(CHECK, MSG) \
  YAC_ASSERT_F((CHECK), \
    "ERROR(%s): " MSG " in YAML configuration file \"%s\"", \
    routine_name, yaml_filename)
#define YAML_ASSERT_F(CHECK, MSG, ...) \
  YAC_ASSERT_F((CHECK), \
    "ERROR(%s): " MSG " in YAML configuration file \"%s\"", \
    routine_name, __VA_ARGS__, yaml_filename)
#define YAML_UNREACHABLE_DEFAULT(MSG) \
  YAC_UNREACHABLE_DEFAULT_F( \
    "ERROR(%s): " MSG " in YAML configuration file \"%s\"", \
    routine_name, yaml_filename)
#define YAML_UNREACHABLE_DEFAULT_F(MSG, ...) \
  YAC_UNREACHABLE_DEFAULT_F( \
    "ERROR(%s): " MSG " in YAML configuration file \"%s\"", \
    routine_name, __VA_ARGS__, yaml_filename)

#define DEF_INTERP_STACK_ADD(NAME, ...) \
  yac_interp_stack_config_add_ ## NAME ( \
    interp_stack, __VA_ARGS__);

#define DEF_INTERP_STACK_ADD_NO_PARAM(NAME) \
  yac_interp_stack_config_add_ ## NAME ( \
    interp_stack);

#define DEF_INTERP_STACK_GET(NAME, ...) \
  yac_interp_stack_config_entry_get_ ## NAME ( \
    interp_stack_entry, __VA_ARGS__);

#define DEF_INTERP_STACK_GET_NO_PARAM(NAME) \
  {}

#define DEF_INTERP_METHOD(YAML_NAME, FUNC_ADD, FUNC_GET, ...) \
  DEF_INTERP_METHOD_ADD_FUNC(YAML_NAME, FUNC_ADD) \
  DEF_INTERP_METHOD_GET_FUNC(YAML_NAME, FUNC_GET) \
  DEF_MAP_SUBPARAMS(YAML_NAME, __VA_ARGS__)

#define DEF_INTERP_METHOD_NO_PARAM(YAML_NAME, UI_NAME) \
  DEF_INTERP_METHOD_ADD_FUNC( \
    YAML_NAME, DEF_INTERP_STACK_ADD_NO_PARAM(UI_NAME)) \
  DEF_INTERP_METHOD_GET_FUNC( \
    YAML_NAME, DEF_INTERP_STACK_GET_NO_PARAM(UI_NAME))

// interpolation method average
DEF_INTERP_METHOD(average,
  DEF_INTERP_STACK_ADD(average,
    (enum yac_interp_avg_weight_type)
      parameter_value.data.map_values[0].data.enum_value,
    parameter_value.data.map_values[1].data.bool_value),
  enum yac_interp_avg_weight_type reduction_type;
  DEF_INTERP_STACK_GET(average,
    &reduction_type,
    &parameter_value->data.map_values[1].data.bool_value)
  parameter_value->data.map_values[0].data.enum_value =
    (int)reduction_type;,
  // paramters(map):
  //   weighted(enum)
  //   partial_coverage(bool)
  DEF_ENUM_PARAM(
    weighted, YAC_INTERP_AVG_WEIGHT_TYPE_DEFAULT,
    DEF_NAME_TYPE_PAIR(distance_weighted, YAC_INTERP_AVG_DIST),
    DEF_NAME_TYPE_PAIR(arithmetic_average, YAC_INTERP_AVG_ARITHMETIC),
    DEF_NAME_TYPE_PAIR(barycentric_coordinate, YAC_INTERP_AVG_BARY)),
  DEF_BOOL_PARAM(
    partial_coverage, YAC_INTERP_AVG_PARTIAL_COVERAGE_DEFAULT))

// interpolation method nearest corner cells
DEF_INTERP_METHOD(ncc,
  DEF_INTERP_STACK_ADD(ncc,
    (enum yac_interp_ncc_weight_type)
      parameter_value.data.map_values[0].data.enum_value,
    parameter_value.data.map_values[1].data.bool_value),
  enum yac_interp_ncc_weight_type weight_type;
  DEF_INTERP_STACK_GET(ncc,
    &weight_type,
    &parameter_value->data.map_values[1].data.bool_value)
  parameter_value->data.map_values[0].data.enum_value =
    (int)weight_type;,
  // parameters(map):
  //   weighted(enum)
  //   partial_coverage(bool)
  DEF_ENUM_PARAM(
    weighted, YAC_INTERP_NCC_WEIGHT_TYPE_DEFAULT,
    DEF_NAME_TYPE_PAIR(arithmetic_average, YAC_INTERP_NCC_AVG),
    DEF_NAME_TYPE_PAIR(distance_weighted, YAC_INTERP_NCC_DIST)),
  DEF_BOOL_PARAM(
    partial_coverage, YAC_INTERP_NCC_PARTIAL_COVERAGE_DEFAULT))

// interpolation method n-nearest neighbor
DEF_INTERP_METHOD(nnn,
  DEF_INTERP_STACK_ADD(nnn,
    (enum yac_interp_nnn_weight_type)
      parameter_value.data.map_values[0].data.enum_value,
    (size_t)parameter_value.data.map_values[1].data.int_value,
    parameter_value.data.map_values[2].data.dble_value,
    parameter_value.data.map_values[3].data.dble_value),
  enum yac_interp_nnn_weight_type type;
  size_t n;
  DEF_INTERP_STACK_GET(nnn,
    &type, &n, &parameter_value->data.map_values[2].data.dble_value,
    &parameter_value->data.map_values[3].data.dble_value)
  parameter_value->data.map_values[0].data.enum_value = (int)type;
  parameter_value->data.map_values[1].data.int_value = (int)n;
  parameter_value->data.map_values[2].data.dble_value /= YAC_RAD;,
  DEF_ENUM_PARAM(
    weighted, YAC_INTERP_NNN_WEIGHTED_DEFAULT,
    DEF_NAME_TYPE_PAIR(distance_weighted,  YAC_INTERP_NNN_DIST),
    DEF_NAME_TYPE_PAIR(gauss_weighted,     YAC_INTERP_NNN_GAUSS),
    DEF_NAME_TYPE_PAIR(arithmetic_average, YAC_INTERP_NNN_AVG),
    DEF_NAME_TYPE_PAIR(zero, YAC_INTERP_NNN_ZERO)),
  DEF_INT_PARAM(n, YAC_INTERP_NNN_N_DEFAULT, 1, INT_MAX),
  // paramters:
  //   weighted(enum)
  //   n(int)
  //   max_search_distance(deg)
  //   gauss_scale(dble)
  DEF_DEG_PARAM(
    max_search_distance, YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT,
    0.0, 179.9999),
  DEF_DBLE_PARAM(
    gauss_scale, YAC_INTERP_NNN_GAUSS_SCALE_DEFAULT, -DBL_MAX, DBL_MAX))

// interpolation method conservative
DEF_INTERP_METHOD(conservative,
  DEF_INTERP_STACK_ADD(conservative,
    parameter_value.data.map_values[0].data.int_value,
    parameter_value.data.map_values[1].data.bool_value,
    parameter_value.data.map_values[2].data.bool_value,
    (enum yac_interp_method_conserv_normalisation)
      parameter_value.data.map_values[3].data.enum_value),
  enum yac_interp_method_conserv_normalisation normalisation;
  DEF_INTERP_STACK_GET(conservative,
    &parameter_value->data.map_values[0].data.int_value,
    &parameter_value->data.map_values[1].data.bool_value,
    &parameter_value->data.map_values[2].data.bool_value,
    &normalisation)
  parameter_value->data.map_values[3].data.enum_value =
    (int)normalisation;,
  // parameter:
  //   order(int)
  //   enforced_conservation(bool)
  //   partial_coverage(bool)
  //   normalisation(enum)
  DEF_INT_PARAM(order, YAC_INTERP_CONSERV_ORDER_DEFAULT, 1, 2),
  DEF_BOOL_PARAM(
    enforced_conservation, YAC_INTERP_CONSERV_ENFORCED_CONSERV_DEFAULT),
  DEF_BOOL_PARAM(
    partial_coverage, YAC_INTERP_CONSERV_PARTIAL_COVERAGE_DEFAULT),
  DEF_ENUM_PARAM(
    normalisation, YAC_INTERP_CONSERV_NORMALISATION_DEFAULT,
    DEF_NAME_TYPE_PAIR(fracarea, YAC_INTERP_CONSERV_FRACAREA),
    DEF_NAME_TYPE_PAIR(destarea, YAC_INTERP_CONSERV_DESTAREA)))

static struct yac_spmap_cell_area_config * parse_cell_area_config_parameters(
  interp_method_parameter_value cell_area_config_param_values,
  char const * routine_name, char const * yaml_filename, char const * type) {

  interp_method_parameter_value yac_values =
    cell_area_config_param_values.data.map_values[0];
  interp_method_parameter_value file_values =
    cell_area_config_param_values.data.map_values[1];

  static struct yac_spmap_cell_area_config * cell_area_config;

  // if no cell area configuration is provided
  if ((!yac_values.is_set) && (!file_values.is_set)) {
      cell_area_config = YAC_INTERP_SPMAP_CELL_AREA_CONFIG_DEFAULT;
  } else {

    // ensure that not two cell area provider types are given
    YAML_ASSERT_F(
      (yac_values.is_set ^ file_values.is_set),
      "configurations for the computation and reading of %s cell areas "
      "were provided (only one is allowed)", type);

    if (yac_values.is_set) {
      cell_area_config =
        yac_spmap_cell_area_config_yac_new(
          yac_values.data.map_values[0].data.dble_value);
    } else {
      cell_area_config =
        yac_spmap_cell_area_config_file_new(
          file_values.data.map_values[0].data.str_value,
          file_values.data.map_values[1].data.str_value,
          file_values.data.map_values[2].data.int_value);
    }
  }

  return cell_area_config;
}

static void get_cell_area_config_parameters(
  struct yac_spmap_cell_area_config const * cell_area_config,
  double * sphere_radius, char const ** filename, char const ** varname,
  int * min_global_id, char const * desc) {

  enum yac_interp_spmap_cell_area_provider type =
    yac_spmap_cell_area_config_get_type(cell_area_config);

  switch(type) {
    YAC_UNREACHABLE_DEFAULT_F(
      "ERROR(get_cell_area_config_parameters): "
      "invalid cell area configuration type (%s)", desc);
    case(YAC_INTERP_SPMAP_CELL_AREA_FILE): {
      *sphere_radius = YAC_INTERP_SPMAP_SPHERE_RADIUS_DEFAULT;
      *filename = yac_spmap_cell_area_config_get_filename(cell_area_config);
      *varname = yac_spmap_cell_area_config_get_varname(cell_area_config);
      *min_global_id =
        (int)yac_spmap_cell_area_config_get_min_global_id(cell_area_config);
      break;
    }
    case(YAC_INTERP_SPMAP_CELL_AREA_YAC): {
      *sphere_radius =
        yac_spmap_cell_area_config_get_sphere_radius(cell_area_config);
      *filename = YAC_INTERP_SPMAP_FILENAME_DEFAULT;
      *varname = YAC_INTERP_SPMAP_VARNAME_DEFAULT;
      *min_global_id = YAC_INTERP_SPMAP_MIN_GLOBAL_ID_DEFAULT;
      break;
    }
  }
}

static void parse_point_coordinates(
  interp_method_parameter_value const * parameters,
  double * lon, double * lat) {

  YAC_ASSERT(
    parameters[0].is_set,
    "ERROR(parse_point_coordinates): longitude coordinate is not defined");
  YAC_ASSERT(
    parameters[1].is_set,
    "ERROR(parse_point_coordinates): longitude coordinate is not defined");

  *lon = parameters[0].data.dble_value;
  *lat = parameters[1].data.dble_value;
}

static struct yac_point_selection * parse_point_selection_parameter_bnd_circle(
  interp_method_parameter_value const parameter) {

  YAC_ASSERT(
    parameter.data.map_values[0].is_set,
    "ERROR(parse_point_selection_parameter_bnd_circle): "
    "bounding circle center is not defined");
  YAC_ASSERT(
    parameter.data.map_values[1].is_set,
    "ERROR(parse_point_selection_parameter_bnd_circle): "
    "bounding circle radius is not defined");

  double center_lon, center_lat, inc_angle;
  parse_point_coordinates(
    parameter.data.map_values[0].data.map_values, &center_lon, &center_lat);
  inc_angle = parameter.data.map_values[1].data.dble_value;

  return
    yac_point_selection_bnd_circle_new(center_lon, center_lat, inc_angle);
}

DEF_MAP_SUBPARAMS(
  bounding_circle_center,
  DEF_DEG_PARAM(lon, DBL_MAX, -720.0, 720.0),
  DEF_DEG_PARAM(lat, DBL_MAX, -90.0, 90.0))
DEF_MAP_SUBPARAMS(
  bounding_circle,
  DEF_MAP_PARAM(center, bounding_circle_center),
  DEF_DEG_PARAM(radius, DBL_MAX, 0.0, 179.9999))
DEF_MAP_SUBPARAMS(
  source_to_target_map_overwrite_condition,
  DEF_MAP_PARAM(bounding_circle, bounding_circle))
static struct yac_point_selection * parse_point_selection_parameter(
  interp_method_parameter_value const parameter) {

  size_t point_selection_type_idx = SIZE_MAX;
  for (size_t i = 0;
       i < source_to_target_map_overwrite_condition_subparams_count; ++i) {
    if (parameter.data.map_values[i].is_set) {
      YAC_ASSERT(
        point_selection_type_idx == SIZE_MAX,
        "ERROR(parse_point_selection_parameter): "
        "more than one point selection type was provided")
      point_selection_type_idx = i;
    }
  }

  YAC_ASSERT(
    point_selection_type_idx != SIZE_MAX,
    "ERROR(parse_point_selection_parameter): point selection is defined")

  struct yac_point_selection * src_point_selection;
  switch (point_selection_type_idx) {

    YAC_UNREACHABLE_DEFAULT("unsupported point selection type");

    // bounding circle
    case (0): {
      src_point_selection =
        parse_point_selection_parameter_bnd_circle(
          parameter.data.map_values[0]);
      break;
    }
  }

  return src_point_selection;
}

static struct yac_spmap_overwrite_config * parse_overwrite_config_parameter(
  interp_method_parameter_value const * parameter_value,
  struct yac_spmap_scale_config * scale_config) {

  /**
   * parameter_value[0]          condition
   * parameter_value[0][0]         bnd_circle
   * parameter_value[0][0][0]        center
   * parameter_value[0][0][0][0]       lon
   * parameter_value[0][0][0][1]       lat
   * parameter_value[0][0][1]        radius
   */

  // extract data from parameter
  struct yac_point_selection * src_point_selection =
    parse_point_selection_parameter(
      parameter_value->data.map_values[0]);
  double spread_distance =
    parameter_value->data.map_values[1].data.dble_value;
  double max_search_distance =
    parameter_value->data.map_values[2].data.dble_value;

  enum yac_interp_spmap_weight_type weight_type =
    (enum yac_interp_spmap_weight_type)
      parameter_value->data.map_values[3].data.enum_value;
  struct yac_interp_spmap_config * spmap_config =
    yac_interp_spmap_config_new(
      spread_distance, max_search_distance, weight_type, scale_config);

  struct yac_spmap_overwrite_config * overwrite_config =
    yac_spmap_overwrite_config_new(src_point_selection, spmap_config);

  yac_interp_spmap_config_delete(spmap_config);
  yac_point_selection_delete(src_point_selection);

  return overwrite_config;
}

static void get_point_selection_parameters(
  struct yac_point_selection const * point_selection,
  interp_method_parameter_value * parameter_value) {

  /**
   * parameter_value[0]          condition
   * parameter_value[0][0]         bnd_circle
   * parameter_value[0][0][0]        center
   * parameter_value[0][0][0][0]       lon
   * parameter_value[0][0][0][1]       lat
   * parameter_value[0][0][1]        radius
   */

  switch (yac_point_selection_get_type(point_selection)) {

    YAC_UNREACHABLE_DEFAULT_F(
      "ERROR(get_point_selection_parameters): %s point selection type",
      (yac_point_selection_get_type(point_selection) ==
       YAC_POINT_SELECTION_TYPE_EMPTY)?"empty":"invalid");

    case (YAC_POINT_SELECTION_TYPE_BND_CIRCLE): {
      double center_lon, center_lat, radius;
      yac_point_selection_bnd_circle_get_config(
        point_selection, &center_lon, &center_lat, &radius);
      parameter_value[0].                        // condition
        data.map_values[0].                      //   bnd_circle
          data.map_values[0].                    //     center
            data.map_values[0].data.dble_value = //       lon
              center_lon / YAC_RAD;
      parameter_value[0].                        // condition
        data.map_values[0].                      //   bnd_circle
          data.map_values[0].                    //     center
            data.map_values[1].data.dble_value = //       lat
              center_lat / YAC_RAD;
      parameter_value[0].                        // condition
        data.map_values[0].                      //   bnd_circle
          data.map_values[1].data.dble_value =   //     radius
            radius / YAC_RAD;
      break;
    }
  }
}

static interp_method_parameter_value
  yaml_interp_method_parameter_get_default(
    struct interp_method_parameter parameter);

DEF_MAP_SUBPARAMS(
  source_to_target_map_overwrite,
  DEF_MAP_PARAM(condition, source_to_target_map_overwrite_condition),
  DEF_DEG_PARAM(
    spread_distance, YAC_INTERP_SPMAP_SPREAD_DISTANCE_DEFAULT, 0.0, 89.9999),
  DEF_DEG_PARAM(
    max_search_distance,
    YAC_INTERP_SPMAP_MAX_SEARCH_DISTANCE_DEFAULT, 0.0, 179.9999),
  DEF_ENUM_PARAM(
    weighted, YAC_INTERP_SPMAP_WEIGHTED_DEFAULT,
    DEF_NAME_TYPE_PAIR(distance_weighted, YAC_INTERP_SPMAP_DIST),
    DEF_NAME_TYPE_PAIR(arithmetic_average, YAC_INTERP_SPMAP_AVG)))
static struct interp_method_parameter
  source_to_target_map_overwrite_seq =
    DEF_MAP_PARAM(overwrite, source_to_target_map_overwrite);
static void get_overwrite_config_parameters(
  struct yac_spmap_overwrite_config const * overwrite_config,
  interp_method_parameter_value * parameter_value) {

  /**
   * parameter_value[0] condition
   * parameter_value[1] spread_distance
   * parameter_value[2] max_search_distance
   * parameter_value[3] weighted
   */

  struct yac_point_selection const * src_point_selection =
    yac_spmap_overwrite_config_get_src_point_selection(overwrite_config);
  struct yac_interp_spmap_config const * spmap_config =
    yac_spmap_overwrite_config_get_spmap_config(overwrite_config);

  *parameter_value =
    yaml_interp_method_parameter_get_default(
      source_to_target_map_overwrite_seq);

  get_point_selection_parameters(
    src_point_selection, parameter_value->data.map_values + 0);
  parameter_value->data.map_values[1].data.dble_value =
    yac_interp_spmap_config_get_spread_distance(spmap_config) / YAC_RAD;
  parameter_value->data.map_values[2].data.dble_value =
    yac_interp_spmap_config_get_max_search_distance(spmap_config) / YAC_RAD;
  parameter_value->data.map_values[3].data.enum_value =
    yac_interp_spmap_config_get_weight_type(spmap_config);
}

// interpolation method source to target mapping
DEF_MAP_SUBPARAMS(
  source_to_target_map_cell_area_provider_yac,
  DEF_DBLE_PARAM(
    sphere_radius, YAC_INTERP_SPMAP_SPHERE_RADIUS_DEFAULT, 1e-9, DBL_MAX))
DEF_MAP_SUBPARAMS(
  source_to_target_map_cell_area_provider_file,
  DEF_STR_PARAM(filename, NULL, YAC_MAX_FILE_NAME_LENGTH),
  DEF_STR_PARAM(varname, NULL, 128),
  DEF_INT_PARAM(
    min_global_id, YAC_INTERP_SPMAP_MIN_GLOBAL_ID_DEFAULT, -INT_MAX, INT_MAX))
DEF_MAP_SUBPARAMS(
  source_to_target_map_cell_area_provider,
  DEF_MAP_PARAM(yac, source_to_target_map_cell_area_provider_yac),
  DEF_MAP_PARAM(file, source_to_target_map_cell_area_provider_file))

DEF_INTERP_METHOD(source_to_target_map,

  // extracting configuration from parametersm, which were read from file
  double spread_distance =
    parameter_value.data.map_values[0].data.dble_value;
  double max_search_distance =
    parameter_value.data.map_values[1].data.dble_value;
  enum yac_interp_spmap_weight_type weight_type =
    (enum yac_interp_spmap_weight_type)
      parameter_value.data.map_values[2].data.enum_value;
  enum yac_interp_spmap_scale_type scale_type =
    (enum yac_interp_spmap_scale_type)
      parameter_value.data.map_values[3].data.enum_value;
  struct yac_spmap_cell_area_config * source_cell_area_config =
    parse_cell_area_config_parameters(
      parameter_value.data.map_values[4],
      routine_name, yaml_filename, "source");
  struct yac_spmap_cell_area_config * target_cell_area_config =
    parse_cell_area_config_parameters(
      parameter_value.data.map_values[5],
      routine_name, yaml_filename, "target");
  struct yac_spmap_scale_config * scale_config =
    yac_spmap_scale_config_new(
      scale_type, source_cell_area_config, target_cell_area_config);
  yac_spmap_cell_area_config_delete(source_cell_area_config);
  yac_spmap_cell_area_config_delete(target_cell_area_config);
  struct yac_interp_spmap_config * default_config =
    yac_interp_spmap_config_new(
      spread_distance, max_search_distance, weight_type, scale_config);
  struct yac_spmap_overwrite_config ** overwrite_configs =
    (parameter_value.data.map_values[6].data.seq.count > 0)?
      xcalloc(
        (parameter_value.data.map_values[6].data.seq.count + 1),
        sizeof(*overwrite_configs)):NULL;
  for (size_t i = 0; i < parameter_value.data.map_values[6].data.seq.count;
       ++i)
    overwrite_configs[i] =
      parse_overwrite_config_parameter(
        parameter_value.data.map_values[6].data.seq.values + i, scale_config);
  yac_spmap_scale_config_delete(scale_config);

  // set configuration
  DEF_INTERP_STACK_ADD(spmap_ext, default_config, overwrite_configs)
  yac_spmap_overwrite_configs_delete(overwrite_configs);
  yac_interp_spmap_config_delete(default_config);,

  // generate parameters from configuration in order to be written to file
  struct yac_interp_spmap_config const * default_config;
  struct yac_spmap_overwrite_config const ** overwrite_configs;
  DEF_INTERP_STACK_GET(spmap_ext, &default_config, &overwrite_configs)
  struct yac_spmap_scale_config const * scale_config =
    yac_interp_spmap_config_get_scale_config(default_config);
  parameter_value->data.map_values[0].data.dble_value =
    yac_interp_spmap_config_get_spread_distance(default_config) / YAC_RAD;
  parameter_value->data.map_values[1].data.dble_value =
    yac_interp_spmap_config_get_max_search_distance(default_config) / YAC_RAD;
  parameter_value->data.map_values[2].data.enum_value =
    (int)yac_interp_spmap_config_get_weight_type(default_config);
  parameter_value->data.map_values[3].data.enum_value =
    (int)yac_spmap_scale_config_get_type(scale_config);
  get_cell_area_config_parameters(
     yac_spmap_scale_config_get_src_cell_area_config(scale_config),
    &parameter_value->data.map_values[4].data.map_values[0].data.map_values[0].data.dble_value,
    &parameter_value->data.map_values[4].data.map_values[1].data.map_values[0].data.str_value,
    &parameter_value->data.map_values[4].data.map_values[1].data.map_values[1].data.str_value,
    &parameter_value->data.map_values[4].data.map_values[1].data.map_values[2].data.int_value,
    "source");
  get_cell_area_config_parameters(
     yac_spmap_scale_config_get_tgt_cell_area_config(scale_config),
    &parameter_value->data.map_values[5].data.map_values[0].data.map_values[0].data.dble_value,
    &parameter_value->data.map_values[5].data.map_values[1].data.map_values[0].data.str_value,
    &parameter_value->data.map_values[5].data.map_values[1].data.map_values[1].data.str_value,
    &parameter_value->data.map_values[5].data.map_values[1].data.map_values[2].data.int_value,
    "target");
  size_t overwrite_config_count = 0;
  for (;
       (overwrite_configs != NULL) &&
       (overwrite_configs[overwrite_config_count] != NULL);
       ++overwrite_config_count);
  parameter_value->data.map_values[6].data.seq.values =
    xmalloc(
      overwrite_config_count *
      sizeof(
        *(parameter_value->data.map_values[6].data.seq.values)));
  parameter_value->data.map_values[6].data.seq.count =
    overwrite_config_count;
  for (size_t i = 0; i < overwrite_config_count; ++i)
    get_overwrite_config_parameters(
      overwrite_configs[i],
        parameter_value->data.map_values[6].data.seq.values + i);,

  // define parameters:
  // parameter_values[0]       spread_distance
  // parameter_values[1]       max_search_distance
  // parameter_values[2]       weighted
  // parameter_values[3]       scale
  // parameter_values[4]       src_cell_area
  // parameter_values[4][0]                 .yac
  // parameter_values[4][0][0]                  .sphere_radius
  // parameter_values[4][1]                 .file
  // parameter_values[4][1][0]                   .filename
  // parameter_values[4][1][1]                   .varname
  // parameter_values[4][1][2]                   .min_global_id
  // parameter_values[5]       tgt_cell_area
  // parameter_values[5][0]                 .yac
  // parameter_values[5][0][0]                  .sphere_radius
  // parameter_values[5][1]                 .file
  // parameter_values[5][1][0]                   .filename
  // parameter_values[5][1][1]                   .varname
  // parameter_values[5][1][2]                   .min_global_id
  DEF_DEG_PARAM(
    spread_distance, YAC_INTERP_SPMAP_SPREAD_DISTANCE_DEFAULT, 0.0, 89.9999),
  DEF_DEG_PARAM(
    max_search_distance,
    YAC_INTERP_SPMAP_MAX_SEARCH_DISTANCE_DEFAULT, 0.0, 179.9999),
  DEF_ENUM_PARAM(
    weighted, YAC_INTERP_SPMAP_WEIGHTED_DEFAULT,
    DEF_NAME_TYPE_PAIR(distance_weighted, YAC_INTERP_SPMAP_DIST),
    DEF_NAME_TYPE_PAIR(arithmetic_average, YAC_INTERP_SPMAP_AVG)),
  DEF_ENUM_PARAM(
    scale, YAC_INTERP_SPMAP_SCALE_TYPE_DEFAULT,
    DEF_NAME_TYPE_PAIR(none, YAC_INTERP_SPMAP_NONE),
    DEF_NAME_TYPE_PAIR(srcarea, YAC_INTERP_SPMAP_SRCAREA),
    DEF_NAME_TYPE_PAIR(invtgtarea, YAC_INTERP_SPMAP_INVTGTAREA),
    DEF_NAME_TYPE_PAIR(fracarea, YAC_INTERP_SPMAP_FRACAREA)),
  DEF_MAP_PARAM(
    src_cell_area, source_to_target_map_cell_area_provider),
  DEF_MAP_PARAM(
    tgt_cell_area, source_to_target_map_cell_area_provider),
  DEF_SEQ_PARAM(
    overwrite, &source_to_target_map_overwrite_seq))

// interpolation method fixed
DEF_INTERP_METHOD(fixed,
  YAML_ASSERT(
    parameter_value.data.map_values[0].data.dble_value != DBL_MAX,
    "parameter 'user_value' of interpolation method 'fixed' is unset");
  DEF_INTERP_STACK_ADD(
    fixed, parameter_value.data.map_values[0].data.dble_value),
  DEF_INTERP_STACK_GET(
    fixed, &parameter_value->data.map_values[0].data.dble_value),
  // paramters:
  //   fixed(dble)
  DEF_DBLE_PARAM(
    user_value, YAC_INTERP_FIXED_VALUE_DEFAULT, -DBL_MAX, DBL_MAX))

// interpolation method user file
DEF_INTERP_METHOD(user_file,
  YAML_ASSERT(
    parameter_value.data.map_values[0].data.str_value,
    "parameter \"filename\" of interpolation method \"user file\" is unset");
  DEF_INTERP_STACK_ADD(user_file,
    (char*)(parameter_value.data.map_values[0].data.str_value),
    (enum yac_interp_file_on_missing_file)
      parameter_value.data.map_values[1].data.enum_value,
    (enum yac_interp_file_on_success)
      parameter_value.data.map_values[2].data.enum_value),
  enum yac_interp_file_on_missing_file on_missing_file;
  enum yac_interp_file_on_success on_success;
  DEF_INTERP_STACK_GET(user_file,
    &parameter_value->data.map_values[0].data.str_value,
    &on_missing_file, &on_success)
  parameter_value->data.map_values[1].data.enum_value =
    (int)on_missing_file;
  parameter_value->data.map_values[2].data.enum_value =
    (int)on_success;,
  // paramter:
  //   filename(str)
  //   on_missing_file(enum)
  //   on_success(enum)
  DEF_STR_PARAM(
    filename, YAC_INTERP_FILE_WEIGHT_FILE_NAME_DEFAULT,
    YAC_MAX_FILE_NAME_LENGTH),
  DEF_ENUM_PARAM(
    on_missing_file, YAC_INTERP_FILE_ON_MISSING_FILE_DEFAULT,
    DEF_NAME_TYPE_PAIR(error, YAC_INTERP_FILE_MISSING_ERROR),
    DEF_NAME_TYPE_PAIR(cont, YAC_INTERP_FILE_MISSING_CONT),
    DEF_NAME_TYPE_PAIR(continue, YAC_INTERP_FILE_MISSING_CONT)),
  DEF_ENUM_PARAM(
    on_success, YAC_INTERP_FILE_ON_SUCCESS_DEFAULT,
    DEF_NAME_TYPE_PAIR(stop, YAC_INTERP_FILE_SUCCESS_STOP),
    DEF_NAME_TYPE_PAIR(cont, YAC_INTERP_FILE_SUCCESS_CONT),
    DEF_NAME_TYPE_PAIR(continue, YAC_INTERP_FILE_SUCCESS_CONT)))

// interpolation method check
DEF_INTERP_METHOD(check,
  DEF_INTERP_STACK_ADD(check,
    (char*)(parameter_value.data.map_values[0].data.str_value),
    (char*)(parameter_value.data.map_values[1].data.str_value)),
  DEF_INTERP_STACK_GET(check,
    &parameter_value->data.map_values[0].data.str_value,
    &parameter_value->data.map_values[1].data.str_value),
  // paramters:
  //   constructor_key(str)
  //   do_search_key(str)
  DEF_STR_PARAM(
    constructor_key, YAC_INTERP_CHECK_CONSTRUCTOR_KEY_DEFAULT,
    YAC_MAX_ROUTINE_NAME_LENGTH),
  DEF_STR_PARAM(
    do_search_key, YAC_INTERP_CHECK_DO_SEARCH_KEY_DEFAULT,
    YAC_MAX_ROUTINE_NAME_LENGTH))

// interpolation method Bernstein Bezier
DEF_INTERP_METHOD_NO_PARAM(bernstein_bezier, hcsbb)

// interpolation method radial basis function
DEF_INTERP_METHOD(rbf,
  DEF_INTERP_STACK_ADD(rbf,
    (size_t)parameter_value.data.map_values[0].data.int_value,
    parameter_value.data.map_values[1].data.dble_value,
    parameter_value.data.map_values[2].data.dble_value),
  size_t n;
  DEF_INTERP_STACK_GET(rbf,
    &n, &parameter_value->data.map_values[1].data.dble_value,
    &parameter_value->data.map_values[2].data.dble_value)
  parameter_value->data.map_values[0].data.int_value = (int)n;
  parameter_value->data.map_values[1].data.dble_value /= YAC_RAD;
  parameter_value->data.map_values[3].data.enum_value = (int)0;,
  // parameters:
  //   n(int)
  //   max_search_distance(dble)
  //   rbf_scale(dble)
  //   rbf_kernel(enum)
  DEF_INT_PARAM(n, YAC_INTERP_RBF_N_DEFAULT, 1, INT_MAX),
  DEF_DEG_PARAM(
    max_search_distance, YAC_INTERP_RBF_MAX_SEARCH_DISTANCE_DEFAULT,
    0.0, 179.9999),
  DEF_DBLE_PARAM(rbf_scale, YAC_INTERP_RBF_SCALE_DEFAULT, -DBL_MAX, DBL_MAX),
  DEF_ENUM_PARAM(
    rbf_kernel, YAC_INTERP_RBF_KERNEL_DEFAULT,
    DEF_NAME_TYPE_PAIR(gauss_kernel, 0)))

// interpolation method creep
DEF_INTERP_METHOD(creep,
  DEF_INTERP_STACK_ADD(
    creep, parameter_value.data.map_values[0].data.int_value),
  DEF_INTERP_STACK_GET(
    creep, &parameter_value->data.map_values[0].data.int_value),
  // parameters:
  //   creep_distance(int)
  DEF_INT_PARAM(
    creep_distance, YAC_INTERP_CREEP_DISTANCE_DEFAULT, -1, INT_MAX))

// interpolation method user_callback
DEF_INTERP_METHOD(user_callback,
  YAML_ASSERT(
    parameter_value.data.map_values[0].data.str_value,
    "parameter \"func_compute_weights\" "
    "of interpolation method \"user callback\" is unset")
  DEF_INTERP_STACK_ADD(
    user_callback,
    (char*)(parameter_value.data.map_values[0].data.str_value)),
  DEF_INTERP_STACK_GET(
    user_callback,
    &parameter_value->data.map_values[0].data.str_value),
  // parameters:
  //   func_compute_weights(str)
  DEF_STR_PARAM(
    func_compute_weights, YAC_INTERP_CALLBACK_COMPUTE_WEIGHTS_KEY_DEFAULT,
    YAC_MAX_ROUTINE_NAME_LENGTH))

#define ADD_INTERPOLATION(NAME, TYPE) \
  {.name = #NAME , \
   .type = TYPE , \
   .add_interpolation = add_interp_method_ ## NAME , \
   .get_interpolation = get_interp_method_ ## NAME , \
   .parameter = DEF_MAP_PARAM(NAME, NAME)}
#define ADD_INTERPOLATION_NO_PARAM(NAME, TYPE) \
  {.name = #NAME , \
   .type = TYPE , \
   .add_interpolation = add_interp_method_ ## NAME , \
   .get_interpolation = get_interp_method_ ## NAME , \
   .parameter = \
     {.name = #NAME, \
      .type = MAP_PARAM, \
      .data.map_param.sub_params = NULL, \
      .data.map_param.num_sub_params = 0, \
      .default_value.data.map_values = NULL}}

// data structure containing all available interpolation methods
struct yac_interpolation_method {
  char const * name;
  enum yac_interpolation_list type;
  void(*add_interpolation)(
    struct yac_interp_stack_config * interp_stack,
    interp_method_parameter_value parameter_value,
    char const * yaml_filename);
  void(*get_interpolation)(
    union yac_interp_stack_config_entry const * interp_stack_entry,
    interp_method_parameter_value * parameter_value);
  struct interp_method_parameter parameter;
} const interpolation_methods[] =
  {ADD_INTERPOLATION(average, YAC_AVERAGE),
   ADD_INTERPOLATION(ncc, YAC_NEAREST_CORNER_CELLS),
   ADD_INTERPOLATION(nnn, YAC_N_NEAREST_NEIGHBOR),
   ADD_INTERPOLATION(conservative, YAC_CONSERVATIVE),
   ADD_INTERPOLATION(source_to_target_map, YAC_SOURCE_TO_TARGET_MAP),
   ADD_INTERPOLATION(fixed, YAC_FIXED_VALUE),
   ADD_INTERPOLATION(user_file, YAC_USER_FILE),
   ADD_INTERPOLATION(check, YAC_CHECK),
   ADD_INTERPOLATION_NO_PARAM(bernstein_bezier, YAC_BERNSTEIN_BEZIER),
   ADD_INTERPOLATION(rbf, YAC_RADIAL_BASIS_FUNCTION),
   ADD_INTERPOLATION(creep, YAC_CREEP),
   ADD_INTERPOLATION(user_callback, YAC_USER_CALLBACK)};
enum {
  NUM_INTERPOLATION_METHODS =
    sizeof(interpolation_methods)/sizeof(interpolation_methods[0]),
};

static char const * yaml_parse_string_value(
  fy_node_t value_node, char const * name, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_string_value";

  YAML_ASSERT_F(
    value_node && fy_node_is_scalar(value_node),
    "unsupported node type for \"%s\" (the node is expected to be scalar)",
    name);

  return fy_node_get_scalar0(value_node);
}

static calendarType yaml_parse_calendar_value(
  fy_node_t value_node, char const * key_name, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_calendar_value";

  char const * calendar_name =
    yaml_parse_string_value(value_node, key_name, yaml_filename);

  int calendar_type =
    yac_name_type_pair_get_type(
      calendar_types, num_calendar_types, calendar_name);

  YAML_ASSERT_F(
    calendar_type != INT_MAX,
    "\"%s\" is not a valid calendar name", calendar_name);

  return (calendarType)calendar_type;
}

static char const * yaml_parse_timestep_value(
  fy_node_t value_node, char const * key_name, char const * yaml_filename,
  enum yac_time_unit_type time_unit) {
  char const * routine_name = "yaml_parse_timestep_value";

  YAML_ASSERT(
    time_unit != TIME_UNIT_UNDEFINED, "time unit is not yet defined");
  YAML_ASSERT_F(
    value_node && fy_node_is_scalar(value_node),
    "unsupported node type for \"%s\" (the node is expected to be scalar)",
    key_name);
  char const * timestep =
    yaml_parse_string_value(value_node, key_name, yaml_filename);
  char const * timestep_iso =
    yac_time_to_ISO(timestep, time_unit);

  YAML_ASSERT_F(
    timestep_iso, "valid to convert timestep \"%s\" to ISO 8601 format",
    timestep);

  return strdup(timestep_iso);
}

static enum yac_time_unit_type yaml_parse_timestep_unit_value(
  fy_node_t value_node, char const * key_name, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_timestep_unit_value";

  char const * timestep_unit_str =
    yaml_parse_string_value(value_node, key_name, yaml_filename);

  int timestep_unit =
    yac_name_type_pair_get_type(
      timestep_units, num_timestep_units, timestep_unit_str);

  YAML_ASSERT_F(
    timestep_unit != INT_MAX,
    "\"%s\" is not a valid time step unit", timestep_unit_str);

  return (enum yac_time_unit_type)timestep_unit;
}

static enum yac_reduction_type yaml_parse_time_reduction_value(
  fy_node_t value_node, char const * key_name, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_time_reduction_value";

  char const * time_reduction_str =
    yaml_parse_string_value(value_node, key_name, yaml_filename);

  int time_reduction =
    yac_name_type_pair_get_type(
      time_operations, num_time_operations, time_reduction_str);

  YAML_ASSERT_F(
    time_reduction != INT_MAX,
    "\"%s\" is not a valid time reduction type in", time_reduction_str);

  return (enum yac_reduction_type)time_reduction;
}

static int yaml_parse_integer_value(
  fy_node_t value_node, char const * key_name, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_integer_value";

  char const * integer_str =
    yaml_parse_string_value(value_node, key_name, yaml_filename);

  char * endptr;
  long int long_value = strtol(integer_str, &endptr, 10);

  YAML_ASSERT_F(
    (endptr != integer_str) && (*endptr == '\0') &&
    (long_value >= INT_MIN) && (long_value <= INT_MAX),
    "\"%s\" is not a valid integer value", integer_str);

  return (int)long_value;
}

static double yaml_parse_double_value(
  fy_node_t value_node, char const * key_name, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_double_value";

  char const * double_str =
    yaml_parse_string_value(value_node, key_name, yaml_filename);

  char * endptr;
  double dble_value = strtod(double_str, &endptr);

  YAML_ASSERT_F(
    (endptr != double_str) && (*endptr == '\0'),
    "\"%s\" is not a valid double value for \"%s\"", double_str, key_name);

  return dble_value;
}

static int yaml_parse_enum_value(
  struct yac_name_type_pair const * valid_values, size_t num_valid_values,
  fy_node_t value_node, char const * key_name, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_enum_value";

  char const * value_str =
    yaml_parse_string_value(value_node, key_name, yaml_filename);

  int value =
    yac_name_type_pair_get_type(
      valid_values, num_valid_values, value_str);

  YAML_ASSERT_F(
    value != INT_MAX,
    "\"%s\" is not a valid enum value for \"%s\" ", value_str, key_name);

  return value;
}

/// @brief reads interpolation name and user-provided parameters from YAML node
/// @param[out] interpolation_type_str name of the interpolation
/// @param[out] parameter_node         YAML node containing user-provided
///                                    interpolation method configuration
/// @param[in]  interp_method_node     YAML node containing interpolation
/// @param[in]  yaml_filename          name of the YAML file
///                                    (used for error messages) 
/// @remark interp_method_node can be either a scalar not (only interpolation
///         name and no parameters) or a map (interpolation name and parameters)
static void yaml_parse_base_interp_method_node(
  char const ** interpolation_type_str, fy_node_t * parameter_node,
  fy_node_t interp_method_node, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_base_interp_method_node";

  fy_node_t interpolation_name_node;

  // determine type of interpolation node
  switch(fy_node_get_type(interp_method_node)) {
    YAML_UNREACHABLE_DEFAULT(
      "unsupported interpolation method node type "
      "(interpolation methods are expected to be defined as either scalar "
      "or maps)");

    // interpolation node only contains the name and not parameters
    case (FYNT_SCALAR):
      interpolation_name_node = interp_method_node;
      *parameter_node = NULL;
      break;

    // interpolation node contains name and parameters
    case (FYNT_MAPPING):

      YAML_ASSERT(
        fy_node_mapping_item_count(interp_method_node) == 1,
        "base interpolation method node is only allowed to have one pair ");

      fy_node_pair_t base_interp_method_pair =
        fy_node_mapping_get_by_index(interp_method_node, 0);

      interpolation_name_node = fy_node_pair_key(base_interp_method_pair);

      fy_node_t base_interp_method_value_node =
        fy_node_pair_value(base_interp_method_pair);

      YAML_ASSERT_F(
        fy_node_is_mapping(base_interp_method_value_node),
        "unsupported base interpolation method value node type "
        "for interpolation method \"%s\" "
        "(interpolation method parameters are expected to be "
        "defined as maps)",
        yaml_parse_string_value(
          interpolation_name_node, "interpolation method name", yaml_filename));

      *parameter_node = base_interp_method_value_node;
      break;
  }

  *interpolation_type_str =
    yaml_parse_string_value(
      interpolation_name_node, "interpolation method name", yaml_filename);
}

/// @brief allocates parameter values and assigns default values based on
///        provided interpolation parameter definition
/// @param[in] parameter parameter definition
/// @return parameter value with default values
static interp_method_parameter_value
  yaml_interp_method_parameter_get_default(
    struct interp_method_parameter parameter) {

  interp_method_parameter_value parameter_value = {.is_set = 0};

  switch(parameter.type) {

    YAC_UNREACHABLE_DEFAULT(
      "ERROR(yaml_interp_method_parameter_get_default): "
      "unsupported parameter type");

    // parameter is a sequence
    case (SEQ_PARAM): {

      // by default all sequences are empty
      parameter_value.data.seq.values = NULL;
      parameter_value.data.seq.count = 0;
      break;
    }

    // parameter is a map
    case (MAP_PARAM): {

      parameter_value.data.map_values =
        xmalloc(
          parameter.data.map_param.num_sub_params *
          sizeof(*(parameter_value.data.map_values)));
      for (size_t i = 0; i < parameter.data.map_param.num_sub_params; ++i)
        parameter_value.data.map_values[i] =
          yaml_interp_method_parameter_get_default(
            parameter.data.map_param.sub_params[i]);
      break;
    }

    // parameter is a scalar
    case (ENUM_PARAM):
    case (INT_PARAM):
    case (DBLE_PARAM):
    case (BOOL_PARAM):
    case (STR_PARAM):
    case (DEG_PARAM):
    {
      parameter_value = parameter.default_value;
      break;
    }
  }

  return parameter_value;
}

/// @brief free memory associated with 
/// @param parameter_value parameter values to be freed
/// @param parameter       parameter definition associated with the values
static void yaml_interp_method_parameter_free(
  interp_method_parameter_value parameter_value,
  struct interp_method_parameter parameter) {

  switch(parameter.type) {

    YAC_UNREACHABLE_DEFAULT(
      "ERROR(yaml_interp_method_parameter_free): "
      "unsupported parameter type");

    // parameter is a sequence
    case (SEQ_PARAM): {

      for (size_t seq_idx = 0;
            seq_idx < parameter_value.data.seq.count; ++seq_idx) {
        yaml_interp_method_parameter_free(
          parameter_value.data.seq.values[seq_idx],
          *parameter.data.seq_param.sub_param);
      }
      free(parameter_value.data.seq.values);

      break;
    }

    // parameter is a map
    case (MAP_PARAM): {

      for (size_t param_idx = 0;
           param_idx < parameter.data.map_param.num_sub_params; ++param_idx)
        yaml_interp_method_parameter_free(
          parameter_value.data.map_values[param_idx],
          parameter.data.map_param.sub_params[param_idx]);
      free(parameter_value.data.map_values);
      break;
    }

    // parameter is a scalar
    case (ENUM_PARAM):
    case (INT_PARAM):
    case (DBLE_PARAM):
    case (BOOL_PARAM):
    case (STR_PARAM):
    case (DEG_PARAM):
      // nothing to be done
      break;
  }
}

static void yaml_parse_interp_method_parameter_value(
  struct interp_method_parameter parameter,
  interp_method_parameter_value * parameter_value,
  fy_node_t value_node, char const * interpolation_name,
  char const * yaml_filename);

/**
 * Parses a YAML map parameter node
 * @param[out]    map_parameter        parameter map configuration
 * @param[in,out] map_parameter_values values of map parameters
 * @param[in]     value_node           YAML containing the map data
 * @param[in]     interpolation_name   name of the interpolation
 *                                     (used for error messages)
 * @param[in]     yaml_filename        name of the YAML file
 *                                     (used for error messages)
 */
static void yaml_parse_interp_method_map_parameter(
  struct interp_method_parameter map_parameter,
  interp_method_parameter_value * map_parameter_values,
  fy_node_t value_node, char const * interpolation_name,
  char const * yaml_filename) {
  char const * routine_name = "yaml_parse_interp_method_map_parameter";

  // the parameter value node has to be a map
  YAML_ASSERT_F(
    !value_node || fy_node_is_mapping(value_node),
    "invalid node type for parameter \"%s\" of interpolation method "
    "\"%s\" (has to be a map)", map_parameter.name, interpolation_name);

  struct interp_method_parameter const * sub_parameters =
    map_parameter.data.map_param.sub_params;
  size_t num_sub_parameters = map_parameter.data.map_param.num_sub_params;

  // for all key value pairs of the map
  void * iter = NULL;
  fy_node_pair_t parameter_pair;
  while ((parameter_pair = fy_node_mapping_iterate(value_node, &iter))) {

    // get the key of the pair and match it with the valid key of the map
    char const * sub_parameter_name =
      yaml_parse_string_value(
        fy_node_pair_key(parameter_pair),
        "interpolation method parameter name", yaml_filename);
    size_t sub_param_idx = SIZE_MAX;
    for (size_t i = 0;
         (i < num_sub_parameters) && (sub_param_idx == SIZE_MAX); ++i)
      if (!strcmp(sub_parameter_name, sub_parameters[i].name))
        sub_param_idx = i;
    YAML_ASSERT_F(
      sub_param_idx != SIZE_MAX,
      "\"%s\" is not a valid parameter for interpolation method \"%s\"",
      sub_parameter_name, interpolation_name);

    fy_node_t parameter_value_node = fy_node_pair_value(parameter_pair);

    yaml_parse_interp_method_parameter_value(
      sub_parameters[sub_param_idx], &map_parameter_values[sub_param_idx],
      parameter_value_node, interpolation_name, yaml_filename);
  }
}

/**
 * Parses a YAML sequence parameter node
 * @param[in]  seq_parameter        sequence parameter configuration
 * @param[out] seq_parameter_values value of the sequence
 * @param[out] seq_parameter_count  number of entries in the sequence
 * @param[in]  value_node           YAML node containing the sequence
 * @param[in]  interpolation_name   name of the interpolation
 *                                  (used for error messages)
 * @param[in]  yaml_filename        name of the YAML file
 *                                  (used for error messages)
 */
static void yaml_parse_interp_method_seq_parameter(
  struct interp_method_parameter seq_parameter,
  interp_method_parameter_value ** seq_parameter_values,
  size_t * seq_parameter_count,
  fy_node_t value_node, char const * interpolation_name,
  char const * yaml_filename) {
  char const * routine_name = "yaml_parse_interp_method_seq_parameter";

  // the parameter value node has to be a sequence
  YAML_ASSERT_F(
    fy_node_is_sequence(value_node),
    "invalid node type for parameter \"%s\" of interpolation method "
    "\"%s\" (has to be a sequence)", seq_parameter.name, interpolation_name);

  // get the number of item in the sequence and allocate sub parameter
  // values accordingly
  *seq_parameter_count =
    (size_t)fy_node_sequence_item_count(value_node);
  *seq_parameter_values =
    xmalloc(*seq_parameter_count * sizeof(**seq_parameter_values));

  // parse the sequence
  for (size_t seq_idx = 0; seq_idx < *seq_parameter_count; ++seq_idx) {
    (*seq_parameter_values)[seq_idx] =
      yaml_interp_method_parameter_get_default(
        *seq_parameter.data.seq_param.sub_param);
    yaml_parse_interp_method_parameter_value(
      *seq_parameter.data.seq_param.sub_param,
      (*seq_parameter_values) + seq_idx,
      fy_node_sequence_get_by_index(
        value_node, (int)seq_idx),
      interpolation_name, yaml_filename);
  }
}

/**
 * reads a user-provided parameter value. Afterwards this value is converted
 * to the associated data-type and it is check for validity.
 * @param[in]  parameter       configuration of the parameter
 * @param[out] parameter_value parsed user-provided parameter value
 * @param[in]  value_node      YAML node containing the user-provoided paramter
 *                             value
 * @param[in]  interpolation_name name of the interpolation
 *                                (used for error messages)
 * @param[in]  yaml_filename      name of the user-provided yaml file
 *                                (used for error messages)
 */
static void yaml_parse_interp_method_parameter_value(
  struct interp_method_parameter parameter,
  interp_method_parameter_value * parameter_value,
  fy_node_t value_node, char const * interpolation_name,
  char const * yaml_filename) {
  char const * routine_name = "yaml_parse_interp_method_parameter_value";

  // ensures that the parameter is not yet set
  YAML_ASSERT_F(
    !parameter_value->is_set,
    "\"%s\" parameter of interpolation method \"%s\" has already been set",
    parameter.name, interpolation_name);

  // check the type of parameter
  switch(parameter.type) {
    YAML_UNREACHABLE_DEFAULT_F(
      "unsupported type parameter \"%s\" of interpolation method \"%s\"",
      parameter.name, interpolation_name);
    case (ENUM_PARAM):
      parameter_value->data.enum_value =
        yaml_parse_enum_value(
          parameter.data.enum_param.valid_values,
          parameter.data.enum_param.num_valid_values,
          value_node, "interpolation method enum parameter value",
          yaml_filename);
      break;
    case (INT_PARAM):
      parameter_value->data.int_value =
        yaml_parse_integer_value(
          value_node, "interpolation method integer parameter value",
          yaml_filename);
      YAML_ASSERT_F(
        (parameter_value->data.int_value >= parameter.data.int_param.valid_min) &&
        (parameter_value->data.int_value <= parameter.data.int_param.valid_max),
        "\"%d\" is not a valid integer parameter value for parameter \"%s\" "
        "of interpolation method \"%s\" "
        "(valid range: %d <= value <= %d)",
        parameter_value->data.int_value, parameter.name, interpolation_name,
        parameter.data.int_param.valid_min,
        parameter.data.int_param.valid_max);
      break;
    case (DBLE_PARAM):
      parameter_value->data.dble_value =
        yaml_parse_double_value(
          value_node, "interpolation method double parameter value",
          yaml_filename);
      YAML_ASSERT_F(
        (parameter_value->data.dble_value >=
         parameter.data.dble_param.valid_min) &&
        (parameter_value->data.dble_value <=
         parameter.data.dble_param.valid_max),
        "\"%lf\" is not a valid double parameter value for parameter \"%s\" "
        "of interpolation method \"%s\" "
        "(valid range: %e <= value <= %e)",
        parameter_value->data.dble_value, parameter.name, interpolation_name,
        parameter.data.dble_param.valid_min,
        parameter.data.dble_param.valid_max);
      break;
    case (DEG_PARAM):
      parameter_value->data.dble_value =
        yaml_parse_double_value(
          value_node, "interpolation method degree parameter value",
          yaml_filename);
      YAML_ASSERT_F(
        (parameter_value->data.dble_value >=
         parameter.data.dble_param.valid_min) &&
        (parameter_value->data.dble_value <=
         parameter.data.dble_param.valid_max),
        "\"%lf\" is not a valid degree parameter value for parameter \"%s\" "
        "of interpolation method \"%s\" "
        "(valid range: %e <= value <= %e)",
        parameter_value->data.dble_value, parameter.name, interpolation_name,
        parameter.data.dble_param.valid_min,
        parameter.data.dble_param.valid_max);
      parameter_value->data.dble_value *= YAC_RAD;
      break;
    case (BOOL_PARAM):
      parameter_value->data.bool_value =
        yaml_parse_enum_value(
          bool_names, num_bool_names, value_node,
          "interpolation method bool parameter value", yaml_filename);
      break;
    case (STR_PARAM):
      parameter_value->data.str_value =
        yaml_parse_string_value(
          value_node, "interpolation method string parameter value",
          yaml_filename);
      YAML_ASSERT_F(
        strlen(parameter_value->data.str_value) <
          parameter.data.str_param.max_str_len,
        "\"%s\" is not a valid string parameter value for parameter \"%s\" "
        "of interpolation method \"%s\" "
        "(maximum string length: %d)",
        parameter_value->data.str_value, parameter.name, interpolation_name,
        (int)(parameter.data.str_param.max_str_len - 1));
      break;
    case (MAP_PARAM): {
      yaml_parse_interp_method_map_parameter(
        parameter, parameter_value->data.map_values,
        value_node, interpolation_name, yaml_filename);
      break;
    }
    case (SEQ_PARAM): {
      yaml_parse_interp_method_seq_parameter(
        parameter, &parameter_value->data.seq.values,
        &parameter_value->data.seq.count,
        value_node, interpolation_name, yaml_filename);
    }
  };

  parameter_value->is_set = 1;
}

static void yaml_parse_interp_method(
  struct yac_interp_stack_config * interp_stack,
  fy_node_t interp_method_node, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_interp_method";

  char const * interpolation_type_str;
  fy_node_t parameter_node;

  // get name of interpolations
  yaml_parse_base_interp_method_node(
    &interpolation_type_str, &parameter_node,
    interp_method_node, yaml_filename);

  // match interpolation name with list of valid interpolation methods
  struct yac_interpolation_method const * interp_method = NULL;
  for (int i = 0; (i < NUM_INTERPOLATION_METHODS) && (!interp_method); ++i)
    if (!strcmp(interpolation_type_str, interpolation_methods[i].name))
      interp_method = &interpolation_methods[i];

  YAML_ASSERT_F(
    interp_method,
    "\"%s\" is not a valid interpolation method",
    interpolation_type_str);

  // allocate and initialise interpolation method parameter with default values
  interp_method_parameter_value interp_method_parameter_value =
    yaml_interp_method_parameter_get_default(interp_method->parameter);
  // parse user-provided interpolation method configuration
  yaml_parse_interp_method_parameter_value(
    interp_method->parameter, &interp_method_parameter_value,
    parameter_node, interp_method->name, yaml_filename);
  // add interpolation method to stack using user-provied interpolation
  // method configuration
  interp_method->add_interpolation(
    interp_stack, interp_method_parameter_value, yaml_filename);
  // free interpolation method parameters
  yaml_interp_method_parameter_free(
    interp_method_parameter_value, interp_method->parameter);
}

static void yaml_parse_interp_stack_value(
  struct yac_interp_stack_config * interp_stack,
  fy_node_t interp_stack_node, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_interp_stack_value";

  YAML_ASSERT(
    fy_node_is_sequence(interp_stack_node),
    "unsupported interpolation stack node type"
    "(interpolation stacks are expected to be defined as a sequence)");

  // parse couplings
  void * iter = NULL;
  fy_node_t interp_stack_item;
  while ((interp_stack_item =
            fy_node_sequence_iterate(interp_stack_node, &iter)))
    yaml_parse_interp_method(
      interp_stack, interp_stack_item, yaml_filename);
}

static struct field_couple_field_names yaml_parse_field_name(
  fy_node_t field_node, const char * yaml_filename) {
  char const * routine_name = "yaml_parse_field_name";

  struct field_couple_field_names field_name;

  switch(fy_node_get_type(field_node)) {

    YAML_UNREACHABLE_DEFAULT(
      "unsupported field name node type "
      "(field name is either scalars or a map)");

    // the node contains one name for both source and target field
    case (FYNT_SCALAR): {

      field_name.src =
        ((field_name.tgt =
            yaml_parse_string_value(field_node, "field name", yaml_filename)));

      break;
    }

    // the node contains different names for the source and target field
    case (FYNT_MAPPING): {

      field_name.src =
        fy_node_mapping_lookup_scalar0_by_simple_key(
          field_node, "src", (size_t)-1);
      field_name.tgt =
        fy_node_mapping_lookup_scalar0_by_simple_key(
          field_node, "tgt", (size_t)-1);

      YAML_ASSERT(
        field_name.src && field_name.tgt &&
        (fy_node_mapping_item_count(field_node) == 2),
        "invalid field name mapping node "
        "(field name mapping node has to contain two maps "
        "with the keys \"src\" and \"tgt\")")

      break;
    }
  }

  return field_name;
}

static void yaml_parse_string_sequence(
  char const *** values, size_t * num_values,
  fy_node_t values_node, char const * sequence_name,
  const char * yaml_filename) {
  char const * routine_name = "yaml_parse_string_sequence";

  YAML_ASSERT_F(
    (*values == NULL) && (*num_values == 0),
    "values have already been set for sequence \"%s\"",
    sequence_name);

  // if the field node contains multiple fields
  if (fy_node_is_sequence(values_node)) {

    *num_values = (size_t)fy_node_sequence_item_count(values_node);
    *values = xmalloc(*num_values * sizeof(**values));
    for (size_t value_idx = 0; value_idx < *num_values; ++value_idx)
      (*values)[value_idx] =
        yaml_parse_string_value(
          fy_node_sequence_get_by_index(values_node, value_idx),
          sequence_name, yaml_filename);
  } else {
    *num_values = 1;
    *values = xmalloc(sizeof(**values));
    **values =
      yaml_parse_string_value(
        values_node, sequence_name, yaml_filename);
  }
}

static void yaml_parse_field_names(
  struct field_couple_field_names ** field_names,
  size_t * num_field_names, fy_node_t fields_node,
  const char * yaml_filename) {

  // if the field node contains multiple fields
  if (fy_node_is_sequence(fields_node)) {

    size_t start_idx = *num_field_names;
    *num_field_names += (size_t)fy_node_sequence_item_count(fields_node);
    *field_names =
      xrealloc(*field_names, *num_field_names * sizeof(**field_names));
    for (size_t i = start_idx; i < *num_field_names; ++i)
      (*field_names)[i] =
        yaml_parse_field_name(
          fy_node_sequence_get_by_index(fields_node, i), yaml_filename);
  } else {
    ++*num_field_names;
    *field_names =
      xrealloc(*field_names, *num_field_names * sizeof(**field_names));
    (*field_names)[*num_field_names-1] =
      yaml_parse_field_name(fields_node, yaml_filename);
  }
}

static void yaml_parse_comp_grid_names(
  char const ** comp_name, char const *** grid_names, size_t * num_grid_names,
  fy_node_t values_node, char const * type_name, const char * yaml_filename) {
  char const * routine_name = "yaml_parse_comp_grids_names";

  YAML_ASSERT_F(
    *comp_name == NULL, "%s component name already set", type_name);
  YAML_ASSERT_F(
    (*grid_names == NULL) && (*num_grid_names == 0),
    "%s grid names already set", type_name);

  YAML_ASSERT(
    fy_node_is_mapping(values_node),
    "unsupported component/grid names node type (has to be a map)");

  // parse component/grid names
  void * iter = NULL;
  fy_node_pair_t pair;
  while((pair = fy_node_mapping_iterate(values_node, &iter))) {

    enum yaml_comp_grid_names_key_types comp_grid_names_key_type =
      (enum yaml_comp_grid_names_key_types)
        yaml_parse_enum_value(
            yaml_comp_grid_names_keys, num_yaml_comp_grid_names_keys,
            fy_node_pair_key(pair),
            "component/grid names parameter name", yaml_filename);
    char const * comp_grid_names_key_name =
      yac_name_type_pair_get_name(
        yaml_comp_grid_names_keys, num_yaml_comp_grid_names_keys,
        comp_grid_names_key_type);

    fy_node_t name_node = fy_node_pair_value(pair);

    switch(comp_grid_names_key_type) {
      YAML_UNREACHABLE_DEFAULT("invalid component/grid name key type");
      case(COMP_NAME):
        *comp_name =
          yaml_parse_string_value(
            name_node, comp_grid_names_key_name, yaml_filename);
        break;
      case(GRID_NAMES):
        yaml_parse_string_sequence(
          grid_names, num_grid_names,
          name_node, comp_grid_names_key_name, yaml_filename);
        break;
    }
  }
}

static void yaml_parse_weight_file_data(
  char const ** weight_file_name,
  enum yac_weight_file_on_existing * on_existing,
  fy_node_t values_node, char const * type_name, const char * yaml_filename) {
  char const * routine_name = "yaml_parse_weight_file";

  YAML_ASSERT_F(
    *weight_file_name == NULL, "%s weight file name already set", type_name);
  YAML_ASSERT_F(
    *on_existing == YAC_WEIGHT_FILE_UNDEFINED,
    "%s \"on existing\" already set ", type_name);

  YAML_ASSERT(
    fy_node_is_mapping(values_node),
    "unsupported weight file data node type (has to be a map)");

  // parse component/grid names
  void * iter = NULL;
  fy_node_pair_t pair;
  while((pair = fy_node_mapping_iterate(values_node, &iter))) {

    enum yaml_weight_file_data_key_types weight_file_data_key_type =
      (enum yaml_weight_file_data_key_types)
        yaml_parse_enum_value(
            yaml_weight_file_data_keys, num_yaml_weight_file_data_keys,
            fy_node_pair_key(pair),
            "weight file data parameter name", yaml_filename);
    char const * weight_file_data_key_name =
      yac_name_type_pair_get_name(
        yaml_weight_file_data_keys, num_yaml_weight_file_data_keys,
        weight_file_data_key_type);

    fy_node_t weight_file_data_node = fy_node_pair_value(pair);

    switch(weight_file_data_key_type) {
      YAML_UNREACHABLE_DEFAULT("invalid weight file data key type");
      case(WEIGHT_FILE_DATA_NAME):
        *weight_file_name =
          yaml_parse_string_value(
            weight_file_data_node, weight_file_data_key_name, yaml_filename);
        break;
      case(WEIGHT_FILE_DATA_ON_EXISTING):
        *on_existing =
          (enum yac_weight_file_on_existing)
            yaml_parse_enum_value(
                weight_file_on_existing_types,
                num_weight_file_on_existing_types,
                weight_file_data_node, weight_file_data_key_name,
                yaml_filename);
        break;
    }
  }
}

static void yaml_parse_couple_map_pair(
  struct field_couple_buffer * field_buffer,
  fy_node_pair_t couple_pair, const char * yaml_filename,
  enum yac_time_unit_type time_unit) {
  char const * routine_name = "yaml_parse_couple_map_pair";

  enum yaml_couple_key_types couple_key_type =
    (enum yaml_couple_key_types)
      yaml_parse_enum_value(
        yaml_couple_keys, num_yaml_couple_keys,
        fy_node_pair_key(couple_pair),
        "couple configuration parameter name", yaml_filename);
  char const * couple_key_name =
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, couple_key_type);

  fy_node_t value_node = fy_node_pair_value(couple_pair);

  switch (couple_key_type) {
    YAML_UNREACHABLE_DEFAULT("invalid couple key type");
    case (SOURCE_NAMES):
      yaml_parse_comp_grid_names(
        &(field_buffer->src.comp_name),
        &(field_buffer->src.grid.name),
        &(field_buffer->src.grid.count),
        value_node, couple_key_name, yaml_filename);
      break;
    case (SOURCE_COMPONENT):
      field_buffer->src.comp_name =
        yaml_parse_string_value(value_node, couple_key_name, yaml_filename);
      break;
    case (SOURCE_GRID):
      yaml_parse_string_sequence(
        &field_buffer->src.grid.name, &field_buffer->src.grid.count,
        value_node, couple_key_name, yaml_filename);
      break;
    case (TARGET_NAMES):
      yaml_parse_comp_grid_names(
        &(field_buffer->tgt.comp_name),
        &(field_buffer->tgt.grid.name),
        &(field_buffer->tgt.grid.count),
        value_node, couple_key_name, yaml_filename);
      break;
    case (TARGET_COMPONENT):
      field_buffer->tgt.comp_name =
        yaml_parse_string_value(value_node, couple_key_name, yaml_filename);
      break;
    case (TARGET_GRID):
      yaml_parse_string_sequence(
        &field_buffer->tgt.grid.name, &field_buffer->tgt.grid.count,
        value_node, couple_key_name, yaml_filename);
      break;
    case (FIELD):
      yaml_parse_field_names(
        &(field_buffer->field_names),
        &(field_buffer->num_field_names),
        value_node, yaml_filename);
      break;
    case (COUPLING_PERIOD):
      field_buffer->coupling_period =
        yaml_parse_timestep_value(
          value_node, couple_key_name, yaml_filename, time_unit);
      break;
    case (TIME_REDUCTION):
      field_buffer->time_reduction =
        yaml_parse_time_reduction_value(
          value_node, couple_key_name, yaml_filename);
      break;
    case (SOURCE_LAG):
      field_buffer->src.lag =
        yaml_parse_integer_value(
          value_node, couple_key_name, yaml_filename);
      break;
    case (TARGET_LAG):
      field_buffer->tgt.lag =
        yaml_parse_integer_value(
          value_node, couple_key_name, yaml_filename);
      break;
    case (WEIGHT_FILE_NAME):
      field_buffer->weight_file.name =
        yaml_parse_string_value(
          value_node, couple_key_name, yaml_filename);
      break;
    case (WEIGHT_FILE_ON_EXISTING):
      field_buffer->weight_file.on_existing =
        (enum yac_weight_file_on_existing)
          yaml_parse_enum_value(
            weight_file_on_existing_types,
            num_weight_file_on_existing_types,
            value_node, couple_key_name, yaml_filename);
      break;
    case (WEIGHT_FILE_DATA):
      yaml_parse_weight_file_data(
        &(field_buffer->weight_file.name),
        &(field_buffer->weight_file.on_existing),
        value_node, couple_key_name, yaml_filename);
      break;
    case (MAPPING_SIDE):
      field_buffer->mapping_on_source =
        yaml_parse_enum_value(
          mapping_sides, num_mapping_sides,
          value_node, couple_key_name, yaml_filename);
      break;
    case (SCALE_FACTOR):
      field_buffer->scale_factor =
        yaml_parse_double_value(
          value_node, couple_key_name, yaml_filename);
      break;
    case (SCALE_SUMMAND):
      field_buffer->scale_summand =
        yaml_parse_double_value(
          value_node, couple_key_name, yaml_filename);
      break;
    case (INTERPOLATION):
      yaml_parse_interp_stack_value(
        field_buffer->interp_stack, value_node, yaml_filename);
      break;
    case (SOURCE_MASK_NAME):
    case (SOURCE_MASK_NAMES):
      yaml_parse_string_sequence(
        &(field_buffer->src_mask_names),
        &(field_buffer->num_src_mask_names),
        value_node, couple_key_name, yaml_filename);
      break;
    case (TARGET_MASK_NAME):
      field_buffer->tgt_mask_name =
        yaml_parse_string_value(
          value_node, couple_key_name, yaml_filename);
      break;
    case (YAXT_EXCHANGER_NAME):
      field_buffer->yaxt_exchanger_name =
        yaml_parse_string_value(
          value_node, couple_key_name, yaml_filename);
      break;
    case (USE_RAW_EXCHANGE):
      field_buffer->use_raw_exchange =
        yaml_parse_enum_value(
          bool_names, num_bool_names,
          value_node, couple_key_name, yaml_filename);
      break;
  }
}

static void yaml_parse_couple(
  struct yac_couple_config * couple_config, fy_node_t couple_node,
  char const * yaml_filename, enum yac_time_unit_type time_unit) {
  char const * routine_name = "yaml_parse_couple";

  YAML_ASSERT(
    fy_node_is_mapping(couple_node),
    "unsupported couple node type "
    "(couples are expected to be defined as a mapping)");

  // initialise field configuration buffer with default values
  struct field_couple_buffer field_buffer = {
    .src.comp_name = NULL,
    .src.grid.name = NULL,
    .src.grid.count = 0,
    .tgt.comp_name = NULL,
    .tgt.grid.name = NULL,
    .tgt.grid.count = 0,
    .field_names = NULL,
    .num_field_names = 0,
    .coupling_period = NULL,
    .time_reduction = TIME_NONE,
    .interp_stack = yac_interp_stack_config_new(),
    .src.lag = 0,
    .tgt.lag = 0,
    .weight_file.name = NULL,
    .weight_file.on_existing = YAC_WEIGHT_FILE_UNDEFINED,
    .mapping_on_source = 1,
    .scale_factor = 1.0,
    .scale_summand = 0.0,
    .src_mask_names = NULL,
    .num_src_mask_names = 0,
    .tgt_mask_name = NULL,
    .yaxt_exchanger_name = NULL,
    .use_raw_exchange = 0};

  // parse couple
  void * iter = NULL;
  fy_node_pair_t pair;
  while ((pair = fy_node_mapping_iterate(couple_node, &iter)))
    yaml_parse_couple_map_pair(
      &field_buffer, pair, yaml_filename, time_unit);

  YAML_ASSERT(
    field_buffer.src.comp_name, "missing source component name");
  YAML_ASSERT_F(
    field_buffer.src.grid.count > 0,
    "missing source grid name (component \"%s\")",
    field_buffer.src.comp_name);
  for (size_t i = 0; i < field_buffer.src.grid.count; ++i)
    YAML_ASSERT_F(
      (field_buffer.src.grid.name[i] != NULL) &&
      (field_buffer.src.grid.name[i][0] != '\0'),
      "invalid source grid name (component \"%s\" grid idx %zu)",
      field_buffer.src.comp_name, i);
  YAML_ASSERT(
    field_buffer.tgt.comp_name,
    "missing target component name");
  YAML_ASSERT_F(
    field_buffer.tgt.grid.count > 0,
    "missing target grid name (component \"%s\")",
    field_buffer.tgt.comp_name);
  for (size_t i = 0; i < field_buffer.tgt.grid.count; ++i)
    YAML_ASSERT_F(
      (field_buffer.tgt.grid.name[i] != NULL) &&
      (field_buffer.tgt.grid.name[i][0] != '\0'),
      "invalid target grid name (component \"%s\" grid idx %zu)",
      field_buffer.tgt.comp_name, i);
  YAML_ASSERT_F(
    field_buffer.num_field_names > 0,
    "missing field names "
    "(source component \"%s\" source grid \"%s\" "
    "target component \"%s\" target grid \"%s\")",
    field_buffer.src.comp_name, field_buffer.src.grid.name[0],
    field_buffer.tgt.comp_name, field_buffer.tgt.grid.name[0]);
  YAML_ASSERT_F(
    field_buffer.coupling_period,
    "missing coupling period "
    "(source component \"%s\" source grid \"%s\" "
    "target component \"%s\" target grid \"%s\")",
    field_buffer.src.comp_name, field_buffer.src.grid.name[0],
    field_buffer.tgt.comp_name, field_buffer.tgt.grid.name[0]);

  for (size_t i = 0; i < field_buffer.num_field_names; ++i)
    for (size_t j = 0; j < field_buffer.src.grid.count; ++j)
      for (size_t k = 0; k < field_buffer.tgt.grid.count; ++k)
        yac_couple_config_def_couple(
          couple_config,
          field_buffer.src.comp_name,
          field_buffer.src.grid.name[j],
          field_buffer.field_names[i].src,
          field_buffer.tgt.comp_name,
          field_buffer.tgt.grid.name[k],
          field_buffer.field_names[i].tgt,
          field_buffer.coupling_period,
          field_buffer.time_reduction,
          field_buffer.interp_stack,
          field_buffer.src.lag,
          field_buffer.tgt.lag,
          field_buffer.weight_file.name,
          (field_buffer.weight_file.on_existing == YAC_WEIGHT_FILE_UNDEFINED)?
            YAC_WEIGHT_FILE_ON_EXISTING_DEFAULT_VALUE:
            field_buffer.weight_file.on_existing,
          field_buffer.mapping_on_source,
          field_buffer.scale_factor,
          field_buffer.scale_summand,
          field_buffer.num_src_mask_names,
          field_buffer.src_mask_names,
          field_buffer.tgt_mask_name,
          field_buffer.yaxt_exchanger_name,
          field_buffer.use_raw_exchange);

  // cleanup
  free((void*)field_buffer.src.grid.name);
  free((void*)field_buffer.tgt.grid.name);
  free((void*)field_buffer.field_names);
  free((void*)field_buffer.coupling_period);
  yac_interp_stack_config_delete(field_buffer.interp_stack);
  free((void*)field_buffer.src_mask_names);
}

static void yaml_parse_coupling(
  struct yac_couple_config * couple_config,
  fy_node_t coupling_node, char const * yaml_filename,
  enum yac_time_unit_type time_unit) {
  char const * routine_name = "yaml_parse_coupling";

  // check if the coupling node is empty -> nothing to be read
  if (!coupling_node) return;

  YAML_ASSERT(
    fy_node_is_sequence(coupling_node),
    "unsupported coupling node type "
    "(couplings are expected to be defined as a sequence)");

  // parse couplings
  void * iter = NULL;
  fy_node_t couple_node;
  while ((couple_node = fy_node_sequence_iterate(coupling_node, &iter)))
    yaml_parse_couple(couple_config, couple_node, yaml_filename, time_unit);
}

static struct debug_config_file yaml_parse_config_file_value(
  fy_node_t config_file_node, char const * file_type_name,
  const char * yaml_filename) {
  char const * routine_name = "yaml_parse_config_file_value";

  struct debug_config_file config_file =
    {.name = NULL, .type = YAC_TEXT_FILETYPE_YAML, .include_definitions = 0};

  char * str_buffer = xmalloc(strlen(file_type_name) + 32);

  switch(fy_node_get_type(config_file_node)) {

    YAML_UNREACHABLE_DEFAULT_F(
    "unsupported config file node type "
    "(%s is either scalar or a map)", file_type_name);

    // the node contains only the filename
    case (FYNT_SCALAR): {

      config_file.name =
        yaml_parse_string_value(
          config_file_node,
          strcat(strcpy(str_buffer, file_type_name), " name"),
          yaml_filename);

      break;
    }

    // the node contains the name and the type
    case (FYNT_MAPPING): {

      fy_node_t filename_node =
        fy_node_mapping_lookup_by_string(
          config_file_node, "filename", (size_t)-1);
      fy_node_t filetype_node =
        fy_node_mapping_lookup_by_string(
          config_file_node, "filetype", (size_t)-1);
      fy_node_t include_definitions_node =
        fy_node_mapping_lookup_by_string(
          config_file_node, "include_definitions", (size_t)-1);

      YAML_ASSERT_F(
        filename_node,
        "invalid %s mapping node "
        "(global config file mapping node has to include a map "
        "with the keys \"filename\")", file_type_name)

      config_file.name =
        yaml_parse_string_value(
          filename_node, strcat(strcpy(str_buffer, file_type_name), " name"),
          yaml_filename);
      config_file.type =
        (filetype_node)?
          (enum yac_text_filetype)
            yaml_parse_enum_value(
              config_filetypes, num_config_filetypes,
              filetype_node, strcat(strcpy(str_buffer, file_type_name), " type"),
              yaml_filename):YAC_TEXT_FILETYPE_YAML;
      config_file.include_definitions =
        (include_definitions_node)?
          yaml_parse_enum_value(
            bool_names, num_bool_names,
            include_definitions_node,
            strcat(strcpy(str_buffer, file_type_name), " include definitions"),
            yaml_filename):0;
      break;
    }
  }

  YAML_ASSERT_F(
    config_file.name, "missing filename for %s", file_type_name);

  free(str_buffer);

  return config_file;
}

static void yaml_parse_debug_config_file_map_pair(
  struct debug_config_file_buffer * config_file_buffer,
  fy_node_pair_t config_file_pair, char const * config_file_type_name,
  const char * yaml_filename) {

  enum yaml_debug_sync_loc_key_types
    debug_sync_loc_key_type =
      (enum yaml_debug_sync_loc_key_types)
        yaml_parse_enum_value(
          yaml_debug_sync_loc_keys,
          num_yaml_debug_sync_loc_keys,
          fy_node_pair_key(config_file_pair),
          "config synchronisation location parameter name",
          yaml_filename);
  char const * debug_sync_loc_key_name =
    yac_name_type_pair_get_name(
      yaml_debug_sync_loc_keys, num_yaml_debug_sync_loc_keys,
      debug_sync_loc_key_type);

  fy_node_t value_node = fy_node_pair_value(config_file_pair);

  char config_file_type_name_sync[
    strlen(config_file_type_name) + strlen(debug_sync_loc_key_name) + 8];
  sprintf(
    config_file_type_name_sync, "%s (%s)",
    config_file_type_name, debug_sync_loc_key_name);

  config_file_buffer->config_file[debug_sync_loc_key_type] =
    yaml_parse_config_file_value(
      value_node, config_file_type_name_sync, yaml_filename);
}

static struct debug_config_file_buffer yaml_parse_debug_config_file_buffer(
  fy_node_t config_file_node, char const * config_file_type_name,
  char const * yaml_filename) {
  char const * routine_name = "yaml_parse_debug_config_file_buffer";

  YAML_ASSERT_F(
    fy_node_is_mapping(config_file_node),
    "unsupported %s node type "
    "(%s is expected to be defined as a mapping)",
    config_file_type_name, config_file_type_name);

  struct debug_config_file_buffer config_file_buffer;
  config_file_buffer.sync_loc_ref[SYNC_LOC_DEF_COMP] =
    YAC_INSTANCE_CONFIG_OUTPUT_REF_COMP;
  config_file_buffer.sync_loc_ref[SYNC_LOC_SYNC_DEF] =
    YAC_INSTANCE_CONFIG_OUTPUT_REF_SYNC;
  config_file_buffer.sync_loc_ref[SYNC_LOC_ENDDEF] =
    YAC_INSTANCE_CONFIG_OUTPUT_REF_ENDDEF;
  for (int i = 0; i < SYNC_LOC_COUNT; ++i) {
    config_file_buffer.config_file[i].name = NULL;
    config_file_buffer.config_file[i].type = YAC_TEXT_FILETYPE_YAML;
    config_file_buffer.config_file[i].include_definitions = 0;
    YAML_ASSERT_F(
      config_file_buffer.sync_loc_ref[i] != NULL,
      "invalid unsupported synchronisation location (%d) for %s",
      i, config_file_type_name);
  }

  // parse couplings
  void * iter = NULL;
  fy_node_pair_t pair;
  while ((pair = fy_node_mapping_iterate(config_file_node, &iter)))
    yaml_parse_debug_config_file_map_pair(
      &config_file_buffer, pair, config_file_type_name, yaml_filename);

  return config_file_buffer;
}

static void yaml_parse_debug_global_config(
  struct yac_couple_config * couple_config, fy_node_t global_config_node,
  char const * yaml_filename) {
  char const * routine_name = "yaml_parse_debug_global_config";

  char const * config_file_type_name = "debug global config file";

  YAML_ASSERT_F(
    fy_node_is_mapping(global_config_node),
    "unsupported %s node type "
    "(%s is expected to be defined as a mapping)",
    config_file_type_name, config_file_type_name);

  struct debug_config_file_buffer global_config_buffer =
    yaml_parse_debug_config_file_buffer(
      global_config_node, config_file_type_name, yaml_filename);

  for (int i = 0; i < SYNC_LOC_COUNT; ++i)
    if (global_config_buffer.config_file[i].name != NULL)
      yac_couple_config_set_config_output_filename(
        couple_config, global_config_buffer.config_file[i].name,
        global_config_buffer.config_file[i].type,
        global_config_buffer.sync_loc_ref[i],
        global_config_buffer.config_file[i].include_definitions);
}

static void yaml_parse_output_grid_pair(
  char const ** grid_name, char const ** file_name,
  fy_node_pair_t output_grid_pair, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_output_grid_pair";

  enum yaml_debug_output_grid_key_types output_grid_key_type =
    (enum yaml_debug_output_grid_key_types)
      yaml_parse_enum_value(
        yaml_debug_output_grid_keys, num_yaml_debug_output_grid_keys,
        fy_node_pair_key(output_grid_pair),
        "output grid parameter name", yaml_filename);
  char const * debug_output_key_name =
    yac_name_type_pair_get_name(
      yaml_debug_output_grid_keys, num_yaml_debug_output_grid_keys,
      output_grid_key_type);

  fy_node_t value_node = fy_node_pair_value(output_grid_pair);

  switch(output_grid_key_type) {
    YAML_UNREACHABLE_DEFAULT("invalid output grid key type");
    case(OUTPUT_GRID_GRID_NAME): {
      *grid_name =
        yaml_parse_string_value(value_node, debug_output_key_name, yaml_filename);
      break;
    }
    case(OUTPUT_GRID_FILE_NAME): {
      *file_name =
        yaml_parse_string_value(value_node, debug_output_key_name, yaml_filename);
      break;
    }
  };
}

static void yaml_parse_output_grid(
  struct yac_couple_config * couple_config, fy_node_t output_grid_node,
  char const * yaml_filename) {
  char const * routine_name = "yaml_parse_output_grid";

  YAML_ASSERT(
    fy_node_is_mapping(output_grid_node),
    "unsupported output grid node type "
    "(output grids are expected to be defined as a mapping)");

  char const * grid_name = NULL;
  char const * file_name = NULL;

  // parse output grid
  void * iter = NULL;
  fy_node_pair_t pair;
  while ((pair = fy_node_mapping_iterate(output_grid_node, &iter)))
    yaml_parse_output_grid_pair(&grid_name, &file_name, pair, yaml_filename);

  YAML_ASSERT(grid_name, "missing grid name");
  YAML_ASSERT_F(file_name, "missing file name for grid \"%s\"", grid_name);

  yac_couple_config_add_grid(couple_config, grid_name);
  yac_couple_config_grid_set_output_filename(
    couple_config, grid_name, file_name);
}

static void yaml_parse_debug_output_grids(
  struct yac_couple_config * couple_config, fy_node_t output_grids_node,
  char const * yaml_filename) {
  char const * routine_name = "yaml_parse_debug_output_grids";

  YAML_ASSERT(
    fy_node_is_sequence(output_grids_node),
    "unsupported debug output grids node type "
    "(debug output grids is expected to be defined as a sequence)");

  // parse output grids
  void * iter = NULL;
  fy_node_t output_grid_node;
  while ((output_grid_node =
            fy_node_sequence_iterate(output_grids_node, &iter)))
    yaml_parse_output_grid(
      couple_config, output_grid_node, yaml_filename);
}

static void yaml_parse_debug_map_pair(
  struct yac_couple_config * couple_config, fy_node_pair_t debug_pair,
  const char * yaml_filename) {
  char const * routine_name = "yaml_parse_debug_map_pair";

  enum yaml_debug_key_types debug_key_type =
    (enum yaml_debug_key_types)
      yaml_parse_enum_value(
        yaml_debug_keys, num_yaml_debug_keys,
        fy_node_pair_key(debug_pair),
        "debug configuration parameter name", yaml_filename);

  fy_node_t value_node = fy_node_pair_value(debug_pair);

  switch (debug_key_type) {
    YAML_UNREACHABLE_DEFAULT("invalid debug_key_type");
    case(GLOBAL_CONFIG):
      yaml_parse_debug_global_config(
        couple_config, value_node, yaml_filename);
      break;
    case(OUTPUT_GRIDS):
      yaml_parse_debug_output_grids(
        couple_config, value_node, yaml_filename);
      break;
    case(MISSING_DEF):
      yac_couple_config_set_missing_definition_is_fatal(
        couple_config,
        yaml_parse_enum_value(
          bool_names, num_bool_names, value_node,
          "\"missing definition is fatal\" bool value", yaml_filename));
    break;
  };
}

static void yaml_parse_debug(
  struct yac_couple_config * couple_config,
  fy_node_t debug_node, char const * yaml_filename) {
  char const * routine_name = "yaml_parse_debug";

  // check if the debug node is empty -> nothing to be read
  if (!debug_node) return;

  YAML_ASSERT(
    fy_node_is_mapping(debug_node),
    "unsupported debug node type "
    "(debug is expected to be defined as a mapping)");

  // parse couplings
  void * iter = NULL;
  fy_node_pair_t pair;
  while ((pair = fy_node_mapping_iterate(debug_node, &iter)))
    yaml_parse_debug_map_pair(couple_config, pair, yaml_filename);
}

static void yaml_parse_base_map_pair(
  struct yac_couple_config * couple_config, fy_node_pair_t base_pair,
  const char * yaml_filename, enum yac_time_unit_type * time_unit,
  char const ** start_datetime, char const ** end_datetime) {
  char const * routine_name = "yaml_parse_base_map_pair";

  fy_node_t key_node = fy_node_pair_key(base_pair);

  YAML_ASSERT(
    fy_node_is_scalar(key_node),
    "unsupported key node type "
    "(key nodes are expected to be scalar nodes)");

  char const * base_key_name = fy_node_get_scalar0(key_node);
  int base_key_type =
    yac_name_type_pair_get_type(
      yaml_base_keys, num_yaml_base_keys, base_key_name);

  fy_node_t value_node = fy_node_pair_value(base_pair);

  switch (base_key_type) {
    case (START_DATE):
      *start_datetime =
        yaml_parse_string_value(
          value_node, base_key_name, yaml_filename);
      break;
    case (END_DATE):
      *end_datetime =
        yaml_parse_string_value(
          value_node, base_key_name, yaml_filename);
      break;
    case (CALENDAR):
      yac_cdef_calendar(
        (int)yaml_parse_calendar_value(
          value_node, base_key_name, yaml_filename));
      break;
    case (TIMESTEP_UNIT): {
      enum yac_time_unit_type time_unit_ =
        yaml_parse_timestep_unit_value(
          value_node, base_key_name, yaml_filename);
      YAML_ASSERT(
        (*time_unit == TIME_UNIT_UNDEFINED) || (*time_unit == time_unit_),
        "inconsistent redefinition of time unit")
      *time_unit = time_unit_;
      break;
    }
    case (COUPLING):
      YAML_ASSERT(
        (*time_unit != TIME_UNIT_UNDEFINED),
        "time unit has to be defined before the couplings")
      yaml_parse_coupling(
        couple_config, value_node, yaml_filename, *time_unit);
      break;
    case (DEBUG):
      yaml_parse_debug(
        couple_config, value_node, yaml_filename);
    default:
      // nothing to be done
      break;
  }
}

static void yaml_parse_document(
  struct yac_couple_config * couple_config, fy_document_t document,
  const char * yaml_filename) {
  char const * routine_name = "yaml_parse_document";

  // get root node of document
  fy_node_t root_node = fy_document_root(document);

  // if the configuration file is empty
  if (!root_node) return;

  YAML_ASSERT(
    fy_node_is_mapping(root_node),
    "unsupported root node type (root node is expected to be a mapping node)");

  char const * start_datetime = NULL;
  char const * end_datetime = NULL;

  // parse base root mappings
  enum yac_time_unit_type time_unit = TIME_UNIT_UNDEFINED;
  void * iter = NULL;
  fy_node_pair_t pair;
  while ((pair = fy_node_mapping_iterate(root_node, &iter)))
    yaml_parse_base_map_pair(
      couple_config, pair, yaml_filename, &time_unit,
      &start_datetime, &end_datetime);

  if ((start_datetime != NULL) || (end_datetime != NULL)) {

    YAC_ASSERT(
      yac_cget_calendar() != YAC_CALENDAR_NOT_SET,
      "ERROR(yaml_parse_document): "
      "cannot set start/end datetime because calendar has not yet been set");

    yac_couple_config_set_datetime(
      couple_config, start_datetime, end_datetime);
  }
}

void yac_yaml_read_coupling(
  struct yac_couple_config * couple_config, const char * yaml_filename,
  int parse_flags) {
  char const * routine_name = "yac_yaml_read_coupling";

  // check whether the yaml configuration file exists
  YAC_ASSERT_F(
    yac_file_exists(yaml_filename),
    "ERROR(%s): YAML configuration file could not be found \"%s\"",
    routine_name, yaml_filename);

  // open yaml configuration file
  FILE * config_file = xfopen(yaml_filename, "r");
  YAC_ASSERT_F(
    config_file,
    "ERROR(%s): could not open YAML configuration file \"%s\"",
    routine_name, yaml_filename);

  // parse yaml configuration file into document
  struct fy_parse_cfg parse_config =
    {.search_path = ".",
     .userdata = NULL,
     .diag = NULL};
  parse_config.flags =
    (enum fy_parse_cfg_flags)parse_flags;

  fy_document_t document =
    fy_document_build_from_fp(&parse_config, config_file);
  YAC_ASSERT_F(
    document,
    "ERROR(%s): could not parse YAML configuration file \"%s\"",
    routine_name, yaml_filename);

  // resolve anchors and merge keys
  YAML_ASSERT(
    !fy_document_resolve(document),
    "could not resolve anchors and merge keys");

  // parse document into couple configuration
  yaml_parse_document(couple_config, document, yaml_filename);

  // cleanup
  fy_document_destroy(document);
  xfclose(config_file);

  return;
}

struct yac_interp_stack_config *
  yac_yaml_parse_interp_stack_config_string(
    char const * str_interp_stack_config, int parse_flags) {
  char const * routine_name =
    "yac_yaml_parse_interp_stack_config_string";
  char const * yaml_filename =
    "user-provided interp stack config string";

  YAML_ASSERT(
    str_interp_stack_config != NULL, "interpolation stack string is NULL");

  // parse string into document
  struct fy_parse_cfg parse_config =
    {.search_path = ".",
     .userdata = NULL,
     .diag = NULL};
  parse_config.flags =
    (enum fy_parse_cfg_flags)parse_flags;

  fy_document_t document =
    fy_document_build_from_string(
      &parse_config, str_interp_stack_config, (size_t)-1);
  YAML_ASSERT(document, "failed parsing");

  // resolve anchors and merge keys
  YAML_ASSERT(
    !fy_document_resolve(document),
    "could not resolve anchors and merge keys");

  fy_node_t interp_stack_config_node = fy_document_root(document);
  YAML_ASSERT(interp_stack_config_node, "invalid root node");

  struct yac_interp_stack_config * interp_stack_config =
    yac_interp_stack_config_new();

  yaml_parse_interp_stack_value(
    interp_stack_config, interp_stack_config_node,
    yaml_filename);

  // cleanup
  fy_document_destroy(document);

  return interp_stack_config;
}

static fy_node_t yac_yaml_create_scalar(
  fy_document_t document, char const * value) {

  // return NULL if value is empty
  if (!value) return (fy_node_t)NULL;

  fy_node_t scalar_node =
    fy_node_create_scalar_copy(document, value, strlen(value));
  YAC_ASSERT(
    scalar_node,
    "ERROR(yac_yaml_create_scalar): failed to create scalar node");

  return scalar_node;
}

static fy_node_t yac_yaml_create_scalar_int(
  fy_document_t document, int value) {

  char str_value[16];
  int value_size = snprintf(str_value, sizeof(str_value), "%d", value);
  YAC_ASSERT_F(
    (value_size >= 0) && ((size_t)value_size < sizeof(str_value)),
    "ERROR(yac_yaml_create_scalar_int): "
    "could not write \"%d\" to string buffer of size %zu",
    value, sizeof(str_value));

  return yac_yaml_create_scalar(document, str_value);
}

static fy_node_t yac_yaml_create_scalar_dble(
  fy_document_t document, double value) {

  char str_value[32];
  int value_size = snprintf(str_value, sizeof(str_value), "%g", value);
  YAC_ASSERT_F(
    (value_size >= 0) && ((size_t)value_size < sizeof(str_value)),
    "ERROR(yac_yaml_create_scalar_dble): "
    "could not write \"%g\" to string buffer of size %zu",
    value, sizeof(str_value));

  return yac_yaml_create_scalar(document, str_value);
}

static fy_node_t yac_yaml_create_sequence_scalar(
  fy_document_t document, char const * const * values, size_t num_values) {

  // return NULL if sequence is empty
  if (num_values == 0) return (fy_node_t)NULL;

  YAC_ASSERT(
    values, "ERROR(yac_yaml_create_sequence_scalar): no values provided");

  // create sequence node
  fy_node_t sequence_node = fy_node_create_sequence(document);
  YAC_ASSERT(
    sequence_node, "ERROR(yac_yaml_create_sequence_scalar): "
    "failed to create sequence node");

  for (size_t value_idx = 0; value_idx < num_values; ++value_idx) {
    char const * value = values[value_idx];
    YAC_ASSERT_F(
      value, "ERROR(yac_yaml_create_sequence_scalar): "
      "invalid value at idx %zu", value_idx);
    int appending_failed =
      fy_node_sequence_append(
        sequence_node, yac_yaml_create_scalar(document, value));
    YAC_ASSERT(
      !appending_failed, "ERROR(yac_yaml_create_sequence_scalar): "
      "failed to append interpolation node");
  }

  return sequence_node;
}

static void yac_yaml_map_append(
  fy_node_t map, char const * key, fy_node_t value) {

  // if the value node is empty, return
  if (!value) return;

  YAC_ASSERT(
    key, "ERROR(yac_yaml_map_append): NULL key is not supported");

  fy_document_t document = fy_node_document(map);
  YAC_ASSERT(
    document,
    "ERROR(yac_yaml_map_append): failed to get document from node");

  // set key and value for root node
  int appending_failed =
    fy_node_mapping_append(
      map, yac_yaml_create_scalar(document, key), value);
  YAC_ASSERT(
    !appending_failed,
    "ERROR(yac_yaml_map_append): failed to append mapping node pair");
}

static void yac_yaml_map_append_scalar(
  fy_node_t map, char const * key, char const * value) {

  fy_document_t document = fy_node_document(map);
  YAC_ASSERT(
    document, "ERROR(yac_yaml_map_append_scalar): "
    "failed to get document from node");

  yac_yaml_map_append(map, key, yac_yaml_create_scalar(document, value));
}

static void yac_yaml_map_append_scalar_int(
  fy_node_t map, char const * key, int value) {

  fy_document_t document = fy_node_document(map);
  YAC_ASSERT(
    document, "ERROR(yac_yaml_map_append_scalar_int): "
    "failed to get document from node");

  yac_yaml_map_append(
    map, key, yac_yaml_create_scalar_int(document, value));
}

static void yac_yaml_map_append_scalar_dble(
  fy_node_t map, char const * key, double value) {

  fy_document_t document = fy_node_document(map);
  YAC_ASSERT(
    document, "ERROR(yac_yaml_map_append_scalar_dble): "
    "failed to get document from node");

  yac_yaml_map_append(
    map, key, yac_yaml_create_scalar_dble(document, value));
}

static fy_node_t yac_yaml_create_field_name_node(
  fy_document_t document, struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx) {

  const char * src_field_name;
  const char * tgt_field_name;
  yac_couple_config_get_field_names(
    couple_config, couple_idx, field_couple_idx,
    &src_field_name, &tgt_field_name);

  // if both names are identical
  if (!strcmp(src_field_name, tgt_field_name))
    return yac_yaml_create_scalar(document, src_field_name);

  // create field name node
  fy_node_t field_name_node = fy_node_create_mapping(document);
  YAC_ASSERT(
    field_name_node, "ERROR(yac_yaml_create_field_name_node): "
    "failed to create mapping node");

  // add source field name
  yac_yaml_map_append_scalar(field_name_node, "src", src_field_name);

  // add target field name
  yac_yaml_map_append_scalar(field_name_node, "tgt", tgt_field_name);

  return field_name_node;
}

static fy_node_t yac_yaml_create_scalar_parameter_node(
  fy_document_t document, struct interp_method_parameter parameter,
  interp_method_parameter_value parameter_value) {

  fy_node_t parameter_node = NULL;

  switch (parameter.type) {
    YAC_UNREACHABLE_DEFAULT(
      "ERROR(yac_yaml_create_scalar_parameter_node): "
      "unsupported parameter type");
    case (ENUM_PARAM): {
      parameter_node =
        yac_yaml_create_scalar(
          document,
          yac_name_type_pair_get_name(
            parameter.data.enum_param.valid_values,
            parameter.data.enum_param.num_valid_values,
            parameter_value.data.enum_value));
      break;
    }
    case (INT_PARAM): {
      parameter_node =
        yac_yaml_create_scalar_int(document, parameter_value.data.int_value);
      break;
    }
    case (DEG_PARAM):
    case (DBLE_PARAM): {
      parameter_node =
        yac_yaml_create_scalar_dble(document, parameter_value.data.dble_value);
      break;
    }
    case (BOOL_PARAM): {
      parameter_node =
        yac_yaml_create_scalar(
          document,
          yac_name_type_pair_get_name(
            bool_names, num_bool_names, parameter_value.data.bool_value));
      break;
    }
    case (STR_PARAM): {
      parameter_node =
        yac_yaml_create_scalar(document, parameter_value.data.str_value);
      break;
    }
  };

  return parameter_node;
}

static void yac_yaml_map_append_parameter(
  fy_document_t document, fy_node_t map,
  struct interp_method_parameter parameter,
  interp_method_parameter_value parameter_value) {

  yac_yaml_map_append(
    map, parameter.name,
    yac_yaml_create_scalar_parameter_node(
      document, parameter, parameter_value));
}

static int compare_parameter_values(
  interp_method_parameter_value const * a,
  interp_method_parameter_value const * b,
  struct interp_method_parameter parameter) {

  int ret = 0;
  switch (parameter.type) {
    YAC_UNREACHABLE_DEFAULT(
      "ERROR(compare_parameter_values): "
      "invalid interpolation method parameter value type");
    case (ENUM_PARAM):
      ret = (a->data.enum_value > b->data.enum_value) -
            (a->data.enum_value < b->data.enum_value);
      break;
    case (INT_PARAM):
      ret = (a->data.int_value > b->data.int_value) -
            (a->data.int_value < b->data.int_value);
      break;
    case (DEG_PARAM):
    case (DBLE_PARAM):
      ret = (a->data.dble_value > b->data.dble_value) -
            (a->data.dble_value < b->data.dble_value);
      break;
    case (BOOL_PARAM):
      ret = (a->data.bool_value > b->data.bool_value) -
            (a->data.bool_value < b->data.bool_value);
      break;
    case (STR_PARAM):
      if ((a->data.str_value != NULL) && (b->data.str_value != NULL))
        ret = strcmp(a->data.str_value, b->data.str_value);
      else
        ret = (a->data.str_value != NULL) - (b->data.str_value != NULL);
      break;
  }
  return ret;
}

static fy_node_t yac_yaml_create_non_default_map_parameter_node(
  fy_document_t document,
  interp_method_parameter_value const * map_parameter_values,
  struct interp_method_parameter parameter);

static fy_node_t yac_yaml_create_seq_parameter_node(
  fy_document_t document,
  interp_method_parameter_value const * seq_parameter_values,
  size_t seq_parameter_values_count,
  struct interp_method_parameter parameter) {

  YAC_ASSERT_F(
    parameter.type == SEQ_PARAM,
    "ERROR(yac_yaml_create_seq_parameter_node): "
    "parameter \"%s\" is not a map", parameter.name);

  // create parameter node
  fy_node_t sequence_node =
    (seq_parameter_values_count > 0)?fy_node_create_sequence(document):NULL;
  YAC_ASSERT(
    sequence_node || (seq_parameter_values_count == 0),
    "ERROR(yac_yaml_create_seq_parameter_node): "
    "failed to create sequence node");

  struct interp_method_parameter sub_parameter =
    *parameter.data.seq_param.sub_param;

  for (size_t seq_idx = 0; seq_idx < seq_parameter_values_count; ++seq_idx) {

    fy_node_t sub_parameter_node = NULL;
    interp_method_parameter_value sub_paramter_value =
      seq_parameter_values[seq_idx];

    switch (sub_parameter.type) {

      YAC_UNREACHABLE_DEFAULT_F(
        "ERROR(yac_yaml_create_seq_parameter_node): "
        "parameter \"%s\" has unsupported parameter type",
        sub_parameter.name);

      // parameter is a map
      case (MAP_PARAM): {

        sub_parameter_node =
          yac_yaml_create_non_default_map_parameter_node(
            document, sub_paramter_value.data.map_values, sub_parameter);
        break;
      }

      /*
       * currently, there is no use-case for this branch, therefore this would
       * be untested -> activate once it is being used

      // parameter is a scalar
      case (ENUM_PARAM):
      case (INT_PARAM):
      case (DBLE_PARAM):
      case (BOOL_PARAM):
      case (STR_PARAM):
      case (DEG_PARAM):
      {
        if (compare_parameter_values(
              &sub_parameter.default_value,
              &sub_paramter_value, sub_parameter))
          sub_parameter_node =
            yac_yaml_create_scalar_parameter_node(
              document, sub_parameter, sub_paramter_value);
        break;
      }
      */
    }

    // if the current parameter value is equal to the default values,
    // the parameter node is NULL, which is fine
    int appending_failed =
      fy_node_sequence_append(sequence_node, sub_parameter_node);
    YAC_ASSERT_F(
      !appending_failed,
      "ERROR(yac_yaml_create_seq_parameter_node): "
      "failed to append parameter \"%s\" to sequence", parameter.name);
  }

  return sequence_node;
}

static fy_node_t yac_yaml_create_non_default_map_parameter_node(
  fy_document_t document,
  interp_method_parameter_value const * map_parameter_values,
  struct interp_method_parameter parameter) {

  YAC_ASSERT_F(
    parameter.type == MAP_PARAM,
    "ERROR(yac_yaml_create_non_default_map_parameter_node): "
    "parameter \"%s\" is not a map", parameter.name);

  // create parameter node
  fy_node_t parameter_node = fy_node_create_mapping(document);
  YAC_ASSERT(
    parameter_node, "ERROR(yac_yaml_create_non_default_map_parameter_node): "
    "failed to create mapping node");

  int contains_non_default_values = 0;
  for (size_t sub_param_idx = 0;
       sub_param_idx < parameter.data.map_param.num_sub_params;
       ++sub_param_idx) {

    struct interp_method_parameter sub_parameter =
      parameter.data.map_param.sub_params[sub_param_idx];
    interp_method_parameter_value sub_parameter_value =
      map_parameter_values[sub_param_idx];

    switch(sub_parameter.type) {

      YAC_UNREACHABLE_DEFAULT_F(
        "ERROR(yac_yaml_create_non_default_map_parameter_node): "
        "parameter \"%s\" has unsupported parameter type",
        sub_parameter.name)

      // parameter is a map
      case (MAP_PARAM): {

        fy_node_t sub_parameter_node =
          yac_yaml_create_non_default_map_parameter_node(
            document, sub_parameter_value.data.map_values, sub_parameter);

        // if the map contained non-default parameter values
        if (sub_parameter_node) {
          yac_yaml_map_append(
            parameter_node, sub_parameter.name, sub_parameter_node);
          contains_non_default_values = 1;
        }
        break;
      }

      // parameter is a sequence
      case (SEQ_PARAM): {

        fy_node_t sub_parameter_node =
          yac_yaml_create_seq_parameter_node(
            document, sub_parameter_value.data.seq.values,
            sub_parameter_value.data.seq.count, sub_parameter);

        // if the sequence contained parameter values
        if (sub_parameter_node) {
          yac_yaml_map_append(
            parameter_node, sub_parameter.name, sub_parameter_node);
          contains_non_default_values = 1;
        }
        break;
      }

      // parameter is a scalar
      case (ENUM_PARAM):
      case (INT_PARAM):
      case (DBLE_PARAM):
      case (BOOL_PARAM):
      case (STR_PARAM):
      case (DEG_PARAM):
      {
        if (compare_parameter_values(
              &sub_parameter.default_value,
              &sub_parameter_value, sub_parameter)) {
          yac_yaml_map_append_parameter(
            document, parameter_node, sub_parameter, sub_parameter_value);
          contains_non_default_values = 1;
        }
        break;
      }
    }
  }

  // if there were no non-default parameters
  if (!contains_non_default_values) {
    fy_node_free(parameter_node);
    parameter_node = NULL;
  }
  return parameter_node;
}

static fy_node_t yac_yaml_create_interpolation_node(
  fy_document_t document,
  union yac_interp_stack_config_entry const * interp_stack_entry) {

  // get the interpolation data matching the interpolation stack entry
  enum yac_interpolation_list interp_type =
    yac_interp_stack_config_entry_get_type(interp_stack_entry);
  struct yac_interpolation_method const * interp_method = NULL;
  for (size_t i = 0;
       (i < NUM_INTERPOLATION_METHODS) && !interp_method; ++i)
    if (interpolation_methods[i].type == interp_type)
      interp_method = interpolation_methods + i;

  // get the interpolation parameter configuration and their default values
  struct interp_method_parameter parameter =
    interp_method->parameter;

  fy_node_t parameter_node;

  interp_method_parameter_value parameter_value =
    yaml_interp_method_parameter_get_default(parameter);

  // get parameter values set in the interpolation stack entry
  interp_method->get_interpolation(interp_stack_entry, &parameter_value);

  // create parameter node containing all non-default values
  parameter_node =
    yac_yaml_create_non_default_map_parameter_node(
      document, parameter_value.data.map_values, parameter);

  yaml_interp_method_parameter_free(parameter_value, parameter);

  fy_node_t interpolation_node;

  // if the interpolation contains non-default parameter values
  if (parameter_node) {

    // create interpolation node
    interpolation_node = fy_node_create_mapping(document);
    YAC_ASSERT(
      interpolation_node, "ERROR(yac_yaml_create_interpolation_node): "
      "failed to create mapping node");

    yac_yaml_map_append(
      interpolation_node, interp_method->name, parameter_node);
  } else {

    interpolation_node =
      yac_yaml_create_scalar(document, interp_method->name);
  }

  return interpolation_node;
}

static fy_node_t yac_yaml_create_interpolation_stack_node(
  fy_document_t document, struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx) {

  struct yac_interp_stack_config * interp_stack =
    yac_couple_config_get_interp_stack(
      couple_config, couple_idx, field_couple_idx);

  YAC_ASSERT(
    interp_stack, "ERROR(yac_yaml_create_interpolation_stack_node): "
    "invalid interpolation stack");

  // create interpolation stack node
  fy_node_t interp_stack_node = fy_node_create_sequence(document);
  YAC_ASSERT(
    interp_stack_node, "ERROR(yac_yaml_create_interpolation_stack_node): "
    "failed to create sequence node");

  size_t interp_stack_size =
    yac_interp_stack_config_get_size(interp_stack);
  YAC_ASSERT(
    interp_stack_size, "ERROR(yac_yaml_create_interpolation_stack_node): "
    "invalid interpolation stack size");

  for (size_t interp_stack_idx = 0; interp_stack_idx < interp_stack_size;
       ++interp_stack_idx) {
    int appending_failed =
      fy_node_sequence_append(
        interp_stack_node,
        yac_yaml_create_interpolation_node(
          document,
          yac_interp_stack_config_get_entry(
            interp_stack, interp_stack_idx)));
    YAC_ASSERT(
      !appending_failed,
      "ERROR(yac_yaml_create_interpolation_stack_node): "
      "failed to append interpolation node");
  }

  return interp_stack_node;
}

static void yac_yaml_append_couple_field_nodes(
  fy_node_t coupling_node, struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx) {

  fy_document_t document = fy_node_document(coupling_node);
  YAC_ASSERT(
    document, "ERROR(yac_yaml_append_couple_field_nodes): "
    "failed to get document from node");

  // create couple node
  fy_node_t field_couple_node = fy_node_create_mapping(document);
  YAC_ASSERT(
    coupling_node, "ERROR(yac_yaml_append_couple_field_nodes): "
    "failed to create mapping node");

  // get component names
  char const * src_component_name;
  char const * tgt_component_name;
  yac_couple_config_get_field_couple_component_names(
    couple_config, couple_idx, field_couple_idx,
    &src_component_name, &tgt_component_name);

  // add source component name
  yac_yaml_map_append_scalar(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, SOURCE_COMPONENT),
    src_component_name);

  // add target component name
  yac_yaml_map_append_scalar(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, TARGET_COMPONENT),
    tgt_component_name);

  // get grid names
  char const * src_grid_name;
  char const * tgt_grid_name;
  yac_couple_config_get_field_grid_names(
    couple_config, couple_idx, field_couple_idx,
    &src_grid_name, &tgt_grid_name);

  // add source grid name
  yac_yaml_map_append_scalar(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, SOURCE_GRID),
    src_grid_name);

  // add target grid name
  yac_yaml_map_append_scalar(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, TARGET_GRID),
    tgt_grid_name);

  // add field names
  yac_yaml_map_append(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, FIELD),
    yac_yaml_create_field_name_node(
      document, couple_config, couple_idx, field_couple_idx));

  // add coupling period
  yac_yaml_map_append_scalar(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, COUPLING_PERIOD),
    yac_couple_config_get_coupling_period(
      couple_config, couple_idx, field_couple_idx));

  // add time reduction
  yac_yaml_map_append_scalar(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, TIME_REDUCTION),
    yac_name_type_pair_get_name(
      time_operations, num_time_operations,
       yac_couple_config_get_coupling_period_operation(
         couple_config, couple_idx, field_couple_idx)));

  // add source lag
  yac_yaml_map_append_scalar_int(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, SOURCE_LAG),
    yac_couple_config_get_source_lag(
      couple_config, couple_idx, field_couple_idx));

  // add target lag
  yac_yaml_map_append_scalar_int(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, TARGET_LAG),
    yac_couple_config_get_target_lag(
      couple_config, couple_idx, field_couple_idx));

  // add weight file name
  if (yac_couple_config_enforce_write_weight_file(
        couple_config, couple_idx, field_couple_idx)) {
    yac_yaml_map_append_scalar(
      field_couple_node,
      yac_name_type_pair_get_name(
        yaml_couple_keys, num_yaml_couple_keys, WEIGHT_FILE_NAME),
      yac_couple_config_get_weight_file_name(
        couple_config, couple_idx, field_couple_idx));
    yac_yaml_map_append_scalar(
      field_couple_node,
      yac_name_type_pair_get_name(
        yaml_couple_keys, num_yaml_couple_keys, WEIGHT_FILE_ON_EXISTING),
      yac_name_type_pair_get_name(
        weight_file_on_existing_types, num_weight_file_on_existing_types,
        yac_couple_config_get_weight_file_on_existing(
          couple_config, couple_idx, field_couple_idx)));
  }

  // add mapping side
  yac_yaml_map_append_scalar(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, MAPPING_SIDE),
    yac_name_type_pair_get_name(
      mapping_sides, num_mapping_sides,
       yac_couple_config_mapping_on_source(
         couple_config, couple_idx, field_couple_idx)));

  // add scale factor
  yac_yaml_map_append_scalar_dble(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, SCALE_FACTOR),
    yac_couple_config_get_scale_factor(
      couple_config, couple_idx, field_couple_idx));

  // add scale summand
  yac_yaml_map_append_scalar_dble(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, SCALE_SUMMAND),
    yac_couple_config_get_scale_summand(
      couple_config, couple_idx, field_couple_idx));

  // add interpolation
  yac_yaml_map_append(
    field_couple_node,
    yac_name_type_pair_get_name(
      yaml_couple_keys, num_yaml_couple_keys, INTERPOLATION),
    yac_yaml_create_interpolation_stack_node(
      document, couple_config, couple_idx, field_couple_idx));

  // add source mask names
  char const * const * src_mask_names;
  size_t num_src_mask_names;
  yac_couple_config_get_src_mask_names(
    couple_config, couple_idx, field_couple_idx,
    &src_mask_names, &num_src_mask_names);
  if (num_src_mask_names == 1)
    yac_yaml_map_append_scalar(
      field_couple_node,
      yac_name_type_pair_get_name(
        yaml_couple_keys, num_yaml_couple_keys, SOURCE_MASK_NAME),
        src_mask_names[0]);
    else if (num_src_mask_names > 1)
      yac_yaml_map_append(
        field_couple_node,
        yac_name_type_pair_get_name(
          yaml_couple_keys, num_yaml_couple_keys, SOURCE_MASK_NAMES),
          yac_yaml_create_sequence_scalar(
            document, src_mask_names, num_src_mask_names));

  // add target mask name
  char const * tgt_mask_name =
    yac_couple_config_get_tgt_mask_name(
      couple_config, couple_idx, field_couple_idx);
  if (tgt_mask_name)
    yac_yaml_map_append_scalar(
      field_couple_node,
      yac_name_type_pair_get_name(
        yaml_couple_keys, num_yaml_couple_keys, TARGET_MASK_NAME),
      tgt_mask_name);

  // add yaxt exchanger name
  char const * yaxt_exchanger_name =
    yac_couple_config_get_yaxt_exchanger_name(
      couple_config, couple_idx, field_couple_idx);
  if (yaxt_exchanger_name)
    yac_yaml_map_append_scalar(
      field_couple_node,
      yac_name_type_pair_get_name(
        yaml_couple_keys, num_yaml_couple_keys, YAXT_EXCHANGER_NAME),
      yaxt_exchanger_name);

  // add use raw exchange
  yac_yaml_map_append_scalar(
    field_couple_node,
      yac_name_type_pair_get_name(
        yaml_couple_keys, num_yaml_couple_keys, USE_RAW_EXCHANGE),
    yac_name_type_pair_get_name(
      bool_names, num_bool_names,
       yac_couple_config_get_use_raw_exchange(
         couple_config, couple_idx, field_couple_idx)));

  int appending_failed =
    fy_node_sequence_append(coupling_node, field_couple_node);
  YAC_ASSERT(
    !appending_failed,
    "ERROR(yac_yaml_append_couple_field_nodes): "
    "failed to append field couple node");
}

static void yac_yaml_append_couple_nodes(
  fy_node_t coupling_node, struct yac_couple_config * couple_config,
  size_t couple_idx) {

  size_t num_couple_fields =
    yac_couple_config_get_num_couple_fields(couple_config, couple_idx);

  for (size_t field_couple_idx = 0;
       field_couple_idx < num_couple_fields; ++field_couple_idx)
    yac_yaml_append_couple_field_nodes(
      coupling_node, couple_config, couple_idx, field_couple_idx);
}

static fy_node_t yac_yaml_create_output_grid_node(
  fy_document_t document, char const * grid_name, char const * file_name) {

  fy_node_t output_grid_node = fy_node_create_mapping(document);
  YAC_ASSERT(
    output_grid_node, "ERROR(yac_yaml_create_output_grid_node): "
    "failed to create mapping node");

  yac_yaml_map_append_scalar(
    output_grid_node,
    yac_name_type_pair_get_name(
      yaml_debug_output_grid_keys, num_yaml_debug_output_grid_keys,
      OUTPUT_GRID_GRID_NAME),
    grid_name);
  yac_yaml_map_append_scalar(
    output_grid_node,
    yac_name_type_pair_get_name(
      yaml_debug_output_grid_keys, num_yaml_debug_output_grid_keys,
      OUTPUT_GRID_FILE_NAME),
    file_name);

  return output_grid_node;
}

static fy_node_t yac_yaml_create_output_grids_node(
  fy_document_t document, struct yac_couple_config * couple_config) {

  // count the number of output grids
  size_t num_output_grids = 0;
  size_t num_grids =
    yac_couple_config_get_num_grids(couple_config);
  for (size_t grid_idx = 0; grid_idx < num_grids; ++grid_idx)
    if (yac_couple_config_grid_get_output_filename(
          couple_config,
          yac_couple_config_get_grid_name(couple_config, grid_idx)) != NULL)
      ++num_output_grids;

  fy_node_t output_grids_node = NULL;

  if (num_output_grids > 0) {

    // create output grids node
    output_grids_node = fy_node_create_sequence(document);
    YAC_ASSERT(
      output_grids_node, "ERROR(yac_yaml_create_output_grids_node): "
      "failed to create sequence node");

    // for all output grids
    for (size_t grid_idx = 0; grid_idx < num_grids; ++grid_idx) {

      char const * grid_name =
        yac_couple_config_get_grid_name(couple_config, grid_idx);
      char const * file_name =
        yac_couple_config_grid_get_output_filename(couple_config, grid_name);

      if (file_name != NULL) {
        int appending_failed =
          fy_node_sequence_append(
            output_grids_node,
            yac_yaml_create_output_grid_node(document, grid_name, file_name));
        YAC_ASSERT(
          !appending_failed, "ERROR(yac_yaml_create_output_grids_node): "
          "failed to append output grid node");
      }
    }
  }

  return output_grids_node;
}

static fy_node_t yac_yaml_create_debug_node(
  fy_document_t document, struct yac_couple_config * couple_config) {

  // create debug node
  fy_node_t debug_node = fy_node_create_mapping(document);
  YAC_ASSERT(
    debug_node, "ERROR(yac_yaml_create_debug_node): "
    "failed to create mapping node");

  // add output grids node
  yac_yaml_map_append(
    debug_node,
    yac_name_type_pair_get_name(
      yaml_debug_keys, num_yaml_debug_keys, OUTPUT_GRIDS),
    yac_yaml_create_output_grids_node(document, couple_config));

  // add "missing_definition_is_fatal" node
  yac_yaml_map_append_scalar(
    debug_node,
    yac_name_type_pair_get_name(
      yaml_debug_keys, num_yaml_debug_keys, MISSING_DEF),
    yac_name_type_pair_get_name(
      bool_names, num_bool_names,
      yac_couple_config_get_missing_definition_is_fatal(couple_config)));

  return debug_node;
}

static fy_node_t yac_yaml_create_coupling_node(
  fy_document_t document, struct yac_couple_config * couple_config) {

  size_t num_couples =
    yac_couple_config_get_num_couples(couple_config);

  if (!num_couples) return NULL;

  // create coupling node
  fy_node_t coupling_node = fy_node_create_sequence(document);
  YAC_ASSERT(
    coupling_node, "ERROR(yac_yaml_create_coupling_node): "
    "failed to create sequence node");

  // for all couples
  for (size_t couple_idx = 0; couple_idx < num_couples; ++couple_idx)
    yac_yaml_append_couple_nodes(
      coupling_node, couple_config, couple_idx);

  return coupling_node;
}

static fy_node_t yac_yaml_create_field_node(
  fy_document_t document, struct yac_couple_config * couple_config,
  size_t component_idx, size_t field_idx) {

  fy_node_t field_node = fy_node_create_mapping(document);
  YAC_ASSERT(
    field_node, "ERROR(yac_yaml_create_field_node): "
    "failed to create mapping node");

  // get field parameters
  char const * component_name =
    yac_couple_config_get_component_name(
      couple_config, component_idx);
  char const * field_name =
    yac_couple_config_get_field_name(
      couple_config, component_idx, field_idx);
  char const * grid_name =
    yac_couple_config_get_field_grid_name(
      couple_config, component_idx, field_idx);
  char const * metadata =
    yac_couple_config_field_get_metadata(
      couple_config, component_name, grid_name, field_name);
  size_t collection_size =
    yac_couple_config_get_field_collection_size(
      couple_config, component_name, grid_name, field_name);
  char const * timestep =
    yac_couple_config_get_field_timestep(
      couple_config, component_name, grid_name, field_name);
  int role =
    yac_couple_config_get_field_role(
      couple_config, component_name, grid_name, field_name);
  double frac_mask_fallback_value =
    yac_couple_config_get_frac_mask_fallback_value(
      couple_config, component_name, grid_name, field_name);

  // add field name
  yac_yaml_map_append_scalar(field_node, "name", field_name);

  // add grid name
  yac_yaml_map_append_scalar(field_node, "grid_name", grid_name);

  // add metadata
  if (metadata)
    yac_yaml_map_append_scalar(field_node, "metadata", metadata);

  // add collection_size
  if (collection_size != SIZE_MAX)
    yac_yaml_map_append_scalar_int(
      field_node, "collection_size", (int)collection_size);

  // add timestep
  if (timestep)
    yac_yaml_map_append_scalar(field_node, "timestep", timestep);

  // add role
  yac_yaml_map_append_scalar(
    field_node, "role",
    yac_name_type_pair_get_name(
      role_types, num_role_types, (enum yac_field_exchange_type)role));

  // add fractional fallback value
  if (YAC_FRAC_MASK_VALUE_IS_VALID(frac_mask_fallback_value))
    yac_yaml_map_append_scalar_dble(
      field_node, "frac_mask_fallback_value", frac_mask_fallback_value);

  return field_node;
}

static fy_node_t yac_yaml_create_fields_node(
  fy_document_t document, struct yac_couple_config * couple_config,
  size_t component_idx) {

  size_t num_fields =
    yac_couple_config_get_num_fields(couple_config, component_idx);

  if (num_fields == 0) return (fy_node_t)NULL;

  fy_node_t fields_node = fy_node_create_sequence(document);
  YAC_ASSERT(
    fields_node, "ERROR(yac_yaml_create_fields_node): "
    "failed to create sequence node");

  for (size_t field_idx = 0; field_idx < num_fields; ++field_idx) {

    int appending_failed =
      fy_node_sequence_append(
        fields_node,
        yac_yaml_create_field_node(
          document, couple_config, component_idx, field_idx));
    YAC_ASSERT(
      !appending_failed, "ERROR(yac_yaml_create_fields_node): "
      "failed to append field node");
  }

  return fields_node;
}

static fy_node_t yac_yaml_create_component_node(
  fy_document_t document, struct yac_couple_config * couple_config,
  size_t component_idx) {

  fy_node_t component_node = fy_node_create_mapping(document);
  YAC_ASSERT(
    component_node, "ERROR(yac_yaml_create_component_node): "
    "failed to create mapping node");

  // get component name, component metadata, and fields
  char const * component_name =
    yac_couple_config_get_component_name(couple_config, component_idx);
  char const * metadata =
    yac_couple_config_component_get_metadata(couple_config, component_name);
  fy_node_t fields_node =
    yac_yaml_create_fields_node(
      document, couple_config, component_idx);

  // add component name
  yac_yaml_map_append_scalar(component_node, "name", component_name);

  // add metadata
  if (metadata)
    yac_yaml_map_append_scalar(component_node, "metadata", metadata);

  // add fields
  yac_yaml_map_append(component_node, "fields", fields_node);

  return component_node;
}

static fy_node_t yac_yaml_create_components_node(
  fy_document_t document, struct yac_couple_config * couple_config) {

  // get number of components in coupling configuration
  size_t num_components = yac_couple_config_get_num_components(couple_config);

  if (num_components == 0) return (fy_node_t)NULL;

  // create sequence node
  fy_node_t components_node = fy_node_create_sequence(document);
  YAC_ASSERT(
    components_node, "ERROR(yac_yaml_create_components_node): "
    "failed to create sequence node");

  for (size_t component_idx = 0; component_idx < num_components;
      ++component_idx) {

    int appending_failed =
      fy_node_sequence_append(
        components_node,
        yac_yaml_create_component_node(document, couple_config, component_idx));
    YAC_ASSERT(
      !appending_failed, "ERROR(yac_yaml_create_components_node): "
      "failed to append component node");
  }

  return components_node;
}

static fy_node_t yac_yaml_create_grid_node(
  fy_document_t document, struct yac_couple_config * couple_config,
  size_t grid_idx) {

  fy_node_t grid_node = fy_node_create_mapping(document);
  YAC_ASSERT(
    grid_node, "ERROR(yac_yaml_create_grid_node): "
    "failed to create mapping node");

  // get grid name and metadata
  char const * grid_name =
    yac_couple_config_get_grid_name(couple_config, grid_idx);
  char const * metadata =
    yac_couple_config_grid_get_metadata(couple_config, grid_name);

  // add grid name
  yac_yaml_map_append_scalar(grid_node, "name", grid_name);

  // add metadata
  if (metadata)
    yac_yaml_map_append_scalar(grid_node, "metadata", metadata);

  return grid_node;
}

static fy_node_t yac_yaml_create_grids_node(
  fy_document_t document, struct yac_couple_config * couple_config) {

  // get number of grids in coupling configuration
  size_t num_grids = yac_couple_config_get_num_grids(couple_config);

  if (num_grids == 0) return (fy_node_t)NULL;

  // create sequence node
  fy_node_t grids_node = fy_node_create_sequence(document);
  YAC_ASSERT(
    grids_node, "ERROR(yac_yaml_create_grids_node): "
    "failed to create sequence node");

  for (size_t grids_idx = 0; grids_idx < num_grids; ++grids_idx) {

    int appending_failed =
      fy_node_sequence_append(
        grids_node,
        yac_yaml_create_grid_node(document, couple_config, grids_idx));
    YAC_ASSERT(
      !appending_failed, "ERROR(yac_yaml_create_grids_node): "
      "failed to append grid node");
  }

  return grids_node;
}

static fy_node_t yac_yaml_create_definitions_node(
  fy_document_t document, struct yac_couple_config * couple_config) {

  // create definition
  fy_node_t definition_node = fy_node_create_mapping(document);
  YAC_ASSERT(
    definition_node,
    "ERROR(yac_yaml_create_definitions_node): "
    "failed to create mapping node");

  // add components
  yac_yaml_map_append(
    definition_node,
    "components", yac_yaml_create_components_node(document, couple_config));

  // add grids
  yac_yaml_map_append(
    definition_node, "grids", yac_yaml_create_grids_node(document, couple_config));

  return definition_node;
}

static fy_node_t yac_yaml_create_couple_config_nodes(
  fy_document_t document, struct yac_couple_config * couple_config,
  int include_definitions) {

  // create root node
  fy_node_t root_node = fy_node_create_mapping(document);
  YAC_ASSERT(
    root_node,
    "ERROR(yac_yaml_create_couple_config_nodes): "
    "failed to create mapping node");

  // add debug
  yac_yaml_map_append(
    root_node,
    yac_name_type_pair_get_name(
      yaml_base_keys, num_yaml_base_keys, DEBUG),
    yac_yaml_create_debug_node(document, couple_config));

  // add user definitions (components, grids, and fields)
  if (include_definitions)
    yac_yaml_map_append(
      root_node, "definitions",
      yac_yaml_create_definitions_node(document, couple_config));

  // add start datetime
  char * start_datetime = yac_couple_config_get_start_datetime(couple_config);
  yac_yaml_map_append_scalar(
    root_node,
    yac_name_type_pair_get_name(
      yaml_base_keys, num_yaml_base_keys, START_DATE), start_datetime);
  free(start_datetime);

  // add end datetime
  char * end_datetime = yac_couple_config_get_end_datetime(couple_config);
  yac_yaml_map_append_scalar(
    root_node,
    yac_name_type_pair_get_name(
      yaml_base_keys, num_yaml_base_keys, END_DATE), end_datetime);
  free(end_datetime);

  // add calendar
  yac_yaml_map_append_scalar(
    root_node,
    yac_name_type_pair_get_name(
      yaml_base_keys, num_yaml_base_keys, CALENDAR),
    yac_name_type_pair_get_name(
      calendar_types, num_calendar_types, getCalendarType()));

  // add timestep unit
  yac_yaml_map_append_scalar(
    root_node,
    yac_name_type_pair_get_name(
      yaml_base_keys, num_yaml_base_keys, TIMESTEP_UNIT),
    yac_name_type_pair_get_name(
      timestep_units, num_timestep_units, C_ISO_FORMAT));

  // add couplings
  yac_yaml_map_append(
    root_node,
    yac_name_type_pair_get_name(
      yaml_base_keys, num_yaml_base_keys, COUPLING),
    yac_yaml_create_coupling_node(document, couple_config));

  return root_node;
}

char * yac_yaml_emit_coupling(
  struct yac_couple_config * couple_config, int emit_flags,
  int include_definitions) {

  // create an empty document
  fy_document_t document = fy_document_create(NULL);
  YAC_ASSERT(
    document, "ERROR(yac_yaml_emit): failed to create document");

  // create nodes from coupling configuration
  fy_node_t root_node =
    yac_yaml_create_couple_config_nodes(
      document, couple_config, include_definitions);

  // set root node of the document
  int setting_root_failed =
    fy_document_set_root(document, root_node);
  YAC_ASSERT(
    !setting_root_failed,
    "ERROR(yac_yaml_emit): failed to add root node to document");

  // emit document to string
  char * str_document =
    fy_emit_document_to_string(
      document, (enum fy_emitter_cfg_flags)emit_flags);

  YAC_ASSERT(
    str_document, "ERROR(yac_yaml_emit): failed to emit document to string");

  // destroy document
  fy_document_destroy(document);

  return str_document;
}
