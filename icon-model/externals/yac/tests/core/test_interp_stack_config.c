// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "tests.h"
#include "dist_grid_utils.h"
#include "interp_stack_config.h"
#include "interp_method_spmap.h"
#include "geometry.h"

/** \file test_interp_stack_config.c
 *  \test
 * Tests for interpolation stack interface routines.
 */

#define ARGS(...) __VA_ARGS__
#define _GET_NTH_ARG(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, N, ...) N
#define EXPAND(x) x
#define FOREACH(name, ...) \
  { \
    enum {NUM_ ## name = sizeof( name ) / sizeof( name [0])}; \
    int name ## _idx[2]; \
    for (name ## _idx[0] = 0; name ## _idx[0] < NUM_ ## name; \
         ++ name ## _idx[0]) { \
      for (name ## _idx[1] = 0; name ## _idx[1] < NUM_ ## name; \
           ++ name ## _idx[1]) { \
        configs_differ += (name ## _idx[0]) != (name ## _idx[1]); \
        {__VA_ARGS__} \
        configs_differ -= (name ## _idx[0]) != (name ## _idx[1]); \
      } \
    } \
  }
#define FOREACH_ENUM(name, values, ...) \
  { \
    enum yac_ ## name name [] = {values}; \
    FOREACH(name, __VA_ARGS__) \
  }
#define FOREACH_TYPE(name, type, values, ...) \
  { \
    type name[] = {values}; \
    FOREACH(name, __VA_ARGS__) \
  }
#define FOREACH_INT(name, values, ...) \
  FOREACH_TYPE(name, int, ARGS(values), __VA_ARGS__)
#define FOREACH_DBLE(name, values, ...) \
  FOREACH_TYPE(name, double, ARGS(values), __VA_ARGS__)
#define FOREACH_BOOL(name, ...) FOREACH_INT(name, ARGS(0, 1), __VA_ARGS__)
#define FOREACH_STRING(name, values, ...) \
  FOREACH_TYPE(name, ARGS(char const *), ARGS(values), __VA_ARGS__)
#define FOREACH_STRUCT(name, struct_name, values, ...) \
  FOREACH_TYPE(name, ARGS(struct struct_name), ARGS(values), __VA_ARGS__)
#define _CHECK_STACKS(interp_name, config) \
  { \
    int config_idx; \
    struct yac_interp_stack_config * a = yac_interp_stack_config_new(); \
    struct yac_interp_stack_config * b = yac_interp_stack_config_new(); \
    config_idx = 0, yac_interp_stack_config_add_ ## interp_name ( a, config ); \
    config_idx = 1, yac_interp_stack_config_add_ ## interp_name ( b, config ); \
    utest_check_compare_stacks(a, b, configs_differ); \
  }
#define _CONFIG_ARGS1(arg_name) arg_name[arg_name ## _idx[config_idx]]
#define _CONFIG_ARGS2(arg_name, ...) \
  _CONFIG_ARGS1(arg_name), _CONFIG_ARGS1(__VA_ARGS__)
#define _CONFIG_ARGS3(arg_name, ...) \
  _CONFIG_ARGS1(arg_name), _CONFIG_ARGS2(__VA_ARGS__)
#define _CONFIG_ARGS4(arg_name, ...) \
  _CONFIG_ARGS1(arg_name), _CONFIG_ARGS3(__VA_ARGS__)
#define _CONFIG_ARGS5(arg_name, ...) \
  _CONFIG_ARGS1(arg_name), _CONFIG_ARGS4(__VA_ARGS__)
#define _CONFIG_ARGS6(arg_name, ...) \
  _CONFIG_ARGS1(arg_name), _CONFIG_ARGS5(__VA_ARGS__)
#define _CONFIG_ARGS7(arg_name, ...) \
  _CONFIG_ARGS1(arg_name), _CONFIG_ARGS6(__VA_ARGS__)
#define _CONFIG_ARGS8(arg_name, ...) \
  _CONFIG_ARGS1(arg_name), _CONFIG_ARGS7(__VA_ARGS__)
#define _CONFIG_ARGS9(arg_name, ...) \
  _CONFIG_ARGS1(arg_name), _CONFIG_ARGS8(__VA_ARGS__)
#define _CONFIG_ARGS10(arg_name, ...) \
  _CONFIG_ARGS1(arg_name), _CONFIG_ARGS9(__VA_ARGS__)
#define _CONFIG_ARGS11(arg_name, ...) \
  _CONFIG_ARGS1(arg_name), _CONFIG_ARGS10(__VA_ARGS__)
#define _CONFIG_ARGS12(arg_name, ...) \
  _CONFIG_ARGS1(arg_name), _CONFIG_ARGS11(__VA_ARGS__)
#define CHECK_STACKS(interp_name, ... ) \
  _CHECK_STACKS(interp_name, \
    EXPAND(_GET_NTH_ARG(__VA_ARGS__, _CONFIG_ARGS12, \
                                     _CONFIG_ARGS11, \
                                     _CONFIG_ARGS10, \
                                     _CONFIG_ARGS9, \
                                     _CONFIG_ARGS8, \
                                     _CONFIG_ARGS7, \
                                     _CONFIG_ARGS6, \
                                     _CONFIG_ARGS5, \
                                     _CONFIG_ARGS4, \
                                     _CONFIG_ARGS3, \
                                     _CONFIG_ARGS2, \
                                     _CONFIG_ARGS1)(__VA_ARGS__)))

static void utest_check_compare_stacks(
  struct yac_interp_stack_config * a, struct yac_interp_stack_config * b,
  int configs_differ);

static void yac_interp_stack_config_add_spmap_(
  struct yac_interp_stack_config * interp_stack_config,
  double spread_distance, double max_search_distance,
  enum yac_interp_spmap_weight_type weight_type,
  enum yac_interp_spmap_scale_type scale_type,
  struct yac_spmap_cell_area_config * src_cell_area_config,
  struct yac_spmap_cell_area_config * tgt_cell_area_config);

static void yac_interp_stack_config_add_spmap_ext_(
  struct yac_interp_stack_config * interp_stack_config,
  struct yac_interp_spmap_config * default_config,
  struct yac_spmap_overwrite_config ** overwrite_configs);

int main (void) {

  MPI_Init(NULL, NULL);
  xt_initialize(MPI_COMM_WORLD);

  int configs_differ = 0;

  { // stack with different sizes
    struct yac_interp_stack_config * a = yac_interp_stack_config_new();
    struct yac_interp_stack_config * b = yac_interp_stack_config_new();
    yac_interp_stack_config_add_average(a, YAC_INTERP_AVG_ARITHMETIC, 1);
    yac_interp_stack_config_add_fixed(a, -1.0);
    yac_interp_stack_config_add_average(a, YAC_INTERP_AVG_ARITHMETIC, 1);
    utest_check_compare_stacks(a, b, 1);
  }

  { // compare empty config
    struct yac_interp_stack_config * a = yac_interp_stack_config_new();
    struct yac_interp_stack_config * b = yac_interp_stack_config_new();
    utest_check_compare_stacks(a, b, 0);
  }

  // compare average config
  FOREACH_ENUM(
    interp_avg_weight_type,
    ARGS(YAC_INTERP_AVG_ARITHMETIC, YAC_INTERP_AVG_DIST, YAC_INTERP_AVG_BARY),
    FOREACH_BOOL(
      partial_coverage,
      CHECK_STACKS(average, interp_avg_weight_type, partial_coverage)))

  // compare ncc config
  FOREACH_ENUM(
    interp_ncc_weight_type,
    ARGS(YAC_INTERP_NCC_AVG, YAC_INTERP_NCC_DIST),
    FOREACH_BOOL(
      partial_coverage,
      CHECK_STACKS(ncc, interp_ncc_weight_type, partial_coverage)))

  // compare nnn config
  //  for YAC_INTERP_NNN_AVG, YAC_INTERP_NNN_DIST, and YAC_INTERP_NNN_ZERO
  //  the scale parameter is being ignored
  FOREACH_ENUM(
    interp_nnn_weight_type,
    ARGS(YAC_INTERP_NNN_AVG, YAC_INTERP_NNN_DIST, YAC_INTERP_NNN_ZERO),
    FOREACH_INT(
      counts, ARGS(1,3,9),
      FOREACH_DBLE(
        max_search_distance, ARGS(0.0, M_PI_2),
        FOREACH_DBLE(
          scales, -1.0,
          CHECK_STACKS(
            nnn, interp_nnn_weight_type, counts,
            max_search_distance, scales)))))

  // compare nnn config
  //   for YAC_INTERP_NNN_GAUSS and YAC_INTERP_NNN_RBF the scale
  //   parameter is being interpreted
  FOREACH_ENUM(
    interp_nnn_weight_type,
    ARGS(YAC_INTERP_NNN_GAUSS, YAC_INTERP_NNN_RBF),
    FOREACH_INT(
      counts, ARGS(1,3,9),
      FOREACH_DBLE(
        max_search_distance, ARGS(0.0, M_PI_2),
        FOREACH_DBLE(
          scales, ARGS(0.5, 1.0),
          CHECK_STACKS(
            nnn, interp_nnn_weight_type, counts,
            max_search_distance, scales)))))

  // compare conservative config
  FOREACH_INT(
    order, ARGS(1,2),
    FOREACH_BOOL(
      enforced_conserv,
      FOREACH_BOOL(
        partial_coverage,
        FOREACH_ENUM(
          interp_method_conserv_normalisation,
          ARGS(YAC_INTERP_CONSERV_DESTAREA, YAC_INTERP_CONSERV_FRACAREA),
          CHECK_STACKS(conservative,
            order, enforced_conserv, partial_coverage,
            interp_method_conserv_normalisation)))))

  // compare source point mapping
  {
    struct yac_spmap_cell_area_config * cell_area_configs[] = {
        yac_spmap_cell_area_config_yac_new(1.0),
        yac_spmap_cell_area_config_yac_new(2.0),
        yac_spmap_cell_area_config_file_new("area.nc", "cell_area", 0),
        yac_spmap_cell_area_config_file_new("area.nc_", "cell_area", 0),
        yac_spmap_cell_area_config_file_new("area.nc", "cell_area_", 0),
        yac_spmap_cell_area_config_file_new("area.nc", "cell_area", 1),
        yac_spmap_cell_area_config_file_new("area.nc_", "cell_area_", 1),
      };
    enum {
      NUM_CELL_AREA_CONFIGS =
        sizeof(cell_area_configs) / sizeof(cell_area_configs[0])};

    FOREACH_TYPE(
      src_cell_area_config, ARGS(struct yac_spmap_cell_area_config *),
      ARGS(
        cell_area_configs[0], cell_area_configs[1], cell_area_configs[2],
        cell_area_configs[3], cell_area_configs[4], cell_area_configs[5],
        cell_area_configs[6]),
      FOREACH_TYPE(
        tgt_cell_area_config, ARGS(struct yac_spmap_cell_area_config *),
        ARGS(
          cell_area_configs[0], cell_area_configs[1], cell_area_configs[2],
          cell_area_configs[3], cell_area_configs[4], cell_area_configs[5],
          cell_area_configs[6]),
        FOREACH_DBLE(
          spread_distance, ARGS(0.0, 0.1),
          FOREACH_DBLE(
            max_search_distance, ARGS(0.0, 0.4),
            FOREACH_ENUM(
              interp_spmap_weight_type,
              ARGS(YAC_INTERP_SPMAP_AVG, YAC_INTERP_SPMAP_DIST),
              FOREACH_ENUM(
                interp_spmap_scale_type,
                ARGS(
                  YAC_INTERP_SPMAP_NONE, YAC_INTERP_SPMAP_SRCAREA,
                  YAC_INTERP_SPMAP_INVTGTAREA, YAC_INTERP_SPMAP_FRACAREA),
                CHECK_STACKS(spmap_,
                  spread_distance, max_search_distance,
                  interp_spmap_weight_type, interp_spmap_scale_type,
                  src_cell_area_config, tgt_cell_area_config)))))))

    for (size_t i = 0; i < NUM_CELL_AREA_CONFIGS; ++i)
      yac_spmap_cell_area_config_delete(cell_area_configs[i]);
  }

  // compare source point mapping extended
  {
    struct yac_spmap_scale_config * scale_config_custom =
      yac_spmap_scale_config_new(YAC_INTERP_SPMAP_FRACAREA, NULL, NULL);
    struct yac_interp_spmap_config * default_configs[] = {
        YAC_INTERP_SPMAP_DEFAULT_CONFIG,
        yac_interp_spmap_config_new(
          1.0,
          YAC_INTERP_SPMAP_MAX_SEARCH_DISTANCE_DEFAULT,
          YAC_INTERP_SPMAP_WEIGHTED_DEFAULT,
          YAC_INTERP_SPMAP_SCALE_CONFIG_DEFAULT),
        yac_interp_spmap_config_new(
          YAC_INTERP_SPMAP_SPREAD_DISTANCE_DEFAULT,
          1.0,
          YAC_INTERP_SPMAP_WEIGHTED_DEFAULT,
          YAC_INTERP_SPMAP_SCALE_CONFIG_DEFAULT),
        yac_interp_spmap_config_new(
          YAC_INTERP_SPMAP_SPREAD_DISTANCE_DEFAULT,
          YAC_INTERP_SPMAP_MAX_SEARCH_DISTANCE_DEFAULT,
          YAC_INTERP_SPMAP_DIST,
          YAC_INTERP_SPMAP_SCALE_CONFIG_DEFAULT),
        yac_interp_spmap_config_new(
          YAC_INTERP_SPMAP_SPREAD_DISTANCE_DEFAULT,
          YAC_INTERP_SPMAP_MAX_SEARCH_DISTANCE_DEFAULT,
          YAC_INTERP_SPMAP_WEIGHTED_DEFAULT,
          scale_config_custom),
      };
    enum {
      NUM_DEFAULT_CONFIGS =
        sizeof(default_configs)/sizeof(default_configs[0])};

    struct yac_point_selection * bnd_point_selection_a =
      yac_point_selection_bnd_circle_new(0.1, 0.1, 1.0);
    struct yac_point_selection * bnd_point_selection_b =
      yac_point_selection_bnd_circle_new(0.1, 0.1, 2.0);
    struct yac_interp_spmap_config * spmap_config_custom =
      yac_interp_spmap_config_new(
          1.0,
          YAC_INTERP_SPMAP_MAX_SEARCH_DISTANCE_DEFAULT,
          YAC_INTERP_SPMAP_DIST,
          YAC_INTERP_SPMAP_SCALE_CONFIG_DEFAULT);

    struct yac_spmap_overwrite_config *** overwrite_configs =
      (struct yac_spmap_overwrite_config **[])
        {NULL,
        (struct yac_spmap_overwrite_config *[]) {
          yac_spmap_overwrite_config_new(bnd_point_selection_a, NULL),
          NULL},
        (struct yac_spmap_overwrite_config *[]) {
          yac_spmap_overwrite_config_new(bnd_point_selection_b, NULL),
          NULL},
        (struct yac_spmap_overwrite_config *[]) {
          yac_spmap_overwrite_config_new(
            bnd_point_selection_b, spmap_config_custom),
          NULL},
        (struct yac_spmap_overwrite_config *[]) {
          yac_spmap_overwrite_config_new(
            bnd_point_selection_a, NULL),
          yac_spmap_overwrite_config_new(
            bnd_point_selection_b, NULL),
          NULL},
        (struct yac_spmap_overwrite_config *[]) {
          yac_spmap_overwrite_config_new(bnd_point_selection_a, NULL),
          yac_spmap_overwrite_config_new(
            bnd_point_selection_b, spmap_config_custom),
          NULL}};
  enum {NUM_OVERWRITE_CONFIGS = 6};

  FOREACH_TYPE(
    default_config, ARGS(struct yac_interp_spmap_config *),
    ARGS(
      default_configs[0], default_configs[1],
      default_configs[2], default_configs[3]),
    FOREACH_TYPE(
      overwrite_config, ARGS(struct yac_spmap_overwrite_config **),
      ARGS(
        overwrite_configs[0], overwrite_configs[1], overwrite_configs[2],
        overwrite_configs[3], overwrite_configs[4], overwrite_configs[5]),
      CHECK_STACKS(spmap_ext_, default_config, overwrite_config)))

    for (size_t i = 0; i < NUM_OVERWRITE_CONFIGS; ++i)
      for (size_t j = 0;
           (overwrite_configs[i] != NULL) && (overwrite_configs[i][j] != NULL);
           ++j)
        yac_spmap_overwrite_config_delete(overwrite_configs[i][j]);
    yac_interp_spmap_config_delete(spmap_config_custom);
    yac_point_selection_delete(bnd_point_selection_a);
    yac_point_selection_delete(bnd_point_selection_b);
    for (size_t i = 0; i < NUM_DEFAULT_CONFIGS; ++i)
      yac_interp_spmap_config_delete(default_configs[i]);
    yac_spmap_scale_config_delete(scale_config_custom);
  }

  // compare user file
  FOREACH_STRING(
    filename, ARGS(
      "test_interp_stack_config_file_a.nc",
      "test_interp_stack_config_file_b.nc"),
    FOREACH_ENUM(
      interp_file_on_missing_file,
      ARGS(YAC_INTERP_FILE_MISSING_ERROR, YAC_INTERP_FILE_MISSING_CONT),
      FOREACH_ENUM(
        interp_file_on_success,
        ARGS(YAC_INTERP_FILE_SUCCESS_STOP, YAC_INTERP_FILE_SUCCESS_CONT),
        CHECK_STACKS(
          user_file, filename, interp_file_on_missing_file,
          interp_file_on_success))))

  // compare fixed
  FOREACH_DBLE(
    fixed_value, ARGS(-1.0, 0.0, 1.0),
    CHECK_STACKS(fixed, fixed_value))

  // compare check
  FOREACH_STRING(
    constructor_key, ARGS(NULL, "constructor_a", "constructor_b"),
    FOREACH_STRING(
      do_search_key, ARGS(NULL, "do_search_key_a", "do_search_key_b"),
      CHECK_STACKS(check, constructor_key, do_search_key)))

  // compare creep
  FOREACH_INT(
    creep_distance, ARGS(-1, 0, 1),
    CHECK_STACKS(creep, creep_distance))

  // compare user callback
  FOREACH_STRING(
    compute_weights_key, ARGS("compute_weights_a", "compute_weights_b"),
    CHECK_STACKS(user_callback, compute_weights_key))

  { // testing spmap interpolation generated from an interp_stack

    // trivial 2x2 grid
    double * coords = (double[]){0.0,0.1,0.2};
    struct yac_basic_grid_data grid_data[2] =
      {yac_generate_basic_grid_data_reg2d(
         coords, coords, (size_t[]){2,2},
         (size_t[]){0,0}, (size_t[]){2,2}, 1),
       yac_generate_basic_grid_data_reg2d(
         coords, coords, (size_t[]){2,2},
         (size_t[]){0,0}, (size_t[]){2,2}, 1)};
    yac_coordinate_pointer cell_center_coords =
      malloc(4 * sizeof(*cell_center_coords));
    for (int i = 0, k = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j, ++k)
        LLtoXYZ(
          (coords[j] + coords[j+1])*0.5, (coords[i] + coords[i+1])*0.5,
          cell_center_coords[k]);

    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new("src_grid", grid_data[0]),
       yac_basic_grid_new("tgt_grid", grid_data[1])};

    size_t src_cell_coord_idx =
      yac_basic_grid_add_coordinates(
        grids[0], YAC_LOC_CELL, cell_center_coords, 4);
    size_t tgt_cell_coord_idx =
      yac_basic_grid_add_coordinates(
        grids[1], YAC_LOC_CELL, cell_center_coords, 4);
    free(cell_center_coords);

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CELL, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CELL, .masks_idx = SIZE_MAX};
    src_fields[0].coordinates_idx = src_cell_coord_idx;
    tgt_field.coordinates_idx = tgt_cell_coord_idx;

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, "src_grid", "tgt_grid",
                          num_src_fields, src_fields, tgt_field);

    enum {OVERWRITE_CONFIG_COUNT = 1};
    struct yac_spmap_overwrite_config *
      overwrite_configs[OVERWRITE_CONFIG_COUNT+1];
    overwrite_configs[OVERWRITE_CONFIG_COUNT] = NULL;

    {
      struct yac_point_selection * src_point_selection =
        yac_point_selection_bnd_circle_new(0.05, 0.05, 0.01);
      struct yac_interp_spmap_config * spmap_config =
        yac_interp_spmap_config_new(
          0.11,
          YAC_INTERP_SPMAP_MAX_SEARCH_DISTANCE_DEFAULT,
          YAC_INTERP_SPMAP_WEIGHTED_DEFAULT,
          YAC_INTERP_SPMAP_SCALE_CONFIG_DEFAULT);
      overwrite_configs[0] =
        yac_spmap_overwrite_config_new(src_point_selection, spmap_config);
      yac_interp_spmap_config_delete(spmap_config);
      yac_point_selection_delete(src_point_selection);
    }

    struct yac_interp_stack_config * interp_stack_config =
      yac_interp_stack_config_new();

    yac_interp_stack_config_add_spmap_ext(
      interp_stack_config, YAC_INTERP_SPMAP_DEFAULT_CONFIG, overwrite_configs);

    struct interp_method ** method_stack =
      yac_interp_stack_config_generate(interp_stack_config);

    yac_interp_stack_config_delete(interp_stack_config);

    struct yac_interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    struct yac_interpolation * interpolation =
      yac_interp_weights_get_interpolation(
        weights, YAC_MAPPING_ON_SRC, 1,
        YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0, NULL, 1, 1);

    {
      double * src_field = (double[]){1.0,2.0,3.0,4.0};
      double ** src_fields = &src_field;
      double * tgt_field = (double[]){0.0,0.0,0.0,0.0};
      double const * ref_tgt_field =
        (double[]){1.0/3.0, 2.0+1.0/3.0, 3.0+1.0/3.0, 4.0};

      yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

      for (int i = 0; i < 4; ++i)
        if (fabs(tgt_field[i] - ref_tgt_field[i]) > 1e-6)
          PUT_ERR("ERROR in yac_interp_stack_config_add_spmap_ext");
    }

    yac_interpolation_delete(interpolation);

    yac_interp_weights_delete(weights);
    yac_interp_method_delete(method_stack);
    for (size_t i = 0; i < OVERWRITE_CONFIG_COUNT; ++i)
      yac_spmap_overwrite_config_delete(overwrite_configs[i]);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(grids[1]);
    yac_basic_grid_delete(grids[0]);
  }

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void utest_check_compare_stacks_(
  struct yac_interp_stack_config * a, struct yac_interp_stack_config * b,
  int configs_differ) {

  configs_differ = configs_differ != 0;

  if (yac_interp_stack_config_compare(a, a))
    PUT_ERR("error in yac_interp_stack_config_compare (a != a)")
  if (yac_interp_stack_config_compare(b, b))
    PUT_ERR("error in yac_interp_stack_config_compare (b != b)")

  int cmp_a = yac_interp_stack_config_compare(a, b);
  int cmp_b = yac_interp_stack_config_compare(b, a);

  if ((cmp_a != cmp_b) ^ configs_differ)
    PUT_ERR("error in yac_interp_stack_config_compare ((a > b) == (a < b))")
  if ((cmp_a != 0) ^ configs_differ)
    PUT_ERR("error in yac_interp_stack_config_compare ((a > b) == 0)")
  if ((cmp_b != 0) ^ configs_differ)
    PUT_ERR("error in yac_interp_stack_config_compare ((a > b) == 0)")

  yac_interp_stack_config_delete(b);
  yac_interp_stack_config_delete(a);
}

static void utest_check_compare_stacks(
  struct yac_interp_stack_config * a, struct yac_interp_stack_config * b,
  int configs_differ) {

  utest_check_compare_stacks_(
    yac_interp_stack_config_copy(a),
    yac_interp_stack_config_copy(b), configs_differ);
  utest_check_compare_stacks_(a, b, configs_differ);
}

static void yac_interp_stack_config_add_spmap_(
  struct yac_interp_stack_config * interp_stack_config,
  double spread_distance, double max_search_distance,
  enum yac_interp_spmap_weight_type weight_type,
  enum yac_interp_spmap_scale_type scale_type,
  struct yac_spmap_cell_area_config * src_cell_area_config,
  struct yac_spmap_cell_area_config * tgt_cell_area_config) {

  enum yac_interp_spmap_cell_area_provider src_cell_area_config_type =
    yac_spmap_cell_area_config_get_type(src_cell_area_config);
  double src_sphere_radius =
    (src_cell_area_config_type == YAC_INTERP_SPMAP_CELL_AREA_YAC)?
      yac_spmap_cell_area_config_get_sphere_radius(src_cell_area_config):0.0;
  char const * src_filename =
    (src_cell_area_config_type == YAC_INTERP_SPMAP_CELL_AREA_FILE)?
      yac_spmap_cell_area_config_get_filename(src_cell_area_config):NULL;
  char const * src_varname =
    (src_cell_area_config_type == YAC_INTERP_SPMAP_CELL_AREA_FILE)?
      yac_spmap_cell_area_config_get_varname(src_cell_area_config):NULL;
  yac_int src_min_global_id =
    (src_cell_area_config_type == YAC_INTERP_SPMAP_CELL_AREA_FILE)?
      yac_spmap_cell_area_config_get_min_global_id(src_cell_area_config):0;

  enum yac_interp_spmap_cell_area_provider tgt_cell_area_config_type =
    yac_spmap_cell_area_config_get_type(tgt_cell_area_config);
  double tgt_sphere_radius =
    (tgt_cell_area_config_type == YAC_INTERP_SPMAP_CELL_AREA_YAC)?
      yac_spmap_cell_area_config_get_sphere_radius(tgt_cell_area_config):0.0;
  char const * tgt_filename =
    (tgt_cell_area_config_type == YAC_INTERP_SPMAP_CELL_AREA_FILE)?
      yac_spmap_cell_area_config_get_filename(tgt_cell_area_config):NULL;
  char const * tgt_varname =
    (tgt_cell_area_config_type == YAC_INTERP_SPMAP_CELL_AREA_FILE)?
      yac_spmap_cell_area_config_get_varname(tgt_cell_area_config):NULL;
  yac_int tgt_min_global_id =
    (tgt_cell_area_config_type == YAC_INTERP_SPMAP_CELL_AREA_FILE)?
      yac_spmap_cell_area_config_get_min_global_id(tgt_cell_area_config):0;

  yac_interp_stack_config_add_spmap(
    interp_stack_config,
    spread_distance, max_search_distance, weight_type, scale_type,
    src_sphere_radius,
    src_filename,
    src_varname,
    src_min_global_id,
    tgt_sphere_radius,
    tgt_filename,
    tgt_varname,
    tgt_min_global_id);

  // test yac_interp_stack_config_entry_get_spmap
  union yac_interp_stack_config_entry const * interp_stack_entry =
    yac_interp_stack_config_get_entry(
      interp_stack_config,
      yac_interp_stack_config_get_size(interp_stack_config) - 1);
  double spread_distance_;
  double max_search_distance_;
  enum yac_interp_spmap_weight_type weight_type_;
  enum yac_interp_spmap_scale_type scale_type_;
  double src_sphere_radius_;
  char const * src_filename_;
  char const * src_varname_;
  int src_min_global_id_;
  double tgt_sphere_radius_;
  char const * tgt_filename_;
  char const * tgt_varname_;
  int tgt_min_global_id_;
  yac_interp_stack_config_entry_get_spmap(
    interp_stack_entry,
    &spread_distance_, &max_search_distance_, &weight_type_, &scale_type_,
    &src_sphere_radius_, &src_filename_, &src_varname_, &src_min_global_id_,
    &tgt_sphere_radius_, &tgt_filename_, &tgt_varname_, &tgt_min_global_id_);

  if (spread_distance_ != spread_distance)
    PUT_ERR("ERROR in yac_interp_stack_config_entry_get_spmap");
  if (max_search_distance_ != max_search_distance)
    PUT_ERR("ERROR in yac_interp_stack_config_entry_get_spmap");
  if (weight_type_ != weight_type)
    PUT_ERR("ERROR in yac_interp_stack_config_entry_get_spmap");
  if (scale_type_ != scale_type)
    PUT_ERR("ERROR in yac_interp_stack_config_entry_get_spmap");
}

static void yac_interp_stack_config_add_spmap_ext_(
  struct yac_interp_stack_config * interp_stack_config,
  struct yac_interp_spmap_config * default_config,
  struct yac_spmap_overwrite_config ** overwrite_configs) {

  yac_interp_stack_config_add_spmap_ext(
    interp_stack_config, default_config, overwrite_configs);
}
