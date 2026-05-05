// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "utils_core.h"
#include "yac.h"
#include "instance.h"
#include "event.h"
#include "yac_mpi_common.h"
#include "fields.h"
#include "component.h"
#include "config_yaml.h"
#include "interpolation_exchange.h"

enum yac_instance_phase {
  INSTANCE_DEFINITION = 0, // after yac_cinit
  INSTANCE_DEFINITION_COMP = 1, // after yac_cdef_comp
  INSTANCE_DEFINITION_SYNC = 2, // after yac_csync_def
  INSTANCE_EXCHANGE = 3, // after yac_cenddef
  INSTANCE_UNKNOWN = 4,
};

// bit mask values to determine which field is available
enum field_availability_type {
  SOURCE_FLAG = 1 << 0, // first bit (source is available)
  TARGET_FLAG = 1 << 1, // second bit (target is available)
  SOURCE_TARGET_FLAG = (1 << 0) | (1 << 1), // first two bits
                                            // (both source and target are
                                            //  available)
};

static const char * yac_instance_phase_str[] =
  {"definition phase",
   "definition phase (after component definition)",
   "definition phase (after synchronisation)",
   "exchange phase",
   "unknown phase"};

#define CHECK_PHASE(FUNC_NAME, REF_PHASE, NEW_PHASE) \
  { \
    enum yac_instance_phase ref_phase_ = (REF_PHASE); \
    YAC_ASSERT_F( \
      instance->phase == (ref_phase_), \
      "ERROR(%s): Invalid phase " \
      "(current phase: \"%s\" expected phase: \"%s\")", \
      #FUNC_NAME, yac_instance_phase_str[instance->phase], \
      yac_instance_phase_str[(ref_phase_)]); \
    instance->phase = (NEW_PHASE); \
  }
#define CHECK_MIN_PHASE(FUNC_NAME, MIN_REF_PHASE) \
  { \
    enum yac_instance_phase ref_min_phase_ = (MIN_REF_PHASE); \
    YAC_ASSERT_F( \
      instance->phase >= (ref_min_phase_), \
      "ERROR(%s): Invalid phase " \
      "(current phase: \"%s\" minimum expected phase: \"%s\")", \
      #FUNC_NAME, yac_instance_phase_str[instance->phase], \
      yac_instance_phase_str[(ref_min_phase_)]); \
  }
#define CHECK_MAX_PHASE(FUNC_NAME, MAX_REF_PHASE) \
  { \
    enum yac_instance_phase ref_max_phase_ = (MAX_REF_PHASE); \
    YAC_ASSERT_F( \
      instance->phase <= (ref_max_phase_), \
      "ERROR(%s): Invalid phase " \
      "(current phase: \"%s\" maximum expected phase: \"%s\")", \
      #FUNC_NAME, \
      yac_instance_phase_str[ \
        MIN(instance->phase,INSTANCE_UNKNOWN)], \
      yac_instance_phase_str[(ref_max_phase_)]); \
  }

struct yac_instance {
  /// Coupling configuration data of this YAC instance
  struct yac_couple_config * couple_config;

  /// Component configuration data of this YAC instance
  struct yac_component_config * comp_config;

  /// Coupling fields added to this YAC instance via yac_instance_add_field
  struct coupling_field ** cpl_fields;

  /// Number of elements in coupling_field ** cpl_fields
  size_t num_cpl_fields;

  /// MPI communicator that contains the processes of this YAC instance
  MPI_Comm comm;

  /// Current phase of this YAC instance
  enum yac_instance_phase phase;
};

enum field_type {
  SRC = 1,
  TGT = 2
};

struct comp_grid_config {
  const char * grid_name;
  const char * comp_name;
};

struct comp_grid_pair_config {
  struct comp_grid_config config[2];
};

static struct field_config_event_data {
  char const * timestep;
  char const * coupling_period;
  int timelag;
  enum yac_reduction_type reduction_operation;
  char const * start_datetime;
  char const * end_datetime;
} empty_event_data =
  {.timestep = NULL,
   .coupling_period = NULL,
   .timelag = 0,
   .reduction_operation = TIME_NONE,
   .start_datetime = NULL,
   .end_datetime = NULL};

struct src_field_config {
  struct coupling_field * field;
  char const * name;

  struct yac_interp_field * interp_fields;
  size_t num_interp_fields;

  struct yac_interp_stack_config * interp_stack;
  const char * weight_file_name;
  enum yac_weight_file_on_existing weight_file_on_existing;
  struct field_config_event_data event_data;

  int avail_config_id;
};

struct tgt_field_config {
  struct coupling_field * field;
  char const * name;

  struct yac_interp_field * interp_field;

  struct field_config_event_data event_data;

  int avail_config_id;
};

struct field_config {

  struct comp_grid_pair_config comp_grid_pair;
  int src_comp_idx;

  struct src_field_config src_interp_config;
  struct tgt_field_config tgt_interp_config;

  enum yac_interp_weights_reorder_type reorder_type;

  double frac_mask_fallback_value;

  double scale_factor, scale_summand;

  char const * yaxt_exchanger_name;

  int use_raw_exchange;

  size_t collection_size;
};

struct output_grid {
  char const * grid_name;
  char const * filename;
  struct yac_basic_grid * grid;
};

struct dist_flag {
  uint64_t * data;
  size_t count;
  size_t idx;
};

static struct yac_basic_grid * get_basic_grid(
  const char * grid_name, struct yac_basic_grid ** grids, size_t num_grids,
  int * delete_flag) {

  struct yac_basic_grid * grid = NULL;
  for (size_t i = 0; (i < num_grids) && (grid == NULL); ++i)
    if (!strcmp(grid_name, yac_basic_grid_get_name(grids[i])))
      grid = grids[i];

  *delete_flag = grid == NULL;
  return (grid == NULL)?yac_basic_grid_empty_new(grid_name):grid;
}

static
int compare_comp_grid_config(const void * a, const void * b) {

  struct comp_grid_config * a_ = (struct comp_grid_config *)a,
                          * b_ = (struct comp_grid_config *)b;

  int ret;
  if ((ret = strcmp(a_->comp_name, b_->comp_name))) return ret;
  else return strcmp(a_->grid_name, b_->grid_name);
}

static struct coupling_field * get_coupling_field(
  char const * component_name, const char * field_name,
  const char * grid_name, size_t num_fields,
  struct coupling_field ** coupling_fields) {

  for (size_t i = 0; i < num_fields; ++i) {
    struct coupling_field * curr_field = coupling_fields[i];
    if (!strcmp(component_name, yac_get_coupling_field_comp_name(curr_field)) &&
        !strcmp(field_name, yac_get_coupling_field_name(curr_field)) &&
        !strcmp(grid_name,
                yac_basic_grid_get_name(
                  yac_coupling_field_get_basic_grid(curr_field))))
      return curr_field;
  }

  return NULL;
}

static int compare_field_config_interpolation_build_config(
  struct field_config * a, struct field_config * b) {

  int ret;

  if ((ret = (a->reorder_type > b->reorder_type) -
             (a->reorder_type < b->reorder_type))) return ret;
  if ((ret = (a->collection_size > b->collection_size) -
             (a->collection_size < b->collection_size))) return ret;
  if ((ret = (a->scale_factor > b->scale_factor) -
             (a->scale_factor < b->scale_factor))) return ret;
  if ((ret = (a->scale_summand > b->scale_summand) -
             (a->scale_summand < b->scale_summand))) return ret;
  if ((ret = (int)(a->yaxt_exchanger_name == NULL) -
             (int)(b->yaxt_exchanger_name == NULL))) return ret;
  if ((ret = (a->yaxt_exchanger_name != NULL) &&
             strcmp(a->yaxt_exchanger_name, b->yaxt_exchanger_name)))
    return ret;
  if ((ret = (a->use_raw_exchange > b->use_raw_exchange) -
             (a->use_raw_exchange < b->use_raw_exchange))) return ret;
  if ((ret = memcmp(
               &(a->frac_mask_fallback_value), &(b->frac_mask_fallback_value),
               sizeof(double)))) return ret;
  if ((ret = (a->src_interp_config.avail_config_id -
              b->src_interp_config.avail_config_id))) return ret;
  if ((ret = (a->tgt_interp_config.avail_config_id -
              b->tgt_interp_config.avail_config_id))) return ret;

  return 0;
}

struct field_config_event_data get_event_data(
  struct yac_instance * instance, int couple_idx, int field_couple_idx,
  enum field_type field_type) {

  struct yac_couple_config * couple_config = instance->couple_config;

  enum yac_reduction_type reduction_operation =
    yac_couple_config_get_coupling_period_operation(
      couple_config, couple_idx, field_couple_idx);

  char const * coupling_period =
    yac_couple_config_get_coupling_period(
      couple_config, couple_idx, field_couple_idx);
  char const * timestep;
  int timelag;
  if ( field_type == SRC ){
    timestep =
      yac_couple_config_get_source_timestep(
        couple_config, couple_idx, field_couple_idx);
    timelag =
      yac_couple_config_get_source_lag(
        couple_config, couple_idx, field_couple_idx);
  }else{
    reduction_operation = TIME_NONE;
    timestep =
      yac_couple_config_get_target_timestep(
        couple_config, couple_idx, field_couple_idx);
    timelag =
      yac_couple_config_get_target_lag(
        couple_config, couple_idx, field_couple_idx);
  }

  return
    (struct field_config_event_data)
      {.timestep = timestep,
       .coupling_period = coupling_period,
       .timelag = timelag,
       .reduction_operation = reduction_operation,
       .start_datetime = yac_couple_config_get_start_datetime(couple_config),
       .end_datetime = yac_couple_config_get_end_datetime(couple_config)};
}

static struct event * generate_event(
  struct field_config_event_data event_data) {

  struct event * event = yac_event_new();
  yac_event_add(
    event, event_data.timestep, event_data.coupling_period,
    event_data.timelag, event_data.reduction_operation,
    event_data.start_datetime, event_data.end_datetime);
  return event;
}


// checks for all coupling configurations whether the associated source and
// target field is available on any process
static int * determine_valid_field_configurations(
  struct yac_couple_config * couple_config, MPI_Comm comm,
  struct coupling_field ** coupling_fields, size_t num_fields) {

  size_t num_couples = yac_couple_config_get_num_couples(couple_config);
  size_t total_num_fields = 0;
  for (size_t couple_idx = 0; couple_idx < num_couples; ++couple_idx)
    total_num_fields +=
      yac_couple_config_get_num_couple_fields(couple_config, couple_idx);

  int * is_valid_field_configuration =
    xcalloc(total_num_fields, sizeof(*is_valid_field_configuration));

  // for all coupling configurations
  for (size_t couple_idx = 0, i = 0; couple_idx < num_couples;
       ++couple_idx) {

    // for all field configuration of the current coupling configuration
    size_t curr_num_fields =
      yac_couple_config_get_num_couple_fields(couple_config, couple_idx);
    for (size_t field_couple_idx = 0; field_couple_idx < curr_num_fields;
         ++field_couple_idx, ++i) {

      // get component, grid, and field names for the current
      // field configuration
      char const * src_component_name;
      char const * tgt_component_name;
      char const * src_grid_name;
      char const * tgt_grid_name;
      const char * src_field_name;
      const char * tgt_field_name;
      yac_couple_config_get_field_couple_component_names(
        couple_config, couple_idx, field_couple_idx,
        &src_component_name, &tgt_component_name);
      yac_couple_config_get_field_grid_names(
        couple_config, couple_idx, field_couple_idx,
        &src_grid_name, &tgt_grid_name);
      yac_couple_config_get_field_names(
        couple_config, couple_idx, field_couple_idx,
          &src_field_name, &tgt_field_name);

      // check availability of respective source and target fields
      if (get_coupling_field(
            src_component_name, src_field_name, src_grid_name,
            num_fields, coupling_fields) != NULL)
        is_valid_field_configuration[i] |= SOURCE_FLAG;
      if (get_coupling_field(
            tgt_component_name, tgt_field_name, tgt_grid_name,
            num_fields, coupling_fields) != NULL)
        is_valid_field_configuration[i] |= TARGET_FLAG;
    }
  }

  // check availability of source and target fields across all processes
  // (using "bit-wise or" reduction operation)
  yac_mpi_call(
    MPI_Allreduce(
      MPI_IN_PLACE, is_valid_field_configuration, total_num_fields,
      MPI_INT, MPI_BOR, comm), comm);

  int comm_rank;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  int const is_root = comm_rank == 0;

  // check for missing fields print warning or abort
  // (depending on configuration of missing_definition_is_fatal)
  int missing_definition_is_fatal =
    yac_couple_config_get_missing_definition_is_fatal(couple_config);
  for (size_t couple_idx = 0, i = 0; couple_idx < num_couples;
       ++couple_idx) {
    size_t curr_num_fields =
      yac_couple_config_get_num_couple_fields(couple_config, couple_idx);
    for (size_t field_couple_idx = 0; field_couple_idx < curr_num_fields;
         ++field_couple_idx, ++i) {
      if ((is_valid_field_configuration[i] != SOURCE_TARGET_FLAG) &&
          is_root) {
        char const * src_component_name;
        char const * tgt_component_name;
        char const * src_grid_name;
        char const * tgt_grid_name;
        const char * src_field_name;
        const char * tgt_field_name;
        yac_couple_config_get_field_couple_component_names(
          couple_config, couple_idx, field_couple_idx,
          &src_component_name, &tgt_component_name);
        yac_couple_config_get_field_grid_names(
          couple_config, couple_idx, field_couple_idx,
          &src_grid_name, &tgt_grid_name);
        yac_couple_config_get_field_names(
          couple_config, couple_idx, field_couple_idx,
            &src_field_name, &tgt_field_name);
        fprintf(stderr, "%s: couple defined for field: \n"
                        "  source (%s):\n"
                        "    component name: \"%s\"\n"
                        "    grid name:      \"%s\"\n"
                        "    field name:     \"%s\"\n"
                        "  target(%s):\n"
                        "    component name: \"%s\"\n"
                        "    grid name:      \"%s\"\n"
                        "    field name:     \"%s\"\n",
                missing_definition_is_fatal?"ERROR":"WARNING",
                (is_valid_field_configuration[i] & SOURCE_FLAG)?
                  "defined":"not defined",
                src_component_name, src_grid_name, src_field_name,
                (is_valid_field_configuration[i] & TARGET_FLAG)?
                  "defined":"not defined",
                tgt_component_name, tgt_grid_name, tgt_field_name);
        YAC_ASSERT(
          !missing_definition_is_fatal,
          "ERROR(get_field_configuration): missing definition")
      }
      is_valid_field_configuration[i] =
        (SOURCE_TARGET_FLAG == is_valid_field_configuration[i]);
    }
  }

  return is_valid_field_configuration;
}

static void get_interp_fields_from_coupling_field(
  struct coupling_field * field, char const * const * mask_names,
  size_t num_mask_names, struct yac_interp_field ** interp_fields,
  size_t * num_fields, MPI_Comm comm) {

  struct yac_interp_field * interp_fields_;
  size_t num_fields_;

  if (field != NULL) {

    num_fields_ = yac_coupling_field_get_num_interp_fields(field);

    YAC_ASSERT_F(
      (num_mask_names == 0) || (num_fields_ == num_mask_names),
      "ERROR(get_interp_fields_from_coupling_field): "
      "missmatch in number of interpolation fields of coupling field \"%s\" "
      "and number of provided mask names (%zu != %zu)",
      yac_get_coupling_field_name(field), num_fields_, num_mask_names)

    uint64_t local_counts[2] =
      {(uint64_t)num_fields_, (uint64_t)num_mask_names};
    uint64_t global_counts[2];

    yac_mpi_call(
      MPI_Allreduce(
        local_counts, global_counts, 2, MPI_UINT64_T, MPI_MAX, comm), comm);

    YAC_ASSERT_F(
      local_counts[0] == global_counts[0],
      "ERROR(get_interp_fields_from_coupling_field): missmatch in number of"
      "local interpolation fields for coupling field \"%s\" and global "
      "(%zu != %zu)",
      yac_get_coupling_field_name(field), num_fields_, (size_t)(global_counts[0]))

    YAC_ASSERT_F(
      (num_mask_names != 0) || (global_counts[1] == 0),
      "ERROR(get_interp_fields_from_coupling_field): local process did "
      "not provide mask names for coupling field \"%s\" while others did",
      yac_get_coupling_field_name(field));

    // make a copy of the interpolation fields of the coupling field
    interp_fields_ = xmalloc(num_fields_ * sizeof(*interp_fields_));
    memcpy(
      interp_fields_, yac_coupling_field_get_interp_fields(field),
      num_fields_ * sizeof(*interp_fields_));

    // if mask names are provided, overwrite already existing masks
    struct yac_basic_grid * grid = yac_coupling_field_get_basic_grid(field);
    for (size_t i = 0; i < num_mask_names; ++i) {
      char const * mask_name = mask_names[i];
      YAC_ASSERT_F(
        mask_name != NULL,
        "ERROR(get_interp_fields_from_coupling_field): "
        "make_names[%zu] is NULL", i);
      interp_fields_[i].masks_idx =
        yac_basic_grid_get_named_mask_idx(
          grid, interp_fields_[i].location, mask_name);
    }

    uint64_t data[num_fields_][3];

    for (size_t i = 0; i < num_fields_; ++i) {
      data[i][0] = (uint64_t)interp_fields_[i].location;
      data[i][1] = (uint64_t)interp_fields_[i].coordinates_idx;
      data[i][2] = (uint64_t)interp_fields_[i].masks_idx;
    }

    yac_mpi_call(
      MPI_Allreduce(
        MPI_IN_PLACE, data, 3 * (int)num_fields_,
        MPI_UINT64_T, MPI_MIN, comm), comm);

    for (size_t i = 0; i < num_fields_; ++i) {
      YAC_ASSERT(
        data[i][0] == (uint64_t)(interp_fields_[i].location),
        "ERROR(get_interp_fields_from_coupling_field): location mismatch")
      YAC_ASSERT(
        data[i][1] == (uint64_t)(interp_fields_[i].coordinates_idx),
        "ERROR(get_interp_fields_from_coupling_field): "
        "coordinates index mismatch")
      YAC_ASSERT(
        data[i][2] == (uint64_t)(interp_fields_[i].masks_idx),
        "ERROR(get_interp_fields_from_coupling_field): "
        "masks index mismatch")
    }

  } else {

    uint64_t zero_counts[2] = {0,0};
    uint64_t counts[2];

    yac_mpi_call(
      MPI_Allreduce(
        zero_counts, counts, 2,
        MPI_UINT64_T, MPI_MAX, comm), comm);

    num_fields_ = (size_t)(counts[0]);
    interp_fields_ = xmalloc(num_fields_ * sizeof(*interp_fields_));

    uint64_t data[num_fields_][3];

    for (size_t i = 0; i < num_fields_; ++i) {
      data[i][0] = (uint64_t)YAC_LOC_UNDEFINED;
      data[i][1] = (uint64_t)UINT64_MAX;
      data[i][2] = (uint64_t)UINT64_MAX;
    }

    yac_mpi_call(
      MPI_Allreduce(
        MPI_IN_PLACE, data, 3 * (int)num_fields_,
        MPI_UINT64_T, MPI_MIN, comm), comm);

    for (size_t i = 0; i < num_fields_; ++i) {
      interp_fields_[i].location = (enum yac_location)data[i][0];
      interp_fields_[i].coordinates_idx = (size_t)data[i][1];
      interp_fields_[i].masks_idx = (size_t)data[i][2];
    }
  }

  *interp_fields = interp_fields_;
  *num_fields = num_fields_;
}

static struct src_field_config get_src_interp_config(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx,
  struct coupling_field * field, MPI_Comm comm) {

  char const * const * src_mask_names;
  size_t num_src_mask_names;
  yac_couple_config_get_src_mask_names(
    couple_config, couple_idx, field_couple_idx,
    &src_mask_names, &num_src_mask_names);

  struct yac_interp_field * interp_fields;
  size_t num_interp_fields;
  get_interp_fields_from_coupling_field(
    field, src_mask_names, num_src_mask_names,
    &interp_fields, &num_interp_fields, comm);

  return
    (struct src_field_config) {
    .field = field,
    .name = NULL,
    .interp_fields = interp_fields,
    .num_interp_fields = num_interp_fields,
    .interp_stack =
      yac_couple_config_get_interp_stack(
        couple_config, couple_idx, field_couple_idx),
    .weight_file_name =
      (yac_couple_config_enforce_write_weight_file(
         couple_config, couple_idx, field_couple_idx))?
        yac_couple_config_get_weight_file_name(
          couple_config, couple_idx, field_couple_idx):NULL,
    .weight_file_on_existing =
      (yac_couple_config_enforce_write_weight_file(
         couple_config, couple_idx, field_couple_idx))?
        yac_couple_config_get_weight_file_on_existing(
          couple_config, couple_idx, field_couple_idx):
        YAC_WEIGHT_FILE_ON_EXISTING_DEFAULT_VALUE,
    .event_data = empty_event_data,
  };
}

static struct tgt_field_config get_tgt_interp_config(
  struct yac_couple_config * couple_config,
  size_t couple_idx, size_t field_couple_idx,
  struct coupling_field * field, MPI_Comm comm) {

  char const * mask_name =
    yac_couple_config_get_tgt_mask_name(
      couple_config, couple_idx, field_couple_idx);

  struct yac_interp_field * interp_field;
  size_t num_interp_field;
  get_interp_fields_from_coupling_field(
    field, &mask_name, mask_name != NULL,
    &interp_field, &num_interp_field, comm);

  YAC_ASSERT(
    num_interp_field == 1,
    "ERROR(get_tgt_interp_config): "
    "only one point set per target field supported")

  return
    (struct tgt_field_config) {
    .field = field,
    .name = NULL,
    .interp_field = interp_field,
    .event_data = empty_event_data,
  };
}

static int compare_dist_flags(const void * a_, const void * b_) {

  struct dist_flag const * a = a_;
  struct dist_flag const * b = b_;

  return memcmp(a->data, b->data, a->count * sizeof(a->data[0]));
}

static void generate_dist_flag_config_ids(
  uint64_t * local_flags, size_t flag_count,
  uint64_t * global_flags_buffer, struct dist_flag * dist_flags,
  int * dist_flag_ids, MPI_Comm comm) {

  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  // number of 64-bit blocks required to store all flags
  size_t flags_num_blocks = (flag_count + 63) / 64;
  // number of 64-bit blocks required to store the value of all ranks for
  // one flag
  size_t ranks_num_blocks = (size_t)((comm_size + 63) / 64);

  // all-gather the flags from all ranks
  yac_mpi_call(
    MPI_Allgather(
      local_flags, (int)flags_num_blocks, MPI_UINT64_T,
      global_flags_buffer, (int)flags_num_blocks, MPI_UINT64_T,
      comm), comm);

  // initialise per-flag storage of all data
  for (size_t flag_idx = 0; flag_idx < flag_count; ++flag_idx) {
    memset(
      dist_flags[flag_idx].data, 0,
      ranks_num_blocks * sizeof(dist_flags[flag_idx].data[0]));
    dist_flags[flag_idx].idx = flag_idx;
    dist_flags[flag_idx].count = ranks_num_blocks;
  }

  // converte from storing flags on per-rank basis to per-flag
  for (int rank = 0; rank < comm_size; ++rank) {

    // get flags for current rank
    uint64_t * rank_flags =
      global_flags_buffer + (size_t)rank * flags_num_blocks;

    size_t rank_block_idx = rank / 64;
    uint64_t rank_flag_mask = ((uint64_t)1) << (rank % 64);

    for (size_t flag_idx = 0; flag_idx < flag_count; ++flag_idx)
      if (rank_flags[flag_idx/64] & (((uint64_t)1) << (flag_idx % 64)))
        dist_flags[flag_idx].data[rank_block_idx] |= rank_flag_mask;
  }

  // sort flags based on availability of on all ranks
  qsort(dist_flags, flag_count, sizeof(*dist_flags), compare_dist_flags);

  // determine unique configurations
  int id = 0;
  struct dist_flag * prev_dist_flag = dist_flags;
  for (size_t flag_idx = 0; flag_idx < flag_count; ++flag_idx) {

    struct dist_flag * curr_dist_flag = dist_flags + flag_idx;
    if (compare_dist_flags(prev_dist_flag, curr_dist_flag)) ++id;

    dist_flag_ids[curr_dist_flag->idx] = id;
    prev_dist_flag = curr_dist_flag;
  }
}

static void generate_coupling_field_avail_config_ids(
  struct field_config * field_configs, size_t num_fields,
  int * is_valid_field_configuration, MPI_Comm comm) {

  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  // number of 64-bit blocks required to store the flags of all fields
  size_t flags_num_blocks = (num_fields + 63) / 64;
  // number of 64-bit blocks required to store the flags of all ranks for
  // one field
  size_t ranks_num_blocks = (size_t)((comm_size + 63) / 64);

  // flags for storing the availability of source/target fields of the local
  // process
  uint64_t * field_is_available_flags =
    xmalloc(flags_num_blocks * sizeof(*field_is_available_flags));

  // flags for storing the availability of source/target fields af all processes
  // (per rank)
  uint64_t * all_field_is_available_flags =
    xmalloc(
      (size_t)comm_size * flags_num_blocks * sizeof(*all_field_is_available_flags));

  // flags for storing the availability of source/target fields af all processes
  // (per field)
  struct dist_flag * dist_flags = xmalloc(num_fields * sizeof(*dist_flags));
  uint64_t * dist_flags_data_buffer =
    xmalloc(num_fields * ranks_num_blocks * sizeof(*dist_flags_data_buffer));
  for (size_t i = 0; i < num_fields; ++i) {
    dist_flags[i].data = dist_flags_data_buffer + i * ranks_num_blocks;
  }

  // individual ids for all flag configurations
  int * field_avail_config_ids =
    xmalloc(num_fields * sizeof(*field_avail_config_ids));

  { // source field availability

    // determine available of a field on the local process
    memset(
      field_is_available_flags, 0,
      flags_num_blocks * sizeof(*field_is_available_flags));
    for (size_t field_idx = 0; field_idx < num_fields; ++field_idx)
      if (is_valid_field_configuration[field_idx] &&
          (field_configs[field_idx].src_interp_config.field != NULL))
        field_is_available_flags[field_idx/64] |= (1 << (field_idx%64));

    generate_dist_flag_config_ids(
      field_is_available_flags, num_fields, all_field_is_available_flags,
      dist_flags, field_avail_config_ids, comm);

    for (size_t field_idx = 0; field_idx < num_fields; ++field_idx)
      field_configs[field_idx].src_interp_config.avail_config_id =
        field_avail_config_ids[field_idx];
  }

  { // target field availability

    // determine available of a field on the local process
    memset(
      field_is_available_flags, 0,
      flags_num_blocks * sizeof(*field_is_available_flags));
    for (size_t field_idx = 0; field_idx < num_fields; ++field_idx)
      if (is_valid_field_configuration[field_idx] &&
          (field_configs[field_idx].tgt_interp_config.field != NULL))
        field_is_available_flags[field_idx/64] |= (1 << (field_idx%64));

    generate_dist_flag_config_ids(
      field_is_available_flags, num_fields, all_field_is_available_flags,
      dist_flags, field_avail_config_ids, comm);

    for (size_t field_idx = 0; field_idx < num_fields; ++field_idx)
      field_configs[field_idx].tgt_interp_config.avail_config_id =
        field_avail_config_ids[field_idx];
  }

  free(field_avail_config_ids);
  free(dist_flags_data_buffer);
  free(dist_flags);
  free(all_field_is_available_flags);
  free(field_is_available_flags);
}

static void get_field_configuration(
  struct yac_instance * instance,
  struct field_config ** field_configs_, size_t * count) {

  struct yac_couple_config * couple_config = instance->couple_config;
  MPI_Comm comm = instance->comm;
  size_t num_fields = instance->num_cpl_fields;
  struct coupling_field ** coupling_fields = instance->cpl_fields;

  size_t num_couples = yac_couple_config_get_num_couples(couple_config);

  // determines for which coupling configurations the source and target fields
  // are available
  int * is_valid_field_configuration =
    determine_valid_field_configurations(
      couple_config, comm, coupling_fields, num_fields);

  size_t total_num_fields = 0;
  for (size_t couple_idx = 0, i = 0; couple_idx < num_couples; ++couple_idx) {
    size_t curr_num_fields =
      yac_couple_config_get_num_couple_fields(couple_config, couple_idx);
    for (size_t field_couple_idx = 0; field_couple_idx < curr_num_fields;
         ++field_couple_idx, ++i)
      if (is_valid_field_configuration[i]) ++total_num_fields;
  }

  // get all field coupling configurations
  // * due to the synchronisation of couple_config beforehand the ordering
  //   of the field coupling configurations on all processes is identical
  // * all information extracted from the coupling configuration is identical
  //   on all processes
  // * only valid configurations (source and target field are defined somewhere)
  size_t field_config_idx = 0;
  struct field_config * field_configs =
    xmalloc(total_num_fields * sizeof(*field_configs));
  for (size_t couple_idx = 0, i = 0; couple_idx < num_couples; ++couple_idx) {

    size_t curr_num_fields =
      yac_couple_config_get_num_couple_fields(couple_config, couple_idx);

    for (size_t field_couple_idx = 0; field_couple_idx < curr_num_fields;
         ++field_couple_idx, ++i) {

      if (!is_valid_field_configuration[i]) continue;

      struct comp_grid_config src_config, tgt_config;

      yac_couple_config_get_field_couple_component_names(
        couple_config, couple_idx, field_couple_idx,
        &(src_config.comp_name), &(tgt_config.comp_name));

      yac_couple_config_get_field_grid_names(
        couple_config, couple_idx, field_couple_idx,
        &(src_config.grid_name), &(tgt_config.grid_name));

      int src_comp_idx =
        compare_comp_grid_config(&src_config, &tgt_config) > 0;

      const char * src_field_name;
      const char * tgt_field_name;
      yac_couple_config_get_field_names(
        couple_config, couple_idx, field_couple_idx,
          &src_field_name, &tgt_field_name);
      struct coupling_field * src_field =
        get_coupling_field(
          src_config.comp_name, src_field_name, src_config.grid_name,
          num_fields, coupling_fields);
      struct coupling_field * tgt_field =
        get_coupling_field(
          tgt_config.comp_name, tgt_field_name, tgt_config.grid_name,
          num_fields, coupling_fields);

      double frac_mask_fallback_value =
        yac_couple_config_get_frac_mask_fallback_value(
          couple_config, src_config.comp_name, src_config.grid_name,
          src_field_name);
      double scale_factor =
        yac_couple_config_get_scale_factor(
          couple_config, couple_idx, field_couple_idx);
      double scale_summand =
        yac_couple_config_get_scale_summand(
          couple_config, couple_idx, field_couple_idx);
      char const * yaxt_exchanger_name =
        yac_couple_config_get_yaxt_exchanger_name(
          couple_config, couple_idx, field_couple_idx);
      int use_raw_exchange =
        yac_couple_config_get_use_raw_exchange(
          couple_config, couple_idx, field_couple_idx);
      size_t collection_size =
        yac_couple_config_get_field_collection_size(
          couple_config, src_config.comp_name, src_config.grid_name,
          src_field_name);

      YAC_ASSERT_F(
        collection_size ==
        yac_couple_config_get_field_collection_size(
          couple_config, tgt_config.comp_name, tgt_config.grid_name,
          tgt_field_name),
        "ERROR: collection sizes do not match for coupled fields (%zu != %zu): \n"
        "  source:\n"
        "    component name: \"%s\"\n"
        "    grid name:      \"%s\"\n"
        "    field name:     \"%s\"\n"
        "  target:\n"
        "    component name: \"%s\"\n"
        "    grid name:      \"%s\"\n"
        "    field name:     \"%s\"\n",
        collection_size,
        yac_couple_config_get_field_collection_size(
          couple_config, tgt_config.comp_name, tgt_config.grid_name,
          tgt_field_name),
        src_config.comp_name, src_config.grid_name, src_field_name,
        tgt_config.comp_name, tgt_config.grid_name, tgt_field_name);

      field_configs[field_config_idx].comp_grid_pair.config[src_comp_idx] =
        src_config;
      field_configs[field_config_idx].comp_grid_pair.config[src_comp_idx^1] =
        tgt_config;
      field_configs[field_config_idx].src_comp_idx = src_comp_idx;
      field_configs[field_config_idx].src_interp_config =
        get_src_interp_config(
          couple_config, couple_idx, field_couple_idx, src_field, comm);
      field_configs[field_config_idx].tgt_interp_config =
        get_tgt_interp_config(
          couple_config, couple_idx, field_couple_idx, tgt_field, comm);
      field_configs[field_config_idx].reorder_type =
        (yac_couple_config_mapping_on_source(
           couple_config, couple_idx, field_couple_idx))?
        (YAC_MAPPING_ON_SRC):(YAC_MAPPING_ON_TGT);
      field_configs[field_config_idx].src_interp_config.event_data =
        get_event_data(instance, couple_idx, field_couple_idx, SRC);
      field_configs[field_config_idx].tgt_interp_config.event_data =
        get_event_data(instance, couple_idx, field_couple_idx, TGT);
      field_configs[field_config_idx].frac_mask_fallback_value =
        frac_mask_fallback_value;
      field_configs[field_config_idx].scale_factor = scale_factor;
      field_configs[field_config_idx].scale_summand = scale_summand;
      field_configs[field_config_idx].yaxt_exchanger_name = yaxt_exchanger_name;
      field_configs[field_config_idx].use_raw_exchange = use_raw_exchange;
      field_configs[field_config_idx].collection_size = collection_size;
      field_configs[field_config_idx].src_interp_config.name = src_field_name;
      field_configs[field_config_idx].tgt_interp_config.name = tgt_field_name;
      ++field_config_idx;
    }
  }

  // computes unique ids for availability configuration of source and target
  // coupling fields
  generate_coupling_field_avail_config_ids(
    field_configs, total_num_fields, is_valid_field_configuration, comm);

  free(is_valid_field_configuration);

  *field_configs_ = field_configs;
  *count = total_num_fields;
}

static int compare_interp_field(
  struct yac_interp_field * a, struct yac_interp_field * b) {

  int ret;
  if ((ret = (a->location > b->location) - (a->location < b->location)))
    return ret;
  if ((ret =
         (a->coordinates_idx > b->coordinates_idx) -
         (a->coordinates_idx < b->coordinates_idx)))
    return ret;
  return (a->masks_idx > b->masks_idx) - (a->masks_idx < b->masks_idx);
}

static int compare_field_config_interp_fields(
  struct field_config * a, struct field_config * b) {

  int ret;
  if ((ret = strcmp(a->comp_grid_pair.config[a->src_comp_idx].grid_name,
                    b->comp_grid_pair.config[b->src_comp_idx].grid_name)))
    return ret;
  if ((ret = strcmp(a->comp_grid_pair.config[a->src_comp_idx^1].grid_name,
                    b->comp_grid_pair.config[b->src_comp_idx^1].grid_name)))
    return ret;
  if ((ret =
         (a->src_interp_config.num_interp_fields >
          b->src_interp_config.num_interp_fields) -
         (a->src_interp_config.num_interp_fields <
          b->src_interp_config.num_interp_fields))) return ret;
  for (size_t i = 0; i < a->src_interp_config.num_interp_fields; ++i)
    if ((ret =
           compare_interp_field(
             a->src_interp_config.interp_fields + i,
             b->src_interp_config.interp_fields + i))) return ret;
  return
    compare_interp_field(
      a->tgt_interp_config.interp_field, b->tgt_interp_config.interp_field);
}

static
int compare_field_config(const void * a, const void * b) {

  struct field_config * a_ = (struct field_config *)a,
                      * b_ = (struct field_config *)b;

  int ret;
  if ((ret = strcmp(a_->comp_grid_pair.config[0].grid_name,
                    b_->comp_grid_pair.config[0].grid_name))) return ret;
  if ((ret = strcmp(a_->comp_grid_pair.config[1].grid_name,
                    b_->comp_grid_pair.config[1].grid_name))) return ret;
  if ((ret = compare_field_config_interp_fields(a_, b_))) return ret;
  if ((ret = yac_interp_stack_config_compare(
               a_->src_interp_config.interp_stack,
               b_->src_interp_config.interp_stack))) return ret;
  if ((ret = (int)(a_->src_interp_config.weight_file_name == NULL) -
             (int)(b_->src_interp_config.weight_file_name == NULL))) return ret;
  if ((a_->src_interp_config.weight_file_name != NULL) &&
      (ret = strcmp(a_->src_interp_config.weight_file_name,
                    b_->src_interp_config.weight_file_name))) return ret;
  if ((a_->src_interp_config.weight_file_name != NULL) &&
      (ret = (int)(a_->src_interp_config.weight_file_on_existing) -
             (int)(b_->src_interp_config.weight_file_on_existing))) return ret;
  if ((ret = compare_field_config_interpolation_build_config(a_, b_))) return ret;
  if ((ret = strcmp(a_->comp_grid_pair.config[0].comp_name,
                    b_->comp_grid_pair.config[0].comp_name))) return ret;
  if ((ret = strcmp(a_->src_interp_config.name,
                    b_->src_interp_config.name))) return ret;
  if ((ret = strcmp(a_->tgt_interp_config.name,
                    b_->tgt_interp_config.name))) return ret;
  YAC_ASSERT_F(
    ret, "ERROR(compare_field_config): "
    "duplicated coupling field configuration detected:\n"
    "\tcomponent name: \n"
    "\t        source:\"%s\"\n"
    "\t        target:\"%s\"\n"
    "\t     grid name:\n"
    "\t        source:\"%s\"\n"
    "\t        target:\"%s\"\n"
    "\t    field name:\n"
    "\t        source:\"%s\"\n"
    "\t        target:\"%s\"\n",
    a_->comp_grid_pair.config[a_->src_comp_idx  ].comp_name,
    a_->comp_grid_pair.config[a_->src_comp_idx^1].comp_name,
    a_->comp_grid_pair.config[a_->src_comp_idx  ].grid_name,
    a_->comp_grid_pair.config[a_->src_comp_idx^1].grid_name,
    a_->src_interp_config.name,
    a_->tgt_interp_config.name);

  return 0;
}

static struct yac_interp_weights * generate_interp_weights(
  struct src_field_config src_interp_config,
  struct yac_interp_grid * interp_grid) {

  struct interp_method ** method_stack =
    yac_interp_stack_config_generate(src_interp_config.interp_stack);
  struct yac_interp_weights * weights =
    yac_interp_method_do_search(method_stack, interp_grid);

  yac_interp_method_delete(method_stack);
  free(method_stack);

  return weights;
}

static void get_output_grids(
  struct yac_instance * instance, struct yac_basic_grid ** local_grids,
  size_t num_local_grids, struct output_grid ** output_grids,
  size_t * output_grid_count) {

  struct yac_couple_config * couple_config = instance->couple_config;

  // count number of output grids
  size_t num_grids = yac_couple_config_get_num_grids(couple_config);
  *output_grid_count = 0;
  for (size_t grid_idx = 0; grid_idx < num_grids; ++grid_idx)
    if (yac_couple_config_grid_get_output_filename(
          couple_config,
          yac_couple_config_get_grid_name(couple_config, grid_idx)) != NULL)
      ++*output_grid_count;

  *output_grids = xmalloc(*output_grid_count * sizeof(**output_grids));

  // extract output grids and check whether the respective grids are
  // locally available
  for (size_t grid_idx = 0, output_grid_idx = 0; grid_idx < num_grids; ++grid_idx) {
    char const * grid_name =
      yac_couple_config_get_grid_name(couple_config, grid_idx);
    char const * filename =
      yac_couple_config_grid_get_output_filename(couple_config, grid_name);
    if (filename != NULL) {
      struct yac_basic_grid * local_grid = NULL;
      for (size_t i = 0; (i < num_local_grids) && (local_grid == NULL); ++i)
        if (!strcmp(grid_name, yac_basic_grid_get_name(local_grids[i])))
          local_grid = local_grids[i];
      (*output_grids)[output_grid_idx].grid_name = grid_name;
      (*output_grids)[output_grid_idx].filename = filename;
      (*output_grids)[output_grid_idx].grid = local_grid;
      ++output_grid_idx;
    }
  }
}

static int compare_output_grids(const void * a, const void * b) {

  return
    strcmp(
      ((struct output_grid*)a)->grid_name,
      ((struct output_grid*)b)->grid_name);
}

static void write_grids_to_file(
  struct yac_instance * instance, struct yac_basic_grid ** grids, size_t num_grids) {

  MPI_Comm comm = instance->comm;

  // get information about all grids that have to be written to file
  struct output_grid * output_grids;
  size_t output_grid_count;
  get_output_grids(
    instance, grids, num_grids, &output_grids, &output_grid_count);

  // sort output grids
  qsort(
    output_grids, output_grid_count, sizeof(*output_grids),
    compare_output_grids);

  // for all grids that have to be written to file
  for (size_t i = 0; i < output_grid_count; ++i) {

    struct yac_basic_grid * grid = output_grids[i].grid;
    int split_key = (grid != NULL)?1:MPI_UNDEFINED;

    // generate a communicator containing all processes that
    // have parts of the current grid locally available
    MPI_Comm output_comm;
    yac_mpi_call(MPI_Comm_split(comm, split_key, 0, &output_comm), comm);

    // if the local process has some data of the grid locally available
    if (grid != NULL) {

      // write grid to file in parallel
      yac_basic_grid_to_file_parallel(
        grid, output_grids[i].filename, output_comm);

      yac_mpi_call(MPI_Comm_free(&output_comm), comm);
    }
  }

  free(output_grids);

  // wait until all grids have been written
  yac_mpi_call(MPI_Barrier(comm), comm);
}

static void generate_interpolations(
  struct yac_instance * instance, struct yac_basic_grid ** grids,
  size_t num_grids) {

  MPI_Comm comm = instance->comm;

  // get information about all fields
  struct field_config * field_configs;
  size_t field_count;
  get_field_configuration(
    instance, &field_configs, &field_count);

  // sort field configurations
  qsort(
    field_configs, field_count, sizeof(*field_configs), compare_field_config);

  struct yac_dist_grid_pair * dist_grid_pair = NULL;
  struct yac_interp_grid * interp_grid = NULL;
  struct yac_interp_weights * interp_weights = NULL;
  struct yac_interpolation * interp = NULL;
  struct yac_interpolation_exchange * interp_exch = NULL;
  struct yac_interp_weights_data interp_weights_data;
  struct comp_grid_pair_config * prev_comp_grid_pair = NULL;
  struct field_config * prev_field_config = NULL;

  yac_interp_weights_data_init(&interp_weights_data);

  // loop over all fields to build interpolations
  for (size_t i = 0; i < field_count; ++i) {

    struct field_config * curr_field_config = field_configs + i;
    struct comp_grid_pair_config * curr_comp_grid_pair =
      &(curr_field_config->comp_grid_pair);

    int is_source = curr_field_config->src_interp_config.field != NULL;
    int is_target = curr_field_config->tgt_interp_config.field != NULL;

    int build_flag = 0;

    // if the current configuration differs from the previous one and the local
    // process is involved in this configuration
    if ((prev_comp_grid_pair == NULL) ||
        (build_flag ||
         strcmp(prev_comp_grid_pair->config[0].grid_name,
                curr_comp_grid_pair->config[0].grid_name) ||
         strcmp(prev_comp_grid_pair->config[1].grid_name,
                curr_comp_grid_pair->config[1].grid_name))) {

      build_flag = 1;

      if (dist_grid_pair != NULL) yac_dist_grid_pair_delete(dist_grid_pair);

      char const * grid_names[2] =
        {curr_comp_grid_pair->config[0].grid_name,
         curr_comp_grid_pair->config[1].grid_name};

      int delete_flags[2];
      struct yac_basic_grid * basic_grid[2] =
        {get_basic_grid(grid_names[0], grids, num_grids, &delete_flags[0]),
         get_basic_grid(grid_names[1], grids, num_grids, &delete_flags[1])};

      dist_grid_pair =
        yac_dist_grid_pair_new(basic_grid[0], basic_grid[1], comm);

      for (int i = 0; i < 2; ++i)
        if (delete_flags[i]) yac_basic_grid_delete(basic_grid[i]);
    }

    // if the current source or target field data differes from the previous
    // one
    if (build_flag ||
        compare_field_config_interp_fields(
          prev_field_config, curr_field_config)) {

      build_flag = 1;

      struct yac_interp_field * src_fields =
        curr_field_config->src_interp_config.interp_fields;
      size_t num_src_fields =
        curr_field_config->src_interp_config.num_interp_fields;
      struct yac_interp_field * tgt_fields =
        curr_field_config->tgt_interp_config.interp_field;

      if (interp_grid != NULL) yac_interp_grid_delete(interp_grid);

      int src_comp_idx = field_configs[i].src_comp_idx;
      interp_grid = yac_interp_grid_new(
        dist_grid_pair,
        curr_comp_grid_pair->config[src_comp_idx].grid_name,
        curr_comp_grid_pair->config[src_comp_idx^1].grid_name,
        num_src_fields, src_fields, *tgt_fields);
    }

    // if the current interpolation method stack differes from the previous
    // configuration
    if (build_flag ||
        yac_interp_stack_config_compare(
          prev_field_config->src_interp_config.interp_stack,
          curr_field_config->src_interp_config.interp_stack)) {

      build_flag = 1;

      if (interp_weights != NULL) yac_interp_weights_delete(interp_weights);

      // generate interp weights
      interp_weights = generate_interp_weights(
        curr_field_config->src_interp_config, interp_grid);
    }

    if (curr_field_config->src_interp_config.weight_file_name != NULL) {

      int src_comp_idx = field_configs[i].src_comp_idx;

      yac_interp_weights_write_to_file(
        interp_weights,
        curr_field_config->src_interp_config.weight_file_name,
        curr_comp_grid_pair->config[src_comp_idx].grid_name,
        curr_comp_grid_pair->config[src_comp_idx^1].grid_name,
        0, 0, curr_field_config->src_interp_config.weight_file_on_existing);
    }

    // if the current weight reorder method differs from the previous
    // configuration
    // (use memcmp to compare frac_mask_fallback_value, because they can be nan)
    if (build_flag ||
        compare_field_config_interpolation_build_config(
          prev_field_config, curr_field_config)) {

      yac_interpolation_delete(interp);
      interp = NULL;
      yac_interpolation_exchange_delete(interp_exch, "generate_interpolations");
      interp_exch = NULL;
      yac_interp_weights_data_free(interp_weights_data);
      yac_interp_weights_data_init(&interp_weights_data);

      // generate interpolation
      if (curr_field_config->use_raw_exchange) {
        yac_interp_weights_get_interpolation_raw(
          interp_weights,
          curr_field_config->collection_size,
          curr_field_config->frac_mask_fallback_value,
          curr_field_config->scale_factor,
          curr_field_config->scale_summand,
          curr_field_config->yaxt_exchanger_name,
          &interp_exch, &interp_weights_data,
          is_source, is_target);
      } else {
        interp =
          yac_interp_weights_get_interpolation(
            interp_weights, curr_field_config->reorder_type,
            curr_field_config->collection_size,
            curr_field_config->frac_mask_fallback_value,
            curr_field_config->scale_factor,
            curr_field_config->scale_summand,
            curr_field_config->yaxt_exchanger_name,
            is_source, is_target);
      }
    }

    if (curr_field_config->use_raw_exchange) {

      struct yac_interpolation_exchange * interp_exch_copy =
        yac_interpolation_exchange_copy(interp_exch);

      if (is_source) {
        yac_set_coupling_field_put_op_raw(
          curr_field_config->src_interp_config.field,
          generate_event(
            curr_field_config->src_interp_config.event_data),
          interp_exch_copy);
        yac_interpolation_exchange_inc_ref_count(interp_exch_copy);
      }

      if (is_target) {
        struct yac_interp_weights_data interp_weights_data_copy =
          yac_interp_weights_data_copy(interp_weights_data);
        yac_set_coupling_field_get_op_raw(
          curr_field_config->tgt_interp_config.field,
          generate_event(
            curr_field_config->tgt_interp_config.event_data),
          interp_exch_copy, interp_weights_data_copy);
        yac_interpolation_exchange_inc_ref_count(interp_exch_copy);
      }

      yac_interpolation_exchange_delete(
        interp_exch_copy, "generate_interpolations");

    } else {

      struct yac_interpolation * interp_copy = yac_interpolation_copy(interp);

      if (is_source) {
        yac_set_coupling_field_put_op(
          curr_field_config->src_interp_config.field,
          generate_event(
            curr_field_config->src_interp_config.event_data),
          interp_copy);
        yac_interpolation_inc_ref_count(interp_copy);
      }

      if (is_target) {
        yac_set_coupling_field_get_op(
          curr_field_config->tgt_interp_config.field,
          generate_event(
            curr_field_config->tgt_interp_config.event_data),
          interp_copy);
        yac_interpolation_inc_ref_count(interp_copy);
      }

      yac_interpolation_delete(interp_copy);
    }

    prev_comp_grid_pair = curr_comp_grid_pair;
    prev_field_config = curr_field_config;
  }

  for (size_t i = 0; i < field_count; ++i) {
    free(field_configs[i].src_interp_config.interp_fields);
    free(field_configs[i].tgt_interp_config.interp_field);
  }

  yac_interpolation_delete(interp);
  yac_interpolation_exchange_delete(interp_exch, "generate_interpolations");
  yac_interp_weights_data_free(interp_weights_data);
  yac_interp_weights_delete(interp_weights);
  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(dist_grid_pair);
  for (size_t i = 0; i < field_count; ++i) {
    free((void*)(field_configs[i].src_interp_config.event_data.start_datetime));
    free((void*)(field_configs[i].src_interp_config.event_data.end_datetime));
    free((void*)(field_configs[i].tgt_interp_config.event_data.start_datetime));
    free((void*)(field_configs[i].tgt_interp_config.event_data.end_datetime));
  }
  free(field_configs);
}

void yac_instance_sync_def(struct yac_instance * instance) {
  CHECK_PHASE(
    yac_instance_sync_def, INSTANCE_DEFINITION_COMP, INSTANCE_DEFINITION_SYNC);
  YAC_ASSERT(instance->comp_config,
    "ERROR(yac_instance_sync_def): no components have been defined");
  yac_couple_config_sync(
    instance->couple_config, instance->comm,
    YAC_INSTANCE_CONFIG_OUTPUT_REF_SYNC);
}

void yac_instance_setup(
  struct yac_instance * instance, struct yac_basic_grid ** grids, size_t num_grids) {
  // if definitions have not yet been synced
  int requires_def_sync = (instance->phase == INSTANCE_DEFINITION_COMP);
  if (requires_def_sync)
    yac_instance_sync_def(instance);
  CHECK_PHASE(yac_instance_setup, INSTANCE_DEFINITION_SYNC, INSTANCE_EXCHANGE);

  YAC_ASSERT(
    instance->comp_config,
    "ERROR(yac_instance_setup): no components have been defined");

  // sync again, in case a process has done additional definitions
  // after the yac_instance_sync_def call
  yac_couple_config_sync(
    instance->couple_config, instance->comm,
    YAC_INSTANCE_CONFIG_OUTPUT_REF_ENDDEF);

  // write grids to file (if enabled in coupling configuration)
  write_grids_to_file(instance, grids, num_grids);

  generate_interpolations(instance, grids, num_grids);
}

char * yac_instance_setup_and_emit_config(
  struct yac_instance * instance, struct yac_basic_grid ** grids,
  size_t num_grids, int emit_flags) {

  yac_instance_setup(instance, grids, num_grids);

  int include_definitions = 0;
  return
    yac_yaml_emit_coupling(
      instance->couple_config, emit_flags, include_definitions);
}

MPI_Comm yac_instance_get_comps_comm(
  struct yac_instance * instance,
  char const ** comp_names, size_t num_comp_names) {
  CHECK_MIN_PHASE(yac_instance_get_comps_comm, INSTANCE_DEFINITION_COMP);
  return
    yac_component_config_get_comps_comm(
      instance->comp_config, comp_names, num_comp_names);
}

int yac_instance_get_nbr_comps(struct yac_instance * instance) {
  CHECK_MIN_PHASE(yac_instance_get_nbr_comps, INSTANCE_DEFINITION_COMP);
  struct yac_couple_config * couple_config =
    yac_instance_get_couple_config(instance);
  return yac_couple_config_get_num_components(couple_config);
}

int yac_instance_get_comp_size(
  struct yac_instance * instance,
  const char* comp_name){
  CHECK_MIN_PHASE(yac_instance_get_comp_size, INSTANCE_DEFINITION_COMP);
  return
    yac_component_config_comp_size(
      instance->comp_config, comp_name);
}

int yac_instance_get_comp_rank(
  struct yac_instance * instance,
  const char* comp_name){
  CHECK_MIN_PHASE(yac_instance_get_comp_rank, INSTANCE_DEFINITION_COMP);
  return
    yac_component_config_comp_rank(
      instance->comp_config, comp_name);
}

struct yac_instance * yac_instance_new(MPI_Comm comm) {

  struct yac_instance * instance = xmalloc(1 * sizeof(*instance));

  instance->couple_config = yac_couple_config_new();

  instance->comp_config = NULL;

  instance->cpl_fields = NULL;
  instance->num_cpl_fields = 0;

  yac_mpi_call(MPI_Comm_split(comm, 0, 0, &(instance->comm)), comm);

  instance->phase = INSTANCE_DEFINITION;

  return instance;
}


void yac_instance_dummy_new(MPI_Comm comm) {

  MPI_Comm dummy_comm;
  yac_mpi_call(MPI_Comm_split(comm, MPI_UNDEFINED, 0, &dummy_comm), comm);
}

void yac_instance_delete(struct yac_instance * instance) {

  if (instance == NULL) return;

  yac_component_config_delete(instance->comp_config);

  for (size_t i = 0; i < instance->num_cpl_fields; ++i)
    yac_coupling_field_delete(instance->cpl_fields[i]);
  free(instance->cpl_fields);

  yac_couple_config_delete(instance->couple_config);

  yac_mpi_call(MPI_Comm_free(&(instance->comm)), MPI_COMM_WORLD);

  free(instance);
}

MPI_Comm yac_instance_get_comm(struct yac_instance * instance) {

  if (instance == NULL) return MPI_COMM_NULL;

  return instance->comm;
}

struct yac_couple_config * yac_instance_get_couple_config(
  struct yac_instance * instance) {

  return instance->couple_config;
}

void yac_instance_set_couple_config(
  struct yac_instance * instance,
    struct yac_couple_config * couple_config) {
  CHECK_MAX_PHASE("yac_instance_set_couple_config", INSTANCE_DEFINITION);
  if (instance->couple_config == couple_config) return;
  yac_couple_config_delete(instance->couple_config);
  instance->couple_config = couple_config;
}

void yac_instance_def_datetime(
  struct yac_instance * instance, const char * start_datetime,
  const char * end_datetime ) {
  CHECK_MAX_PHASE("yac_instance_def_datetime", INSTANCE_DEFINITION_COMP);
  yac_couple_config_set_datetime(
    instance->couple_config, start_datetime, end_datetime);
}

char * yac_instance_get_start_datetime(struct yac_instance * instance) {
  CHECK_MIN_PHASE("yac_instance_get_start_datetime", INSTANCE_DEFINITION_COMP);
  return (char*)yac_couple_config_get_start_datetime(instance->couple_config);
}

char * yac_instance_get_end_datetime(struct yac_instance * instance) {
  CHECK_MIN_PHASE("yac_instance_get_end_datetime", INSTANCE_DEFINITION_COMP);
  return (char*)yac_couple_config_get_end_datetime(instance->couple_config);
}

void yac_instance_def_components(
  struct yac_instance * instance,
  char const ** comp_names, size_t num_comps) {
  CHECK_PHASE(
    yac_instance_def_components, INSTANCE_DEFINITION, INSTANCE_DEFINITION_COMP);

  YAC_ASSERT(
    !instance->comp_config,
    "ERROR(yac_instance_def_components): components have already been defined")

  // add components to coupling configuration
  for (size_t i = 0; i < num_comps; ++i)
    yac_couple_config_add_component(instance->couple_config, comp_names[i]);

  // synchronise coupling configuration
  yac_couple_config_sync(
    instance->couple_config, instance->comm,
    YAC_INSTANCE_CONFIG_OUTPUT_REF_COMP);

  instance->comp_config =
    yac_component_config_new(
      instance->couple_config, comp_names, num_comps, instance->comm);
}

int yac_instance_components_are_defined(
  struct yac_instance * instance) {
  return instance->phase > INSTANCE_DEFINITION_COMP;
}

struct coupling_field * yac_instance_add_field(
  struct yac_instance * instance, char const * field_name,
  char const * comp_name, struct yac_basic_grid * grid,
    struct yac_interp_field * interp_fields, size_t num_interp_fields,
    int collection_size, char const * timestep) {
  CHECK_MIN_PHASE(yac_instance_add_field, INSTANCE_DEFINITION_COMP);
  CHECK_MAX_PHASE(yac_instance_add_field, INSTANCE_DEFINITION_SYNC);

  struct yac_couple_config * couple_config = instance->couple_config;
  char const * grid_name = yac_basic_grid_get_name(grid);
  if(!yac_couple_config_contains_grid_name(couple_config, grid_name))
    yac_couple_config_add_grid(couple_config, grid_name);

  YAC_ASSERT(
    field_name, "ERROR(yac_instance_add_field): "
    "\"NULL\" is not a valid field name")
  YAC_ASSERT(
    strlen(field_name) <= YAC_MAX_CHARLEN,
    "ERROR(yac_instance_add_field): field name is too long "
    "(maximum is YAC_MAX_CHARLEN)")
  YAC_ASSERT_F(
    (collection_size > 0) && (collection_size < INT_MAX),
    "ERROR(yac_instance_add_field): \"%d\" is not a valid collection size "
    "(component \"%s\" grid \"%s\" field \"%s\")",
    collection_size, comp_name, grid_name, field_name)

  // add field to coupling configuration
  yac_couple_config_component_add_field(
    couple_config, comp_name, grid_name, field_name,
    timestep, collection_size);

  // check whether the field is already defined
  for (size_t i = 0; i < instance->num_cpl_fields; ++i) {
    struct coupling_field * cpl_field = instance->cpl_fields[i];
    YAC_ASSERT_F(
      strcmp(yac_get_coupling_field_name(cpl_field), field_name) ||
      (strcmp(
        yac_get_coupling_field_comp_name(cpl_field), comp_name)) ||
      (yac_coupling_field_get_basic_grid(cpl_field) != grid),
      "ERROR(yac_instance_add_field): "
      "field with the name \"%s\" has already been defined",
      field_name);
  }

  struct coupling_field * cpl_field =
    yac_coupling_field_new(
      field_name, comp_name, grid, interp_fields, num_interp_fields,
        collection_size, timestep);

  instance->cpl_fields =
    xrealloc(
      instance->cpl_fields,
      (instance->num_cpl_fields + 1) * sizeof(*(instance->cpl_fields)));
  instance->cpl_fields[instance->num_cpl_fields] = cpl_field;
  instance->num_cpl_fields++;

  return cpl_field;
}

void yac_instance_def_couple(
  struct yac_instance * instance,
  char const * src_comp_name, char const * src_grid_name, char const * src_field_name,
  char const * tgt_comp_name, char const * tgt_grid_name, char const * tgt_field_name,
  char const * coupling_period, int time_reduction,
  struct yac_interp_stack_config * interp_stack_config, int src_lag, int tgt_lag,
  const char* weight_file_name, int weight_file_on_existing,
  int mapping_on_source, double scale_factor, double scale_summand,
  size_t num_src_mask_names,
  char const * const * src_mask_names, char const * tgt_mask_name,
  char const * yaxt_exchanger_name, int use_raw_exchange) {

  CHECK_MIN_PHASE(yac_instance_def_couple, INSTANCE_DEFINITION);
  CHECK_MAX_PHASE(yac_instance_def_couple, INSTANCE_DEFINITION_SYNC);

  yac_couple_config_def_couple(
    instance->couple_config, src_comp_name, src_grid_name, src_field_name,
    tgt_comp_name, tgt_grid_name, tgt_field_name, coupling_period,
    time_reduction, interp_stack_config, src_lag, tgt_lag, weight_file_name,
    weight_file_on_existing, mapping_on_source, scale_factor, scale_summand,
    num_src_mask_names, src_mask_names, tgt_mask_name,
    yaxt_exchanger_name, use_raw_exchange);
}

struct coupling_field* yac_instance_get_field(struct yac_instance * instance,
  const char * comp_name, const char* grid_name, const char * field_name){
  CHECK_MIN_PHASE("yac_instance_get_field", INSTANCE_DEFINITION_COMP);
  return get_coupling_field(comp_name, field_name,
    grid_name, instance->num_cpl_fields, instance->cpl_fields);
}
