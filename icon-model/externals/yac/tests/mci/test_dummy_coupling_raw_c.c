// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <mpi.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "tests.h"
#include "test_common.h"
#include "yac.h"

/** \file test_dummy_coupling_raw_c.c
 *  \test
 * This test checks the raw data exchange feature.
 */

#define RAD (0.01745329251994329576923690768489) // M_PI / 180
#define FRAC_MASK_TOL (1e-12)
#define RESULT_TOL (1e-2)

struct interp_weights_data {
  double frac_mask_fallback_value;
  double scaling_factor;
  double scaling_summand;
  size_t num_fixed_values;
  double * fixed_values;
  size_t * num_tgt_per_fixed_value;
  size_t * tgt_idx_fixed;
  size_t num_wgt_tgt;
  size_t * wgt_tgt_idx;
  size_t * num_src_per_tgt;
  size_t * src_indptr;
  double * weights;
  size_t * src_field_idx;
  size_t * src_idx;
  size_t num_src_fields;
  size_t * src_field_buffer_sizes;
};

static void compute_tgt_field(
  double *** src_field_buffer, double *** src_frac_mask_buffer,
  double ** tgt_field, size_t collection_size,
  struct interp_weights_data interp_weights_data,
  size_t num_tgt_points, int use_csr_format);

static void interp_weights_data_free(
  struct interp_weights_data * interp_weights_data);

static void multi_compute_weights(
  double const tgt_coords[3], int src_cell_id, size_t src_cell_idx,
  int const ** global_results_points, double ** result_weights,
  size_t * result_count, void * user_data);

static struct interp_weights_data get_interp_weights_data(
  int field_id, int use_csr_format);

static void check_results(
  int is_tgt, int info, size_t collection_size, size_t num_tgt_points,
  double *** src_field_buffer, double *** src_frac_mask_buffer,
  struct interp_weights_data interp_weights_data, int * tgt_global_ids,
  double * ref_tgt_field_data, int use_csr_format);

/**
 * This test checks the user interface of the raw data exchange. Various
 * configurations are checked:
 * - with/without fractional masking (two different fractional masks)
 * - with/without asynchronous get
 * - process roles (only source, only target, both source and target, neither)
 * - single source field
 * - multiple source fields
 */

int main(void) {

  // initialise YAC
  yac_cinit();
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  yac_cdef_datetime("2000-01-01T00:00:00", "2000-01-02T00:00:00");

  int size, rank;
  MPI_Comm_rank ( MPI_COMM_WORLD, &rank );
  MPI_Comm_size ( MPI_COMM_WORLD, &size );

  if (size != 4) {
    PUT_ERR("ERROR: wrong number of processes");
    yac_cfinalize();
    return TEST_EXIT_CODE;
  }

  int is_src = rank <= 1;
  int is_tgt = (rank >= 1) && (rank < 3);
  int is_dummy = !is_src && !is_tgt;
  int src_rank = rank;
  int tgt_rank = rank - 1;

  // define local components
  int comp_ids[2], src_comp_id, tgt_comp_id, dummy_comp_id;
  int num_comps = 0;
  char const * src_comp_name = "source_component";
  char const * tgt_comp_name = "target_component";
  char const * dummy_comp_name = "dummy_component";
  char const * comp_names[2];
  if (is_src) comp_names[num_comps++] = src_comp_name;
  if (is_tgt) comp_names[num_comps++] = tgt_comp_name;
  if (is_dummy) comp_names[num_comps++] = dummy_comp_name;
  yac_cdef_comps(comp_names, num_comps, comp_ids);
  dummy_comp_id = is_dummy?comp_ids[--num_comps]:-1;
  tgt_comp_id = is_tgt?comp_ids[--num_comps]:-1;
  src_comp_id = is_src?comp_ids[--num_comps]:-1;

  // define local grids
  int src_grid_id, tgt_grid_id, dummy_grid_id;
  // quarter of the source grid is masked out
  int src_cell_mask_id_quarter;
  // half of the source grid is masked out (one complete process)
  int src_cell_mask_id_half_a;
  int src_cell_mask_id_half_b;
  // half of the target grid is masked out (one complete process)
  int tgt_cell_mask_id_half_a;
  int tgt_cell_mask_id_half_b;
  char const * src_grid_name = "source_grid";
  char const * tgt_grid_name = "target_grid";
  char const * dummy_grid_name = "dummy_grid";
  int src_cell_global_ids[2][2] = {{0,2}, {1,3}};
  int tgt_cell_global_ids[2][6] = {{0,1,2, 3,4,5}, {3,4,5, 6,7,8}};
  int tgt_vertex_global_ids[2][12] = {{0,1,2,3, 4,5,6,7, 8,9,10,11},
                                      {4,5,6,7, 8,9,10,11, 12,13,14,15}};
  if (is_src) {
    int nbr_vertices[2] = {2,3};
    int cyclic[2] = {0, 0};
    double x_vertices[2][2] = {{-1.0,0.0}, {0.0,1.0}};
    double y_vertices[3] = {-1.0,0.0,1.0};
    int src_cell_is_valid[3][2][2] =
      {{{0,1},{1,1}},{{1,1},{0,0}},{{0,0},{1,1}}};
    for (int i = 0; i < nbr_vertices[0]; ++i) x_vertices[src_rank][i] *= RAD;
    for (int i = 0; i < nbr_vertices[1]; ++i) y_vertices[i] *= RAD;
    yac_cdef_grid_reg2d(
      src_grid_name, nbr_vertices, cyclic,
      x_vertices[src_rank], y_vertices, &src_grid_id);
    yac_cset_global_index(
      src_cell_global_ids[src_rank], YAC_LOCATION_CELL, src_grid_id);
    yac_cdef_mask(
      src_grid_id, (nbr_vertices[0] - 1) * (nbr_vertices[1] - 1),
      YAC_LOCATION_CELL, src_cell_is_valid[0][src_rank],
      &src_cell_mask_id_quarter);
    yac_cdef_mask(
      src_grid_id, (nbr_vertices[0] - 1) * (nbr_vertices[1] - 1),
      YAC_LOCATION_CELL, src_cell_is_valid[1][src_rank],
      &src_cell_mask_id_half_a);
    yac_cdef_mask(
      src_grid_id, (nbr_vertices[0] - 1) * (nbr_vertices[1] - 1),
      YAC_LOCATION_CELL, src_cell_is_valid[2][src_rank],
      &src_cell_mask_id_half_b);
  } else {
    src_grid_id = -1;
    src_cell_mask_id_quarter = -1;
    src_cell_mask_id_half_a = -1;
    src_cell_mask_id_half_b = -1;
  }
  if (is_tgt) {
    int nbr_vertices[2] = {4,3};
    int cyclic[2] = {0, 0};
    double x_vertices[4] = {-1.5,-0.5,0.5,1.5};
    double y_vertices[2][3] = {{-1.5,-0.5,0.5}, {-0.5,0.5,1.5}};
    int tgt_cell_is_valid_a[2][6] = {{0,0,0, 0,0,0},{0,0,0, 1,1,1}};
    int tgt_cell_is_valid_b[2][6] = {{1,1,1, 0,0,0},{0,0,0, 0,0,0}};
    for (int i = 0; i < nbr_vertices[0]; ++i) x_vertices[i] *= RAD;
    for (int i = 0; i < nbr_vertices[1]; ++i) y_vertices[tgt_rank][i] *= RAD;
    yac_cdef_grid_reg2d(
      tgt_grid_name, nbr_vertices, cyclic,
      x_vertices, y_vertices[tgt_rank], &tgt_grid_id);
    yac_cset_global_index(
      tgt_cell_global_ids[tgt_rank], YAC_LOCATION_CELL, tgt_grid_id);
    yac_cset_global_index(
      tgt_vertex_global_ids[tgt_rank], YAC_LOCATION_CORNER, tgt_grid_id);
    yac_cdef_mask(
      tgt_grid_id, (nbr_vertices[0] - 1) * (nbr_vertices[1] - 1),
      YAC_LOCATION_CELL, tgt_cell_is_valid_a[tgt_rank],
      &tgt_cell_mask_id_half_a);
    yac_cdef_mask(
      tgt_grid_id, (nbr_vertices[0] - 1) * (nbr_vertices[1] - 1),
      YAC_LOCATION_CELL, tgt_cell_is_valid_b[tgt_rank],
      &tgt_cell_mask_id_half_b);
  } else {
    tgt_grid_id = -1;
    tgt_cell_mask_id_half_a = -1;
    tgt_cell_mask_id_half_b = -1;
  }
  yac_cdef_grid_reg2d(
    dummy_grid_name, (int[]){2,2}, (int[]){0,0},
    (double[]){-1,1}, (double[]){-1,1}, &dummy_grid_id);

  // define cell points
  int src_cell_point_id, tgt_cell_point_id, tgt_vertex_point_id,
      dummy_cell_point_id;
  if (is_src) {
    int num_cells[2] = {1,2};
    double x_cells[2][2] = {{-0.5}, {0.5}};
    double y_cells[2] = {-0.5,0.5};
    for (int i = 0; i < num_cells[0]; ++i) x_cells[src_rank][i] *= RAD;
    for (int i = 0; i < num_cells[1]; ++i) y_cells[i] *= RAD;
    yac_cdef_points_reg2d(
      src_grid_id, num_cells, YAC_LOCATION_CELL,
      x_cells[src_rank], y_cells, &src_cell_point_id);
  } else {
    src_cell_point_id = -1;
  }
  if (is_tgt) {
    int num_cells[2] = {3,2};
    int num_vertices[2] = {4,3};
    double x_cells[3] = {-1.0,0.0,1.0};
    double y_cells[2][2] = {{-1.0,0.0}, {0.0,1.0}};
    double x_vertices[4] = {-1.5,-0.5,0.5,1.5};
    double y_vertices[2][3] = {{-1.5,-0.5,0.5}, {-0.5,0.5,1.5}};
    for (int i = 0; i < num_cells[0]; ++i) x_cells[i] *= RAD;
    for (int i = 0; i < num_cells[1]; ++i) y_cells[tgt_rank][i] *= RAD;
    for (int i = 0; i < num_vertices[0]; ++i) x_vertices[i] *= RAD;
    for (int i = 0; i < num_vertices[1]; ++i) y_vertices[tgt_rank][i] *= RAD;
    yac_cdef_points_reg2d(
      tgt_grid_id, num_cells, YAC_LOCATION_CELL,
      x_cells, y_cells[tgt_rank], &tgt_cell_point_id);
    yac_cdef_points_reg2d(
      tgt_grid_id, num_vertices, YAC_LOCATION_CORNER,
      x_vertices, y_vertices[tgt_rank], &tgt_vertex_point_id);
  } else {
    tgt_cell_point_id = -1;
    tgt_vertex_point_id = -1;
  }
  yac_cdef_points_reg2d(
    dummy_grid_id, (int[]){1,1}, YAC_LOCATION_CELL,
    (double[]){0}, (double[]){0}, &dummy_cell_point_id);

  // additional definitions
  size_t collection_size = 1;
  if (is_src)
    yac_cadd_compute_weights_callback(
      multi_compute_weights, NULL, "multi_compute_weights");
  double fixed_value = -2.0;
  double frac_mask_fallback_value = 3.0;
  int interp_stack_conserv, interp_stack_callback;
  yac_cget_interp_stack_config(&interp_stack_conserv);
  yac_cadd_interp_stack_config_conservative(
    interp_stack_conserv, 1, 0, 1, YAC_CONSERV_DESTAREA);
  yac_cadd_interp_stack_config_fixed(interp_stack_conserv, fixed_value);
  yac_cget_interp_stack_config(&interp_stack_callback);
  yac_cadd_interp_stack_config_user_callback(
    interp_stack_callback, "multi_compute_weights");
  yac_cadd_interp_stack_config_fixed(interp_stack_callback, fixed_value);
  int ext_couple_config;
  yac_cget_ext_couple_config(&ext_couple_config);
  yac_cset_ext_couple_config_use_raw_exchange(ext_couple_config, 1);

  // all test configurations
  enum {SRC = 0, TGT = 1};
  struct {
    struct {
      char const * name;
      int * point_ids;
      int * mask_ids;
      int ** global_ids;
      int num_points;
      int id;
    } field[2];
    char const * coupling_period;
    int interp_stack_config;
    double frac_mask_fallback_value;
    int use_csr_format;
    struct interp_weights_data interp_weights_data;
    double *** src_field_buffer;
    double *** src_frac_mask_buffer;
  } test_configs[] =
  {
    { // TEST_IDX = 0
      .field =
        {{.name = "source_field",
          .point_ids = (int[]){src_cell_point_id},
          .mask_ids = (int[]){src_cell_mask_id_quarter},
          .global_ids = (int*[]){is_src?src_cell_global_ids[src_rank]:NULL},
          .num_points = 1},
         {.name = "target_field",
          .point_ids = (int[]){tgt_cell_point_id},
          .global_ids = (int*[]){is_tgt?tgt_cell_global_ids[tgt_rank]:NULL},
          .mask_ids = (int[]){-1},
          .num_points = 1}},
      .coupling_period = "1",
      .interp_stack_config = interp_stack_conserv,
      .frac_mask_fallback_value = YAC_FRAC_MASK_NO_VALUE,
      .use_csr_format = 0
    },
    { // TEST_IDX = 1
      .field =
        {{.name = "source_field_frac",
          .point_ids = (int[]){src_cell_point_id},
          .mask_ids = (int[]){src_cell_mask_id_quarter},
          .global_ids = (int*[]){is_src?src_cell_global_ids[src_rank]:NULL},
          .num_points = 1},
         {.name = "target_field_frac",
          .point_ids = (int[]){tgt_cell_point_id},
          .global_ids = (int*[]){is_tgt?tgt_cell_global_ids[tgt_rank]:NULL},
          .mask_ids = (int[]){-1},
          .num_points = 1}},
      .coupling_period = "1",
      .interp_stack_config = interp_stack_conserv,
      .frac_mask_fallback_value = frac_mask_fallback_value,
      .use_csr_format = 0
    },
    { // TEST_IDX = 2
      // test in which the process that has the role of source and target
      // does not receive any data
      // in addition this test check time reduction capabilities
      .field =
        {{.name = "source_field_full_a",
          .point_ids = (int[]){src_cell_point_id},
          .mask_ids = (int[]){-1},
          .global_ids = (int*[]){is_src?src_cell_global_ids[src_rank]:NULL},
          .num_points = 1},
         {.name = "target_field_half_a",
          .point_ids = (int[]){tgt_cell_point_id},
          .global_ids = (int*[]){is_tgt?tgt_cell_global_ids[tgt_rank]:NULL},
          .mask_ids = (int[]){tgt_cell_mask_id_half_a},
          .num_points = 1}},
      .coupling_period = "2",
      .interp_stack_config = interp_stack_conserv,
      .frac_mask_fallback_value = YAC_FRAC_MASK_NO_VALUE,
      .use_csr_format = 0
    },
    { // TEST_IDX = 3
      // test in which the process that has only the role of target
      // does not receive any data
      // in addition this test check time reduction capabilities
      .field =
        {{.name = "source_field_full_b",
          .point_ids = (int[]){src_cell_point_id},
          .mask_ids = (int[]){-1},
          .global_ids = (int*[]){is_src?src_cell_global_ids[src_rank]:NULL},
          .num_points = 1},
         {.name = "target_field_half_b",
          .point_ids = (int[]){tgt_cell_point_id},
          .global_ids = (int*[]){is_tgt?tgt_cell_global_ids[tgt_rank]:NULL},
          .mask_ids = (int[]){tgt_cell_mask_id_half_b},
          .num_points = 1}},
      .coupling_period = "2",
      .interp_stack_config = interp_stack_conserv,
      .frac_mask_fallback_value = YAC_FRAC_MASK_NO_VALUE,
      .use_csr_format = 0
    },
    { // TEST_IDX = 4
      // test in which the process that has the role of source and target
      // does not send any data
      // in addition this test check time reduction capabilities
      .field =
        {{.name = "source_field_half_a",
          .point_ids = (int[]){src_cell_point_id},
          .mask_ids = (int[]){src_cell_mask_id_half_a},
          .global_ids = (int*[]){is_src?src_cell_global_ids[src_rank]:NULL},
          .num_points = 1},
         {.name = "target_field_full_a",
          .point_ids = (int[]){tgt_cell_point_id},
          .global_ids = (int*[]){is_tgt?tgt_cell_global_ids[tgt_rank]:NULL},
          .mask_ids = (int[]){-1},
          .num_points = 1}},
      .coupling_period = "2",
      .interp_stack_config = interp_stack_conserv,
      .frac_mask_fallback_value = YAC_FRAC_MASK_NO_VALUE,
      .use_csr_format = 0
    },
    { // TEST_IDX = 5
      // test in which the process that has the only role of source
      // does not send any data
      // in addition this test check time reduction capabilities
      .field =
        {{.name = "source_field_half_b",
          .point_ids = (int[]){src_cell_point_id},
          .mask_ids = (int[]){src_cell_mask_id_half_b},
          .global_ids = (int*[]){is_src?src_cell_global_ids[src_rank]:NULL},
          .num_points = 1},
         {.name = "target_field_full_b",
          .point_ids = (int[]){tgt_cell_point_id},
          .global_ids = (int*[]){is_tgt?tgt_cell_global_ids[tgt_rank]:NULL},
          .mask_ids = (int[]){-1},
          .num_points = 1}},
      .coupling_period = "2",
      .interp_stack_config = interp_stack_conserv,
      .frac_mask_fallback_value = YAC_FRAC_MASK_NO_VALUE,
      .use_csr_format = 0
    },
    { // TEST_IDX = 6
      // test in which the process that has the role of source and target
      // does not receive any data (with fractional masking)
      // in addition this test check time reduction capabilities
      .field =
        {{.name = "source_field_frac_full_a",
          .point_ids = (int[]){src_cell_point_id},
          .mask_ids = (int[]){-1},
          .global_ids = (int*[]){is_src?src_cell_global_ids[src_rank]:NULL},
          .num_points = 1},
         {.name = "target_field_frac_half_a",
          .point_ids = (int[]){tgt_cell_point_id},
          .global_ids = (int*[]){is_tgt?tgt_cell_global_ids[tgt_rank]:NULL},
          .mask_ids = (int[]){tgt_cell_mask_id_half_a},
          .num_points = 1}},
      .coupling_period = "2",
      .interp_stack_config = interp_stack_conserv,
      .frac_mask_fallback_value = frac_mask_fallback_value,
      .use_csr_format = 0
    },
    { // TEST_IDX = 7
      // test in which the process that has only the role of target
      // does not receive any data (with fractional masking)
      // in addition this test check time reduction capabilities
      .field =
        {{.name = "source_field_frac_full_b",
          .point_ids = (int[]){src_cell_point_id},
          .mask_ids = (int[]){-1},
          .global_ids = (int*[]){is_src?src_cell_global_ids[src_rank]:NULL},
          .num_points = 1},
         {.name = "target_field_frac_half_b",
          .point_ids = (int[]){tgt_cell_point_id},
          .global_ids = (int*[]){is_tgt?tgt_cell_global_ids[tgt_rank]:NULL},
          .mask_ids = (int[]){tgt_cell_mask_id_half_b},
          .num_points = 1}},
      .coupling_period = "2",
      .interp_stack_config = interp_stack_conserv,
      .frac_mask_fallback_value = frac_mask_fallback_value,
      .use_csr_format = 0
    },
    { // TEST_IDX = 8
      // test in which the process that has the role of source and target
      // does not send any data (with fractional masking)
      // in addition this test check time reduction capabilities
      .field =
        {{.name = "source_field_frac_half_a",
          .point_ids = (int[]){src_cell_point_id},
          .mask_ids = (int[]){src_cell_mask_id_half_a},
          .global_ids = (int*[]){is_src?src_cell_global_ids[src_rank]:NULL},
          .num_points = 1},
         {.name = "target_field_frac_full_a",
          .point_ids = (int[]){tgt_cell_point_id},
          .global_ids = (int*[]){is_tgt?tgt_cell_global_ids[tgt_rank]:NULL},
          .mask_ids = (int[]){-1},
          .num_points = 1}},
      .coupling_period = "2",
      .interp_stack_config = interp_stack_conserv,
      .frac_mask_fallback_value = frac_mask_fallback_value,
      .use_csr_format = 0
    },
    { // TEST_IDX = 9
      // test in which the process that has the only role of source
      // does not send any data (with fractional masking)
      // in addition this test check time reduction capabilities
      .field =
        {{.name = "source_field_frac_half_b",
          .point_ids = (int[]){src_cell_point_id},
          .mask_ids = (int[]){src_cell_mask_id_half_b},
          .global_ids = (int*[]){is_src?src_cell_global_ids[src_rank]:NULL},
          .num_points = 1},
         {.name = "target_field_frac_full_b",
          .point_ids = (int[]){tgt_cell_point_id},
          .global_ids = (int*[]){is_tgt?tgt_cell_global_ids[tgt_rank]:NULL},
          .mask_ids = (int[]){-1},
          .num_points = 1}},
      .coupling_period = "2",
      .interp_stack_config = interp_stack_conserv,
      .frac_mask_fallback_value = frac_mask_fallback_value,
      .use_csr_format = 0
    },
    { // TEST_IDX = 10
      // test in the target field is generated from multiple source point
      // sets
      .field =
        {{.name = "source_field_multi",
          .point_ids = (int[]){src_cell_point_id, src_cell_point_id},
          .mask_ids = (int[]){-1,-1},
          .global_ids = (int*[]){is_src?src_cell_global_ids[src_rank]:NULL},
          .num_points = 2},
         {.name = "target_field_multi",
          .point_ids = (int[]){tgt_vertex_point_id},
          .global_ids = (int*[]){is_tgt?tgt_vertex_global_ids[tgt_rank]:NULL},
          .mask_ids = (int[]){-1},
          .num_points = 1}},
      .coupling_period = "1",
      .interp_stack_config = interp_stack_callback,
      .frac_mask_fallback_value = YAC_FRAC_MASK_NO_VALUE,
      .use_csr_format = 0
    },
    { // TEST_IDX = 11
      // test in which a source field is sent to two target fields
      // (one using raw data exchange and the other does not)
      .field =
        {{.name = "source_field_comp",
          .point_ids = (int[]){src_cell_point_id},
          .mask_ids = (int[]){src_cell_mask_id_quarter},
          .global_ids = (int*[]){is_src?src_cell_global_ids[src_rank]:NULL},
          .num_points = 1},
         {.name = "target_field_comp_raw",
          .point_ids = (int[]){tgt_cell_point_id},
          .global_ids = (int*[]){is_tgt?tgt_cell_global_ids[tgt_rank]:NULL},
          .mask_ids = (int[]){-1},
          .num_points = 1}},
      .coupling_period = "1",
      .interp_stack_config = interp_stack_conserv,
      .frac_mask_fallback_value = YAC_FRAC_MASK_NO_VALUE,
      .use_csr_format = 0
    },
    { // TEST_IDX = 12
      .field =
        {{.name = "source_field_csr",
          .point_ids = (int[]){src_cell_point_id},
          .mask_ids = (int[]){src_cell_mask_id_quarter},
          .global_ids = (int*[]){is_src?src_cell_global_ids[src_rank]:NULL},
          .num_points = 1},
         {.name = "target_field_csr",
          .point_ids = (int[]){tgt_cell_point_id},
          .global_ids = (int*[]){is_tgt?tgt_cell_global_ids[tgt_rank]:NULL},
          .mask_ids = (int[]){-1},
          .num_points = 1}},
      .coupling_period = "1",
      .interp_stack_config = interp_stack_conserv,
      .frac_mask_fallback_value = YAC_FRAC_MASK_NO_VALUE,
      .use_csr_format = 1
    }
  };

  // define configurations
  for (size_t i = 0; i < sizeof(test_configs) / sizeof(test_configs[0]); ++i) {

    // define source and target fields
    for (int j = 0; j < 2; ++j) {
      if (test_configs[i].field[j].point_ids[0] != -1) {
        int comp_id = (j == SRC)?src_comp_id:tgt_comp_id;
        if (test_configs[i].field[j].mask_ids[0] != -1)
          yac_cdef_field_mask(
            test_configs[i].field[j].name, comp_id,
            test_configs[i].field[j].point_ids,
            test_configs[i].field[j].mask_ids,
            test_configs[i].field[j].num_points,
            1, "1", YAC_TIME_UNIT_SECOND,
            &test_configs[i].field[j].id);
        else
          yac_cdef_field(
            test_configs[i].field[j].name, comp_id,
            test_configs[i].field[j].point_ids,
            test_configs[i].field[j].num_points,
            1, "1", YAC_TIME_UNIT_SECOND,
            &test_configs[i].field[j].id);
      } else {
        test_configs[i].field[j].id = -1;
      }
    }

    // enable fractional masking, if required
    if (is_src && (test_configs[i].frac_mask_fallback_value !=
                   YAC_FRAC_MASK_NO_VALUE))
      yac_cenable_field_frac_mask(
        src_comp_name, src_grid_name, test_configs[i].field[SRC].name,
        test_configs[i].frac_mask_fallback_value);

    // define couple
    yac_cdef_couple_custom(
      src_comp_name, src_grid_name, test_configs[i].field[SRC].name,
      tgt_comp_name, tgt_grid_name, test_configs[i].field[TGT].name,
      test_configs[i].coupling_period,
      YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_AVERAGE,
      test_configs[i].interp_stack_config, 0, 0, ext_couple_config);
  }

  int tgt_field_non_raw_id;
  if (is_tgt) {
    // define a file that will not use raw data exchange
    // it will receive data from a source field, which also sends data
    // to a target that uses raw data exchange
    yac_cdef_field(
      "target_field_comb_non_raw", tgt_comp_id, &tgt_cell_point_id,
      1, 1, "1", YAC_TIME_UNIT_SECOND, &tgt_field_non_raw_id);
    yac_cdef_couple(
      src_comp_name, src_grid_name, "source_field_comp",
      tgt_comp_name, tgt_grid_name, "target_field_comb_non_raw",
      "1", YAC_TIME_UNIT_SECOND, YAC_REDUCTION_TIME_AVERAGE,
      interp_stack_conserv, 0, 0);
  } else {
    tgt_field_non_raw_id = -1;
  }

  // uncoupled dummy fields
  int src_dummy_field_id = -1, tgt_dummy_field_id = -1;
  switch (1 * is_src + 2 * is_tgt) {
    case (0): { // neither source nor target
      yac_cdef_field(
        "target_field_dummy", dummy_comp_id, &dummy_cell_point_id,
        1, 1, "1", YAC_TIME_UNIT_SECOND, &tgt_dummy_field_id);
      yac_cdef_field(
        "source_field_dummy", dummy_comp_id, &dummy_cell_point_id,
        1, 1, "1", YAC_TIME_UNIT_SECOND, &src_dummy_field_id);
      break;
    }
    case (1): { // only source
      yac_cdef_field(
        "target_field_dummy", src_comp_id, &dummy_cell_point_id,
        1, 1, "1", YAC_TIME_UNIT_SECOND, &tgt_dummy_field_id);
      break;
    }
    case (2): { // only target
      yac_cdef_field(
        "source_field_dummy", tgt_comp_id, &dummy_cell_point_id,
        1, 1, "1", YAC_TIME_UNIT_SECOND, &src_dummy_field_id);
      break;
    }
    case (3): { // source and target
      break;
    }
    default: {
      PUT_ERR("ERROR: invalid process type");
      yac_cfinalize();
      return TEST_EXIT_CODE;
    }
  } // role

  // end definition phase
  yac_cenddef();

  // extract interpolation weight data and allocate source field buffer
  for (size_t i = 0; i < sizeof(test_configs) / sizeof(test_configs[0]); ++i) {
    test_configs[i].interp_weights_data =
      get_interp_weights_data(
        test_configs[i].field[TGT].id, test_configs[i].use_csr_format);
    size_t num_src_fields = test_configs[i].interp_weights_data.num_src_fields;
    size_t * src_field_buffer_sizes =
      test_configs[i].interp_weights_data.src_field_buffer_sizes;
    size_t sum_src_field_buffer_sizes = 0;
    for (size_t i = 0; i < num_src_fields; ++i)
      sum_src_field_buffer_sizes += src_field_buffer_sizes[i];
    double * src_field_buffer_1d =
      (sum_src_field_buffer_sizes > 0)?
        malloc(
          collection_size * sum_src_field_buffer_sizes *
          sizeof(*src_field_buffer_1d)):NULL;
    double *** src_field_buffer =
      malloc(collection_size * sizeof(*src_field_buffer));
    for (size_t i = 0, offset = 0; i < collection_size; ++i) {
      src_field_buffer[i] = malloc(num_src_fields * sizeof(**src_field_buffer));
      for (size_t j = 0; j < num_src_fields; ++j) {
        src_field_buffer[i][j] = src_field_buffer_1d + offset;
        offset += src_field_buffer_sizes[j];
      }
    }
    test_configs[i].src_field_buffer = src_field_buffer;
    if (test_configs[i].frac_mask_fallback_value != YAC_FRAC_MASK_NO_VALUE) {
      double * src_frac_mask_buffer_1d =
      (sum_src_field_buffer_sizes > 0)?
          malloc(
            collection_size * sum_src_field_buffer_sizes *
            sizeof(*src_frac_mask_buffer_1d)):NULL;
      double *** src_frac_mask_buffer =
        malloc(collection_size * sizeof(*src_frac_mask_buffer));
      for (size_t i = 0, offset = 0; i < collection_size; ++i) {
        src_frac_mask_buffer[i] =
          malloc(num_src_fields * sizeof(**src_frac_mask_buffer));
        for (size_t j = 0; j < num_src_fields; ++j) {
          src_frac_mask_buffer[i][j] = src_frac_mask_buffer_1d + offset;
          offset += src_field_buffer_sizes[j];
        }
      }
      test_configs[i].src_frac_mask_buffer = src_frac_mask_buffer;
    } else {
      test_configs[i].src_frac_mask_buffer = NULL;
    }
  }

  //-----------------------------------------------------
  // do tests
  //-----------------------------------------------------

  enum {
    MAX_COLLECTION_SIZE = 1,
    MAX_NUM_SRC_FIELDS = 2,
    GLOBAL_NUM_SRC_CELLS = 4,
    GLOBAL_NUM_TGT_CELLS = 9,
    GLOBAL_NUM_TGT_VERTICES = 16,
    NUM_SRC_CELLS = 2,
    NUM_TGT_CELLS = 6,
    NUM_TGT_VERTICES = 12,
  };
  double global_src_field_data[MAX_NUM_SRC_FIELDS][GLOBAL_NUM_SRC_CELLS] =
    {{1,2,3,4}, {-1,-2,-3,-4}};
  double src_field_data[MAX_COLLECTION_SIZE][MAX_NUM_SRC_FIELDS][NUM_SRC_CELLS];
  double * src_field_[MAX_COLLECTION_SIZE][MAX_NUM_SRC_FIELDS];
  double ** src_field[MAX_COLLECTION_SIZE];
  if (is_src) {
    for (size_t i = 0; i < MAX_COLLECTION_SIZE; ++i) {
      for (size_t j = 0; j < MAX_NUM_SRC_FIELDS; ++j) {
        for (size_t k = 0; k < NUM_SRC_CELLS; ++k)
          src_field_data[i][j][k] =
            global_src_field_data[j][src_cell_global_ids[src_rank][k]];
        src_field_[i][j] = src_field_data[i][j];
      }
      src_field[i] = src_field_[i];
    }
  }

  int info, send_info, recv_info, ierror;

  { // some basic tests

    enum {TEST_IDX = 0, NUM_TESTS_BASIC = 3, COLLECTION_SIZE = 1};
    int src_field_id = test_configs[TEST_IDX].field[SRC].id;
    int tgt_field_id = test_configs[TEST_IDX].field[TGT].id;
    double *** src_field_buffer = test_configs[TEST_IDX].src_field_buffer;
    double ref_tgt_field_data[1][GLOBAL_NUM_TGT_CELLS] =
      {{fixed_value,0.25*(2),0.25*(2),
        0.25*(3),0.25*(2+3+4),0.25*(2+4),
        0.25*(3),0.25*(3+4),0.25*(4)}};

    for (size_t t = 0; t < NUM_TESTS_BASIC; ++t) {

      switch (1 * is_src + 2 * is_tgt) {
        case (0): { // neither source nor target
          break;
        }
        case (1): { // only source

          // send source data
          yac_cput(
            src_field_id, COLLECTION_SIZE, src_field, &info, &ierror);
          yac_cwait(src_field_id);
          break;
        }
        case (2): { // only target

          // receive source field buffer
          switch (t) {
            default:
            case(0): // blocking get
              yac_cget_raw(
                tgt_field_id, COLLECTION_SIZE, src_field_buffer,
                &info, &ierror);
              break;
            case(1): // async get using source field buffer pointers
              yac_cget_raw_async(
                tgt_field_id, COLLECTION_SIZE, src_field_buffer,
                &info, &ierror);
              yac_cwait(tgt_field_id);
              break;
            case(2): // async get using 1d source field buffer
              yac_cget_raw_async_(
                tgt_field_id, COLLECTION_SIZE, src_field_buffer[0][0],
                &info, &ierror);
              yac_cwait(tgt_field_id);
              break;
          }
          break;
        }
        case (3): { // source and target

          yac_cexchange_raw(
            src_field_id, tgt_field_id, COLLECTION_SIZE,
            src_field, src_field_buffer, &send_info, &recv_info, &ierror);
          info = recv_info;
          break;
        }
        default: {
          PUT_ERR("ERROR: invalid process type");
          yac_cfinalize();
          return TEST_EXIT_CODE;
        }
      } // role

      check_results(
        is_tgt, info, COLLECTION_SIZE, NUM_TGT_CELLS,
        src_field_buffer, NULL, test_configs[TEST_IDX].interp_weights_data,
        test_configs[TEST_IDX].field[TGT].global_ids[0],
        ref_tgt_field_data[0], test_configs[TEST_IDX].use_csr_format);
    } // basic test idx
  }

  { // some basic tests using fractional masking

    enum {TEST_IDX = 1, NUM_FRAC_MASKS = 3,
          NUM_FRAC_TESTS = 4, COLLECTION_SIZE = 1, NUM_SRC_FIELDS = 1};
    int src_field_id = test_configs[TEST_IDX].field[SRC].id;
    int tgt_field_id = test_configs[TEST_IDX].field[TGT].id;
    double *** src_field_buffer = test_configs[TEST_IDX].src_field_buffer;
    double *** src_frac_mask_buffer =
      test_configs[TEST_IDX].src_frac_mask_buffer;
    double ref_tgt_field_data[2][GLOBAL_NUM_TGT_CELLS] =
      {{fixed_value,0.25*(2),0.25*(2),
        0.25*(3),0.25*(2+3+4),0.25*(2+4),
        0.25*(3),0.25*(3+4),0.25*(4)},
        {fixed_value, frac_mask_fallback_value, frac_mask_fallback_value,
        (0.25*(1.0*3))/((0.25*1.0)/(0.25)),
        (0.25*(0.0*2+0.5*3+1.0*4))/((0.25*(0.0+0.5+1.0))/(0.25+0.25+0.25)),
        (0.25*(0.0*2+1.0*4))/((0.25*(0.0+1.0))/(0.25+0.25)),
        (0.25*(0.5*3))/((0.25*0.5)/(0.25)),
        (0.25*(0.5*3+1.0*4))/((0.25*(0.5+1.0))/(0.25+0.25)),
        (0.25*(1.0*4))/((0.25*1.0)/(0.25))}};

    double global_src_frac_mask_data
      [NUM_FRAC_MASKS][NUM_SRC_FIELDS][GLOBAL_NUM_SRC_CELLS] =
        {{{1,1,1,1}},{{0.5,0.0,0.5,1.0}}};
    double src_frac_mask_data
      [NUM_FRAC_MASKS][MAX_COLLECTION_SIZE][NUM_SRC_FIELDS][NUM_SRC_CELLS];
    double * src_frac_mask_
      [NUM_FRAC_MASKS][MAX_COLLECTION_SIZE][NUM_SRC_FIELDS];
    double ** src_frac_mask[NUM_FRAC_MASKS][MAX_COLLECTION_SIZE];
    if (is_src) {
      for (size_t t = 0; t < NUM_FRAC_MASKS; ++t) {
        for (size_t i = 0; i < MAX_COLLECTION_SIZE; ++i) {
          for (size_t j = 0; j < NUM_SRC_FIELDS; ++j) {
            for (size_t k = 0; k < NUM_SRC_CELLS; ++k)
              src_frac_mask_data[t][i][j][k] =
                global_src_frac_mask_data
                  [t][j][src_cell_global_ids[src_rank][k]];
            src_frac_mask_[t][i][j] = src_frac_mask_data[t][i][j];
          }
          src_frac_mask[t][i] = src_frac_mask_[t][i];
        }
      }
    }

    for (size_t t = 0; t < NUM_FRAC_TESTS; ++t) {

      switch (1 * is_src + 2 * is_tgt) {
        case (0): { // neither source nor target
          break;
        }
        case (1): { // only source

          // send source data and fractional mask
          yac_cput_frac(
            src_field_id, COLLECTION_SIZE,
            src_field, src_frac_mask[t&1], &info, &ierror);
          yac_cwait(src_field_id);
          break;
        }
        case (2): { // only target

          // receive source field and source fractional mask
          switch (t) {
            default:
            case(0):
              yac_cget_raw_frac_async(
                tgt_field_id, COLLECTION_SIZE, src_field_buffer,
                src_frac_mask_buffer, &info, &ierror);
              yac_cwait(tgt_field_id);
              break;
            case(1):
              yac_cget_raw_frac(
                tgt_field_id, COLLECTION_SIZE, src_field_buffer,
                src_frac_mask_buffer, &info, &ierror);
              break;
            case(2):
              yac_cget_raw_frac_async_(
                tgt_field_id, COLLECTION_SIZE, src_field_buffer[0][0],
                src_frac_mask_buffer[0][0], &info, &ierror);
              yac_cwait(tgt_field_id);
              break;
            case(3):
              yac_cget_raw_frac_(
                tgt_field_id, COLLECTION_SIZE, src_field_buffer[0][0],
                src_frac_mask_buffer[0][0], &info, &ierror);
              break;
          }
          break;
        }
        case (3): { // source and target

          // exchange source field data and source fractional mask
          yac_cexchange_raw_frac(
            src_field_id, tgt_field_id, COLLECTION_SIZE,
            src_field, src_frac_mask[t&1],
            src_field_buffer, src_frac_mask_buffer,
            &send_info, &recv_info, &ierror);
          info = recv_info;
          break;
        }
        default: {
          PUT_ERR("ERROR: invalid process type");
          yac_cfinalize();
          return TEST_EXIT_CODE;
        }
      } // role

      check_results(
        is_tgt, info, COLLECTION_SIZE, NUM_TGT_CELLS,
        src_field_buffer, src_frac_mask_buffer,
        test_configs[TEST_IDX].interp_weights_data,
        test_configs[TEST_IDX].field[TGT].global_ids[0],
        ref_tgt_field_data[t&1], test_configs[TEST_IDX].use_csr_format);
    } // frac test idx
  }

  { // a couple of tests in which one process does not send/receive any data
    // in addition this test check time reduction capabilities

    enum {MIN_TEST_IDX = 2, MAX_TEST_IDX = 9,
          COLLECTION_SIZE = 1, NUM_TIMESTEPS = 6, NUM_SRC_FIELDS = 1};
    double ref_tgt_field_data[8][1][GLOBAL_NUM_TGT_CELLS] =
      {{{-1.0,-1.0,-1.0,
          -1.0,-1.0,-1.0,
          0.25*(3),0.25*(3+4),0.25*(4)}},
        {{0.25*(1),0.25*(1+2),0.25*(2),
          -1.0,-1.0,-1.0,
          -1.0,-1.0,-1.0}},
        {{0.25*(1),0.25*(1),fixed_value,
          0.25*(1+3),0.25*(1+3),fixed_value,
          0.25*(3),0.25*(3),fixed_value}},
        {{fixed_value,0.25*(2),0.25*(2),
          fixed_value,0.25*(2+4),0.25*(2+4),
          fixed_value,0.25*(4),0.25*(4)}},
        {{-1.0,-1.0,-1.0,
          -1.0,-1.0,-1.0,
          0.25*(3),0.25*(3+4),0.25*(4)}},
        {{0.25*(1),0.25*(1+2),0.25*(2),
          -1.0,-1.0,-1.0,
          -1.0,-1.0,-1.0}},
        {{0.25*(1),0.25*(1),fixed_value,
          0.25*(1+3),0.25*(1+3),fixed_value,
          0.25*(3),0.25*(3),fixed_value}},
        {{fixed_value,0.25*(2),0.25*(2),
          fixed_value,0.25*(2+4),0.25*(2+4),
          fixed_value,0.25*(4),0.25*(4)}}};
    double global_src_frac_mask_data
      [NUM_SRC_FIELDS][GLOBAL_NUM_SRC_CELLS] = {{1,1,1,1}};
    double src_frac_mask_data
      [MAX_COLLECTION_SIZE][NUM_SRC_FIELDS][NUM_SRC_CELLS];
    double * src_frac_mask_ [MAX_COLLECTION_SIZE][NUM_SRC_FIELDS];
    double ** src_frac_mask[MAX_COLLECTION_SIZE];
    if (is_src) {
      for (size_t i = 0; i < MAX_COLLECTION_SIZE; ++i) {
        for (size_t j = 0; j < NUM_SRC_FIELDS; ++j) {
          for (size_t k = 0; k < NUM_SRC_CELLS; ++k)
            src_frac_mask_data[i][j][k] =
              global_src_frac_mask_data[j][src_cell_global_ids[src_rank][k]];
          src_frac_mask_[i][j] = src_frac_mask_data[i][j];
        }
        src_frac_mask[i] = src_frac_mask_[i];
      }
    }

    for (int use_only_exchange = 0; use_only_exchange <= 1;
         ++use_only_exchange) {
      for (size_t test_idx = MIN_TEST_IDX; test_idx < MAX_TEST_IDX; ++test_idx) {

        int src_field_id =
          is_src?test_configs[test_idx].field[SRC].id:src_dummy_field_id;
        int tgt_field_id =
          is_tgt?test_configs[test_idx].field[TGT].id:tgt_dummy_field_id;
        double *** src_field_buffer = test_configs[test_idx].src_field_buffer;
        double *** src_frac_mask_buffer =
          test_configs[test_idx].src_frac_mask_buffer;
        int with_frac_mask =
          test_configs[test_idx].frac_mask_fallback_value !=
          YAC_FRAC_MASK_NO_VALUE;

        for (size_t t = 0; t < NUM_TIMESTEPS; ++t) {

          info = YAC_ACTION_NONE;

          switch (1 * is_src + 2 * is_tgt) {
            case (0): { // neither source nor target

              if (use_only_exchange) {
                if (with_frac_mask)
                  yac_cexchange_raw_frac(
                    src_field_id, tgt_field_id, COLLECTION_SIZE,
                    NULL, NULL, NULL, NULL,
                    &send_info, &recv_info, &ierror);
                else
                  yac_cexchange_raw(
                    src_field_id, tgt_field_id, COLLECTION_SIZE,
                    NULL, NULL, &send_info, &recv_info, &ierror);
                info = recv_info;
              }
              break;
            }
            case (1): { // only source

              // send source data
              if (use_only_exchange) {
                if (with_frac_mask)
                  yac_cexchange_raw_frac(
                    src_field_id, tgt_field_id, COLLECTION_SIZE,
                    src_field, src_frac_mask, NULL, NULL,
                    &send_info, &recv_info, &ierror);
                else
                  yac_cexchange_raw(
                    src_field_id, tgt_field_id, COLLECTION_SIZE,
                    src_field, NULL, &send_info, &recv_info, &ierror);
                info = recv_info;
              } else {
                // send source data
                if (with_frac_mask)
                  yac_cput_frac(
                    src_field_id, COLLECTION_SIZE, src_field, src_frac_mask,
                    &info, &ierror);
                else
                  yac_cput(
                    src_field_id, COLLECTION_SIZE, src_field, &info, &ierror);
              }
              yac_cwait(src_field_id);
              break;
            }
            case (2): { // only target

              // receive source field buffer
              if (use_only_exchange) {
                if (with_frac_mask)
                  yac_cexchange_raw_frac(
                    src_field_id, tgt_field_id, COLLECTION_SIZE,
                    NULL, NULL, src_field_buffer, src_frac_mask_buffer,
                    &send_info, &recv_info, &ierror);
                else
                  yac_cexchange_raw(
                    src_field_id, tgt_field_id, COLLECTION_SIZE,
                    NULL, src_field_buffer, &send_info, &recv_info, &ierror);
              } else {
                if (with_frac_mask)
                  yac_cget_raw_frac(
                    tgt_field_id, COLLECTION_SIZE,
                    src_field_buffer, src_frac_mask_buffer, &info, &ierror);
                else
                  yac_cget_raw(
                    tgt_field_id, COLLECTION_SIZE, src_field_buffer, &info, &ierror);
              }
              break;
            }
            case (3): { // source and target

              if (with_frac_mask)
                yac_cexchange_raw_frac(
                  src_field_id, tgt_field_id, COLLECTION_SIZE,
                  src_field, src_frac_mask,
                  src_field_buffer, src_frac_mask_buffer,
                  &send_info, &recv_info, &ierror);
              else
                yac_cexchange_raw(
                  src_field_id, tgt_field_id, COLLECTION_SIZE,
                  src_field, src_field_buffer, &send_info, &recv_info, &ierror);
              info = recv_info;
              break;
            }
            default: {
              PUT_ERR("ERROR: invalid process type");
              yac_cfinalize();
              return TEST_EXIT_CODE;
            }
          } // role

          check_results(
            is_tgt, info, COLLECTION_SIZE, NUM_TGT_CELLS,
            src_field_buffer, with_frac_mask?src_frac_mask_buffer:NULL,
            test_configs[test_idx].interp_weights_data,
            test_configs[test_idx].field[TGT].global_ids[0],
            ref_tgt_field_data[test_idx-MIN_TEST_IDX][0],
            test_configs[test_idx].use_csr_format);
        } // time step
      } // test_idx
    } // use_only_exchange
  }

  { // test using multiple source fields

    enum {TEST_IDX = 10, COLLECTION_SIZE = 1};
    int src_field_id = test_configs[TEST_IDX].field[SRC].id;
    int tgt_field_id = test_configs[TEST_IDX].field[TGT].id;
    double *** src_field_buffer = test_configs[TEST_IDX].src_field_buffer;
    double ref_tgt_field_data[1][GLOBAL_NUM_TGT_VERTICES] =
      {{fixed_value, fixed_value, fixed_value, fixed_value,
        fixed_value, 1, 2, fixed_value,
        fixed_value, -3, -4, fixed_value,
        fixed_value, fixed_value, fixed_value, fixed_value}};

    switch (1 * is_src + 2 * is_tgt) {
      case (0): { // neither source nor target
        break;
      }
      case (1): { // only source

        // send source data
        yac_cput(
          src_field_id, COLLECTION_SIZE, src_field, &info, &ierror);
        yac_cwait(src_field_id);
        break;
      }
      case (2): { // only target

        // receive source field buffer
        yac_cget_raw(
          tgt_field_id, COLLECTION_SIZE, src_field_buffer,
          &info, &ierror);
        break;
      }
      case (3): { // source and target

        yac_cexchange_raw(
          src_field_id, tgt_field_id, COLLECTION_SIZE,
          src_field, src_field_buffer, &send_info, &recv_info, &ierror);
        info = recv_info;
        break;
      }
      default: {
        PUT_ERR("ERROR: invalid process type");
        yac_cfinalize();
        return TEST_EXIT_CODE;
      }
    } // role

    check_results(
      is_tgt, info, COLLECTION_SIZE, NUM_TGT_VERTICES,
      src_field_buffer, NULL,
      test_configs[TEST_IDX].interp_weights_data,
      test_configs[TEST_IDX].field[TGT].global_ids[0],
      ref_tgt_field_data[0], test_configs[TEST_IDX].use_csr_format);
  }

  { // test with multiple targets

    enum {TEST_IDX = 11, COLLECTION_SIZE = 1};
    int src_field_id = test_configs[TEST_IDX].field[SRC].id;
    int tgt_field_id = test_configs[TEST_IDX].field[TGT].id;
    double *** src_field_buffer = test_configs[TEST_IDX].src_field_buffer;
    double ref_tgt_field_data[1][GLOBAL_NUM_TGT_CELLS] =
      {{fixed_value,0.25*(2),0.25*(2),
        0.25*(3),0.25*(2+3+4),0.25*(2+4),
        0.25*(3),0.25*(3+4),0.25*(4)}};

    switch (1 * is_src + 2 * is_tgt) {
      case (0): { // neither source nor target
        break;
      }
      case (1): { // only source

        // send source data
        yac_cput(
          src_field_id, COLLECTION_SIZE, src_field, &info, &ierror);
        yac_cwait(src_field_id);
        break;
      }
      case (2): { // only target

        // receive source field buffer
        yac_cget_raw(
          tgt_field_id, COLLECTION_SIZE, src_field_buffer,
          &info, &ierror);
        break;
      }
      case (3): { // source and target

        yac_cput(
          src_field_id, COLLECTION_SIZE, src_field, &send_info, &ierror);
        yac_cget_raw(
          tgt_field_id, COLLECTION_SIZE, src_field_buffer,
          &recv_info, &ierror);
        info = recv_info;
        break;
      }
      default: {
        PUT_ERR("ERROR: invalid process type");
        yac_cfinalize();
        return TEST_EXIT_CODE;
      }
    } // role

    check_results(
      is_tgt, info, COLLECTION_SIZE, NUM_TGT_CELLS,
      src_field_buffer, NULL,
      test_configs[TEST_IDX].interp_weights_data,
      test_configs[TEST_IDX].field[TGT].global_ids[0],
      ref_tgt_field_data[0], test_configs[TEST_IDX].use_csr_format);

    if (is_tgt) {

      double tgt_field_data[COLLECTION_SIZE][NUM_TGT_CELLS];
      double * tgt_field[COLLECTION_SIZE];
      for (size_t i = 0; i < COLLECTION_SIZE; ++i)
        tgt_field[i] = tgt_field_data[i];

      // initialise target field
      for (size_t i = 0; i < COLLECTION_SIZE; ++i)
        for (size_t j = 0; j < NUM_TGT_CELLS; ++j)
          tgt_field[i][j] = -1.0;

      yac_cget(
        tgt_field_non_raw_id, COLLECTION_SIZE, tgt_field, &info, &ierror);

      // if data was received
      if ((info == YAC_ACTION_COUPLING) ||
          (info == YAC_ACTION_GET_FOR_RESTART)) {

        // check results
        int * tgt_global_ids =
          test_configs[TEST_IDX].field[TGT].global_ids[0];
        for (size_t i = 0; i < COLLECTION_SIZE; ++i)
          for (size_t j = 0; j < NUM_TGT_CELLS; ++j)
            if (fabs(
                  tgt_field_data[i][j] -
                  ref_tgt_field_data[0][tgt_global_ids[j]]) >
                RESULT_TOL)
              PUT_ERR("wrong results");
      }
    }
  }

  { // testing of interpolation weights data in csr format

    enum {TEST_IDX = 12, NUM_TESTS_BASIC = 3, COLLECTION_SIZE = 1};
    int src_field_id = test_configs[TEST_IDX].field[SRC].id;
    int tgt_field_id = test_configs[TEST_IDX].field[TGT].id;
    double *** src_field_buffer = test_configs[TEST_IDX].src_field_buffer;
    double ref_tgt_field_data[1][GLOBAL_NUM_TGT_CELLS] =
      {{fixed_value,0.25*(2),0.25*(2),
        0.25*(3),0.25*(2+3+4),0.25*(2+4),
        0.25*(3),0.25*(3+4),0.25*(4)}};

    for (size_t t = 0; t < NUM_TESTS_BASIC; ++t) {

      switch (1 * is_src + 2 * is_tgt) {
        case (0): { // neither source nor target
          break;
        }
        case (1): { // only source

          // send source data
          yac_cput(
            src_field_id, COLLECTION_SIZE, src_field, &info, &ierror);
          yac_cwait(src_field_id);
          break;
        }
        case (2): { // only target

          // receive source field buffer
          switch (t) {
            default:
            case(0): // blocking get
              yac_cget_raw(
                tgt_field_id, COLLECTION_SIZE, src_field_buffer,
                &info, &ierror);
              break;
            case(1): // async get using source field buffer pointers
              yac_cget_raw_async(
                tgt_field_id, COLLECTION_SIZE, src_field_buffer,
                &info, &ierror);
              yac_cwait(tgt_field_id);
              break;
            case(2): // async get using 1d source field buffer
              yac_cget_raw_async_(
                tgt_field_id, COLLECTION_SIZE, src_field_buffer[0][0],
                &info, &ierror);
              yac_cwait(tgt_field_id);
              break;
          }
          break;
        }
        case (3): { // source and target

          yac_cexchange_raw(
            src_field_id, tgt_field_id, COLLECTION_SIZE,
            src_field, src_field_buffer, &send_info, &recv_info, &ierror);
          info = recv_info;
          break;
        }
        default: {
          PUT_ERR("ERROR: invalid process type");
          yac_cfinalize();
          return TEST_EXIT_CODE;
        }
      } // role

      check_results(
        is_tgt, info, COLLECTION_SIZE, NUM_TGT_CELLS,
        src_field_buffer, NULL, test_configs[TEST_IDX].interp_weights_data,
        test_configs[TEST_IDX].field[TGT].global_ids[0],
        ref_tgt_field_data[0], test_configs[TEST_IDX].use_csr_format);
    } // basic test idx
  }

  yac_cfree_ext_couple_config(ext_couple_config);
  yac_cfree_interp_stack_config(interp_stack_conserv);
  yac_cfree_interp_stack_config(interp_stack_callback);
  for (size_t i = 0; i < sizeof(test_configs) / sizeof(test_configs[0]); ++i) {
    if (test_configs[i].interp_weights_data.num_src_fields > 0)
      free(test_configs[i].src_field_buffer[0][0]);
    for (size_t j = 0; j < collection_size; ++j)
      free(test_configs[i].src_field_buffer[j]);
    free(test_configs[i].src_field_buffer);
    if (test_configs[i].frac_mask_fallback_value != YAC_FRAC_MASK_NO_VALUE) {
      if (test_configs[i].interp_weights_data.num_src_fields > 0)
        free(test_configs[i].src_frac_mask_buffer[0][0]);
      for (size_t j = 0; j < collection_size; ++j)
        free(test_configs[i].src_frac_mask_buffer[j]);
      free(test_configs[i].src_frac_mask_buffer);
    }
    if (is_tgt) interp_weights_data_free(&test_configs[i].interp_weights_data);
  }

  yac_cfinalize();
  return TEST_EXIT_CODE;
}

static void compute_tgt_field(
  double *** src_field_buffer, double *** src_frac_mask_buffer,
  double ** tgt_field, size_t collection_size,
  struct interp_weights_data interp_weights_data,
  size_t num_tgt_points, int use_csr_format) {

  // ignore weights below a certain
  double const wgt_tol = 1e-6;

  // we have to use memcpy to compare against YAC_FRAC_MASK_NO_VALUE,
  // because nan is a valid value for frac_mask_fallback_value and depending
  // on compiler optimisation a comparision agains nan can produce varying
  // results
  int with_frac_mask =
    memcmp(&interp_weights_data.frac_mask_fallback_value,
           &YAC_FRAC_MASK_NO_VALUE, sizeof(YAC_FRAC_MASK_NO_VALUE));

  for (size_t collection_idx = 0; collection_idx < collection_size;
       ++collection_idx) {

    // set fixed targets
    for (size_t i = 0, k = 0; i < interp_weights_data.num_fixed_values; ++i) {

      double fixed_value = interp_weights_data.fixed_values[i];
      size_t num_fixed_tgt = interp_weights_data.num_tgt_per_fixed_value[i];

      for (size_t j = 0; j < num_fixed_tgt; ++j, ++k)
        tgt_field[collection_idx][interp_weights_data.tgt_idx_fixed[k]] =
          fixed_value;
    }

    size_t num_wgt_tgt =
      use_csr_format?num_tgt_points:interp_weights_data.num_wgt_tgt;

    if (with_frac_mask) {

      // set weighted targets
      for (size_t i = 0, k = 0; i < num_wgt_tgt; ++i) {

        double tgt_value = 0.0;
        double frac_weight_sum = 0.0;
        double weight_sum = 0.0;
        size_t num_src =
          use_csr_format?
            (interp_weights_data.src_indptr[i+1] -
             interp_weights_data.src_indptr[i]):
            interp_weights_data.num_src_per_tgt[i];
        size_t tgt_idx =
          use_csr_format?i:interp_weights_data.wgt_tgt_idx[i];

        if (num_src == 0) continue;

        for (size_t j = 0; j < num_src; ++j, ++k) {
          if (fabs(interp_weights_data.weights[k]) < wgt_tol) continue;
          tgt_value +=
            interp_weights_data.weights[k] *
            src_field_buffer
              [collection_idx]
              [interp_weights_data.src_field_idx[k]]
              [interp_weights_data.src_idx[k]];
          frac_weight_sum +=
            interp_weights_data.weights[k] *
            src_frac_mask_buffer
              [collection_idx]
              [interp_weights_data.src_field_idx[k]]
              [interp_weights_data.src_idx[k]];
          weight_sum += interp_weights_data.weights[k];
        }

        tgt_field[collection_idx][tgt_idx] =
          (fabs(frac_weight_sum) > FRAC_MASK_TOL)?
          interp_weights_data.scaling_factor *
          (tgt_value / frac_weight_sum) * weight_sum +
          interp_weights_data.scaling_summand:
          interp_weights_data.frac_mask_fallback_value;
      }

    } else {

      // set weighted targets
      for (size_t i = 0, k = 0; i < num_wgt_tgt; ++i) {

        double tgt_value = 0.0;
        size_t num_src =
          use_csr_format?
            (interp_weights_data.src_indptr[i+1] -
             interp_weights_data.src_indptr[i]):
            interp_weights_data.num_src_per_tgt[i];
        size_t tgt_idx =
          use_csr_format?i:interp_weights_data.wgt_tgt_idx[i];

        if (num_src == 0) continue;

        for (size_t j = 0; j < num_src; ++j, ++k) {
          if (fabs(interp_weights_data.weights[k]) < wgt_tol) continue;
          tgt_value +=
            interp_weights_data.weights[k] *
            src_field_buffer
              [collection_idx]
              [interp_weights_data.src_field_idx[k]]
              [interp_weights_data.src_idx[k]];
        }

        tgt_field[collection_idx][tgt_idx] =
          interp_weights_data.scaling_factor * tgt_value +
          interp_weights_data.scaling_summand;
      }
    }
  }
}

static void interp_weights_data_free(
  struct interp_weights_data * interp_weights_data) {

  free(interp_weights_data->fixed_values);
  free(interp_weights_data->num_tgt_per_fixed_value);
  free(interp_weights_data->tgt_idx_fixed);
  free(interp_weights_data->wgt_tgt_idx);
  free(interp_weights_data->num_src_per_tgt);
  free(interp_weights_data->src_indptr);
  free(interp_weights_data->weights);
  free(interp_weights_data->src_field_idx);
  free(interp_weights_data->src_idx);
  free(interp_weights_data->src_field_buffer_sizes);
}

static void multi_compute_weights(
  double const tgt_coords[3], int src_cell_id, size_t src_cell_idx,
  int const ** global_results_points, double ** result_weights,
  size_t * result_count, void * user_data) {

  // unused arguments
  (void)(tgt_coords);
  (void)(src_cell_idx);
  (void)(user_data);

  // only the cells in the upper row of the source grid receive their
  // data from the second source field
  int use_second_src_field = src_cell_id >= 2;

  static int results_point;
  static double result_weight;

  results_point = src_cell_id;
  result_weight = 1.0;

  global_results_points[use_second_src_field] = &results_point;
  global_results_points[use_second_src_field^1] = NULL;
  result_weights[use_second_src_field] = &result_weight;
  result_weights[use_second_src_field^1] = NULL;
  result_count[use_second_src_field] = 1;
  result_count[use_second_src_field^1] = 0;
}

static struct interp_weights_data get_interp_weights_data(
  int field_id, int use_csr_format) {

  struct interp_weights_data data;

  if (field_id > -1) {

    if (use_csr_format) {
      yac_cget_raw_interp_weights_data_csr(
        field_id,
        &data.frac_mask_fallback_value,
        &data.scaling_factor, &data.scaling_summand,
        &data.num_fixed_values, &data.fixed_values,
        &data.num_tgt_per_fixed_value, &data.tgt_idx_fixed,
        &data.src_indptr, &data.weights, &data.src_field_idx, &data.src_idx,
        &data.num_src_fields, &data.src_field_buffer_sizes);
      data.num_wgt_tgt = 0;
      data.wgt_tgt_idx = NULL;
      data.num_src_per_tgt = NULL;
    } else {
      yac_cget_raw_interp_weights_data(
        field_id,
        &data.frac_mask_fallback_value,
        &data.scaling_factor, &data.scaling_summand,
        &data.num_fixed_values, &data.fixed_values,
        &data.num_tgt_per_fixed_value, &data.tgt_idx_fixed,
        &data.num_wgt_tgt, &data.wgt_tgt_idx, &data.num_src_per_tgt,
        &data.weights, &data.src_field_idx, &data.src_idx,
        &data.num_src_fields, &data.src_field_buffer_sizes);
      data.src_indptr = NULL;
    }
  } else {
    data.num_src_fields = 0;
  }

  return data;
}

static void check_results(
  int is_tgt, int info, size_t collection_size, size_t num_tgt_points,
  double *** src_field_buffer, double *** src_frac_mask_buffer,
  struct interp_weights_data interp_weights_data, int * tgt_global_ids,
  double * ref_tgt_field_data, int use_csr_format) {

  // if data was received
  if (is_tgt &&
      ((info == YAC_ACTION_COUPLING) ||
        (info == YAC_ACTION_GET_FOR_RESTART))) {

    double tgt_field_data[collection_size][num_tgt_points];
    double * tgt_field[collection_size];
    for (size_t i = 0; i < collection_size; ++i)
      tgt_field[i] = tgt_field_data[i];

    // initialise target field
    for (size_t i = 0; i < collection_size; ++i)
      for (size_t j = 0; j < num_tgt_points; ++j)
        tgt_field[i][j] = -1.0;

    // compute target field
    compute_tgt_field(
      src_field_buffer, src_frac_mask_buffer,
      (double**)tgt_field, collection_size, interp_weights_data,
      num_tgt_points, use_csr_format);

    // check results
    for (size_t i = 0; i < collection_size; ++i)
      for (size_t j = 0; j < num_tgt_points; ++j)
        if (fabs(
              tgt_field_data[i][j] - ref_tgt_field_data[tgt_global_ids[j]]) >
            RESULT_TOL)
          PUT_ERR("wrong results");
  }
}
