// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <yaxt.h>
#include <netcdf.h>
#include <string.h>

#include "tests.h"
#include "test_common.h"
#include "geometry.h"
#include "read_icon_grid.h"
#include "interp_method.h"
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
#include "weight_file_common.h"
#include "yac_mpi.h"

/** \file test_interp_method_parallel.c
 *  \test
 * A test for general functional of "abstract" interp_method type.
 */

static void utest_compute_weights_callback(
  double const tgt_coords[3], int src_cell_id, size_t src_cell_idx,
  int const ** global_results_points, double ** result_weights,
  size_t * result_count, void * user_data);

int main(int argc, char** argv) {

  MPI_Init(NULL, NULL);

  xt_initialize(MPI_COMM_WORLD);

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  set_even_io_rank_list(MPI_COMM_WORLD);

  if (argc != 2) {
    PUT_ERR("ERROR: missing grid file directory");
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  char * filenames[2];
  char * grid_filenames[] ={"icon_grid_0030_R02B03_G.nc", "icon_grid_0043_R02B04_G.nc"};
  for (int i = 0; i < 2; ++i)
    filenames[i] =
      strcat(
        strcpy(
          malloc(strlen(argv[1]) + strlen(grid_filenames[i]) + 2), argv[1]),
        grid_filenames[i]);

  struct yac_basic_grid_data grid_data[2];

  for (int i = 0; i < 2; ++i)
    grid_data[i] =
      yac_read_icon_basic_grid_data_parallel(
        filenames[i], MPI_COMM_WORLD);

  struct yac_basic_grid * grids[2] =
    {yac_basic_grid_new(filenames[0], grid_data[0]),
     yac_basic_grid_new(filenames[1], grid_data[1])};

  struct yac_dist_grid_pair * grid_pair =
    yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

  struct yac_interp_field src_fields[] =
    {{.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
  size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
  struct yac_interp_field tgt_field =
    {.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

  struct yac_interp_grid * interp_grid =
    yac_interp_grid_new(grid_pair, filenames[0], filenames[1],
                        num_src_fields, src_fields, tgt_field);

  struct interp_method * method_stack[] =
    {yac_interp_method_fixed_new(-1.0), NULL};

  struct yac_interp_weights * weights =
    yac_interp_method_do_search(method_stack, interp_grid);

  struct yac_interpolation * interpolation_src =
    yac_interp_weights_get_interpolation(
      weights, YAC_MAPPING_ON_SRC, 1,
      YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0, NULL, 1, 1);
  struct yac_interpolation * interpolation_tgt =
    yac_interp_weights_get_interpolation(
      weights, YAC_MAPPING_ON_TGT, 1,
      YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0, NULL, 1, 1);

  yac_interpolation_delete(interpolation_src);
  yac_interpolation_delete(interpolation_tgt);

  yac_interp_weights_delete(weights);

  yac_interp_method_delete(method_stack);

  yac_interp_grid_delete(interp_grid);
  yac_dist_grid_pair_delete(grid_pair);

  for (int i = 0; i < 2; ++i) {
    yac_basic_grid_delete(grids[i]);
    free(filenames[i]);
  }

 { // test interpolation_complete argument

    char const * weight_file_empty =
      "test_interp_method_parallel_weights_empty.nc";
    char const * weight_file_fixed =
      "test_interp_method_parallel_weights_fixed.nc";

    char const * src_grid_name = "src_grid";
    char const * tgt_grid_name = "tgt_grid";

    // create weight file (one empty and one with a fixed value)
    if (comm_rank == 0) {

      { // empty file
        int * tgt_indices = NULL;
        int * src_indices = NULL;
        double * weights = NULL;
        size_t num_links = 0;
        enum yac_location src_locations[] = {YAC_LOC_CELL};
        enum yac_location tgt_location = YAC_LOC_CELL;
        enum {
          NUM_SRC_FIELDS = sizeof(src_locations) / sizeof(src_locations[0])};
        int num_links_per_field[NUM_SRC_FIELDS] = {num_links};
        int * tgt_id_fixed = NULL;
        size_t num_fixed_tgt = 0;
        double * fixed_values = NULL;
        int * num_tgt_per_fixed_value = NULL;
        size_t num_fixed_values = 0;

        write_weight_file(
          weight_file_empty, src_indices, tgt_indices, weights, num_links,
          src_locations, NUM_SRC_FIELDS, num_links_per_field, tgt_id_fixed,
          num_fixed_tgt, fixed_values, num_tgt_per_fixed_value,
          num_fixed_values, tgt_location, src_grid_name, tgt_grid_name);
      }

      { // file with fixed value
        int * tgt_indices = NULL;
        int * src_indices = NULL;
        double * weights = NULL;
        size_t num_links = 0;
        enum yac_location src_locations[] = {YAC_LOC_CELL};
        enum yac_location tgt_location = YAC_LOC_CELL;
        enum {
          NUM_SRC_FIELDS = sizeof(src_locations) / sizeof(src_locations[0])};
        int num_links_per_field[NUM_SRC_FIELDS] = {num_links};
        int tgt_id_fixed[] = {0, 1, 2, 3};
        size_t num_fixed_tgt = 4;
        double fixed_values[] = {999.0};
        int num_tgt_per_fixed_value[] = {4};
        size_t num_fixed_values = 1;

        write_weight_file(
          weight_file_fixed, src_indices, tgt_indices, weights, num_links,
          src_locations, NUM_SRC_FIELDS, num_links_per_field, tgt_id_fixed,
          num_fixed_tgt, fixed_values, num_tgt_per_fixed_value,
          num_fixed_values, tgt_location, src_grid_name, tgt_grid_name);
      }
    }
    double coordinates_x[] = {0.0, 1.0, 2.0};
    double coordinates_y[] = {0.0, 1.0, 2.0};
    size_t num_vertices[2] = {3,3};
    int cyclic[2] = {0,0};

    struct yac_basic_grid * src_grid =
      yac_basic_grid_reg_2d_deg_new(
        src_grid_name, num_vertices, cyclic, coordinates_x, coordinates_y);
    struct yac_basic_grid * tgt_grid =
      yac_basic_grid_reg_2d_deg_new(
        tgt_grid_name, num_vertices, cyclic, coordinates_x, coordinates_y);

    double cell_coords[4][3];
    double cell_coords_x[] = {0.5, 1.5};
    double cell_coords_y[] = {0.5, 1.5};
    for (int i = 0, k = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j, ++k)
        LLtoXYZ_deg(cell_coords_x[j], cell_coords_y[i], cell_coords[k]);
    yac_basic_grid_add_coordinates(src_grid, YAC_LOC_CELL, cell_coords, 4);
    yac_basic_grid_add_coordinates(tgt_grid, YAC_LOC_CELL, cell_coords, 4);

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(src_grid, tgt_grid, MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CELL, .coordinates_idx = 0, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CELL, .coordinates_idx = 0, .masks_idx = SIZE_MAX};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, src_grid_name, tgt_grid_name,
                          num_src_fields, src_fields, tgt_field);

    //----------------------------------------
    // test generation of interpolation method
    //----------------------------------------

    struct interp_method * method_stack[] = {
      yac_interp_method_file_new( // this method stops the stack
        weight_file_empty, YAC_INTERP_FILE_MISSING_ERROR,
        YAC_INTERP_FILE_SUCCESS_STOP),
      yac_interp_method_avg_new(
        (enum yac_interp_avg_weight_type)
          YAC_INTERP_AVG_WEIGHT_TYPE_DEFAULT,
        YAC_INTERP_AVG_PARTIAL_COVERAGE_DEFAULT),
      yac_interp_method_callback_new(utest_compute_weights_callback, NULL),
      yac_interp_method_check_new( NULL, NULL, NULL, NULL),
      yac_interp_method_conserv_new(
        YAC_INTERP_CONSERV_ORDER_DEFAULT,
        YAC_INTERP_CONSERV_ENFORCED_CONSERV_DEFAULT,
        YAC_INTERP_CONSERV_PARTIAL_COVERAGE_DEFAULT,
        (enum yac_interp_method_conserv_normalisation)
          YAC_INTERP_CONSERV_NORMALISATION_DEFAULT),
      yac_interp_method_creep_new(YAC_INTERP_CREEP_DISTANCE_DEFAULT),
      yac_interp_method_file_new(
        weight_file_fixed, YAC_INTERP_FILE_MISSING_ERROR,
        YAC_INTERP_FILE_SUCCESS_CONT),
      yac_interp_method_fixed_new(YAC_INTERP_FIXED_VALUE_DEFAULT),
      yac_interp_method_hcsbb_new(),
      yac_interp_method_ncc_new(
        (enum yac_interp_ncc_weight_type)
          YAC_INTERP_NCC_WEIGHT_TYPE_DEFAULT,
        YAC_INTERP_NCC_PARTIAL_COVERAGE_DEFAULT),
      yac_interp_method_nnn_new(
        (struct yac_nnn_config) {
          .type =
            (enum yac_interp_nnn_weight_type)YAC_INTERP_NNN_WEIGHTED_DEFAULT,
          .n = YAC_INTERP_NNN_N_DEFAULT,
          .max_search_distance = YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT,
          .data.gauss_scale = YAC_INTERP_NNN_GAUSS_SCALE_DEFAULT}),
      yac_interp_method_spmap_new(
        YAC_INTERP_SPMAP_DEFAULT_CONFIG, YAC_INTERP_SPMAP_OVERWRITE_DEFAULT),
      NULL};

    //-----------------
    // generate weights
    //-----------------

    struct yac_interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    if (yac_interp_weights_get_interp_count(weights) != 0)
      PUT_ERR(
        "error in handling of interpolation_complete argument by do_search");

    //---------------
    // cleanup
    //---------------

    yac_interp_method_delete(method_stack);
    yac_interp_weights_delete(weights);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(tgt_grid);
    yac_basic_grid_delete(src_grid);

    if (comm_rank == 0) {
      unlink(weight_file_empty);
      unlink(weight_file_fixed);
    }
  }

  xt_finalize();

  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static void utest_compute_weights_callback(
  double const tgt_coords[3], int src_cell_id, size_t src_cell_idx,
  int const ** global_results_points, double ** result_weights,
  size_t * result_count, void * user_data) {

  UNUSED(tgt_coords);
  UNUSED(src_cell_id);
  UNUSED(src_cell_idx);
  UNUSED(global_results_points);
  UNUSED(result_weights);
  UNUSED(result_count);
  UNUSED(user_data);

  PUT_ERR("ERROR: in handling interpolation_complete in do_search_callback");
}
