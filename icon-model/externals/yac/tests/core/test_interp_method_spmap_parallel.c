// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <unistd.h>

#include "tests.h"
#include "test_common.h"
#include "geometry.h"
#include "interp_method.h"
#include "interp_method_spmap.h"
#include "interp_method_fixed.h"
#include "dist_grid_utils.h"
#include "yac_mpi.h"
#include "io_utils.h"
#include "point_selection.h"

#include <mpi.h>
#include <yaxt.h>
#include <netcdf.h>

/** \file test_interp_method_spmap_parallel.c
 *  \test
 * A test for the parallel source point mapping interpolation method.
 */

enum yac_interp_spmap_weight_type weight_types[] = {
  YAC_INTERP_SPMAP_AVG,
  YAC_INTERP_SPMAP_DIST
};
size_t num_weight_types = sizeof(weight_types) / sizeof(weight_types[0]);

enum yac_interp_spmap_scale_type scale_types[] = {
  YAC_INTERP_SPMAP_NONE,
  YAC_INTERP_SPMAP_SRCAREA,
  YAC_INTERP_SPMAP_INVTGTAREA,
  YAC_INTERP_SPMAP_FRACAREA
};
size_t num_scale_types = sizeof(scale_types) / sizeof(scale_types[0]);

enum yac_interp_weights_reorder_type reorder_types[] = {
  YAC_MAPPING_ON_SRC,
  YAC_MAPPING_ON_TGT
};
size_t num_reorder_types = sizeof(reorder_types) / sizeof(reorder_types[0]);

static char const * grid_names[2] = {"grid1", "grid2"};

static double utest_get_sq_sphere_radius(
  struct yac_spmap_cell_area_config const * cell_area_confi);
static double utest_compute_scale(
  enum yac_interp_spmap_scale_type scale_type,
  double src_cell_area, double tgt_cell_area);
static void utest_write_area_file(
  char const * filename, double * cell_areas, size_t dim_x, size_t dim_y);

int main(void) {

  MPI_Init(NULL, NULL);

  xt_initialize(MPI_COMM_WORLD);

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  if (comm_size != 6) {
    PUT_ERR("ERROR: wrong number of processes");
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  MPI_Comm split_comm;
  MPI_Comm_split(MPI_COMM_WORLD, comm_rank < 1, 0, &split_comm);

  int split_comm_rank, split_comm_size;
  MPI_Comm_rank(split_comm, &split_comm_rank);
  MPI_Comm_size(split_comm, &split_comm_size);

  { // parallel interpolation process
    // corner and cell ids for a 7 x 7 grid
    // 56-----57-----58-----59-----60-----61-----62-----63
    //  |      |      |      |      |      |      |      |
    //  |  42  |  43  |  44  |  45  |  46  |  47  |  48  |
    //  |      |      |      |      |      |      |      |
    // 48-----49-----50-----51-----52-----53-----54-----55
    //  |      |      |      |      |      |      |      |
    //  |  35  |  36  |  37  |  38  |  39  |  40  |  41  |
    //  |      |      |      |      |      |      |      |
    // 40-----41-----42-----43-----44-----45-----46-----47
    //  |      |      |      |      |      |      |      |
    //  |  28  |  29  |  30  |  31  |  32  |  33  |  34  |
    //  |      |      |      |      |      |      |      |
    // 32-----33-----34-----35-----36-----37-----38-----39
    //  |      |      |      |      |      |      |      |
    //  |  21  |  22  |  23  |  24  |  25  |  26  |  27  |
    //  |      |      |      |      |      |      |      |
    // 24-----25-----26-----27-----28-----29-----30-----31
    //  |      |      |      |      |      |      |      |
    //  |  14  |  15  |  16  |  17  |  18  |  19  |  20  |
    //  |      |      |      |      |      |      |      |
    // 16-----17-----18-----19-----20-----21-----22-----23
    //  |      |      |      |      |      |      |      |
    //  |  07  |  08  |  09  |  10  |  11  |  12  |  13  |
    //  |      |      |      |      |      |      |      |
    // 08-----09-----10-----11-----12-----13-----14-----15
    //  |      |      |      |      |      |      |      |
    //  |  00  |  01  |  02  |  03  |  04  |  05  |  06  |
    //  |      |      |      |      |      |      |      |
    // 00-----01-----02-----03-----04-----05-----06-----07
    //
    // the grid is distributed among the processes as follows:
    // (index == process)
    //
    // 4---4---4---4---4---4---4---4
    // | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
    // 4---4---4---4---4---4---4---4
    // | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
    // 3---3---3---3---4---4---4---4
    // | 2 | 2 | 2 | 2 | 3 | 3 | 3 |
    // 2---2---2---2---3---3---3---3
    // | 2 | 2 | 2 | 2 | 3 | 3 | 3 |
    // 1---1---1---1---2---2---2---2
    // | 1 | 1 | 1 | 1 | 3 | 3 | 3 |
    // 0---0---0---0---1---1---1---1
    // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
    // 0---0---0---0---0---0---0---0
    // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
    // 0---0---0---0---0---0---0---0
    //
    // mask
    // land = 0
    // coast = 1
    // ocean = 2
    // +---+---+---+---+---+---+---+
    // | 1 | 2 | 2 | 2 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 1 | 2 | 2 | 1 | 1 | 1 | 2 |
    // +---+---+---+---+---+---+---+
    // | 1 | 2 | 2 | 1 | 0 | 1 | 2 |
    // +---+---+---+---+---+---+---+
    // | 1 | 2 | 2 | 1 | 1 | 1 | 2 |
    // +---+---+---+---+---+---+---+
    // | 1 | 1 | 2 | 2 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 0 | 1 | 1 | 2 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 0 | 0 | 1 | 1 | 1 | 1 | 1 |
    // +---+---+---+---+---+---+---+
    //
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[8] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    double coordinates_y[8] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    size_t const num_cells[2] = {7,7};
    size_t local_start[2][5][2] = {{{0,0},{0,2},{0,3},{4,2},{0,5}}, {{0,0}}};
    size_t local_count[2][5][2] = {{{7,2},{4,1},{4,2},{3,3},{7,2}}, {{7,7}}};
    int global_mask[7*7] = {
      0,0,1,1,1,1,1,
      0,1,1,2,2,2,2,
      1,1,2,2,2,2,2,
      1,2,2,1,1,1,2,
      1,2,2,1,0,1,2,
      1,2,2,1,1,1,2,
      1,2,2,2,2,2,2};
    int with_halo = 0;
    for (size_t i = 0; i <= num_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct yac_basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_cells,
        local_start[is_tgt][split_comm_rank],
        local_count[is_tgt][split_comm_rank], with_halo);

    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[is_tgt], grid_data),
       yac_basic_grid_empty_new(grid_names[is_tgt^1])};

    {
      int valid_mask_value = (is_tgt)?2:1;
      yac_coordinate_pointer point_coordinates =
        xmalloc(grid_data.num_cells * sizeof(*point_coordinates));
      int * mask = xmalloc(grid_data.num_cells * sizeof(*mask));
      for (size_t i = 0; i < grid_data.num_cells; ++i) {
        double * middle_point = point_coordinates[i];
        for (size_t k = 0; k < 3; ++k) middle_point[k] = 0.0;
        size_t * curr_vertices =
          grid_data.cell_to_vertex + grid_data.cell_to_vertex_offsets[i];
        size_t curr_num_vertices = grid_data.num_vertices_per_cell[i];
        for (size_t j = 0; j < curr_num_vertices; ++j) {
          double * curr_vertex_coord =
            grid_data.vertex_coordinates[curr_vertices[j]];
          for (size_t k = 0; k < 3; ++k)
            middle_point[k] += curr_vertex_coord[k];
        }
        normalise_vector(middle_point);
        mask[i] = global_mask[grid_data.cell_ids[i]] == valid_mask_value;
      }
      yac_basic_grid_add_coordinates_nocpy(
        grids[0], YAC_LOC_CELL, point_coordinates);
      yac_basic_grid_add_mask_nocpy(
        grids[0], YAC_LOC_CELL, mask, NULL);
    }

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CELL, .coordinates_idx = 0, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CELL, .coordinates_idx = 0, .masks_idx = 0};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                          num_src_fields, src_fields, tgt_field);

    struct interp_method * method_stack[] =
      {yac_interp_method_spmap_new(
         YAC_INTERP_SPMAP_DEFAULT_CONFIG, YAC_INTERP_SPMAP_OVERWRITE_DEFAULT),
       yac_interp_method_fixed_new(-2.0), NULL};

    struct yac_interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    for (size_t i = 0; i < num_reorder_types; ++i) {

      struct yac_interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_types[i], 1,
          YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0, NULL, 1, 1);

      // check generated interpolation
      {
        double * src_field = NULL;
        double ** src_fields = &src_field;
        double * tgt_field = NULL;

        if (is_tgt) {
          tgt_field = xmalloc(grid_data.num_cells * sizeof(*tgt_field));
          for (size_t i = 0; i < grid_data.num_cells; ++i) tgt_field[i] = -1;
        } else {
          src_field = xmalloc(grid_data.num_cells * sizeof(*src_field));
          for (size_t i = 0; i < grid_data.num_cells; ++i)
            src_field[i] = (double)(grid_data.cell_ids[i]);
        }

        yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

        double ref_tgt_field[7*7] = {
          -1,-1,-1,-1,-1,-1,-1,
          -1,-1,-1,2+3+9,4,5,6,
          -1,-1,8+15,-2,25,-2,-2,
          -1,14+21,24,-1,-1,-1,26,
          -1,28,31,-1,-1,-1,33,
          -1,35,38,-1,-1,-1,40,
          -1,42,-2,-2,39,-2,-2};

        if (is_tgt)
          for (size_t i = 0; i < grid_data.num_cells; ++i)
            if (ref_tgt_field[grid_data.cell_ids[i]] != tgt_field[i])
              PUT_ERR("wrong interpolation result");

        free(tgt_field);
        free(src_field);
      }

      yac_interpolation_delete(interpolation);
    }

    yac_interp_weights_delete(weights);
    yac_interp_method_delete(method_stack);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(grids[1]);
    yac_basic_grid_delete(grids[0]);
  }

  { // parallel interpolation process
    // corner and cell ids for a 7 x 7 grid
    // 56-----57-----58-----59-----60-----61-----62-----63
    //  |      |      |      |      |      |      |      |
    //  |  42  |  43  |  44  |  45  |  46  |  47  |  48  |
    //  |      |      |      |      |      |      |      |
    // 48-----49-----50-----51-----52-----53-----54-----55
    //  |      |      |      |      |      |      |      |
    //  |  35  |  36  |  37  |  38  |  39  |  40  |  41  |
    //  |      |      |      |      |      |      |      |
    // 40-----41-----42-----43-----44-----45-----46-----47
    //  |      |      |      |      |      |      |      |
    //  |  28  |  29  |  30  |  31  |  32  |  33  |  34  |
    //  |      |      |      |      |      |      |      |
    // 32-----33-----34-----35-----36-----37-----38-----39
    //  |      |      |      |      |      |      |      |
    //  |  21  |  22  |  23  |  24  |  25  |  26  |  27  |
    //  |      |      |      |      |      |      |      |
    // 24-----25-----26-----27-----28-----29-----30-----31
    //  |      |      |      |      |      |      |      |
    //  |  14  |  15  |  16  |  17  |  18  |  19  |  20  |
    //  |      |      |      |      |      |      |      |
    // 16-----17-----18-----19-----20-----21-----22-----23
    //  |      |      |      |      |      |      |      |
    //  |  07  |  08  |  09  |  10  |  11  |  12  |  13  |
    //  |      |      |      |      |      |      |      |
    // 08-----09-----10-----11-----12-----13-----14-----15
    //  |      |      |      |      |      |      |      |
    //  |  00  |  01  |  02  |  03  |  04  |  05  |  06  |
    //  |      |      |      |      |      |      |      |
    // 00-----01-----02-----03-----04-----05-----06-----07
    //
    // the grid is distributed among the processes as follows:
    // (index == process)
    //
    // 4---4---4---4---4---4---4---4
    // | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
    // 4---4---4---4---4---4---4---4
    // | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
    // 3---3---3---3---4---4---4---4
    // | 2 | 2 | 2 | 2 | 3 | 3 | 3 |
    // 2---2---2---2---3---3---3---3
    // | 2 | 2 | 2 | 2 | 3 | 3 | 3 |
    // 1---1---1---1---2---2---2---2
    // | 1 | 1 | 1 | 1 | 3 | 3 | 3 |
    // 0---0---0---0---1---1---1---1
    // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
    // 0---0---0---0---0---0---0---0
    // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
    // 0---0---0---0---0---0---0---0
    //
    // mask
    // land = 0
    // coast = 1
    // ocean = 2
    // +---+---+---+---+---+---+---+
    // | 1 | 2 | 2 | 2 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 1 | 2 | 2 | 1 | 1 | 1 | 2 |
    // +---+---+---+---+---+---+---+
    // | 1 | 2 | 2 | 1 | 0 | 1 | 2 |
    // +---+---+---+---+---+---+---+
    // | 1 | 2 | 2 | 1 | 1 | 1 | 2 |
    // +---+---+---+---+---+---+---+
    // | 1 | 1 | 2 | 2 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 0 | 1 | 1 | 2 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 0 | 0 | 1 | 1 | 1 | 1 | 1 |
    // +---+---+---+---+---+---+---+
    //
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[8] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    double coordinates_y[8] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    size_t const num_cells[2] = {7,7};
    size_t local_start[2][5][2] = {{{0,0},{0,2},{0,3},{4,2},{0,5}}, {{0,0}}};
    size_t local_count[2][5][2] = {{{7,2},{4,1},{4,2},{3,3},{7,2}}, {{7,7}}};
    int global_mask[7*7] = {
      0,0,1,1,1,1,1,
      0,1,1,2,2,2,2,
      1,1,2,2,2,2,2,
      1,2,2,1,1,1,2,
      1,2,2,1,0,1,2,
      1,2,2,1,1,1,2,
      1,2,2,2,2,2,2};
    int with_halo = 0;
    for (size_t i = 0; i <= num_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct yac_basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_cells,
        local_start[is_tgt][split_comm_rank],
        local_count[is_tgt][split_comm_rank], with_halo);

    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[is_tgt], grid_data),
       yac_basic_grid_empty_new(grid_names[is_tgt^1])};

    {
      int valid_mask_value = (is_tgt)?2:1;
      yac_coordinate_pointer point_coordinates =
        xmalloc(grid_data.num_cells * sizeof(*point_coordinates));
      int * mask = xmalloc(grid_data.num_cells * sizeof(*mask));
      for (size_t i = 0; i < grid_data.num_cells; ++i) {
        double * middle_point = point_coordinates[i];
        for (size_t k = 0; k < 3; ++k) middle_point[k] = 0.0;
        size_t * curr_vertices =
          grid_data.cell_to_vertex + grid_data.cell_to_vertex_offsets[i];
        size_t curr_num_vertices = grid_data.num_vertices_per_cell[i];
        for (size_t j = 0; j < curr_num_vertices; ++j) {
          double * curr_vertex_coord =
            grid_data.vertex_coordinates[curr_vertices[j]];
          for (size_t k = 0; k < 3; ++k)
            middle_point[k] += curr_vertex_coord[k];
        }
        normalise_vector(middle_point);
        mask[i] = global_mask[grid_data.cell_ids[i]] == valid_mask_value;
      }
      yac_basic_grid_add_coordinates_nocpy(
        grids[0], YAC_LOC_CELL, point_coordinates);
      yac_basic_grid_add_mask_nocpy(
        grids[0], YAC_LOC_CELL, mask, NULL);
    }

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CELL, .coordinates_idx = 0, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CELL, .coordinates_idx = 0, .masks_idx = 0};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                          num_src_fields, src_fields, tgt_field);

    struct yac_interp_spmap_config * spmap_config =
      yac_interp_spmap_config_new(
        YAC_RAD * 1.1,
        YAC_INTERP_SPMAP_MAX_SEARCH_DISTANCE_DEFAULT,
        YAC_INTERP_SPMAP_WEIGHTED_DEFAULT,
        YAC_INTERP_SPMAP_SCALE_CONFIG_DEFAULT);
    struct interp_method * method_stack[] =
      {yac_interp_method_spmap_new(
         spmap_config, YAC_INTERP_SPMAP_OVERWRITE_DEFAULT),
       yac_interp_method_fixed_new(-2.0), NULL};
    yac_interp_spmap_config_delete(spmap_config);

    struct yac_interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    for (size_t i = 0; i < num_reorder_types; ++i) {

      struct yac_interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_types[i], 1,
          YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0, NULL, 1, 1);

      // check generated interpolation
      {
        double * src_field = NULL;
        double ** src_fields = &src_field;
        double * tgt_field = NULL;

        if (is_tgt) {
          tgt_field = xmalloc(grid_data.num_cells * sizeof(*tgt_field));
          for (size_t i = 0; i < grid_data.num_cells; ++i) tgt_field[i] = -1;
        } else {
          src_field = xmalloc(grid_data.num_cells * sizeof(*src_field));
          for (size_t i = 0; i < grid_data.num_cells; ++i)
            src_field[i] = (double)(grid_data.cell_ids[i]);
        }

        yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

        double ref_tgt_field[7*7] = {
          -1,-1,-1,-1,-1,-1,-1,
          -1,-1,-1, 0, 0, 0, 0,
          -1,-1, 0, 0, 0, 0, 0,
          -1, 0, 0,-1,-1,-1, 0,
          -1, 0, 0,-1,-1,-1, 0,
          -1, 0, 0,-1,-1,-1, 0,
          -1, 0, 0, 0, 0, 0, 0};
        size_t coast_point[] = {2,3,4,5,6,
                                8,9,
                                14,15,
                                21,24,25,26,
                                28,31,33,
                                35,38,39,40,
                                42};
        size_t num_coast_points = sizeof(coast_point)/sizeof(coast_point[0]);
        size_t num_tgt_per_coast[] = {3,3,4,4,3,
                                      3,3,
                                      3,3,
                                      3,4,4,3,
                                      4,4,3,
                                      4,4,3,3,
                                      3};
        size_t tgts[] = {10,11,17, 10,11,17, 10,11,12,18, 11,12,13,19, 12,13,20,
                         16,17,23, 10,11,17,
                         22,23,29, 16,17,23,
                         22,23,29, 16,22,23,30, 11,17,18,19, 20,27,34,
                         22,29,30,36, 23,29,30,37, 27,34,41,
                         29,36,37,43, 30,36,37,44, 45,46,47, 34,41,48,
                         36,43,44};

        for (size_t i = 0, k = 0; i < num_coast_points; ++i) {
          size_t curr_num_tgt = num_tgt_per_coast[i];
          double curr_data = (double)(coast_point[i]) / (double)curr_num_tgt;
          for (size_t j = 0; j < curr_num_tgt; ++j, ++k)
            ref_tgt_field[tgts[k]] += curr_data;
        }

        if (is_tgt)
          for (size_t i = 0; i < grid_data.num_cells; ++i)
            if (fabs(ref_tgt_field[grid_data.cell_ids[i]] -
                     tgt_field[i]) > 1e-6)
              PUT_ERR("wrong interpolation result");

        free(tgt_field);
        free(src_field);
      }

      yac_interpolation_delete(interpolation);
    }

    yac_interp_weights_delete(weights);
    yac_interp_method_delete(method_stack);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(grids[1]);
    yac_basic_grid_delete(grids[0]);
  }

  { // parallel interpolation process
    // corner and cell ids for a 7 x 9 grid
    // 72-----73-----74-----75-----76-----77-----78-----79
    //  |      |      |      |      |      |      |      |
    //  |  56  |  57  |  58  |  59  |  60  |  61  |  62  |
    //  |      |      |      |      |      |      |      |
    // 64-----65-----66-----67-----68-----69-----70-----71
    //  |      |      |      |      |      |      |      |
    //  |  49  |  50  |  51  |  52  |  53  |  54  |  55  |
    //  |      |      |      |      |      |      |      |
    // 56-----57-----58-----59-----60-----61-----62-----63
    //  |      |      |      |      |      |      |      |
    //  |  42  |  43  |  44  |  45  |  46  |  47  |  48  |
    //  |      |      |      |      |      |      |      |
    // 48-----49-----50-----51-----52-----53-----54-----55
    //  |      |      |      |      |      |      |      |
    //  |  35  |  36  |  37  |  38  |  39  |  40  |  41  |
    //  |      |      |      |      |      |      |      |
    // 40-----41-----42-----43-----44-----45-----46-----47
    //  |      |      |      |      |      |      |      |
    //  |  28  |  29  |  30  |  31  |  32  |  33  |  34  |
    //  |      |      |      |      |      |      |      |
    // 32-----33-----34-----35-----36-----37-----38-----39
    //  |      |      |      |      |      |      |      |
    //  |  21  |  22  |  23  |  24  |  25  |  26  |  27  |
    //  |      |      |      |      |      |      |      |
    // 24-----25-----26-----27-----28-----29-----30-----31
    //  |      |      |      |      |      |      |      |
    //  |  14  |  15  |  16  |  17  |  18  |  19  |  20  |
    //  |      |      |      |      |      |      |      |
    // 16-----17-----18-----19-----20-----21-----22-----23
    //  |      |      |      |      |      |      |      |
    //  |  07  |  08  |  09  |  10  |  11  |  12  |  13  |
    //  |      |      |      |      |      |      |      |
    // 08-----09-----10-----11-----12-----13-----14-----15
    //  |      |      |      |      |      |      |      |
    //  |  00  |  01  |  02  |  03  |  04  |  05  |  06  |
    //  |      |      |      |      |      |      |      |
    // 00-----01-----02-----03-----04-----05-----06-----07
    //
    // the grid is distributed among the processes as follows:
    // (index == process)
    //
    // 4---4---4---4---4---4---4---4
    // | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
    // 4---4---4---4---4---4---4---4
    // | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
    // 4---4---4---4---4---4---4---4
    // | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
    // 3---3---3---3---4---4---4---4
    // | 2 | 2 | 2 | 2 | 3 | 3 | 3 |
    // 2---2---2---2---3---3---3---3
    // | 2 | 2 | 2 | 2 | 3 | 3 | 3 |
    // 1---1---1---1---2---2---2---2
    // | 1 | 1 | 1 | 1 | 3 | 3 | 3 |
    // 0---0---0---0---1---1---1---1
    // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
    // 0---0---0---0---0---0---0---0
    // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
    // 0---0---0---0---0---0---0---0
    // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
    // 0---0---0---0---0---0---0---0
    //
    // mask
    // land = 0
    // coast = 1
    // ocean = 2
    // +---+---+---+---+---+---+---+
    // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
    // +---+---+---+---+---+---+---+
    // | 1 | 0 | 0 | 0 | 0 | 1 | 2 |
    // +---+---+---+---+---+---+---+
    // | 2 | 1 | 0 | 0 | 1 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 2 | 2 | 1 | 1 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 2 | 2 | 1 | 1 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 2 | 1 | 0 | 0 | 1 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 1 | 0 | 0 | 0 | 0 | 1 | 2 |
    // +---+---+---+---+---+---+---+
    // | 0 | 0 | 0 | 0 | 0 | 0 | 1 |
    // +---+---+---+---+---+---+---+
    // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
    // +---+---+---+---+---+---+---+
    //
    //---------------
    // setup
    //---------------

    enum {NUM_X = 7, NUM_Y = 9};
    int is_tgt = split_comm_size == 1;
    double coordinates_x[NUM_X+1] = {0.0,1.01,2.03,3.06,4.1,5.15,6.21,7.28};
    double coordinates_y[NUM_Y+1] =
      {-1.0,0.0,1.01,2.03,3.06,4.1,5.15,6.21,7.28,8.36};
    size_t const num_cells[2] = {NUM_X,NUM_Y};
    size_t local_start[2][5][2] = {{{0,0},{0,3},{0,4},{4,3},{0,6}}, {{0,0}}};
    size_t local_count[2][5][2] = {{{7,3},{4,1},{4,2},{3,3},{7,3}},
                                   {{NUM_X,NUM_Y}}};
    int global_mask[NUM_X * NUM_Y] = {
      0,0,0,0,0,0,0,
      0,0,0,0,0,0,1,
      1,0,0,0,0,1,2,
      2,1,0,0,1,2,2,
      2,2,1,1,2,2,2,
      2,2,1,1,2,2,2,
      2,1,0,0,1,2,2,
      1,0,0,0,0,1,2,
      0,0,0,0,0,0,0};
    int with_halo = 0;
    double const grid_scale = 10.0;
    for (size_t i = 0; i <= NUM_X; ++i)
      coordinates_x[i] *= YAC_RAD * grid_scale;
    for (size_t i = 0; i <= NUM_Y; ++i)
      coordinates_y[i] *= YAC_RAD * grid_scale;

    struct yac_basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_cells,
        local_start[is_tgt][split_comm_rank],
        local_count[is_tgt][split_comm_rank], with_halo);

    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[is_tgt], grid_data),
       yac_basic_grid_empty_new(grid_names[is_tgt^1])};

    yac_coordinate_pointer point_coordinates =
      xmalloc(grid_data.num_cells * sizeof(*point_coordinates));
    {
      int valid_mask_value = (is_tgt)?2:1;
      int * mask = xmalloc(grid_data.num_cells * sizeof(*mask));
      for (size_t i = 0; i < grid_data.num_cells; ++i) {
        double * middle_point = point_coordinates[i];
        for (size_t k = 0; k < 3; ++k) middle_point[k] = 0.0;
        size_t * curr_vertices =
          grid_data.cell_to_vertex + grid_data.cell_to_vertex_offsets[i];
        size_t curr_num_vertices = grid_data.num_vertices_per_cell[i];
        for (size_t j = 0; j < curr_num_vertices; ++j) {
          double * curr_vertex_coord =
            grid_data.vertex_coordinates[curr_vertices[j]];
          for (size_t k = 0; k < 3; ++k)
            middle_point[k] += curr_vertex_coord[k];
        }
        normalise_vector(middle_point);
        mask[i] = global_mask[grid_data.cell_ids[i]] == valid_mask_value;
      }
      yac_basic_grid_add_coordinates_nocpy(
        grids[0], YAC_LOC_CELL, point_coordinates);
      yac_basic_grid_add_mask_nocpy(
        grids[0], YAC_LOC_CELL, mask, NULL);
    }

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    double grid_cell_areas[NUM_X*NUM_Y];
    for (int i = 0, k = 0; i < NUM_Y; ++i)
      for (int j = 0; j < NUM_X; ++j, ++k)
        grid_cell_areas[k] =
          fabs(
            (coordinates_x[j+1] - coordinates_x[j+0]) *
            (sin(coordinates_y[i+0]) - sin(coordinates_y[i+1])));

    char const * area_filename = "test_interp_method_spmap_area.nc";
    if (comm_rank == 0)
      utest_write_area_file(area_filename, grid_cell_areas, NUM_X, NUM_Y);
    MPI_Barrier(MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CELL, .coordinates_idx = 0, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CELL, .coordinates_idx = 0, .masks_idx = 0};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                          num_src_fields, src_fields, tgt_field);

    double const spread_distances[] = {3.7 * grid_scale, 0.0, 0.001};
    enum {NUM_SPREAD_DISTANCES =
      sizeof(spread_distances) / sizeof(spread_distances[0])};

    struct {
      struct yac_spmap_cell_area_config * src;
      struct yac_spmap_cell_area_config * tgt;
    } scale_configs[] =
      {{.src = yac_spmap_cell_area_config_yac_new(1.0),
        .tgt = yac_spmap_cell_area_config_yac_new(2.0)},
       {.src = yac_spmap_cell_area_config_yac_new(1.0),
        .tgt =
          yac_spmap_cell_area_config_file_new(area_filename, "area_1d", 0)},
       {.src =
          yac_spmap_cell_area_config_file_new(area_filename, "area_1d", 0),
        .tgt = yac_spmap_cell_area_config_yac_new(1.0)},
       {.src =
          yac_spmap_cell_area_config_file_new(area_filename, "area_2d", 0),
        .tgt =
          yac_spmap_cell_area_config_file_new(area_filename, "area_2d", 0)},
      };
    enum {NUM_SCALE_CONFIGS =
      sizeof(scale_configs) / sizeof(scale_configs[0])};

    char const * io_ranks[] = {"0", "3,5", "1,3,5"};
    enum {NUM_IO_CONFIGS = sizeof(io_ranks) / sizeof(io_ranks[0])};

    for (size_t spread_distance_idx = 0;
         spread_distance_idx < NUM_SPREAD_DISTANCES; ++spread_distance_idx) {

      for (size_t weight_type_idx = 0; weight_type_idx < num_weight_types;
          ++weight_type_idx) {

        for (size_t scale_type_idx = 0; scale_type_idx < num_scale_types;
            ++scale_type_idx) {

          for (size_t scale_config_idx = 0;
               scale_config_idx < NUM_SCALE_CONFIGS; ++scale_config_idx) {

            struct yac_spmap_scale_config * scale_config =
              yac_spmap_scale_config_new(
                scale_types[scale_type_idx],
                scale_configs[scale_config_idx].src,
                scale_configs[scale_config_idx].tgt);

            for (size_t io_config_idx = 0; io_config_idx < NUM_IO_CONFIGS;
                 ++io_config_idx) {

              // clear environment
              clear_yac_io_env();

              setenv("YAC_IO_RANK_LIST", io_ranks[io_config_idx], 1);
              setenv("YAC_IO_MAX_NUM_RANKS_PER_NODE", "12", 1);

              struct yac_interp_spmap_config * spmap_config =
                yac_interp_spmap_config_new(
                  YAC_RAD * spread_distances[spread_distance_idx],
                  YAC_INTERP_SPMAP_MAX_SEARCH_DISTANCE_DEFAULT,
                  weight_types[weight_type_idx],
                  scale_config);

              struct interp_method * method_stack[] =
                {yac_interp_method_spmap_new(
                   spmap_config, YAC_INTERP_SPMAP_OVERWRITE_DEFAULT),
                yac_interp_method_fixed_new(-2.0), NULL};
              yac_interp_spmap_config_delete(spmap_config);

              struct yac_interp_weights * weights =
                yac_interp_method_do_search(method_stack, interp_grid);

              for (size_t reorder_type_idx = 0;
                    reorder_type_idx < num_reorder_types;
                    ++reorder_type_idx) {

                struct yac_interpolation * interpolation =
                  yac_interp_weights_get_interpolation(
                    weights, reorder_types[reorder_type_idx], 1,
                    YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0, NULL, 1, 1);

                // check generated interpolation
                {
                  double * src_field = NULL;
                  double ** src_fields = &src_field;
                  double * tgt_field = NULL;

                  if (is_tgt) {
                    tgt_field =
                      xmalloc(grid_data.num_cells * sizeof(*tgt_field));
                    for (size_t i = 0; i < grid_data.num_cells; ++i)
                      tgt_field[i] = -1;
                  } else {
                    src_field =
                      xmalloc(grid_data.num_cells * sizeof(*src_field));
                    for (size_t i = 0; i < grid_data.num_cells; ++i)
                      src_field[i] = (double)(grid_data.cell_ids[i]);
                  }

                  yac_interpolation_execute(
                    interpolation, &src_fields, &tgt_field);

                  if (is_tgt) {

                    double ref_tgt_field[NUM_X * NUM_Y] = {
                      -1,-1,-1,-1,-1,-1,-1,
                      -1,-1,-1,-1,-1,-1,-1,
                      -1,-1,-1,-1,-1,-1, 0,
                      0,-1,-1,-1,-1, 0, 0,
                      0, 0,-1,-1, 0, 0, 0,
                      0, 0,-1,-1, 0, 0, 0,
                      0,-1,-1,-1,-1, 0, 0,
                      -1,-1,-1,-1,-1,-1, 0,
                      -1,-1,-1,-1,-1,-1,-1};
                    size_t coast_point[] = {14,22,30,37,43,49,
                                            13,19,25,31,38,46,54};
                    enum {
                      NUM_COAST_POINTS =
                        sizeof(coast_point)/sizeof(coast_point[0]),
                      MAX_NUM_TGT_PER_COAST = 12};
                    size_t num_tgt_per_coast[NUM_SPREAD_DISTANCES]
                                            [NUM_COAST_POINTS] =
                      {{6,6,6,6,6,6,
                        9,9,11,12,12,11,9},
                      {1,1,1,1,1,1,
                        1,1,1,1,1,1,1},
                      {1,1,1,1,1,1,
                        1,1,1,1,1,1,1}};
                    size_t tgts[NUM_SPREAD_DISTANCES]
                              [NUM_COAST_POINTS]
                              [MAX_NUM_TGT_PER_COAST] =
                      {{{21,28,29,35,36,42},
                        {21,28,29,35,36,42},
                        {21,28,29,35,36,42},
                        {21,28,29,35,36,42},
                        {21,28,29,35,36,42},
                        {21,28,29,35,36,42},

                        {20,26,27,32,33,34,39,40,41},
                        {20,26,27,32,33,34,39,40,41},
                        {20,26,27,32,33,34,39,40,41,47,48},
                        {20,26,27,32,33,34,39,40,41,47,48,55},
                        {20,26,27,32,33,34,39,40,41,47,48,55},
                        {26,27,32,33,34,39,40,41,47,48,55},
                        {32,33,34,39,40,41,47,48,55}},
                       {{21},
                        {21},
                        {29},
                        {36},
                        {42},
                        {42},

                        {20},
                        {20},
                        {26},
                        {32},
                        {39},
                        {47},
                        {55}},
                       {{21},
                        {21},
                        {29},
                        {36},
                        {42},
                        {42},

                        {20},
                        {20},
                        {26},
                        {32},
                        {39},
                        {47},
                        {55}}};

                    switch (weight_types[weight_type_idx]) {
                      default:
                      case(YAC_INTERP_SPMAP_AVG): {
                        for (size_t i = 0; i < NUM_COAST_POINTS; ++i) {
                          size_t curr_num_tgt =
                            num_tgt_per_coast[spread_distance_idx][i];
                          double curr_data =
                            (double)(coast_point[i]) / (double)curr_num_tgt;
                          for (size_t j = 0; j < curr_num_tgt; ++j)
                            ref_tgt_field[tgts[spread_distance_idx][i][j]] +=
                              curr_data *
                              utest_compute_scale(
                                scale_types[scale_type_idx],
                                grid_cell_areas[coast_point[i]] *
                                utest_get_sq_sphere_radius(
                                  yac_spmap_scale_config_get_src_cell_area_config(
                                    scale_config)),
                                grid_cell_areas[
                                  tgts[spread_distance_idx][i][j]] *
                                utest_get_sq_sphere_radius(
                                  yac_spmap_scale_config_get_tgt_cell_area_config(
                                    scale_config)));
                        }
                        break;
                      }
                      case(YAC_INTERP_SPMAP_DIST): {
                        for (size_t i = 0, k = 0; i < NUM_COAST_POINTS; ++i) {
                          size_t * curr_tgts = tgts[spread_distance_idx][i];
                          size_t curr_num_tgt =
                            num_tgt_per_coast[spread_distance_idx][i];
                          size_t curr_coast_point = coast_point[i];
                          double curr_src_data = (double)(curr_coast_point);
                          double inv_distances[curr_num_tgt];
                          double inv_distances_sum = 0.0;
                          for (size_t j = 0; j < curr_num_tgt; ++j)
                            inv_distances_sum +=
                              ((inv_distances[j] =
                                  1.0 / get_vector_angle(
                                          point_coordinates[curr_coast_point],
                                          point_coordinates[curr_tgts[j]])));
                          for (size_t j = 0; j < curr_num_tgt; ++j, ++k)
                            ref_tgt_field[curr_tgts[j]] +=
                              curr_src_data
                              * (inv_distances[j] / inv_distances_sum) *
                              utest_compute_scale(
                                scale_types[scale_type_idx],
                                grid_cell_areas[coast_point[i]] *
                                utest_get_sq_sphere_radius(
                                  yac_spmap_scale_config_get_src_cell_area_config(
                                    scale_config)),
                                grid_cell_areas[
                                  tgts[spread_distance_idx][i][j]] *
                                utest_get_sq_sphere_radius(
                                  yac_spmap_scale_config_get_tgt_cell_area_config(
                                    scale_config)));
                        }
                      }
                    }

                    for (size_t i = 0; i < grid_data.num_cells; ++i)
                      if (((ref_tgt_field[grid_data.cell_ids[i]] == 0.0) &&
                          (tgt_field[i] != -2.0)) ||
                          ((ref_tgt_field[grid_data.cell_ids[i]] != 0.0) &&
                          (fabs(ref_tgt_field[grid_data.cell_ids[i]] -
                                tgt_field[i]) > 1.0e-6)))
                        PUT_ERR("wrong interpolation result");

                    double src_sum = 0.0;
                    for (size_t i = 0; i < NUM_COAST_POINTS; ++i) {
                      double curr_src_data = (double)(coast_point[i]);
                      if (
                        (scale_types[scale_type_idx] ==
                        YAC_INTERP_SPMAP_SRCAREA) ||
                        (scale_types[scale_type_idx] ==
                        YAC_INTERP_SPMAP_FRACAREA))
                        curr_src_data *=
                          grid_cell_areas[coast_point[i]] *
                          utest_get_sq_sphere_radius(
                            yac_spmap_scale_config_get_src_cell_area_config(
                              scale_config));
                      src_sum += curr_src_data;
                    }
                    double tgt_sum = 0.0;
                    for (size_t i = 0; i < grid_data.num_cells; ++i) {
                      if ((tgt_field[i] != -1.0) &&
                          (tgt_field[i] != -2.0)) {
                        double curr_tgt_data = tgt_field[i];
                        if (
                          (scale_types[scale_type_idx] ==
                          YAC_INTERP_SPMAP_INVTGTAREA) ||
                          (scale_types[scale_type_idx] ==
                          YAC_INTERP_SPMAP_FRACAREA))
                          curr_tgt_data *=
                            grid_cell_areas[grid_data.cell_ids[i]] *
                            utest_get_sq_sphere_radius(
                              yac_spmap_scale_config_get_tgt_cell_area_config(
                                scale_config));
                        tgt_sum += curr_tgt_data;
                      }
                    }
                    if (fabs(src_sum - tgt_sum) > 1.0e-6)
                      PUT_ERR("wrong interpolation result (not conservative)");
                  }

                  free(tgt_field);
                  free(src_field);
                } // check

                yac_interpolation_delete(interpolation);
              } // reorder_type_idx

              yac_interp_weights_delete(weights);
              yac_interp_method_delete(method_stack);
            } // io_config_idx
            yac_spmap_scale_config_delete(scale_config);
          } // scale_config_idx
        } // scale_type_idx
      } // weight_type_idx
    } // spread_distance_idx

    for (size_t i = 0; i < NUM_SCALE_CONFIGS; ++i) {
      yac_spmap_cell_area_config_delete(scale_configs[i].src);
      yac_spmap_cell_area_config_delete(scale_configs[i].tgt);
    }
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(grids[1]);
    yac_basic_grid_delete(grids[0]);

    MPI_Barrier(MPI_COMM_WORLD);
    if (comm_rank == 0) unlink(area_filename);
  }

  { // parallel interpolation process
    // corner and cell ids for a 7 x 7 grid
    // 56-----57-----58-----59-----60-----61-----62-----63
    //  |      |      |      |      |      |      |      |
    //  |  42  |  43  |  44  |  45  |  46  |  47  |  48  |
    //  |      |      |      |      |      |      |      |
    // 48-----49-----50-----51-----52-----53-----54-----55
    //  |      |      |      |      |      |      |      |
    //  |  35  |  36  |  37  |  38  |  39  |  40  |  41  |
    //  |      |      |      |      |      |      |      |
    // 40-----41-----42-----43-----44-----45-----46-----47
    //  |      |      |      |      |      |      |      |
    //  |  28  |  29  |  30  |  31  |  32  |  33  |  34  |
    //  |      |      |      |      |      |      |      |
    // 32-----33-----34-----35-----36-----37-----38-----39
    //  |      |      |      |      |      |      |      |
    //  |  21  |  22  |  23  |  24  |  25  |  26  |  27  |
    //  |      |      |      |      |      |      |      |
    // 24-----25-----26-----27-----28-----29-----30-----31
    //  |      |      |      |      |      |      |      |
    //  |  14  |  15  |  16  |  17  |  18  |  19  |  20  |
    //  |      |      |      |      |      |      |      |
    // 16-----17-----18-----19-----20-----21-----22-----23
    //  |      |      |      |      |      |      |      |
    //  |  07  |  08  |  09  |  10  |  11  |  12  |  13  |
    //  |      |      |      |      |      |      |      |
    // 08-----09-----10-----11-----12-----13-----14-----15
    //  |      |      |      |      |      |      |      |
    //  |  00  |  01  |  02  |  03  |  04  |  05  |  06  |
    //  |      |      |      |      |      |      |      |
    // 00-----01-----02-----03-----04-----05-----06-----07
    //
    // the grid is distributed among the processes as follows:
    // (index == process)
    //
    // 4---4---4---4---4---4---4---4
    // | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
    // 4---4---4---4---4---4---4---4
    // | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
    // 3---3---3---3---4---4---4---4
    // | 2 | 2 | 2 | 2 | 3 | 3 | 3 |
    // 2---2---2---2---3---3---3---3
    // | 2 | 2 | 2 | 2 | 3 | 3 | 3 |
    // 1---1---1---1---2---2---2---2
    // | 1 | 1 | 1 | 1 | 3 | 3 | 3 |
    // 0---0---0---0---1---1---1---1
    // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
    // 0---0---0---0---0---0---0---0
    // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
    // 0---0---0---0---0---0---0---0
    //
    // mask
    // land = 0
    // coast = 1
    // ocean = 2
    // +---+---+---+---+---+---+---+
    // | 0 | 0 | 0 | 1 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 0 | 0 | 0 | 0 | 1 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 0 | 0 | 0 | 0 | 1 | 1 | 2 |
    // +---+---+---+---+---+---+---+
    // | 0 | 0 | 0 | 1 | 0 | 0 | 1 |
    // +---+---+---+---+---+---+---+
    // | 0 | 0 | 1 | 0 | 0 | 0 | 0 |
    // +---+---+---+---+---+---+---+
    // | 0 | 1 | 0 | 0 | 0 | 0 | 0 |
    // +---+---+---+---+---+---+---+
    // | 1 | 0 | 0 | 0 | 0 | 0 | 0 |
    // +---+---+---+---+---+---+---+
    //
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[8] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    double coordinates_y[8] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    size_t const num_cells[2] = {7,7};
    size_t local_start[2][5][2] = {{{0,0},{0,2},{0,3},{4,2},{0,5}}, {{0,0}}};
    size_t local_count[2][5][2] = {{{7,2},{4,1},{4,2},{3,3},{7,2}}, {{7,7}}};
    int global_mask[7*7] = {
      1,0,0,0,0,0,0,
      0,1,0,0,0,0,0,
      0,0,1,0,0,0,0,
      0,0,0,1,0,0,1,
      0,0,0,0,1,1,2,
      0,0,0,0,1,2,2,
      0,0,0,1,2,2,2};
    int with_halo = 0;
    for (size_t i = 0; i <= num_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct yac_basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_cells,
        local_start[is_tgt][split_comm_rank],
        local_count[is_tgt][split_comm_rank], with_halo);

    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[is_tgt], grid_data),
       yac_basic_grid_empty_new(grid_names[is_tgt^1])};

    yac_coordinate_pointer point_coordinates =
      xmalloc(grid_data.num_cells * sizeof(*point_coordinates));
    {
      int valid_mask_value = (is_tgt)?2:1;
      int * mask = xmalloc(grid_data.num_cells * sizeof(*mask));
      for (size_t i = 0; i < grid_data.num_cells; ++i) {
        double * middle_point = point_coordinates[i];
        for (size_t k = 0; k < 3; ++k) middle_point[k] = 0.0;
        size_t * curr_vertices =
          grid_data.cell_to_vertex + grid_data.cell_to_vertex_offsets[i];
        size_t curr_num_vertices = grid_data.num_vertices_per_cell[i];
        for (size_t j = 0; j < curr_num_vertices; ++j) {
          double * curr_vertex_coord =
            grid_data.vertex_coordinates[curr_vertices[j]];
          for (size_t k = 0; k < 3; ++k)
            middle_point[k] += curr_vertex_coord[k];
        }
        normalise_vector(middle_point);
        mask[i] = global_mask[grid_data.cell_ids[i]] == valid_mask_value;
      }
      yac_basic_grid_add_coordinates_nocpy(
        grids[0], YAC_LOC_CELL, point_coordinates);
      yac_basic_grid_add_mask_nocpy(
        grids[0], YAC_LOC_CELL, mask, NULL);
    }

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CELL, .coordinates_idx = 0, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CELL, .coordinates_idx = 0, .masks_idx = 0};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                          num_src_fields, src_fields, tgt_field);

    struct yac_interp_spmap_config * spmap_config =
      yac_interp_spmap_config_new(
        YAC_INTERP_SPMAP_SPREAD_DISTANCE_DEFAULT,
        YAC_RAD * 3.0,
        YAC_INTERP_SPMAP_WEIGHTED_DEFAULT,
        YAC_INTERP_SPMAP_SCALE_CONFIG_DEFAULT);

    struct interp_method * method_stack[] =
      {yac_interp_method_spmap_new(
         spmap_config, YAC_INTERP_SPMAP_OVERWRITE_DEFAULT),
       yac_interp_method_fixed_new(-2.0), NULL};
    yac_interp_spmap_config_delete(spmap_config);

    struct yac_interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    for (size_t i = 0; i < num_reorder_types; ++i) {

      struct yac_interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_types[i], 1,
          YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0, NULL, 1, 1);

      // check generated interpolation
      {
        double * src_field = NULL;
        double ** src_fields = &src_field;
        double * tgt_field = NULL;

        if (is_tgt) {
          tgt_field = xmalloc(grid_data.num_cells * sizeof(*tgt_field));
          for (size_t i = 0; i < grid_data.num_cells; ++i) tgt_field[i] = -1;
        } else {
          src_field = xmalloc(grid_data.num_cells * sizeof(*src_field));
          for (size_t i = 0; i < grid_data.num_cells; ++i)
            src_field[i] = (double)(grid_data.cell_ids[i]);
        }

        yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

        if (is_tgt) {

          double ref_tgt_field[7*7] = {
             -1,-1,-1,-1,-1,-1,-1,
             -1,-1,-1,-1,-1,-1,-1,
             -1,-1,-1,-1,-1,-1,-1,
             -1,-1,-1,-1,-1,-1,-1,
             -1,-1,-1,-1,-1,-1, 0,
             -1,-1,-1,-1,-1, 0,-2,
             -1,-1,-1,-1, 0,-2,-2};
          size_t coast_point[] = {0,8,16,24,27,32,33,39,45};
          size_t num_coast_points = sizeof(coast_point)/sizeof(coast_point[0]);
          size_t num_tgt_per_coast[] = {0,0,0,1,1,1,1,1,1};
          size_t tgts[] = {40,34,40,34,40,46};

          for (size_t i = 0, k = 0; i < num_coast_points; ++i) {
            size_t curr_num_tgt = num_tgt_per_coast[i];
            if (curr_num_tgt == 0) continue;
            double curr_data = (double)(coast_point[i]) / (double)curr_num_tgt;
            for (size_t j = 0; j < curr_num_tgt; ++j, ++k)
              ref_tgt_field[tgts[k]] += curr_data;
          }

          for (size_t i = 0; i < grid_data.num_cells; ++i)
            if (fabs(ref_tgt_field[grid_data.cell_ids[i]] -
                     tgt_field[i]) > 1e-6)
              PUT_ERR("wrong interpolation result");
        }

        free(tgt_field);
        free(src_field);
      }

      yac_interpolation_delete(interpolation);
    }

    yac_interp_weights_delete(weights);
    yac_interp_method_delete(method_stack);

    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(grids[1]);
    yac_basic_grid_delete(grids[0]);
  }

  { // testing YAC_INTERP_SPMAP_DIST
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[4] = {0.0,1.0,2.0,3.0};
    double coordinates_y[4] = {0.0,1.0,2.0,3.0};
    size_t const num_cells[2] = {3,3};
    size_t local_start[2][5][2] = {{{0,0},{0,0},{0,0},{0,0},{0,0}}, {{0,0}}};
    size_t local_count[2][5][2] = {{{3,3},{3,3},{3,3},{3,3},{3,3}}, {{3,3}}};
    int global_mask[3*3] = {
      1,1,2,
      1,3,2,
      2,2,2};
    int with_halo = 0;

    for (size_t i = 0; i <= num_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct yac_basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_cells,
        local_start[is_tgt][split_comm_rank],
        local_count[is_tgt][split_comm_rank], with_halo);

    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[is_tgt], grid_data),
       yac_basic_grid_empty_new(grid_names[is_tgt^1])};

    yac_coordinate_pointer point_coordinates =
      xmalloc(grid_data.num_cells * sizeof(*point_coordinates));
    {
      int valid_mask_value = (is_tgt)?2:1;
      int * mask = xmalloc(grid_data.num_cells * sizeof(*mask));
      for (size_t i = 0; i < grid_data.num_cells; ++i) {
        double * middle_point = point_coordinates[i];
        for (size_t k = 0; k < 3; ++k) middle_point[k] = 0.0;
        size_t * curr_vertices =
          grid_data.cell_to_vertex + grid_data.cell_to_vertex_offsets[i];
        size_t curr_num_vertices = grid_data.num_vertices_per_cell[i];
        for (size_t j = 0; j < curr_num_vertices; ++j) {
          double * curr_vertex_coord =
            grid_data.vertex_coordinates[curr_vertices[j]];
          for (size_t k = 0; k < 3; ++k)
            middle_point[k] += curr_vertex_coord[k];
        }
        normalise_vector(middle_point);
        mask[i] = (global_mask[grid_data.cell_ids[i]] & valid_mask_value) > 0;
      }
      yac_basic_grid_add_coordinates_nocpy(
        grids[0], YAC_LOC_CELL, point_coordinates);
      yac_basic_grid_add_mask_nocpy(
        grids[0], YAC_LOC_CELL, mask, NULL);
    }

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CELL, .coordinates_idx = 0, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CELL, .coordinates_idx = 0, .masks_idx = 0};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                          num_src_fields, src_fields, tgt_field);

    struct yac_interp_spmap_config * spmap_config =
      yac_interp_spmap_config_new(
        1.1 * YAC_RAD,
        YAC_INTERP_SPMAP_MAX_SEARCH_DISTANCE_DEFAULT,
        YAC_INTERP_SPMAP_DIST,
        YAC_INTERP_SPMAP_SCALE_CONFIG_DEFAULT);

    struct interp_method * method_stack[] =
      {yac_interp_method_spmap_new(
         spmap_config, YAC_INTERP_SPMAP_OVERWRITE_DEFAULT),
       yac_interp_method_fixed_new(-2.0), NULL};
    yac_interp_spmap_config_delete(spmap_config);

    struct yac_interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    struct yac_interpolation * interpolation =
      yac_interp_weights_get_interpolation(
        weights, YAC_MAPPING_ON_SRC, 1,
        YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0, NULL, 1, 1);

    // check generated interpolation
    {
      double * src_field = NULL;
      double ** src_fields = &src_field;
      double * tgt_field = NULL;

      if (is_tgt) {
        tgt_field = xmalloc(grid_data.num_cells * sizeof(*tgt_field));
        for (size_t i = 0; i < grid_data.num_cells; ++i) tgt_field[i] = -1;
      } else {
        src_field = xmalloc(grid_data.num_cells * sizeof(*src_field));
        for (size_t i = 0; i < grid_data.num_cells; ++i)
          src_field[i] = (double)(grid_data.cell_ids[i] + 1);
      }

      yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

      if (is_tgt) {
        double ref_tgt_field[3*3] = {-1,-1,0, -1,0,0, -2,0,-2};
        size_t coast_point[] = {0,1,3,4};
        enum {
          NUM_COAST_POINTS = sizeof(coast_point)/sizeof(coast_point[0]),
          MAX_NUM_TGT_PER_COAST = 3};
        size_t num_tgt_per_coast[NUM_COAST_POINTS] = {3,2,3,1};
        size_t tgts[NUM_COAST_POINTS][MAX_NUM_TGT_PER_COAST] =
          {{4,5,7},{2,5},{4,5,7},{4}};


        for (size_t i = 0, k = 0; i < NUM_COAST_POINTS; ++i) {
          size_t * curr_tgts = tgts[i];
          size_t curr_num_tgt = num_tgt_per_coast[i];
          size_t curr_coast_point = coast_point[i];
          double curr_src_data = (double)(curr_coast_point + 1);
          double inv_distances[curr_num_tgt];
          double inv_distances_sum = 0.0;
          if ( curr_num_tgt == 1) {
            ref_tgt_field[curr_tgts[0]] += curr_src_data;
          } else {
            for (size_t j = 0; j < curr_num_tgt; ++j)
              inv_distances_sum +=
                ((inv_distances[j] =
                    1.0 / get_vector_angle(
                            point_coordinates[curr_coast_point],
                            point_coordinates[curr_tgts[j]])));
            for (size_t j = 0; j < curr_num_tgt; ++j, ++k)
              ref_tgt_field[curr_tgts[j]] +=
                curr_src_data * (inv_distances[j] / inv_distances_sum);
          }
        }

        for (size_t i = 0; i < grid_data.num_cells; ++i)
          if (fabs(ref_tgt_field[grid_data.cell_ids[i]] - tgt_field[i]) > 1e-9)
            PUT_ERR("wrong interpolation result");
      }

      free(tgt_field);
      free(src_field);
    }

    yac_interpolation_delete(interpolation);

    yac_interp_weights_delete(weights);
    yac_interp_method_delete(method_stack);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(grids[1]);
    yac_basic_grid_delete(grids[0]);
  }

  { // parallel interpolation process
    // corner and cell ids for a 7 x 7 grid
    // 56-----57-----58-----59-----60-----61-----62-----63
    //  |      |      |      |      |      |      |      |
    //  |  42  |  43  |  44  |  45  |  46  |  47  |  48  |
    //  |      |      |      |      |      |      |      |
    // 48-----49-----50-----51-----52-----53-----54-----55
    //  |      |      |      |      |      |      |      |
    //  |  35  |  36  |  37  |  38  |  39  |  40  |  41  |
    //  |      |      |      |      |      |      |      |
    // 40-----41-----42-----43-----44-----45-----46-----47
    //  |      |      |      |      |      |      |      |
    //  |  28  |  29  |  30  |  31  |  32  |  33  |  34  |
    //  |      |      |      |      |      |      |      |
    // 32-----33-----34-----35-----36-----37-----38-----39
    //  |      |      |      |      |      |      |      |
    //  |  21  |  22  |  23  |  24  |  25  |  26  |  27  |
    //  |      |      |      |      |      |      |      |
    // 24-----25-----26-----27-----28-----29-----30-----31
    //  |      |      |      |      |      |      |      |
    //  |  14  |  15  |  16  |  17  |  18  |  19  |  20  |
    //  |      |      |      |      |      |      |      |
    // 16-----17-----18-----19-----20-----21-----22-----23
    //  |      |      |      |      |      |      |      |
    //  |  07  |  08  |  09  |  10  |  11  |  12  |  13  |
    //  |      |      |      |      |      |      |      |
    // 08-----09-----10-----11-----12-----13-----14-----15
    //  |      |      |      |      |      |      |      |
    //  |  00  |  01  |  02  |  03  |  04  |  05  |  06  |
    //  |      |      |      |      |      |      |      |
    // 00-----01-----02-----03-----04-----05-----06-----07
    //
    // the grid is distributed among the processes as follows:
    // (index == process)
    //
    // 4---4---4---4---4---4---4---4
    // | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
    // 4---4---4---4---4---4---4---4
    // | 4 | 4 | 4 | 4 | 4 | 4 | 4 |
    // 3---3---3---3---4---4---4---4
    // | 2 | 2 | 2 | 2 | 3 | 3 | 3 |
    // 2---2---2---2---3---3---3---3
    // | 2 | 2 | 2 | 2 | 3 | 3 | 3 |
    // 1---1---1---1---2---2---2---2
    // | 1 | 1 | 1 | 1 | 3 | 3 | 3 |
    // 0---0---0---0---1---1---1---1
    // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
    // 0---0---0---0---0---0---0---0
    // | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
    // 0---0---0---0---0---0---0---0
    //
    // global_mask
    // land = 0
    // coast = 1
    // ocean = 2
    // +---+---+---+---+---+---+---+
    // | 1 | 2 | 2 | 2 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 1 | 2 | 2 | 2 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 1 | 2 | 2 | 2 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 1 | 2 | 2 | 2 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 1 | 1 | 2 | 2 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 0 | 1 | 1 | 2 | 2 | 2 | 2 |
    // +---+---+---+---+---+---+---+
    // | 0 | 0 | 1 | 1 | 1 | 1 | 1 |
    // +---+---+---+---+---+---+---+
    //
    //---------------
    // setup
    //---------------

    int is_tgt = split_comm_size == 1;
    double coordinates_x[8] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    double coordinates_y[8] = {0.0,1.0,2.0,3.0,4.0,5.0,6.0,7.0};
    size_t const num_cells[2] = {7,7};
    size_t local_start[2][5][2] = {{{0,0},{0,2},{0,3},{4,2},{0,5}}, {{0,0}}};
    size_t local_count[2][5][2] = {{{7,2},{4,1},{4,2},{3,3},{7,2}}, {{7,7}}};
    int global_mask[7*7] = {
      0,0,1,1,1,1,1,
      0,1,1,2,2,2,2,
      1,1,2,2,2,2,2,
      1,2,2,2,2,2,2,
      1,2,2,2,2,2,2,
      1,2,2,2,2,2,2,
      1,2,2,2,2,2,2};
    int with_halo = 0;
    for (size_t i = 0; i <= num_cells[0]; ++i)
      coordinates_x[i] *= YAC_RAD;
    for (size_t i = 0; i <= num_cells[1]; ++i)
      coordinates_y[i] *= YAC_RAD;

    struct yac_basic_grid_data grid_data =
      yac_generate_basic_grid_data_reg2d(
        coordinates_x, coordinates_y, num_cells,
        local_start[is_tgt][split_comm_rank],
        local_count[is_tgt][split_comm_rank], with_halo);

    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[is_tgt], grid_data),
       yac_basic_grid_empty_new(grid_names[is_tgt^1])};

    {
      int valid_mask_value = (is_tgt)?2:1;
      yac_coordinate_pointer point_coordinates =
        xmalloc(grid_data.num_cells * sizeof(*point_coordinates));
      int * mask = xmalloc(grid_data.num_cells * sizeof(*mask));
      for (size_t i = 0; i < grid_data.num_cells; ++i) {
        double * middle_point = point_coordinates[i];
        for (size_t k = 0; k < 3; ++k) middle_point[k] = 0.0;
        size_t * curr_vertices =
          grid_data.cell_to_vertex + grid_data.cell_to_vertex_offsets[i];
        size_t curr_num_vertices = grid_data.num_vertices_per_cell[i];
        for (size_t j = 0; j < curr_num_vertices; ++j) {
          double * curr_vertex_coord =
            grid_data.vertex_coordinates[curr_vertices[j]];
          for (size_t k = 0; k < 3; ++k)
            middle_point[k] += curr_vertex_coord[k];
        }
        normalise_vector(middle_point);
        mask[i] = global_mask[grid_data.cell_ids[i]] == valid_mask_value;
      }
      yac_basic_grid_add_coordinates_nocpy(
        grids[0], YAC_LOC_CELL, point_coordinates);
      yac_basic_grid_add_mask_nocpy(
        grids[0], YAC_LOC_CELL, mask, NULL);
    }

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CELL, .coordinates_idx = 0, .masks_idx = 0}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CELL, .coordinates_idx = 0, .masks_idx = 0};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(grid_pair, grid_names[0], grid_names[1],
                          num_src_fields, src_fields, tgt_field);

    enum {OVERWRITE_CONFIG_COUNT = 3};
    struct yac_spmap_overwrite_config *
      overwrite_configs[OVERWRITE_CONFIG_COUNT+1];
    overwrite_configs[OVERWRITE_CONFIG_COUNT] = NULL;

    // select bounding circle containing cell with global id 8 and spread its
    // value to all surrounding cell of the initial target
    {
      struct yac_point_selection * src_point_selection =
        yac_point_selection_bnd_circle_new(
          1.5 * YAC_RAD, 1.5 * YAC_RAD, 0.1 * YAC_RAD);
      struct yac_interp_spmap_config * spmap_config =
        yac_interp_spmap_config_new(
          1.5 * YAC_RAD,
          YAC_INTERP_SPMAP_MAX_SEARCH_DISTANCE_DEFAULT,
          YAC_INTERP_SPMAP_WEIGHTED_DEFAULT,
          YAC_INTERP_SPMAP_SCALE_CONFIG_DEFAULT);
      overwrite_configs[0] =
        yac_spmap_overwrite_config_new(src_point_selection, spmap_config);
      yac_interp_spmap_config_delete(spmap_config);
      yac_point_selection_delete(src_point_selection);
    }

    // select bounding circle containing cells with global id 8, 9, 15 and
    // spread its value to all orthogonally adjacent cell of the initial target
    // (src 8 will be dealt with by the first overwrite config)
    {
      struct yac_point_selection * src_point_selection =
        yac_point_selection_bnd_circle_new(
          1.5 * YAC_RAD, 1.5 * YAC_RAD, 1.1 * YAC_RAD);
      struct yac_interp_spmap_config * spmap_config =
        yac_interp_spmap_config_new(
          1.1 * YAC_RAD,
          YAC_INTERP_SPMAP_MAX_SEARCH_DISTANCE_DEFAULT,
          YAC_INTERP_SPMAP_WEIGHTED_DEFAULT,
          YAC_INTERP_SPMAP_SCALE_CONFIG_DEFAULT);
      overwrite_configs[1] =
        yac_spmap_overwrite_config_new(src_point_selection, spmap_config);
      yac_interp_spmap_config_delete(spmap_config);
      yac_point_selection_delete(src_point_selection);
    }

    // select bounding circle containing cell with global id 6 and
    // and set its maximum search distance so low, that no matching target
    // point is found
    {
      struct yac_point_selection * src_point_selection =
        yac_point_selection_bnd_circle_new(
          6.5 * YAC_RAD, 0.5 * YAC_RAD, 0.1 * YAC_RAD);
      struct yac_interp_spmap_config * spmap_config =
        yac_interp_spmap_config_new(
          YAC_INTERP_SPMAP_SPREAD_DISTANCE_DEFAULT,
          0.1 * YAC_RAD,
          YAC_INTERP_SPMAP_WEIGHTED_DEFAULT,
          YAC_INTERP_SPMAP_SCALE_CONFIG_DEFAULT);
      overwrite_configs[2] =
        yac_spmap_overwrite_config_new(src_point_selection, spmap_config);
      yac_interp_spmap_config_delete(spmap_config);
      yac_point_selection_delete(src_point_selection);
    }

    double const fixed_value = -2.0;
    struct interp_method * method_stack[] =
      {yac_interp_method_spmap_new(
         YAC_INTERP_SPMAP_DEFAULT_CONFIG,
         (struct yac_spmap_overwrite_config const * const *)overwrite_configs),
       yac_interp_method_fixed_new(fixed_value), NULL};

    struct yac_interp_weights * weights =
      yac_interp_method_do_search(method_stack, interp_grid);

    for (size_t reorder_idx = 0; reorder_idx < num_reorder_types;
         ++reorder_idx) {

      struct yac_interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_types[reorder_idx], 1,
          YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0, NULL, 1, 1);

      // check generated interpolation
      {
        double * src_field = NULL;
        double ** src_fields = &src_field;
        double * tgt_field = NULL;

        if (is_tgt) {
          tgt_field = xmalloc(grid_data.num_cells * sizeof(*tgt_field));
          for (size_t i = 0; i < grid_data.num_cells; ++i) tgt_field[i] = -1;
        } else {
          src_field = xmalloc(grid_data.num_cells * sizeof(*src_field));
          for (size_t i = 0; i < grid_data.num_cells; ++i)
            src_field[i] = (double)(grid_data.cell_ids[i]);
        }

        yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

        double ref_tgt_field[7*7];
        for (size_t i = 0; i < 7*7; ++i) ref_tgt_field[i] = 0.0;

        size_t unset_tgt[] = {0,1,2,3,4,5,6, 7,8,9, 14,15, 21, 28, 35, 42};
        for (size_t i = 0; i < sizeof(unset_tgt)/sizeof(unset_tgt[0]); ++i)
          ref_tgt_field[unset_tgt[i]] = -1.0;

        struct {
          yac_int global_id;
          size_t tgt_idx;
        } direct[] =
        {{.global_id =  2, .tgt_idx = 10}, {.global_id =  3, .tgt_idx = 10},
         {.global_id =  4, .tgt_idx = 11}, {.global_id =  5, .tgt_idx = 12},
         {.global_id = 14, .tgt_idx = 22}, {.global_id = 21, .tgt_idx = 22},
         {.global_id = 28, .tgt_idx = 29}, {.global_id = 35, .tgt_idx = 36},
         {.global_id = 42, .tgt_idx = 43}};
        for (size_t i = 0; i < sizeof(direct)/sizeof(direct[0]); ++i)
          ref_tgt_field[direct[i].tgt_idx] += (double)(direct[i].global_id);

        struct {
          yac_int global_id;
          size_t * tgt_idx;
          size_t tgt_count;
        } wsum[] =
        {{.global_id = 8,
          .tgt_idx = (size_t[]){10, 16,17, 22,23,24},
          .tgt_count = 6},
         {.global_id = 9,
          .tgt_idx = (size_t[]){10,11, 17},
          .tgt_count = 3},
         {.global_id = 15,
          .tgt_idx = (size_t[]){16,17, 23},
          .tgt_count = 3}};
        for (size_t i = 0; i < sizeof(wsum)/sizeof(wsum[0]); ++i)
          for (size_t j = 0; j < wsum[i].tgt_count; ++j)
            ref_tgt_field[wsum[i].tgt_idx[j]] +=
              (double)(wsum[i].global_id) / (double)(wsum[i].tgt_count);

        size_t fixed_tgt[] =
          {13, 18,19,20, 25,26,27, 30,31,32,33,34,
           37,38,39,40,41, 44,45,46,47,48};
        for (size_t i = 0; i < sizeof(fixed_tgt)/sizeof(fixed_tgt[0]); ++i)
          ref_tgt_field[fixed_tgt[i]] = fixed_value;

        if (is_tgt)
          for (size_t i = 0; i < grid_data.num_cells; ++i)
            if (fabs(
                  ref_tgt_field[grid_data.cell_ids[i]] -
                  tgt_field[i]) > 1e-9) PUT_ERR("wrong interpolation result");

        free(tgt_field);
        free(src_field);
      }

      yac_interpolation_delete(interpolation);
    }

    yac_interp_weights_delete(weights);
    yac_interp_method_delete(method_stack);
    for (size_t i = 0; i < OVERWRITE_CONFIG_COUNT; ++i)
      yac_spmap_overwrite_config_delete(overwrite_configs[i]);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(grids[1]);
    yac_basic_grid_delete(grids[0]);
  }

  MPI_Comm_free(&split_comm);

  xt_finalize();

  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static double utest_get_sq_sphere_radius(
  struct yac_spmap_cell_area_config const * cell_area_config) {

  double sq_sphere_radius;
  if (
    yac_spmap_cell_area_config_get_type(cell_area_config) ==
    YAC_INTERP_SPMAP_CELL_AREA_YAC) {
    double sphere_radius =
      yac_spmap_cell_area_config_get_sphere_radius(cell_area_config);
    sq_sphere_radius = sphere_radius * sphere_radius;
  } else sq_sphere_radius = 1.0;
  return sq_sphere_radius;
}

static double utest_compute_scale(
  enum yac_interp_spmap_scale_type scale_type,
  double src_cell_area, double tgt_cell_area) {

  double scale = 1.0;

  if ((scale_type == YAC_INTERP_SPMAP_SRCAREA) ||
      (scale_type == YAC_INTERP_SPMAP_FRACAREA))
    scale *= src_cell_area;
  if ((scale_type == YAC_INTERP_SPMAP_INVTGTAREA) ||
      (scale_type == YAC_INTERP_SPMAP_FRACAREA))
    scale /= tgt_cell_area;

  return scale;
}

static void utest_write_area_file(
  char const * filename, double * cell_areas, size_t dim_x, size_t dim_y) {

  // create file
  int ncid;
  yac_nc_create(filename, NC_CLOBBER, &ncid);

  int dim_x_id;
  int dim_y_id;
  int dim_xy_id;

  // define dimensions
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "x", dim_x, &dim_x_id));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "y", dim_y, &dim_y_id));
  YAC_HANDLE_ERROR(nc_def_dim(ncid, "xy", dim_x * dim_y, &dim_xy_id));

  char const * area_2d_varname = "area_2d";
  char const * area_1d_varname = "area_1d";

  int area_2d_dim_ids[2] = {dim_x_id, dim_y_id};
  int area_1d_dim_ids[1] = {dim_xy_id};

  int area_2d_varid;
  int area_1d_varid;

  // define variable
  YAC_HANDLE_ERROR(
    nc_def_var(
      ncid, area_2d_varname, NC_DOUBLE, 2, area_2d_dim_ids, &area_2d_varid));
  YAC_HANDLE_ERROR(
    nc_def_var(
      ncid, area_1d_varname, NC_DOUBLE, 1, area_1d_dim_ids, &area_1d_varid));

  // end definition
  YAC_HANDLE_ERROR(nc_enddef(ncid));

  // write grid data
  YAC_HANDLE_ERROR(nc_put_var_double(ncid, area_2d_varid, cell_areas));
  YAC_HANDLE_ERROR(nc_put_var_double(ncid, area_1d_varid, cell_areas));

  // close file
  YAC_HANDLE_ERROR(nc_close(ncid));
}
