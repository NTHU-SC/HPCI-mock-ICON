// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "tests.h"
#include "test_common.h"
#include "read_icon_grid.h"
#include "geometry.h"
#include "grid2vtk.h"
#include "io_utils.h"
#include "test_read_icon_common.h"

#include <netcdf.h>

#include <mpi.h>

/** \file test_read_icon_parallel_c.c
 *  \test
 * This contains examples for parallel use of yac_read_icon_basic_grid_data.
 */

int main(void) {

  MPI_Init(NULL, NULL);

  int comm_size, comm_rank;

  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  if ((comm_size != 4) && (comm_rank == 0)) {
    fputs("wrong number of processes (has to be four)\n", stderr);
    exit(EXIT_FAILURE);
  }

  char const filename[] = "test_read_icon_parallel_grid_c.nc";



  int ref_nbr_cells[4] = {11, 9, 14, 9};
  int ref_nbr_vertices[4] = {11, 10, 13, 10};
  int ref_cell_to_vertex[14][3];
  int ref_global_cell_ids[4][14] = {{0,1,2,3,4,5,7,8,9,10,11},
                                    {4,5,6,7,2,3,8,9,15},
                                    {8,9,10,11,0,1,2,3,4,5,7,12,14,15},
                                    {12,13,14,15,7,8,9,10,11}};
  int ref_cell_owner[14];
  int ref_global_vertex_ids[4][13] = {{0,1,2,5,6,7,9,10,11,12,13},
                                      {5,6,7,8,9,10,11,12,13,14},
                                      {0,1,2,3,5,6,7,8,9,10,11,12,13},
                                      {1,2,3,4,6,7,8,10,11,13}};
  int ref_vertex_owner[13];

  {
    int all_cell_to_vertex[16][3] =
      {{0,1,5}, {1,5,6}, {5,6,9}, {6,9,10}, {9,10,12}, {10,12,13}, {12,13,14},
      {10,11,13}, {7,10,11}, {6,7,10}, {2,6,7}, {1,2,6}, {2,3,7}, {3,4,8},
      {3,7,8}, {7,8,11}};
    for (int i = 0; i < ref_nbr_cells[comm_rank]; ++i) {
      for (int j = 0; j < 3; ++j) {
        for (int k = 0; k < ref_nbr_vertices[comm_rank]; ++k) {
          if (all_cell_to_vertex[ref_global_cell_ids[comm_rank][i]][j] ==
              ref_global_vertex_ids[comm_rank][k]) {
            ref_cell_to_vertex[i][j] = k;
            break;
          }
        }
      }
    }
    int global_cell_owner[16] = {0,0,0,0, 1,1,1,1, 2,2,2,2, 3,3,3,3};
    for (int i = 0; i < ref_nbr_cells[comm_rank]; ++i) {
      ref_cell_owner[i] =
        global_cell_owner[ref_global_cell_ids[comm_rank][i]];
      if (ref_cell_owner[i] == comm_rank) ref_cell_owner[i] = -1;
    }
    int global_vertex_owner[15] = {0,2,3,3,3,0,2,3,3,1,2,3,1,1,1};
    for (int i = 0; i < ref_nbr_vertices[comm_rank]; ++i) {
      ref_vertex_owner[i] =
        global_vertex_owner[ref_global_vertex_ids[comm_rank][i]];
      if (ref_vertex_owner[i] == comm_rank) ref_vertex_owner[i] = -1;
    }
  }

  enum coord_units coord_unit[] = {RAD, DEG};
  enum {NUM_COORD_UNIT = sizeof(coord_unit) / sizeof(coord_unit[0])};

  for (int coord_unit_idx = 0; coord_unit_idx < NUM_COORD_UNIT;
       ++coord_unit_idx) {

    if (comm_rank == 0)
      write_test_grid_file(filename, coord_unit[coord_unit_idx]);

    // ensure that the grid file exists
    MPI_Barrier(MPI_COMM_WORLD);

    { // test yac_read_icon_grid_information_parallel (incl. halos)
      int nbr_vertices, nbr_cells;
      int * num_vertices_per_cell;
      int * cell_to_vertex;
      int * global_cell_ids;
      int * cell_owner;
      int * global_vertex_ids;
      int * vertex_owner;
      double * x_vertices, * y_vertices;
      double * x_cells, * y_cells;
      int * cell_mask;

      yac_read_icon_grid_information_parallel(
        filename, MPI_COMM_WORLD, &nbr_vertices, &nbr_cells, &num_vertices_per_cell,
        &cell_to_vertex, &global_cell_ids, &cell_owner, &global_vertex_ids,
        &vertex_owner, &x_vertices, &y_vertices, &x_cells, &y_cells, &cell_mask);

      if (nbr_cells != ref_nbr_cells[comm_rank])
        PUT_ERR("wrong number of vertices\n");
      if (nbr_vertices != ref_nbr_vertices[comm_rank])
        PUT_ERR("wrong number of vertices\n");
      for (int i = 0; i < nbr_cells; ++i)
        if (num_vertices_per_cell[i] != 3)
          PUT_ERR("wrong number of vertices per cell\n");
      for (int i = 0; i < nbr_cells; ++i) {
        for (int j = 0; j < 3; ++j)
          if (cell_to_vertex[3 * i + j] != ref_cell_to_vertex[i][j])
            PUT_ERR("error in cell_to_vertex\n");
        if (global_cell_ids[i] != ref_global_cell_ids[comm_rank][i])
          PUT_ERR("wrong global cell id\n");
        if (cell_owner[i] != ref_cell_owner[i])
          PUT_ERR("wrong cell owner\n");
        if (double_are_unequal(
              x_cells[i], clon[ref_global_cell_ids[comm_rank][i]] * YAC_RAD))
          PUT_ERR("wrong cell longitude\n");
        if (double_are_unequal(
              y_cells[i], clat[ref_global_cell_ids[comm_rank][i]] * YAC_RAD))
          PUT_ERR("wrong cell latitude\n");
        if (cell_mask[i] != mask[ref_global_cell_ids[comm_rank][i]])
          PUT_ERR("wrong cell mask\n");
      }
      for (int i = 0; i < nbr_vertices; ++i) {
        if (global_vertex_ids[i] != ref_global_vertex_ids[comm_rank][i])
          PUT_ERR("wrong global vertex id\n");
        if (vertex_owner[i] != ref_vertex_owner[i]) PUT_ERR("wrong vertex owner\n");
        if (double_are_unequal(
              x_vertices[i], vlon[ref_global_vertex_ids[comm_rank][i]] * YAC_RAD))
          PUT_ERR("wrong vertex longitude\n");
        if (double_are_unequal(
              y_vertices[i], vlat[ref_global_vertex_ids[comm_rank][i]] * YAC_RAD))
          PUT_ERR("wrong vertex latitude\n");

      }
      free(num_vertices_per_cell);
      free(cell_to_vertex);
      free(cell_owner);
      free(vertex_owner);
      free(global_cell_ids);
      free(global_vertex_ids);
      free(x_vertices);
      free(y_vertices);
      free(x_cells);
      free(y_cells);
      free(cell_mask);
    }

    {
      set_even_io_rank_list(MPI_COMM_WORLD);
      struct yac_basic_grid_data grid_data =
        yac_read_icon_basic_grid_data_parallel(filename, MPI_COMM_WORLD);
      if (grid_data.num_cells != 4)
        PUT_ERR("error in yac_read_icon_basic_grid_data_parallel");
      if (grid_data.num_vertices != 6)
        PUT_ERR("error in yac_read_icon_basic_grid_data_parallel");
      yac_basic_grid_data_free(grid_data);
      clear_yac_io_env();
    }

    {
      struct yac_basic_grid * icon_grid =
        yac_read_icon_basic_grid_parallel(
          filename, "icon_grid", MPI_COMM_WORLD);
      if (yac_basic_grid_get_data_size(icon_grid, YAC_LOC_CELL) != 4)
        PUT_ERR("error in yac_read_icon_basic_grid_parallel");
      if (yac_basic_grid_get_data_size(icon_grid, YAC_LOC_CORNER) != 6)
        PUT_ERR("error in yac_read_icon_basic_grid_parallel");
      yac_basic_grid_delete(icon_grid);
    }

    {
      struct yac_basic_grid * icon_grid;
      size_t cell_coordinates_idx;
      int * cell_mask;
      yac_read_icon_basic_grid_parallel_2(
        filename, "icon_grid", MPI_COMM_WORLD, &icon_grid,
        &cell_coordinates_idx, &cell_mask);
      if (yac_basic_grid_get_data_size(icon_grid, YAC_LOC_CELL) != 4)
        PUT_ERR("error in yac_read_icon_basic_grid_parallel_2");
      if (yac_basic_grid_get_data_size(icon_grid, YAC_LOC_CORNER) != 6)
        PUT_ERR("error in yac_read_icon_basic_grid_parallel_2");
      if (cell_coordinates_idx != 0)
        PUT_ERR("error in yac_read_icon_basic_grid_parallel_2");
      yac_basic_grid_delete(icon_grid);
      free(cell_mask);
    }

    // ensure that all processes finished reading the file
    MPI_Barrier(MPI_COMM_WORLD);
    if (comm_rank == 0) unlink(filename);
  }

  MPI_Finalize();

  return TEST_EXIT_CODE;
}
