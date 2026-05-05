// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <netcdf.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <mpi.h>

#include "tests.h"
#include "test_common.h"
#include "grid_file_common.h"
#include "read_scrip_grid.h"
#include "grid2vtk.h"
#include "io_utils.h"

/** \file test_read_scrip_parallel_c.c
 *  \test
 * This contains examples for yac_read_scrip_basic_grid_parallel
 */

// #define HAVE_OASIS_FILES
// #define WRITE_VTK_GRID_FILE
#define DUMMY_GRID_NAME ("dummy_grid")

int main(void) {

  MPI_Init(NULL, NULL);

  int comm_size, comm_rank;

  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  if ((comm_size != 4) && (comm_rank == 0)) {
    fputs("wrong number of processes (has to be four)\n", stderr);
    exit(EXIT_FAILURE);
  }

#ifdef HAVE_OASIS_FILES
  char * grid_filename = "../examples/OASIS-grid/grids.nc";
  char * mask_filename = "../examples/OASIS-grid/masks_no_atm.nc";
#else
  char * grid_filename = "test_read_scrip_parallel_c_grids.nc";
  char * mask_filename = "test_read_scrip_parallel_c_masks.nc";
  if (comm_rank == 0)
    write_dummy_scrip_grid_file(
      DUMMY_GRID_NAME, grid_filename, mask_filename, 1,
      360, 180, (double[]){0.0,360.0}, (double[]){-90.0, 90.0});
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  int valid_mask_value = 0;

  if (!yac_file_exists(grid_filename)) return EXIT_FAILURE;
  if (!yac_file_exists(mask_filename)) return EXIT_FAILURE;

  { // read unstructured grid
#ifdef HAVE_OASIS_FILES
    // char * gridname = "bggd";
    // char * gridname = "icoh";
    // char * gridname = "icos";
    // char * gridname = "nogt";
    char * gridname = "sse7";
    // char * gridname = "torc";
    // char * gridname = "ssea";
#else
    char * gridname = DUMMY_GRID_NAME;
#endif

    size_t * duplicated_cell_idx = NULL;
    yac_int * orig_cell_global_ids = NULL;
    size_t nbr_duplicated_cells = 0;

    size_t cell_coord_idx;

    set_even_io_rank_list(MPI_COMM_WORLD);

    for (int reader_idx = 0; reader_idx < 2; ++reader_idx) {

      struct yac_basic_grid * scrip_grid;

      if (reader_idx == 0)
        scrip_grid =
          yac_read_scrip_basic_grid_parallel(
            grid_filename, mask_filename, MPI_COMM_WORLD,
            gridname, valid_mask_value, gridname,
            0, &cell_coord_idx, &duplicated_cell_idx,
            &orig_cell_global_ids, &nbr_duplicated_cells);
      else {
        int point_location;
        scrip_grid =
          yac_read_scrip_generic_basic_grid_parallel(
            grid_filename, mask_filename, MPI_COMM_WORLD,
            gridname, valid_mask_value, gridname,
            0, &cell_coord_idx, &duplicated_cell_idx,
            &orig_cell_global_ids, &nbr_duplicated_cells, &point_location);
        if (point_location != YAC_LOC_CELL)
          PUT_ERR("error in yac_read_scrip_generic_basic_grid_parallel")
      }

      free(duplicated_cell_idx);
      free(orig_cell_global_ids);

#ifdef WRITE_VTK_GRID_FILE
      char vtk_gridname[64];
      snprintf(vtk_gridname, sizeof(vtk_gridname), "%s_%d", gridname, comm_rank);
      yac_write_basic_grid_to_file(scrip_grid, vtk_gridname);
#endif // WRITE_VTK_GRID_FILE

      yac_basic_grid_delete(scrip_grid);
    }
    clear_yac_io_env();
  }

  { // read lon/lat grid
#ifdef HAVE_OASIS_FILES
    // char * gridname = "bggd";
    // char * gridname = "icoh";
    // char * gridname = "icos";
    // char * gridname = "nogt";
    // char * gridname = "sse7";
    char * gridname = "torc";
    // char * gridname = "ssea";
#else
    char * gridname = DUMMY_GRID_NAME;
#endif

    set_even_io_rank_list(MPI_COMM_WORLD);
    struct yac_basic_grid * scrip_grid =
      yac_read_scrip_basic_grid_parallel(
        grid_filename, mask_filename, MPI_COMM_WORLD,
        gridname, valid_mask_value, gridname,
        1, NULL, NULL, NULL, NULL);

#ifdef WRITE_VTK_GRID_FILE
    char vtk_gridname[64];
    snprintf(
      vtk_gridname, sizeof(vtk_gridname), "%s_ll_%d", gridname, comm_rank);
    yac_write_basic_grid_to_file(scrip_grid, vtk_gridname);
#endif // WRITE_VTK_GRID_FILE

    yac_basic_grid_delete(scrip_grid);
    clear_yac_io_env();
  }

  { // read cloud grid (containing only cell centers)

    enum {NUM_LAT = 181, NUM_LON = 380};

    // create grid files
    char * grid_filename = "test_read_scrip_cloud_parallel_c_grids.nc";
    char * mask_filename = "test_read_scrip_cloud_parallel_c_masks.nc";
    if (comm_rank == 0)
      write_dummy_scrip_grid_file(
        DUMMY_GRID_NAME, grid_filename, mask_filename, 0,
        NUM_LON, NUM_LAT, (double[]){0.0, 380.0}, (double[]){-90.0, 90.0});
    MPI_Barrier(MPI_COMM_WORLD);

    setenv("YAC_IO_RANK_LIST", "0,1,2,3", 1);
    setenv("YAC_IO_MAX_NUM_RANKS_PER_NODE", "4", 1);

    for (int reader_idx = 0; reader_idx < 2; ++reader_idx) {

      size_t vertex_coord_idx;
      size_t * duplicated_vertex_idx;
      yac_int * orig_vertex_global_ids;
      size_t nbr_duplicated_vertices;
      struct yac_basic_grid * scrip_cloud_grid;

      if (reader_idx == 0)
        scrip_cloud_grid =
          yac_read_scrip_cloud_basic_grid_parallel(
            grid_filename, mask_filename, MPI_COMM_WORLD, DUMMY_GRID_NAME, 0,
            DUMMY_GRID_NAME, &vertex_coord_idx, &duplicated_vertex_idx,
            &orig_vertex_global_ids, &nbr_duplicated_vertices);
      else {
        int point_location;
        scrip_cloud_grid =
          yac_read_scrip_generic_basic_grid_parallel(
            grid_filename, mask_filename, MPI_COMM_WORLD, DUMMY_GRID_NAME, 0,
            DUMMY_GRID_NAME, 0, &vertex_coord_idx, &duplicated_vertex_idx,
            &orig_vertex_global_ids, &nbr_duplicated_vertices, &point_location);
        if (point_location != YAC_LOC_CORNER)
          PUT_ERR("error in yac_read_scrip_generic_basic_grid_parallel");
      }

#ifdef WRITE_VTK_GRID_FILE
      char vtk_gridname[64];
      snprintf(
        vtk_gridname, sizeof(vtk_gridname), "%s_%d", "scrip_cloud", comm_rank);
      yac_write_basic_grid_to_file(scrip_cloud_grid, vtk_gridname);
#endif // WRITE_VTK_GRID_FILE

      if (strcmp(yac_basic_grid_get_name(scrip_cloud_grid), DUMMY_GRID_NAME))
        PUT_ERR("wrong grid name");

      struct yac_basic_grid_data * scrip_cloud_grid_data =
        yac_basic_grid_get_data(scrip_cloud_grid);

      if (scrip_cloud_grid_data->num_cells != 0)
        PUT_ERR("wrong number of cells");

      int local_num_vertices = (int)scrip_cloud_grid_data->num_vertices;
      int local_num_unique_vertices =
        local_num_vertices - (int)nbr_duplicated_vertices;
      int global_num_vertices, global_num_unique_vertices;
      MPI_Allreduce(
        &local_num_vertices, &global_num_vertices, 1,
        MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(
        &local_num_unique_vertices, &global_num_unique_vertices, 1,
        MPI_INT, MPI_SUM, MPI_COMM_WORLD);
      if (global_num_vertices != NUM_LAT * NUM_LON)
        PUT_ERR("wrong global number of vertices");
      if (global_num_unique_vertices != 181 * 360)
        PUT_ERR("wrong global number of unique vertices");

      if (scrip_cloud_grid_data->num_edges != 0)
        PUT_ERR("wrong number of edges");

      free(duplicated_vertex_idx);
      free(orig_vertex_global_ids);
      yac_basic_grid_delete(scrip_cloud_grid);
    }

    clear_yac_io_env();

    // delete grid files
    MPI_Barrier(MPI_COMM_WORLD);
    if (comm_rank == 0) {
      unlink(grid_filename);
      unlink(mask_filename);
    }
  }

#ifndef HAVE_OASIS_FILES
  MPI_Barrier(MPI_COMM_WORLD);
  if (comm_rank == 0) {
    unlink(grid_filename);
    unlink(mask_filename);
  }
#endif

  MPI_Finalize();

  return TEST_EXIT_CODE;
}

