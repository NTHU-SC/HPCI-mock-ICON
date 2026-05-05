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
#include "read_exodus_grid.h"
#include "grid_file_common.h"

#include <netcdf.h>

#include <mpi.h>

/** \file test_read_exodus_parallel.c
 *  \test
 * This contains examples for reading of exodus-formated netcdf grid files.
 */

int main(void) {

  MPI_Init(NULL, NULL);
  xt_initialize(MPI_COMM_WORLD);

  int comm_size, comm_rank;

  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);

  if ((comm_size != 4) && (comm_rank == 0)) {
    fputs("wrong number of processes (has to be four)\n", stderr);
    exit(EXIT_FAILURE);
  }

  {
    char const filename[] = "test_read_exodus_parallel_grid.nc";
    size_t num_lon = 271, num_lat = 181;
    double lon_range[] = {-135.0, 135.0}, lat_range[] = {-90.0, 90.0};
    if (comm_rank == 0)
      write_dummy_exodus_grid_file(filename, num_lon, num_lat, lon_range, lat_range);

    // ensure that the grid file exists
    MPI_Barrier(MPI_COMM_WORLD);

    {
      int use_ll_edges = 0;
      struct yac_basic_grid * exodus_grid =
        yac_read_exodus_basic_grid_parallel(
          filename, "exodus", use_ll_edges, MPI_COMM_WORLD);

      size_t local_size =
        yac_basic_grid_get_data_size(exodus_grid, YAC_LOC_CELL);

      if (local_size != ((num_lon - 1) * (num_lat - 1)) / 4)
        PUT_ERR("ERROR in yac_read_exodus_basic_grid_parallel");

      yac_basic_grid_delete(exodus_grid);
    }

    // ensure that all tests are completed
    MPI_Barrier(MPI_COMM_WORLD);

    if (comm_rank == 0) unlink(filename);
  }

  {
    char const filename[] = "test_read_exodus_parallel_grid.nc";
    size_t num_lon = 271, num_lat = 181;
    double lon_range[] = {0.0, 10.0}, lat_range[] = {80.0, 90.0};
    if (comm_rank == 0)
      write_dummy_exodus_grid_file(filename, num_lon, num_lat, lon_range, lat_range);

    // ensure that the grid file exists
    MPI_Barrier(MPI_COMM_WORLD);

    {
      int use_ll_edges = 1;
      struct yac_basic_grid * exodus_grid =
        yac_read_exodus_basic_grid_parallel(
          filename, "exodus", use_ll_edges, MPI_COMM_WORLD);

      size_t local_size =
        yac_basic_grid_get_data_size(exodus_grid, YAC_LOC_CELL);

      if (local_size != ((num_lon - 1) * (num_lat - 1)) / 4)
        PUT_ERR("ERROR in yac_read_exodus_basic_grid_parallel");

      yac_basic_grid_delete(exodus_grid);
    }

    // ensure that all tests are completed
    MPI_Barrier(MPI_COMM_WORLD);

    if (comm_rank == 0) unlink(filename);
  }

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}
