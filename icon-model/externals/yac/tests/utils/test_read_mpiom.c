// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <string.h>

#include "tests.h"
#include "test_common.h"
#include "read_mpiom_grid.h"
#include "grid2vtk.h"
#include "io_utils.h"

/** \file test_read_mpiom.c
 *  \test
 * This contains examples for yac_read_mpiom_basic_grid_data.
 */

int main(int argc, char** argv) {

  if (argc != 2) {
    PUT_ERR("ERROR: missing grid file directory");
    return TEST_EXIT_CODE;
  }

  char * grid_filename =
    strcat(
      strcpy(
        malloc(strlen(argv[1]) + 32), argv[1]), "GR30_lsm.nc");

  if (!yac_file_exists(grid_filename)) return EXIT_SKIP_TEST;

  {
    struct yac_basic_grid_data mpiom_grid =
      yac_read_mpiom_basic_grid_data(grid_filename);

    if (mpiom_grid.num_cells != 12120)
      PUT_ERR("wrong number of grid cells");

// #define WRITE_VTK_GRID_FILE
#ifdef WRITE_VTK_GRID_FILE
    yac_write_basic_grid_data_to_file(&mpiom_grid, "mpiom");
#endif // WRITE_VTK_GRID_FILE

    yac_basic_grid_data_free(mpiom_grid);
  }

  {
    struct yac_basic_grid * mpiom_grid =
      yac_read_mpiom_basic_grid(grid_filename, "mpiom");


    if (yac_basic_grid_get_data_size(mpiom_grid, YAC_LOC_CELL) != 12120)
      PUT_ERR("wrong number of grid cells");

    yac_basic_grid_delete(mpiom_grid);
  }

  free(grid_filename);

  return TEST_EXIT_CODE;
}

