// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <netcdf.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "tests.h"
#include "test_common.h"
#include "grid_file_common.h"
#include "read_scrip_grid.h"
#include "grid2vtk.h"
#include "io_utils.h"

/** \file test_read_scrip_c.c
 *  \test
 * This contains examples for yac_read_scrip_basic_grid
 */

// #define HAVE_OASIS_FILES
// #define WRITE_VTK_GRID_FILE
#define DUMMY_GRID_NAME ("dummy_grid")

int main(void) {

  char * grid_filename = "test_read_scrip_c_grids.nc";
  char * mask_filename = "test_read_scrip_c_masks.nc";
  write_dummy_scrip_grid_file(
    DUMMY_GRID_NAME, grid_filename, mask_filename, 1,
    360, 10, (double[]){0.0,360.0}, (double[]){0.0, 10.0});

  size_t ref_num_cells = 360 * 10;
  size_t ref_num_vertices = 360 * 11;

  int valid_mask_value = 0;

  if (!yac_file_exists(grid_filename)) return EXIT_FAILURE;
  if (!yac_file_exists(mask_filename)) return EXIT_FAILURE;

  { // read unstructured grid

    char * gridname = DUMMY_GRID_NAME;

    { // test yac_read_scrip_basic_grid
      size_t * duplicated_cell_idx = NULL;
      yac_int * orig_cell_global_ids = NULL;
      size_t nbr_duplicated_cells = 0;

      size_t cell_coord_idx;

      struct yac_basic_grid * scrip_grid =
        yac_read_scrip_basic_grid(
          grid_filename, mask_filename, gridname, valid_mask_value, gridname,
          0, &cell_coord_idx, &duplicated_cell_idx, &orig_cell_global_ids,
          &nbr_duplicated_cells);

      if (ref_num_vertices !=
          yac_basic_grid_get_data_size(scrip_grid, YAC_LOC_CORNER))
        PUT_ERR("error in yac_read_scrip_grid");
      if (ref_num_cells !=
          yac_basic_grid_get_data_size(scrip_grid, YAC_LOC_CELL))
        PUT_ERR("error in yac_read_scrip_grid");

      free(duplicated_cell_idx);
      free(orig_cell_global_ids);

  #ifdef WRITE_VTK_GRID_FILE
      yac_write_basic_grid_to_file(scrip_grid, gridname);
  #else
      char const * vtk_gridname = "test_read_scrip_grid_c";
      char const * vtk_filename = "test_read_scrip_grid_c.vtk";

      yac_write_basic_grid_to_file(scrip_grid, vtk_gridname);

      if (!yac_file_exists(vtk_filename))
        PUT_ERR("ERROR in yac_write_basic_grid_to_file");

      unlink(vtk_filename);
  #endif // WRITE_VTK_GRID_FILE

      yac_basic_grid_delete(scrip_grid);
    }

    { // test yac_read_scrip_grid_information

      size_t num_vertices;
      size_t num_cells;
      int * num_vertices_per_cell;
      double * x_vertices;
      double * y_vertices;
      double * x_cells;
      double * y_cells;
      int * cell_to_vertex;
      int * cell_core_mask;
      size_t * duplicated_cell_idx;
      size_t * orig_cell_idx;
      size_t nbr_duplicated_cells;

      yac_read_scrip_grid_information(
        grid_filename, mask_filename, gridname, valid_mask_value,
        &num_vertices, &num_cells, &num_vertices_per_cell,
        &x_vertices, &y_vertices, &x_cells, &y_cells,
        &cell_to_vertex, &cell_core_mask,
        &duplicated_cell_idx, &orig_cell_idx, &nbr_duplicated_cells);

      if (ref_num_vertices != num_vertices)
        PUT_ERR("error in yac_read_scrip_grid_information");
      if (ref_num_cells != num_cells)
        PUT_ERR("error in yac_read_scrip_grid_information");

      free(num_vertices_per_cell);
      free(x_vertices);
      free(y_vertices);
      free(x_cells);
      free(y_cells);
      free(cell_to_vertex);
      free(cell_core_mask);
      free(duplicated_cell_idx);
      free(orig_cell_idx);
    }
  }

  { // read lon/lat grid
    char * gridname = DUMMY_GRID_NAME;

    struct yac_basic_grid * scrip_grid =
      yac_read_scrip_basic_grid(
        grid_filename, mask_filename, gridname, valid_mask_value, gridname,
        1, NULL, NULL, NULL, NULL);

#ifdef WRITE_VTK_GRID_FILE
    char vtk_gridname[64];
    yac_write_basic_grid_to_file(
      scrip_grid, strcat(strcpy(vtk_gridname, gridname), "_ll"));
#endif // WRITE_VTK_GRID_FILE

    yac_basic_grid_delete(scrip_grid);
  }

  { // read cloud grid (containing only cell centers)

    // create grid files
    char * grid_filename = "test_read_scrip_cloud_c_grids.nc";
    char * mask_filename = "test_read_scrip_cloud_c_masks.nc";
    write_dummy_scrip_grid_file(
      DUMMY_GRID_NAME, grid_filename, mask_filename, 0,
      380, 180, (double[]){0.0,380.0}, (double[]){-90.0, 90.0});

    size_t vertex_coord_idx;
    size_t * duplicated_vertex_idx;
    yac_int * orig_vertex_global_ids;
    size_t nbr_duplicated_vertices;
    struct yac_basic_grid * scrip_cloud_grid =
      yac_read_scrip_cloud_basic_grid(
        grid_filename, mask_filename, DUMMY_GRID_NAME, 0,
        DUMMY_GRID_NAME, &vertex_coord_idx, &duplicated_vertex_idx,
        &orig_vertex_global_ids, &nbr_duplicated_vertices);

#ifdef WRITE_VTK_GRID_FILE
    yac_write_basic_grid_to_file(scrip_cloud_grid, DUMMY_GRID_NAME);
#endif // WRITE_VTK_GRID_FILE

    if (strcmp(yac_basic_grid_get_name(scrip_cloud_grid), DUMMY_GRID_NAME))
      PUT_ERR("wrong grid name");

    {
      struct yac_basic_grid_data * scrip_cloud_grid_data =
        yac_basic_grid_get_data(scrip_cloud_grid);

      if (scrip_cloud_grid_data->num_cells != 0)
        PUT_ERR("wrong number of cells");
      if (scrip_cloud_grid_data->num_vertices != 380 * 180)
        PUT_ERR("wrong global number of vertices");
      if (nbr_duplicated_vertices != 20 * 180)
        PUT_ERR("wrong global number of duplicated vertices");
      if (scrip_cloud_grid_data->num_edges != 0)
        PUT_ERR("wrong number of edges");
    }

    {
      struct yac_basic_grid_data scrip_cloud_grid_data =
        yac_read_scrip_cloud_basic_grid_data(
          grid_filename, mask_filename, DUMMY_GRID_NAME, 0);

      if (scrip_cloud_grid_data.num_cells != 0)
        PUT_ERR("wrong number of cells");
      if (scrip_cloud_grid_data.num_vertices != 380 * 180)
        PUT_ERR("wrong global number of vertices");
      if (scrip_cloud_grid_data.num_edges != 0)
        PUT_ERR("wrong number of edges");

      yac_basic_grid_data_free(scrip_cloud_grid_data);
    }

    free(duplicated_vertex_idx);
    free(orig_vertex_global_ids);
    yac_basic_grid_delete(scrip_cloud_grid);

    // delete grid files
    unlink(grid_filename);
    unlink(mask_filename);
  }

  unlink(grid_filename);
  unlink(mask_filename);

  return TEST_EXIT_CODE;
}
