// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include "basic_grid.h"

// YAC PUBLIC HEADER START

#include <mpi.h>

/**
 * Reads in an icon grid netcdf file and generates a
 * \ref yac_basic_grid_data from it.
 * @param[in] filename name of the icon grid netcdf file
 * @returns yac_basic_grid_data that contains the icon grid
 */
struct yac_basic_grid_data yac_read_icon_basic_grid_data(
  char const * filename);

/**
 * Reads in an icon grid netcdf file and generates a \ref yac_basic_grid from it.
 * @param[in] filename name of the icon grid netcdf file
 * @param[in] gridname name of the grid
 * @returns yac_basic_grid structure that contains the icon grid
 */
struct yac_basic_grid * yac_read_icon_basic_grid(
  char const * filename, char const * gridname);

/**
 * Reads in an icon grid netcdf file and return the grid information in
 * a format that is supported by the YAC user interface.
 *
 * @param[in]  filename              name of the icon grid netcdf file
 * @param[out] num_vertices          number of vertices in the grid
 * @param[out] num_cells             number of cells in the grid
 * @param[out] num_vertices_per_cell number of vertices per cell
 * @param[out] cell_to_vertex        vertex indices for each cell
 * @param[out] x_vertices            longitudes of vertices
 * @param[out] y_vertices            latitudes of vertices
 * @param[out] x_cells               longitudes of cell center
 * @param[out] y_cells               latitudes of cell center
 * @param[out] cell_mask             mask for cells
 */
void yac_read_icon_grid_information(const char * filename, int * num_vertices,
                                    int * num_cells, int ** num_vertices_per_cell,
                                    int ** cell_to_vertex, double ** x_vertices,
                                    double ** y_vertices, double ** x_cells,
                                    double ** y_cells, int ** cell_mask);

/**
 * Reads in an icon grid netcdf file and return the grid information in
 * a format that is supported by the YAC user interface.
 *
 * @param[in]  filename              name of the icon grid netcdf file
 * @param[out] num_vertices          number of vertices in the grid
 * @param[out] num_cells             number of cells in the grid
 * @param[out] num_vertices_per_cell number of vertices per cell
 * @param[out] cell_to_vertex        vertex indices for each cell
 * @param[out] x_vertices            longitudes of vertices
 * @param[out] y_vertices            latitudes of vertices
 * @param[out] x_cells               longitudes of cell center
 * @param[out] y_cells               latitudes of cell center
 * @param[out] global_cell_id        global cell IDs
 * @param[out] cell_mask             mask for cells
 * @param[out] cell_core_mask        cell core mask
 * @param[out] global_corner_id      global corner IDs
 * @param[out] corner_core_mask      corner core mask
 * @param[out] rank                  local MPI rank
 * @param[out] size                  number of MPI ranks
 * @remark The process reads in the fields of the whole grid and then extracts
 *         parts of it based on the "rank" and "size" parameter.
 */
void yac_read_part_icon_grid_information(
  const char * filename, int * num_vertices, int * num_cells,
  int ** num_vertices_per_cell, int ** cell_to_vertex,
  double ** x_vertices, double ** y_vertices,
  double ** x_cells, double ** y_cells, int ** global_cell_id,
  int ** cell_mask, int ** cell_core_mask,
  int ** global_corner_id, int ** corner_core_mask, int rank, int size);

/**
 * Reads in an icon grid netcdf file and returns the grid information in
 * a format that is supported by the YAC user interface. The reading is done
 * in parallel and a basic domain decomposition is applied.
 *
 * @param[in]  filename              name of the icon grid netcdf file
 * @param[in]  comm                  MPI communicator containing all proceses
 *                                   that will get a part of the grid
 * @param[out] num_vertices          number of vertices in the local part of the
 *                                   grid
 * @param[out] num_cells             number of cells in the local part of the
 *                                   grid
 * @param[out] num_vertices_per_cell number of vertices per cell
 * @param[out] cell_to_vertex        vertex indices for each cell
 * @param[out] global_cell_ids       global ids of local cells (core and halo
 *                                   cells)
 * @param[out] cell_owner            owner of each cell (locally owned cells
 *                                   are marked with -1)
 * @param[out] global_vertex_ids     global ids of local vertices (core and halo
 *                                   vertices)
 * @param[out] vertex_owner          owner of each vertex (locally owned
 *                                   vertices are marked with -1)
 * @param[out] x_vertices            longitudes of vertices
 * @param[out] y_vertices            latitudes of vertices
 * @param[out] x_cells               longitudes of cell center
 * @param[out] y_cells               latitudes of cell center
 * @param[out] cell_mask             mask for cells
 *                                   (variable `cell_sea_land_mask` in
 *                                    icon grid file)
 * @remark The data of each process contains one layer of halo cells.
 * @remark If `(x_cells == NULL) && (y_cells == NULL)` no cell center
 *         coordinates will be read.
 * @remark If `(cell_mask == NULL)` no cell mask information will be read.
 * @see \ref io_config_detail
 */
void yac_read_icon_grid_information_parallel(
  const char * filename, MPI_Comm comm, int * num_vertices, int * num_cells,
  int ** num_vertices_per_cell, int ** cell_to_vertex, int ** global_cell_ids,
  int ** cell_owner, int ** global_vertex_ids, int ** vertex_owner,
  double ** x_vertices, double ** y_vertices, double ** x_cells,
  double ** y_cells, int ** cell_mask);

/**
 * Reads in an icon grid netcdf file and returns the grid information in
 * a format that allows the user to set up a \ref yac_basic_grid_data structure.
 * The reading is done in parallel and a basic domain decomposition is applied.
 *
 * @param[in]  filename              name of the icon grid netcdf file
 * @param[in]  comm                  MPI communicator containing all proceses
 *                                   that will get a part of the grid
 * @param[out] x_vertices            longitudes of vertices
 * @param[out] y_vertices            latitudes of vertices
 * @param[out] cell_ids              global ids of local cells
 * @param[out] vertex_ids            global ids of local vertices
 * @param[out] edge_ids              global ids of local edges
 * @param[out] num_cells             number of cells in the local part of the
 *                                   grid
 * @param[out] num_vertices          number of vertices in the local part of the
 *                                   grid
 * @param[out] num_edges             number of edges in the local part of the
 *                                   grid
 * @param[out] num_vertices_per_cell number of vertices per cell
 * @param[out] num_cells_per_vertex  number of cells per vertex
 * @param[out] cell_to_vertex        local vertex indices for each local cell
 * @param[out] cell_to_edge          local edge indices for each local cell
 * @param[out] vertex_to_cell        local cell indices for each local vertex
 * @param[out] edge_to_vertex        local vertex indices for each local edge
 * @param[out] edge_type             type of each local edge
 * @param[out] x_cells               longitudes of local cell centers
 * @param[out] y_cells               latitudes of local cell centers
 * @param[out] cell_mask             mask for local cells
 *                                   (variable `cell_sea_land_mask` in
 *                                    icon grid file)
 * @remark The data of each process contains no halo cells.
 * @remark If `(x_cells == NULL) && (y_cells == NULL)` no cell center
 *         coordinates will be read.
 * @remark If `(cell_mask == NULL)` no cell mask information will be read.
 * @see \ref io_config_detail
 */
void yac_read_icon_grid_information_parallel_2(
  const char * filename, MPI_Comm comm,
  double ** x_vertices, double ** y_vertices,
  yac_int ** cell_ids, yac_int ** vertex_ids, yac_int ** edge_ids,
  size_t * num_cells, size_t * num_vertices, size_t * num_edges,
  int ** num_vertices_per_cell, int ** num_cells_per_vertex,
  size_t ** cell_to_vertex, size_t ** cell_to_edge, size_t ** vertex_to_cell,
  size_t ** edge_to_vertex, enum yac_edge_type ** edge_type,
  double ** x_cells, double ** y_cells, int ** cell_mask);

/**
 * Reads in an icon grid netcdf file in parallel and returns a
 * \ref yac_basic_grid_data built from it.
 *
 * @param[in] filename name of the icon grid netcdf file
 * @param[in] comm     MPI communicator containing all proceses
 *                     that will get a part of the grid
 * @returns yac_basic_grid_data structure containing part of the icon grid.
 * @remark The data of each process contains no halo cells.
 * @see \ref io_config_detail
 */
struct yac_basic_grid_data yac_read_icon_basic_grid_data_parallel(
  const char * filename, MPI_Comm comm);

/**
 * Reads in an icon grid netcdf file in parallel and returns a
 * \ref yac_basic_grid built from it.
 *
 * @param[in] filename name of the icon grid netcdf file
 * @param[in] gridname name of grid
 * @param[in] comm     MPI communicator containing all proceses
 *                     that will get a part of the grid
 * @returns yac_basic_grid structure containing part of the icon grid.
 * @remark The data of each process contains no halo cells.
 * @see \ref io_config_detail
 */
struct yac_basic_grid * yac_read_icon_basic_grid_parallel(
  char const * filename, char const * gridname, MPI_Comm comm);

/**
 * Reads in an icon grid netcdf file in parallel  and returns a
 * \ref yac_basic_grid built from it.
 *
 * @param[in]  filename            name of the icon grid netcdf file
 * @param[in]  gridname            name of grid
 * @param[in]  comm                MPI communicator containing all proceses
 *                                 that will get a part of the grid
 * @param[out] grid                \ref yac_basic_grid structure containing
 *                                 part of the icon grid.
 * @param[out] cell_coordinate_idx points index at which the cell centers are
 *                                 registerd in the grid
 * @param[out] cell_mask           mask for local cells
 *                                 (variable `cell_sea_land_mask` in
 *                                  icon grid file)
 * @remark If `(cell_coordinate_idx == NULL)` no cell center coordinates
 *         will be read.
 * @remark If `(cell_mask == NULL)` no cell mask information will be read.
 * @remark The data of each process contains no halo cells.
 * @see \ref io_config_detail
 */
void yac_read_icon_basic_grid_parallel_2(
  char const * filename, char const * gridname, MPI_Comm comm,
  struct yac_basic_grid ** grid, size_t * cell_coordinate_idx,
  int ** cell_mask);

/**
 * destroys remaining icon grid data
 *
 * @param[out] cell_mask             mask for cells
 * @param[out] global_cell_id        global cell IDs
 * @param[out] cell_core_mask        cell core mask
 * @param[out] num_vertices_per_cell number of vertices per cell
 * @param[out] global_corner_id      global corner IDs
 * @param[out] corner_core_mask      corner core mask
 * @param[out] cell_to_vertex        vertex indices for each cell
 * @param[out] x_cells               longitudes of cell center
 * @param[out] y_cells               latitudes of cell center
 * @param[out] x_vertices            longitudes of vertices
 * @param[out] y_vertices            latitudes of vertices
 */

void yac_delete_icon_grid_data( int ** cell_mask,
                                int ** global_cell_id,
                                int ** cell_core_mask,
                                int ** num_vertices_per_cell,
                                int ** global_corner_id,
                                int ** corner_core_mask,
                                int ** cell_to_vertex,
                                double ** x_cells,
                                double ** y_cells,
                                double ** x_vertices,
                                double ** y_vertices);

// YAC PUBLIC HEADER STOP
