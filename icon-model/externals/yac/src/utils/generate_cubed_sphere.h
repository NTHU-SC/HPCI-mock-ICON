// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef GENERATE_CUBED_SPHERE_H
#define GENERATE_CUBED_SPHERE_H

#include "basic_grid.h"

// YAC PUBLIC HEADER START

/**
 * Creates a cubed sphere grid with n subdivisions
 *
 * This routine is based on Matlab code provided by Mike Hobson and written by
 * Mike Rezny (both from MetOffice)
 *
 * @param[in] n number of subdivisions of the cubed sphere grid
 * @param[out] num_cells number of cells in the grid
 * @param[out] num_vertices number of vertices in the grid
 * @param[out] x_vertices x coordinate of vertices
 * @param[out] y_vertices y coordinate of vertices
 * @param[out] z_vertices z coordinate of vertices
 * @param[out] vertices_of_cell vertex indices for each cell
 * @param[out] face_id id for orientation of face w.r.t cube faces (unused)
 * @remark Vertex coordinates are provided in Cartesian coordinates on a unit sphere,
 *         with the center of the sphere being the origin of the coordinate system.
 */
void yac_generate_cubed_sphere_grid_information(
  unsigned n, unsigned * num_cells, unsigned * num_vertices,
  double ** x_vertices, double ** y_vertices, double ** z_vertices,
  unsigned ** vertices_of_cell, unsigned ** face_id);

/**
 * Creates a cubed sphere grid with n subdivisions decomposed for a given number of ranks
 *
 * @param[in] n number of subdivisions of the cubed sphere grid
 * @param[out] nbr_vertices number of vertices in the grid
 * @param[out] nbr_cells number of cells in the grid
 * @param[out] num_vertices_per_cell number of vertices per cell
 * @param[out] cell_to_vertex vertex indices for each cell
 * @param[out] x_vertices longitudes of vertices in radians
 * @param[out] y_vertices latitudes of vertices in radians
 * @param[out] x_cells longitudes of cell center in radians
 * @param[out] y_cells latitudes of cell center in radians
 * @param[out] global_cell_id global cell IDs
 * @param[out] cell_core_mask cell core mask
 * @param[out] global_corner_id global corner IDs
 * @param[out] corner_core_mask corner core mask
 * @param[in] rank id of this rank
 * @param[in] size total number of ranks
 *
 * @see \ref yac_generate_cubed_sphere_grid_information for details
 */
void yac_generate_part_cube_grid_information(
  unsigned n, unsigned * nbr_vertices, unsigned * nbr_cells,
  unsigned ** num_vertices_per_cell, unsigned ** cell_to_vertex,
  double ** x_vertices, double ** y_vertices, double ** x_cells,
  double ** y_cells, int ** global_cell_id, int ** cell_core_mask,
  int ** global_corner_id, int ** corner_core_mask, int rank, int size);

/**
 * Creates a cubed sphere grid with n subdivisions
 * @param[in] n number of subdivisions of the cubed sphere grid
 * @returns yac_basic_grid_data containing data for a cubed sphere grid
 *
 * @see \ref yac_generate_cubed_sphere_grid_information for details
 */
struct yac_basic_grid_data yac_generate_cubed_sphere_grid(unsigned n);

/**
 * Creates a cubed sphere grid with n subdivisions and given name
 * @param[in] name name of the grid
 * @param[in] n number of subdivisions of the cubed sphere grid
 * @returns yac_basic_grid holding a cubed sphere grid
 *
 * @see \ref yac_generate_cubed_sphere_grid_information for details
 */
struct yac_basic_grid * yac_generate_cubed_sphere_basic_grid(
  char const * name, size_t n);

/**
 * Creates a cubed sphere grid with n subdivisions and writes it to a netCDF file
 * using an ICON-grid-file like layout.
 * @param[in] n number of subdivisions of the cubed sphere grid
 * @param[in] filename file name
 *
 * @see \ref yac_generate_cubed_sphere_grid_information for details
 */
void yac_write_cubed_sphere_grid(unsigned n, char const * filename);

// YAC PUBLIC HEADER STOP

#endif // GENERATE_CUBED_SPHERE_H

