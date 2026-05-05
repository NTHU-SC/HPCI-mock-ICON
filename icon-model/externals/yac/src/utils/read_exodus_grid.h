// Copyright (c) 2025 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include "basic_grid.h"

// YAC PUBLIC HEADER START

#include <mpi.h>

/**
 * Reads in an EXODUS-formated grid netcdf file and returns the grid information
 * in a format that allows the user to set up a \ref yac_basic_grid_data structure.
 * The reading is done in parallel and a basic domain decomposition is applied.
 *
 * @param[in]  filename           name of the EXODUS-formated grid netcdf file
 * @param[in]  comm               MPI communicator containing all proceses
 *                                that will get a part of the grid
 * @param[out] node_coords        Cartesian coordinates of all nodes/vertices
 * @param[out] elem_ids           global ids of local elements/cells
 * @param[out] node_ids           global ids of local nodes/vertices
 * @param[out] num_elem           number of elements/cells in the local part of
 *                                the grid
 * @param[out] num_nodes          number of nodes/vertices in the local part of the
 *                                the grid
 * @param[out] num_nodes_per_elem number of nodes/vertices per element/cell
 * @param[out] num_elem_per_node  number of elements/cells per node/vertex
 * @param[out] elem_to_node       local node/vertex indices for each local
 *                                element/cell
 * @param[out] node_to_elem       local element/cell indices for each local
 *                                node/vertex
 * @see \ref io_config_detail
 */
void yac_read_exodus_grid_information_parallel(
  const char * filename, MPI_Comm comm, yac_coordinate_pointer * node_coords,
  yac_int ** elem_ids, yac_int ** node_ids,
  size_t * num_elem, size_t * num_nodes,
  int ** num_nodes_per_elem, int ** num_elem_per_node,
  size_t ** elem_to_node, size_t ** node_to_elem);

/**
 * Reads in an EXODUS-formated grid netcdf file in parallel and returns a
 * \ref yac_basic_grid_data built from it.
 *
 * @param[in] filename     name of the EXODUS-formated grid netcdf file
 * @param[in] use_ll_edges if != 0, assume that all edges of the grid follow
 *                         circles of either constant longitude or constant
 *                         latitude and set edge types accordingly. Otherwise
 *                         it is assumed that all edges follow great circles.
 * @param[in] comm         MPI communicator containing all proceses
 *                         that will get a part of the grid
 * @returns yac_basic_grid_data structure containing part of the grid.
 * @see \ref io_config_detail
 */
struct yac_basic_grid_data yac_read_exodus_basic_grid_data_parallel(
  const char * filename, int use_ll_edges, MPI_Comm comm);

/**
 * Reads in an EXODUS-formated grid netcdf file in parallel and returns a
 * \ref yac_basic_grid built from it.
 *
 * @param[in] filename     name of the EXODUS-formated grid netcdf file
 * @param[in] gridname     name of grid
 * @param[in] use_ll_edges if != 0, assume that all edges of the grid follow
 *                         circles of either constant longitude or constant
 *                         latitude and set edge types accordingly. Otherwise
 *                         it is assumed that all edges follow great circles.
 * @param[in] comm         MPI communicator containing all proceses
 *                         that will get a part of the grid
 * @returns yac_basic_grid structure containing part of the grid.
 * @remark The data of each process contains no halo cells.
 * @see \ref io_config_detail
 */
struct yac_basic_grid * yac_read_exodus_basic_grid_parallel(
  char const * filename, char const * gridname, int use_ll_edges,
  MPI_Comm comm);

// YAC PUBLIC HEADER STOP
