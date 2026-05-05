// Copyright (c) 2025 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef POINT_SELECTION_H
#define POINT_SELECTION_H

#include <stdint.h>
#include <mpi.h>

#include "yac_types.h"

// YAC PUBLIC HEADER START

/**
 * Available point selection types
 */
enum yac_point_selection_type {
  YAC_POINT_SELECTION_TYPE_EMPTY,      //!< empty selection
                                       //!< (no points will match)
  YAC_POINT_SELECTION_TYPE_BND_CIRCLE, //!< point selection based on a
                                       //!< bounding circle
};

struct yac_point_selection;

/**
 * Generates a point selection based on a bounding circle
 * @param[in] center_lon longitude coordinate of the center of the
 *                       bounding circle (in radians)
 * @param[in] center_lat latitude coordinate of the center of the
 *                       bounding circle (in radians)
 * @param[in] inc_angle  the angle between center vector and a vector pointing
 *                       to any point on the bounding circle (in radians)
 * @returns point selection
 *
 * @remark use yac_point_selection_delete to delete yac_point_selection
 */
struct yac_point_selection * yac_point_selection_bnd_circle_new(
  double center_lon, double center_lat, double inc_angle);

/**
 * Deletes a given point selection
 */
void yac_point_selection_delete(struct yac_point_selection * point_select);

// YAC PUBLIC HEADER STOP

/**
 * Copies a given point selection
 * @param[in] point_select point selection to be copied
 * @returns copy of the provided point selection
 * @remark use yac_point_selection_delete to delete yac_point_selection
 */
struct yac_point_selection * yac_point_selection_copy(
  struct yac_point_selection const * point_select);

/**
 * Computes the minimum size required for packing a give point selection
 * using `MPI_Pack`
 * @param[in] point_select point selection
 * @param[in] comm         MPI communicator going to be used for the packing
 * @returns minimum packing size
 */
size_t yac_point_selection_get_pack_size(
  struct yac_point_selection const * point_select, MPI_Comm comm);

/**
 * Packs a given point selection into a buffer using `MPI_Pack`
 * @param[in]     point_select point selection
 * @param[in,out] buffer       packing buffer
 * @param[in]     buffer_size  size of packing buffer
 * @param[in,out] position     current packing position
 * @param[in]     comm         MPI communicator used to communicate the buffer
 */
void yac_point_selection_pack(
  struct yac_point_selection const * point_select,
  void * buffer, int buffer_size, int * position, MPI_Comm comm);

/**
 * Unpack a point selection from a given buffer using `MPI_Unpack`
 * @param[in]     buffer      packing buffer
 * @param[in]     buffer_size buffer size
 * @param[in,out] position    current unpacking position
 * @param[in]     comm        MPI communicator used to receive the buffer
 * @returns unpacked point selection
 */
struct yac_point_selection * yac_point_selection_unpack(
  void const * buffer, int buffer_size, int * position, MPI_Comm comm);

/**
 * Compares two point selections
 * @param[in] a point selection a
 * @param[in] b point selection b
 * @returns -1, 0, or 1 depending on a and b
 */
int yac_point_selection_compare(
  struct yac_point_selection const * a, struct yac_point_selection const * b);

/**
 * Generates a list of points that match a given point selection criterion.
 * @param[in]     point_select           point selection criterion
 * @param[in,out] point_coordinates      point coordinates
 * @param[in,out] point_indices          indices of point associated with the
 *                                       provided coordinates
 * @param[in]     num_points             number of entries in point_indices and
 *                                       point_coordinates
 * @param[out]    num_selected_points    number of entries in
 *                                       selected_point_indices
 * @remark coordinates and indices of selected points are moved to the end
 *         of the respective arrays
 */
void yac_point_selection_apply(
  struct yac_point_selection const * point_select,
  yac_coordinate_pointer point_coordinates, size_t * point_indices,
  size_t num_points, size_t * num_selected_points);

/**
 * Gets the type of a provided point selection
 * @param[in] point_select point selection criterion
 * @return point selection type
 */
enum yac_point_selection_type yac_point_selection_get_type(
  struct yac_point_selection const * point_select);

/**
 * Gets the configuration of a bounding circle point selection
 * @param[in]  point_selection point selection criterion
 * @param[out] center_lon      longitude coordinate of the center of the
 *                             bounding circle (in radians)
 * @param[out] center_lat      latitude coordinate of the center of the
 *                             bounding circle (in radians)
 * @param[out] inc_angle       the angle between center vector and a vector pointing
 *                             to any point on the bounding circle (in radians)
 */
void yac_point_selection_bnd_circle_get_config(
  struct yac_point_selection const * point_selection,
  double * center_lon, double * center_lat, double * inc_angle);

#endif // POINT_SELECTION_H
