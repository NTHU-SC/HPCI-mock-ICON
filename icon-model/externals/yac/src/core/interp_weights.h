// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERP_WEIGHTS_H
#define INTERP_WEIGHTS_H

#include "yac_types.h"
#include "location.h"
#include "interpolation.h"

// YAC PUBLIC HEADER START

struct yac_interp_weights_data {

  double frac_mask_fallback_value; // user-defined fractional mask fallback value
  double scaling_factor;           // user-defined scaling factor
  double scaling_summand;          // user-defined scaling summand

  size_t num_fixed_values;          // number of fixed values
  double * fixed_values;            // fixed values
  size_t * num_tgt_per_fixed_value; // number of target points per fixed value
  size_t * tgt_idx_fixed;           // local ids of fixed target points

  size_t num_wgt_tgt;             // number of target points that receive a
                                  // weighted sum of source points
  size_t * wgt_tgt_idx;           // local ids of weighted target points
  size_t * num_src_per_tgt;       // number of source points per target
  double * weights;               // weights
  size_t * src_field_idx;         // source field index for each source point
  size_t * src_idx;               // index of source points in source field
                                  // buffer
  size_t num_src_fields;          // number of source fields
  size_t * src_field_buffer_size; // buffer sizes required for receiving all
                                  // required source data
                                  // (array size is num_src_fields)
};

struct yac_interpolation_exchange;

enum yac_interp_weights_reorder_type {
  YAC_MAPPING_ON_SRC, //!< weights will be applied at source processes
  YAC_MAPPING_ON_TGT, //!< weights will be applied at target processes
};

enum yac_weight_file_on_existing {
  YAC_WEIGHT_FILE_ERROR     = 0, //!< error when weight file existis already
  YAC_WEIGHT_FILE_KEEP      = 1, //!< keep existing weight file
  YAC_WEIGHT_FILE_OVERWRITE = 2, //!< overwrite existing weight file
  YAC_WEIGHT_FILE_UNDEFINED = 3,
};

#define YAC_WEIGHT_FILE_ON_EXISTING_DEFAULT_VALUE (YAC_WEIGHT_FILE_OVERWRITE)

struct yac_interp_weights;

/**
 * Constructor for interpolation weights.
 * @param[in] comm           MPI communicator
 * @param[in] tgt_location   location of target field
 * @param[in] src_locations  locations of source fields
 * @param[in] num_src_fields number of source fields
 * @return interpolation weights
 */
struct yac_interp_weights * yac_interp_weights_new(
  MPI_Comm comm, enum yac_location tgt_location,
  enum yac_location * src_locations, size_t num_src_fields);

/**
 * writes interpolation weights to file
 * @param[in] weights       interpolation weights
 * @param[in] filename      file name
 * @param[in] src_grid_name name of the source grid
 * @param[in] tgt_grid_name name of the target grid
 * @param[in] src_grid_size global size of the source grid
 * @param[in] tgt_grid_size global size of the target grid
 * @param[in] on_existing   specifies how YAC is supposed to handle the case
 *                          of an already existing file with the same name
 * @remark this call is collective
 * @remark Global grid size argument can be either the global grid size or
 *         zero. If a valid global grid size was provided by at least one
 *         process, it will be added as a dimension to the weight file.
 */
void yac_interp_weights_write_to_file(
  struct yac_interp_weights * weights, char const * filename,
  char const * src_grid_name, char const * tgt_grid_name,
  size_t src_grid_size, size_t tgt_grid_size,
  enum yac_weight_file_on_existing on_existing);

/**
 * generates an interpolation from interpolation weights
 * @param[in] weights                  interpolation weights
 * @param[in] reorder                  determines at which processes the
 *                                     weights are
 *                                     to be applied
 * @param[in] collection_size          collection size
 * @param[in] frac_mask_fallback_value fallback value for dynamic
 *                                     fractional masking
 * @param[in] scaling_factor           scaling factor
 * @param[in] scaling_summand          scaling summand
 * @param[in] yaxt_exchanger_name      name of the yaxt exchanger that is to
 *                                     be used in the interpolation
 * @param[in] is_source                defined whether the local process is a
 *                                     source
 * @param[in] is_target                defined whether the local process is a
 *                                     target
 * @return interpolation
 * @remark if frac_mask_fallback_value != YAC_FRAC_MASK_NO_VALUE, dynamic
 *         fractional masking will be used
 * @remark all target field values, whose source points are not masked by
 *         the fractional mask, that receive an interpolation value, which is
 *         not a fixed value will by scaled by the following formula:\n
 *         y = scaling_factor * x + scaling_summand
 * @remark if yaxt_exchanger_name == NULL, the default exchanger will be used
 * @remark The interpolation weights may contain stencils for processes, which
 *         are not actually a target (is_target == 0). In this case the
 *         respective stencils will be ignored.
 */
struct yac_interpolation * yac_interp_weights_get_interpolation(
  struct yac_interp_weights * weights,
  enum yac_interp_weights_reorder_type reorder,
  size_t collection_size, double frac_mask_fallback_value,
  double scaling_factor, double scaling_summand,
  char const * yaxt_exchanger_name, int is_source, int is_target);

/**
 * generates a raw interpolation from interpolation weights
 *
 * In an exchange, the interpolation will receive on all target process the
 * source points required to compute the local target points based on the
 * interpolation data, which is provied by this routine as well.
 * @param[in]  weights                  interpolation weights
 * @param[in]  collection_size          collection size
 * @param[in]  frac_mask_fallback_value fallback value for dynamic
 *                                      fractional masking
 * @param[in]  scaling_factor           scaling factor
 * @param[in]  scaling_summand          scaling summand
 * @param[in]  yaxt_exchanger_name      name of the yaxt exchanger that is to
 *                                      be used in the interpolation
 * @param[out] interpolation_exchange   interpolation exchange structure
 * @param[out] interp_weights_data      interpolation data required to compute
 *                                      local target points from the received
 *                                      source points
 * @param[in] is_source                 defined whether the local process is a
 *                                      source
 * @param[in] is_target                 defined whether the local process is a
 *                                      target
 * @remark if yaxt_exchanger_name == NULL, the default exchanger will be used
 * @remark memory associated with interp_weights_data can be free by a call to
 *         \ref yac_interp_weights_data_free
 * @remark The interpolation weights may contain stencils for processes, which
 *         are not actually a target (is_target == 0). In this case the
 *         respective stencils will be ignored.
 */
void yac_interp_weights_get_interpolation_raw(
  struct yac_interp_weights * weights,
  size_t collection_size, double frac_mask_fallback_value,
  double scaling_factor, double scaling_summand,
  char const * yaxt_exchanger_name,
  struct yac_interpolation_exchange ** interpolation_exchange,
  struct yac_interp_weights_data * interp_weights_data,
  int is_source, int is_target);

/**
 * returns the count of all target for which the weights contain a stencil
 * @param[in] weights interpolation weights
 * @return count of all targets in weights with a stencil
 */
size_t yac_interp_weights_get_interp_count(
  struct yac_interp_weights * weights);

/**
 * returns the global ids of all targets for which the weights contain a
 * stencil
 * @param[in] weights interpolation weights
 * @return global ids of all targets in weights with a stencil
 */
yac_int * yac_interp_weights_get_interp_tgt(
  struct yac_interp_weights * weights);

/**
 * Destructor for interpolation weights.
 * @param[inout] weights interpolation weights
 */
void yac_interp_weights_delete(struct yac_interp_weights * weights);

/**
 * Initialises an instance of type yac_interp_weights_data with empty data
 * @param[inout] interp_weights_data data to be initialised
 */
void yac_interp_weights_data_init(
  struct yac_interp_weights_data * interp_weights_data);

/**
 * Generates a copy of a provided instance of type yac_interp_weights_data
 * @param[in] interp_weights_data data to be copied
 * @return copy of provided instance of type yac_interp_weights_data
 */
struct yac_interp_weights_data yac_interp_weights_data_copy(
  struct yac_interp_weights_data interp_weights_data);

/**
 * Frees data associated with an instance of type yac_interp_weights_data
 * @param[inout] interp_weights_data data to be free
 * @remark an object of type yac_interp_weights_data can be generate by
 *         \ref yac_interp_weights_get_interpolation_raw
 */
void yac_interp_weights_data_free(
  struct yac_interp_weights_data interp_weights_data);

// YAC PUBLIC HEADER STOP

#endif // INTERP_WEIGHTS_H
