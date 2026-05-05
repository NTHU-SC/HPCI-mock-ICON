// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERP_METHOD_SPMAP_H
#define INTERP_METHOD_SPMAP_H

#include "interp_method.h"
#include "point_selection.h"

// YAC PUBLIC HEADER START

enum yac_interp_spmap_weight_type {
  YAC_INTERP_SPMAP_AVG  = 0, // simple average
  YAC_INTERP_SPMAP_DIST = 1, // distance weighted
};

enum yac_interp_spmap_scale_type {
  YAC_INTERP_SPMAP_NONE       = 0, //!< weights are not scaled
  YAC_INTERP_SPMAP_SRCAREA    = 1, //!< weights are multiplied by
                                   //!< the area of the associated source cell
  YAC_INTERP_SPMAP_INVTGTAREA = 2, //!< weights are muliplied by
                                   //!< the inverse of the area of the
                                   //!< associated target cell
  YAC_INTERP_SPMAP_FRACAREA  = 3,  //!< weights are multiplied by
                                   //!< the area of the associated source cell
                                   //!< and the inverse of the area of the
                                   //!< associated target cell
};

enum yac_interp_spmap_cell_area_provider {
  YAC_INTERP_SPMAP_CELL_AREA_FILE = 0, // read cell areas from file
  YAC_INTERP_SPMAP_CELL_AREA_YAC  = 1, // YAC computes the cell areas
};

struct yac_interp_spmap_config;
struct yac_spmap_overwrite_config;
struct yac_spmap_scale_config;
struct yac_spmap_cell_area_config;

#define YAC_INTERP_SPMAP_SPREAD_DISTANCE_DEFAULT (0.0)
#define YAC_INTERP_SPMAP_MAX_SEARCH_DISTANCE_DEFAULT (0.0)
#define YAC_INTERP_SPMAP_WEIGHTED_DEFAULT (YAC_INTERP_SPMAP_AVG)
#define YAC_INTERP_SPMAP_SCALE_TYPE_DEFAULT (YAC_INTERP_SPMAP_NONE)
#define YAC_INTERP_SPMAP_SPHERE_RADIUS_DEFAULT (1.0)
#define YAC_INTERP_SPMAP_FILENAME_DEFAULT (NULL)
#define YAC_INTERP_SPMAP_VARNAME_DEFAULT (NULL)
#define YAC_INTERP_SPMAP_MIN_GLOBAL_ID_DEFAULT (0)
#define YAC_INTERP_SPMAP_CELL_AREA_CONFIG_DEFAULT (NULL)
#define YAC_INTERP_SPMAP_SCALE_CONFIG_DEFAULT (NULL)
#define YAC_INTERP_SPMAP_DEFAULT_CONFIG (NULL)
#define YAC_INTERP_SPMAP_OVERWRITE_DEFAULT (NULL)

/**
 * Constructor for an interpolation method of type interp_method_spmap\n
 * This method searches for each unmasked source point the closest unmasked
 * target point.\n
 * If the maximum search distance is > 0.0, only target points that are within
 * this distance from the source points are being considered.
 * If spread_distance is > 0.0, the method uses the previously found target
 * points as a starting point. Around each starting point, a bounding circle is
 * generated. Afterwards for each starting point all target cells whose bounding
 * circles intersect with the generated one are put into a list. Out of this
 * list all target cells connected directly or indirectly through other cells
 * from this list to the target cell of the starting are selected for the
 * interpolation. Then a weighting method is applied to the selected target
 * cells to generate the weights. Afterwards the weights are scaled. \n
 * The default configuration can be overwritten for selected source cells. An
 * arbitrary number of alternative configuration along with selection criteria
 * can be provided. These configurations will be processed one-by-one
 * starting with the first entry in the array. Any source cell which has already
 * been processed previously will be ignored by any further configuration. Once all alternative
 * configurations are processed the default configuration is applied to the
 * remaining source points.
 * @param[in] default_config         default configuration
 * @param[in] overwrite_configs      alternative configurations along with
 *                                   source cell selection criteria
 * @remark * the unit for the spread and maximum search distance is Radian
 *         * if (overwrite_configs != NULL) the last entry in the array has to
 *           be NULL
 */
struct interp_method * yac_interp_method_spmap_new(
  struct yac_interp_spmap_config const * default_config,
  struct yac_spmap_overwrite_config const * const * overwrite_configs);

/**
 * Constructs a Source to Target mapping configuration object
 * @param[in] spread_distance     spread distance
 * @param[in] max_search_distance maximum search distance
 * @param[in] weight_type         weightening type
 * @param[in] scale_config        scaling configuration
 * @return Source to Target mapping configuration object
 * @remark see \ref interp_method_spmap for more information
 */
struct yac_interp_spmap_config * yac_interp_spmap_config_new(
  double spread_distance, double max_search_distance,
  enum yac_interp_spmap_weight_type weight_type,
  struct yac_spmap_scale_config const * scale_config);

/**
 * Destructor for a Source to Target mapping configuration object
 * @param[in] spmap_config Source to Target mapping configuration object
 */
void yac_interp_spmap_config_delete(
  struct yac_interp_spmap_config * spmap_config);

/**
 * Returns spread distance of Source to Target mapping configuration object
 * @param[in] spmap_config Source to Target mapping configuration object
 * @return spread distance
 */
double yac_interp_spmap_config_get_spread_distance(
  struct yac_interp_spmap_config const * spmap_config);

/**
 * Returns maximum search distance distance of Source to Target mapping
 * configuration object
 * @param[in] spmap_config Source to Target mapping configuration object
 * @return maximum search distance
 */
double yac_interp_spmap_config_get_max_search_distance(
  struct yac_interp_spmap_config const * spmap_config);

/**
 * Returns weight type of Source to Target mapping configuration object
 * @param[in] spmap_config Source to Target mapping configuration object
 * @return spread distance
 */
enum yac_interp_spmap_weight_type yac_interp_spmap_config_get_weight_type(
  struct yac_interp_spmap_config const * spmap_config);

/**
 * Returns scaling configuration of Source to Target mapping configuration
 * object
 * @param[in] spmap_config Source to Target mapping configuration object
 * @return scaling configuration
 */
struct yac_spmap_scale_config const * yac_interp_spmap_config_get_scale_config(
  struct yac_interp_spmap_config const * spmap_config);

/**
 * Constructs a Source to Target mapping scaling configuration object
 * @param[in] scale_type              scaling type
 * @param[in] source_cell_area_config source cell area configuration
 * @param[in] target_cell_area_config target cell area configuration
 * @return Source to Target mapping scaling configuration object
 * @remark see \ref interp_method_spmap for more information
 */
struct yac_spmap_scale_config * yac_spmap_scale_config_new(
  enum yac_interp_spmap_scale_type scale_type,
  struct yac_spmap_cell_area_config const * source_cell_area_config,
  struct yac_spmap_cell_area_config const * target_cell_area_config);

/**
 * Destructor for a Source to Target mapping scaling configuration object
 * @param[in] scale_config scaling configuration
 */
void yac_spmap_scale_config_delete(
  struct yac_spmap_scale_config * scale_config);

/**
 * Gets type of  a Source to Target mapping scaling configuration object
 * @param[in] scale_config scaling configuration
 * @return scaling configuration type
 */
enum yac_interp_spmap_scale_type yac_spmap_scale_config_get_type(
  struct yac_spmap_scale_config const * scale_config);

/**
 * Gets source cell area configuration of a Source to Target mapping scaling
 * configuration object
 * @param[in] scale_config scaling configuration
 * @return source cell area configuration
 */
struct yac_spmap_cell_area_config const *
  yac_spmap_scale_config_get_src_cell_area_config(
    struct yac_spmap_scale_config const * scale_config);

/**
 * Gets target cell area configuration of a Source to Target mapping scaling
 * configuration object
 * @param[in] scale_config scaling configuration
 * @return target cell area configuration
 */
struct yac_spmap_cell_area_config const *
  yac_spmap_scale_config_get_tgt_cell_area_config(
    struct yac_spmap_scale_config const * scale_config);

/**
 * Constructs a Source to Target mapping cell area configuration object
 *
 * This type enable internal computation of cell areas.
 * @param[in] sphere_radius sets sphere radius to be used for area computation
 * @return Source to Target mapping cell area configuration object
 */
struct yac_spmap_cell_area_config * yac_spmap_cell_area_config_yac_new(
  double sphere_radius);

/**
 * Constructs a Source to Target mapping cell area configuration object
 *
 * This type enable reading of cell areas from a netCDF file.
 * @param[in] filename      name of the file containing the cell areas
 * @param[in] varname       name of the variable containing the cell areas
 * @param[in] min_global_id minimum global cell id (usually "0" or "1")
 * @return Source to Target mapping cell area configuration object
 * @remark Cell area for cell with global id X is assumed to be at
 *         (using C indexing; starting at "0"): varname[X-min_global_id]
 */
struct yac_spmap_cell_area_config * yac_spmap_cell_area_config_file_new(
  char const * filename, char const * varname, yac_int min_global_id);

/**
 * Destructor for a Source to Target mapping cell area configuration object
 * @param[in] cell_area_config Source to Target mapping cell area
 *                             configuration object
 */
void yac_spmap_cell_area_config_delete(
  struct yac_spmap_cell_area_config * cell_area_config);

/**
 * Gets the type of a cell area configuration object
 * @param[in] cell_area_config cell area configuration object
 * @return cell area configuration type
 */
enum yac_interp_spmap_cell_area_provider yac_spmap_cell_area_config_get_type(
  struct yac_spmap_cell_area_config const * cell_area_config);

/**
 * Gets the sphere radius of a cell area configuration object
 * @param[in] cell_area_config cell area configuration object
 * @return sphere radius
 * @remark this call is only valid for cell area configuration objects of type
 *         YAC_INTERP_SPMAP_CELL_AREA_YAC
 */
double yac_spmap_cell_area_config_get_sphere_radius(
  struct yac_spmap_cell_area_config const * cell_area_config);

/**
 * Gets the file name of a cell area configuration object
 * @param[in] cell_area_config cell area configuration object
 * @return file name
 * @remark this call is only valid for cell area configuration objects of type
 *         YAC_INTERP_SPMAP_CELL_AREA_FILE
 */
char const * yac_spmap_cell_area_config_get_filename(
  struct yac_spmap_cell_area_config const * cell_area_config);

/**
 * Gets the variable name of a cell area configuration object
 * @param[in] cell_area_config cell area configuration object
 * @return variable name
 * @remark this call is only valid for cell area configuration objects of type
 *         YAC_INTERP_SPMAP_CELL_AREA_FILE
 */
char const * yac_spmap_cell_area_config_get_varname(
  struct yac_spmap_cell_area_config const * cell_area_config);

/**
 * Gets the minimum global id of a cell area configuration object
 * @param[in] cell_area_config cell area configuration object
 * @return minimum global id
 * @remark this call is only valid for cell area configuration objects of type
 *         YAC_INTERP_SPMAP_CELL_AREA_FILE
 */
yac_int yac_spmap_cell_area_config_get_min_global_id(
  struct yac_spmap_cell_area_config const * cell_area_config);

/**
 * Constructs an alternative Source to Target mapping configuration
 * @param[in] src_point_selection Specifies source points the alternative
 *                                configuration is to be appied to
 * @param[in] config              Alternative configuration to be used
 */
struct yac_spmap_overwrite_config * yac_spmap_overwrite_config_new(
  struct yac_point_selection const * src_point_selection,
  struct yac_interp_spmap_config const * config);

/**
 * Destructor for an alternative Source to Target mapping configuration
 * @param[in] overwrite_config Alternative Source to Target mapping
 *                             configuration
 */
void yac_spmap_overwrite_config_delete(
  struct yac_spmap_overwrite_config * overwrite_config);

/**
 * Destructor for an array of alternative Source to Target mapping
 * configurations
 * @param[in] overwrite_configs Array of alternative Source to Target
 *                              mapping configurations
 */
void yac_spmap_overwrite_configs_delete(
  struct yac_spmap_overwrite_config ** overwrite_configs);

/**
 * Gets the source point selection method of an alternative
 * Source to Target mapping configuration
 * @param[in] overwrite_config Alternative Source to Target mapping
 *                             configurations
 * @return source point selection method
 */
struct yac_point_selection const *
  yac_spmap_overwrite_config_get_src_point_selection(
    struct yac_spmap_overwrite_config const * overwrite_config);

/**
 * Gets the configuration of an alternative Source to Target mapping
 * configuration
 * @param[in] overwrite_config Alternative Source to Target mapping
 *                             configurations
 * @return configuration
 */
struct yac_interp_spmap_config const *
  yac_spmap_overwrite_config_get_spmap_config(
    struct yac_spmap_overwrite_config const * overwrite_config);

// YAC PUBLIC HEADER STOP

//------------------------------------------------------------------------------
// some utility rouines for datatypes associated with spmap
//------------------------------------------------------------------------------

struct yac_interp_spmap_config * yac_interp_spmap_config_copy(
  struct yac_interp_spmap_config const * spmap_config);
struct yac_spmap_overwrite_config ** yac_spmap_overwrite_configs_copy(
  struct yac_spmap_overwrite_config const * const * overwrite_configs);

int yac_interp_spmap_config_compare(
  struct yac_interp_spmap_config const * a,
  struct yac_interp_spmap_config const * b);
int yac_spmap_overwrite_config_compare(
  struct yac_spmap_overwrite_config const * a,
  struct yac_spmap_overwrite_config const * b);
size_t yac_interp_spmap_config_get_pack_size(
  struct yac_interp_spmap_config const * spmap_config, MPI_Comm comm);
size_t  yac_spmap_overwrite_configs_get_pack_size(
  struct yac_spmap_overwrite_config const * const * overwrite_configs,
  MPI_Comm comm);
void yac_interp_spmap_config_pack(
  struct yac_interp_spmap_config const * spmap_config,
  void * buffer, int buffer_size, int * position, MPI_Comm comm);
void yac_spmap_overwrite_configs_pack(
  struct yac_spmap_overwrite_config const * const * overwrite_configs,
  void * buffer, int buffer_size, int * position, MPI_Comm comm);
void yac_interp_spmap_config_unpack(
  void * buffer, int buffer_size, int * position,
  struct yac_interp_spmap_config ** spmap_config, MPI_Comm comm);
void yac_spmap_overwrite_configs_unpack(
  void * buffer, int buffer_size, int * position,
  struct yac_spmap_overwrite_config *** overwrite_configs, MPI_Comm comm);

#endif // INTERP_METHOD_SPMAP_H
