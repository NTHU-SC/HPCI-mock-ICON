// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
// Get the definition of the 'restrict' keyword.
#include "config.h"
#endif

#include <float.h>
#include <string.h>

#include "interp_method_internal.h"
#include "interp_method_spmap.h"
#include "ensure_array_size.h"
#include "area.h"
#include "io_utils.h"
#include "point_selection.h"
#include "yac_mpi_internal.h"

static size_t do_search_spmap(struct interp_method * method,
                              struct yac_interp_grid * interp_grid,
                              size_t * tgt_points, size_t count,
                              struct yac_interp_weights * weights,
                              int * interpolation_complete);
static void delete_spmap(struct interp_method * method);

static struct interp_method_vtable
  interp_method_spmap_vtable = {
    .do_search = do_search_spmap,
    .delete = delete_spmap
};

/**
 * Datatype used to store a single link
 * (reference to source and target cell and the associated weight).
 * Using this datatye makes links easily comparable, which produces
 * decomposition independent weight order, if the compare is based on
 * user-provided global ids.
 */
typedef struct interp_link {
  size_t tgt_point;
  struct {
    yac_int global;
    size_t local;
  } src_point;
  double weight;
} InterpLink;

struct yac_spmap_cell_area_config {
  enum yac_interp_spmap_cell_area_provider type;
  union {
    struct {
      double sphere_radius; // used for computing the cell areas
    } yac;
    struct yac_spmap_cell_area_file_config {
      char const * filename;
      char const * varname;
      yac_int min_global_id;
    } file_config;        // used for read the cell areas from file
  };
};

static struct yac_spmap_cell_area_config cell_area_config_default =
  {.type = YAC_INTERP_SPMAP_CELL_AREA_YAC,
   .yac.sphere_radius = YAC_INTERP_SPMAP_SPHERE_RADIUS_DEFAULT};

struct yac_spmap_scale_config {
  enum yac_interp_spmap_scale_type type;
  struct yac_spmap_cell_area_config * src;
  struct yac_spmap_cell_area_config * tgt;
};

static struct yac_spmap_scale_config spmap_scale_config_default =
  {.type = YAC_INTERP_SPMAP_SCALE_TYPE_DEFAULT,
   .src = NULL,
   .tgt = NULL};

struct yac_interp_spmap_config {
  double spread_distance;
  double max_search_distance;
  enum yac_interp_spmap_weight_type weight_type;
  struct yac_spmap_scale_config * scale_config;
};

static struct yac_interp_spmap_config spmap_config_default =
 {.spread_distance = YAC_INTERP_SPMAP_SPREAD_DISTANCE_DEFAULT,
  .max_search_distance = YAC_INTERP_SPMAP_MAX_SEARCH_DISTANCE_DEFAULT,
  .weight_type = YAC_INTERP_SPMAP_WEIGHTED_DEFAULT,
  .scale_config = YAC_INTERP_SPMAP_SCALE_CONFIG_DEFAULT};

struct yac_spmap_overwrite_config {
  struct yac_point_selection * src_point_selection;
  struct yac_interp_spmap_config * config;
};

static struct yac_spmap_overwrite_config overwrite_config_default =
  {.src_point_selection = NULL,
   .config = NULL};

struct interp_method_spmap {

  struct interp_method_vtable * vtable;
  struct yac_interp_spmap_config * default_config;
  struct yac_spmap_overwrite_config ** overwrite_configs;
};

static inline int compare_size_t(const void * a, const void * b) {

  size_t const * a_ = a, * b_ = b;

  return (*a_ > *b_) - (*b_ > *a_);
}

/**
 * Checks initial results from bounding circle search
 * (if spread_distance > 0.0) and removes target points that are not connected
 * with the initial search (single) search result or are too far away from it.
 * @param[in]     interp_grid
 * @param[in]     spread_distance          sin and cos values of spread distance
 * @param[in]     num_src_points           number of source points/initial
 *                                         target points
 * @param[in]     tgt_result_points        initial target result points (the
 *                                         one closest to the source point)
 * @param[in,out] num_tgt_per_src          number of bounding circle search
 *                                         results per source points/initial
 *                                         target point
 * @param[in,out] spread_tgt_result_points results of bounding circle search
 * @remark num_tgt_per_src and spread_tgt_result_points are updated by this
 *         routine
 */
static size_t check_tgt_result_points(
  struct yac_interp_grid * interp_grid, struct sin_cos_angle spread_distance,
  size_t num_src_points, size_t const * const tgt_result_points,
  size_t * num_tgt_per_src, size_t * spread_tgt_result_points) {

  struct yac_const_basic_grid_data * tgt_basic_grid_data =
    yac_interp_grid_get_basic_grid_data_tgt(interp_grid);
  size_t local_num_tgt_cells = tgt_basic_grid_data->count[YAC_LOC_CELL];
  yac_const_coordinate_pointer tgt_field_coords =
    yac_interp_grid_get_tgt_field_coords(interp_grid);
  yac_size_t_2_pointer edge_to_cell =
    yac_interp_grid_generate_tgt_edge_to_cell(interp_grid);

  enum {
    OUTSIDE = 0,   // not part of initial bounding circle search results
                   // --> not to be considered
    INCLUDED = 0,  // already included in the final list of target result points
    CANDIDATE = 1, // potential candidate based on initial bounding
                   // circle search results
  } * tgt_state =
    xmalloc(local_num_tgt_cells * sizeof(*tgt_state));

  size_t * curr_tgt_results = spread_tgt_result_points;
  size_t * curr_bnd_search_tgt_results = spread_tgt_result_points;

  size_t total_num_weights = 0;

  // for all source points
  for (size_t i = 0; i < num_src_points; ++i) {

    // first mark all target cell as being outside and the only set results
    // from the bounding circle search as potential candidates that can be
    // considered for the final result
    memset(
       tgt_state, OUTSIDE, local_num_tgt_cells * sizeof(*tgt_state));
     size_t curr_bnd_search_result_count = num_tgt_per_src[i];
     for (size_t j = 0; j < curr_bnd_search_result_count; ++j) {
       tgt_state[curr_bnd_search_tgt_results[j]] = CANDIDATE;
     }
     curr_bnd_search_tgt_results += curr_bnd_search_result_count;

    // add original results to list of final results
    size_t tgt_start_point = tgt_result_points[i];
    size_t curr_num_tgt_results = 1;
    curr_tgt_results[0] = tgt_start_point;
    tgt_state[tgt_start_point] = INCLUDED;
    double const * start_coord = tgt_field_coords[tgt_start_point];
    size_t * prev_added_tgts = curr_tgt_results;

    size_t prev_num_added_tgt = 0; // number of added targets in prev iteration
    size_t curr_num_added_tgt = 1; // number of added targets in curr iteration
    do {

      prev_added_tgts += prev_num_added_tgt;
      prev_num_added_tgt = curr_num_added_tgt;
      curr_num_added_tgt = 0;

      // for all targets that were added in the previous iteration
      for (size_t j = 0; j < prev_num_added_tgt; ++j) {

        // get target that was added in previous iteration
        size_t curr_tgt_result = prev_added_tgts[j];

        size_t const * curr_cell_to_edge =
          tgt_basic_grid_data->cell_to_edge +
          tgt_basic_grid_data->cell_to_edge_offsets[curr_tgt_result];
        size_t curr_num_edges =
          tgt_basic_grid_data->num_vertices_per_cell[curr_tgt_result];

        // for all neighbour cells
        for (size_t k = 0; k < curr_num_edges; ++k) {

          size_t * curr_edge_to_cell = edge_to_cell[curr_cell_to_edge[k]];
          size_t tgt_candidate =
            curr_edge_to_cell[curr_edge_to_cell[0] == curr_tgt_result];

          // if there actually is a cell (SIZE_MAX -> no cell) and
          // if the target is a potential candiate
          if ((tgt_candidate != SIZE_MAX) &&
              (tgt_state[tgt_candidate] == CANDIDATE)) {

            // add target if within spread distance
            struct sin_cos_angle distance =
              get_vector_angle_2(start_coord, tgt_field_coords[tgt_candidate]);
            if (compare_angles(distance, spread_distance) <= 0) {

              curr_tgt_results[curr_num_tgt_results] = tgt_candidate;
              ++curr_num_tgt_results;
              ++curr_num_added_tgt;

              // mark target as being included in the result list
              tgt_state[tgt_candidate] = INCLUDED;

            } else {

              // mark target as being to far away for the initial search result
              tgt_state[tgt_candidate] = OUTSIDE;
            }
          } // is valid new tgt
        } // num edges
      } // prev num added tgt

    } while (curr_num_added_tgt > 0); // continue until no new targets were added

    curr_tgt_results += curr_num_tgt_results;
    num_tgt_per_src[i] = curr_num_tgt_results;
    total_num_weights += curr_num_tgt_results;
  }

  free(edge_to_cell);
  free(tgt_state);

  return total_num_weights;
}

/**
 * Computes the weights
 * @param[in] interp_grid
 * @param[in] weight_type       determines how the source value is distribted
 *                              to the target points
 * @param[in] num_src_points    number of source points
 * @param[in] src_points        local indices of all interpolated source points
 * @param[in] num_tgt_per_src   number of target points to which each source
 *                              value is to be distributed
 * @param[in] total_num_tgt     total size of tgt_result_points or
 *                              SUM(num_tgt_per_src(:))
 * @param[in] tgt_result_points local indices of all targets
 * @returns weights associated with each entry in tgt_result_points
 */
static double * compute_weights(
  struct yac_interp_grid * interp_grid,
  enum yac_interp_spmap_weight_type weight_type, size_t num_src_points,
  size_t const * const src_points, size_t const * const num_tgt_per_src,
  size_t total_num_tgt, size_t const * const tgt_result_points) {

  double * weights = xmalloc(total_num_tgt * sizeof(*weights));

  switch (weight_type) {

    YAC_UNREACHABLE_DEFAULT("ERROR(do_search_spmap): invalid weight_type");

    // simple average
    case (YAC_INTERP_SPMAP_AVG): {
      for (size_t i = 0, offset = 0; i < num_src_points; ++i) {
        size_t curr_num_tgt = num_tgt_per_src[i];
        if (curr_num_tgt == 0) continue;
        if (curr_num_tgt > 1) {
          double curr_weight_data = 1.0 / (double)(curr_num_tgt);
          for (size_t j = 0; j < curr_num_tgt; ++j, ++offset) {
            weights[offset] = curr_weight_data;
          }
        } else {
          weights[offset] = 1.0;
          ++offset;
        }
      }
      break;
    }

    // distance weighted
    case (YAC_INTERP_SPMAP_DIST): {

      yac_const_coordinate_pointer src_field_coords =
        yac_interp_grid_get_src_field_coords(interp_grid, 0);
      yac_const_coordinate_pointer tgt_field_coords =
        yac_interp_grid_get_tgt_field_coords(interp_grid);

      // for each source point
      for (size_t i = 0, offset = 0; i < num_src_points; ++i) {

        size_t curr_num_tgt = num_tgt_per_src[i];

        if (curr_num_tgt == 0) continue;

        double * curr_weights = weights + offset;

        if (curr_num_tgt > 1) {

          size_t const * const  curr_result_points =
            tgt_result_points + offset;
          double const * curr_src_coord = src_field_coords[src_points[offset]];
          offset += curr_num_tgt;

          int match_flag = 0;

          for (size_t j = 0; j < curr_num_tgt; ++j) {

            double distance =
              get_vector_angle(
                (double*)curr_src_coord,
                (double*)tgt_field_coords[curr_result_points[j]]);

            if (distance < yac_angle_tol) {
              for (size_t k = 0; k < curr_num_tgt; ++k) curr_weights[k] = 0.0;
              curr_weights[j] = 1.0;
              match_flag = 1;
              break;
            }
            curr_weights[j] = 1.0 / distance;
          }

          if (!match_flag) {

            // compute scaling factor for the weights
            double inv_distance_sum = 0.0;
            for (size_t j = 0; j < curr_num_tgt; ++j)
              inv_distance_sum += curr_weights[j];
            double scale = 1.0 / inv_distance_sum;

            for (size_t j = 0; j < curr_num_tgt; ++j) curr_weights[j] *= scale;
          }
        } else {
          *curr_weights = 1.0;
          ++offset;
        }
      }
      break;
    }
  };

  return weights;
}

static void compute_cell_areas(
  struct yac_const_basic_grid_data * basic_grid_data,
  int const * required_cell_areas, double area_scale, char const * type,
  double * cell_areas) {

  struct yac_grid_cell grid_cell;
  yac_init_grid_cell(&grid_cell);

  size_t num_cells = basic_grid_data->count[YAC_LOC_CELL];

  for (size_t i = 0; i < num_cells; ++i) {
    if (required_cell_areas[i]) {
      yac_const_basic_grid_data_get_grid_cell(basic_grid_data, i, &grid_cell);
      double cell_area = yac_huiliers_area(grid_cell) * area_scale;
      YAC_ASSERT_F(
        cell_area > YAC_AREA_TOL,
        "ERROR(get_cell_areas): "
        "area of %s cell (global id %"XT_INT_FMT") is close to zero (%e)",
        type, basic_grid_data->ids[YAC_LOC_CELL][i], cell_area);
      cell_areas[i] = cell_area;
    } else {
      cell_areas[i] = 0.0;
    }
  }

  yac_free_grid_cell(&grid_cell);
}

static void dist_read_cell_areas(
  char const * filename, char const * varname, MPI_Comm comm,
  int * const io_ranks, int io_rank_idx, int num_io_ranks,
  double ** dist_cell_areas, size_t * dist_cell_areas_global_size) {

#ifndef YAC_NETCDF_ENABLED

  UNUSED(filename);
  UNUSED(varname);
  UNUSED(comm);
  UNUSED(io_ranks);
  UNUSED(io_rank_idx);
  UNUSED(num_io_ranks);

  *dist_cell_areas = NULL;
  *dist_cell_areas_global_size = 0;

  die(
    "ERROR(interp_method_spmap::dist_read_cell_areas): "
    "YAC is built without the NetCDF support");
#else

  if ((io_rank_idx >= 0) && (io_rank_idx < num_io_ranks)) {

    // open file
    int ncid;
    yac_nc_open(filename, NC_NOWRITE, &ncid);

    // get variable id
    int varid;
    yac_nc_inq_varid(ncid, varname, &varid);

    // get dimension ids
    int ndims;
    int dimids[NC_MAX_VAR_DIMS];
    YAC_HANDLE_ERROR(
      nc_inq_var(ncid, varid, NULL, NULL, &ndims, dimids, NULL));

    YAC_ASSERT_F(
      (ndims == 1) || (ndims == 2),
      "ERROR(dist_read_cell_areas): "
      "invalid number of dimensions for cell area variable \"%s\" from "
      "file \"%s\" (is %d, but should be either 1 or 2)",
      varname, filename, ndims);

    // get size of dimensions
    size_t dimlen[2];
    *dist_cell_areas_global_size = 1;
    for (int i = 0; i < ndims; ++i) {
      YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dimids[i], &dimlen[i]));
      YAC_ASSERT_F(
        dimlen[i] > 0,
        "ERROR(dist_read_cell_areas): "
        "invalid dimension size for cell area variable \"%s\" from "
        "file \"%s\" (is %zu, should by > 0)",
        varname, filename, dimlen[i]);
      *dist_cell_areas_global_size *= dimlen[i];
    }

    // compute start/count
    // (in 2D case we have to round a little bit and then adjust the
    //  data afterwards)
    size_t start[2], count[2], offset, read_size;
    size_t local_start =
      (size_t)(
        ((long long)*dist_cell_areas_global_size * (long long)io_rank_idx) /
        (long long)num_io_ranks);
    size_t local_count =
      (size_t)(
        ((long long)*dist_cell_areas_global_size *
          (long long)(io_rank_idx+1)) / (long long)num_io_ranks) - local_start;
    if (ndims == 1) {
      start[0] = local_start;
      count[0] = local_count;
      offset = 0;
      read_size = local_count;
    } else {
      start[0] = local_start / dimlen[1];
      count[0] =
        (local_start + local_count + dimlen[1] - 1) / dimlen[1] - start[0];
      start[1] = 0;
      count[1] = dimlen[1];
      offset = local_start - start[0] * dimlen[1];
      read_size = count[0] * count[1];
    }

    // read data
    *dist_cell_areas = xmalloc(read_size * sizeof(**dist_cell_areas));
    YAC_HANDLE_ERROR(
      nc_get_vara_double(ncid, varid, start, count, *dist_cell_areas));

    // adjust data if necessary
    if (ndims == 2)
      memmove(
        *dist_cell_areas, *dist_cell_areas + offset,
        local_count * sizeof(**dist_cell_areas));

    // close file
    YAC_HANDLE_ERROR(nc_close(ncid));

  } else {
    *dist_cell_areas = xmalloc(1*sizeof(**dist_cell_areas));
    *dist_cell_areas_global_size = 0;
  }

  yac_mpi_call(
    MPI_Bcast(
      dist_cell_areas_global_size, 1, YAC_MPI_SIZE_T, io_ranks[0], comm),
    comm);

#endif // YAC_NETCDF_ENABLED
}

static void read_cell_areas(
  struct yac_const_basic_grid_data * basic_grid_data,
  struct yac_spmap_cell_area_file_config file_config, MPI_Comm comm,
  int const * required_cell_areas, char const * type,
  double * cell_areas) {

  char const * routine = "read_cell_areas";

  size_t num_cells = basic_grid_data->count[YAC_LOC_CELL];

  // get io configuration
  int local_is_io, * io_ranks, num_io_ranks;
  yac_get_io_ranks(comm, &local_is_io, &io_ranks, &num_io_ranks);

  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  int io_rank_idx = 0;
  while ((io_rank_idx < num_io_ranks) &&
          (comm_rank != io_ranks[io_rank_idx]))
    ++io_rank_idx;
  YAC_ASSERT_F(
    !local_is_io || (io_rank_idx < num_io_ranks),
    "ERROR(%s): unable to determine io_rank_idx", routine);

  double * dist_cell_areas;
  size_t dist_cell_areas_global_size;

  // read the data on the io processes
  dist_read_cell_areas(
    file_config.filename, file_config.varname, comm,
    io_ranks, io_rank_idx, num_io_ranks,
    &dist_cell_areas, &dist_cell_areas_global_size);

  // count the number of locally required cell areas
  size_t num_required_cell_areas = 0;
  for (size_t i = 0; i < num_cells; ++i)
    if (required_cell_areas[i]) ++num_required_cell_areas;

  size_t * global_idx = xmalloc(num_required_cell_areas * sizeof(*global_idx));
  size_t * reorder_idx =
    xmalloc(num_required_cell_areas * sizeof(*reorder_idx));

  size_t * sendcounts, * recvcounts, * sdispls, * rdispls;
  yac_get_comm_buffers(
    1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);

  // determine global indices of locally required cell areas
  yac_int const * global_cell_ids = basic_grid_data->ids[YAC_LOC_CELL];
  yac_int min_global_id = file_config.min_global_id;
  for (size_t i = 0, j = 0; i < num_cells; ++i) {
    if (required_cell_areas[i]) {
      YAC_ASSERT_F(
        global_cell_ids[i] >= min_global_id,
        "ERROR(%s): %s global id %" XT_INT_FMT " is smaller than "
        "the minimum global id provided by the user (%" XT_INT_FMT ")",
        routine, type, global_cell_ids[i], min_global_id);
      size_t idx = (size_t)(global_cell_ids[i] - min_global_id);
      YAC_ASSERT_F(
        idx < dist_cell_areas_global_size,
        "ERROR(%s): %s global id %" XT_INT_FMT " exceeds "
        "available size of array \"%s\" in file \"%s\" "
        "(min_global_id %" XT_INT_FMT ")", routine, type, global_cell_ids[i],
        file_config.varname, file_config.filename, file_config.min_global_id);
      global_idx[j] = idx;
      reorder_idx[j] = i;
      int dist_rank_idx =
        ((long long)idx * (long long)num_io_ranks +
         (long long)num_io_ranks - 1) / (long long)dist_cell_areas_global_size;
      sendcounts[io_ranks[dist_rank_idx]]++;
      ++j;
    }
  }
  free(io_ranks);

  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);

  // sort required data by their global index
  yac_quicksort_index_size_t_size_t(
    global_idx, num_required_cell_areas, reorder_idx);

  // send required points to io processes
  size_t request_count = recvcounts[comm_size - 1] + rdispls[comm_size - 1];
  size_t * request_global_idx =
    xmalloc(request_count * sizeof(*request_global_idx));
  yac_alltoallv_p2p(
    global_idx, sendcounts, sdispls+1,
    request_global_idx, recvcounts, rdispls,
    sizeof(*global_idx), YAC_MPI_SIZE_T, comm, routine, __LINE__);
  free(global_idx);

  // pack requested cell areas
  double * requested_cell_areas =
    xmalloc(request_count * sizeof(*requested_cell_areas));
  size_t global_idx_offset =
    ((long long)io_rank_idx * (long long)dist_cell_areas_global_size) /
    (long long)num_io_ranks;
  for (size_t i = 0; i < request_count; ++i)
    requested_cell_areas[i] =
      dist_cell_areas[request_global_idx[i] - global_idx_offset];
  free(request_global_idx);
  free(dist_cell_areas);

  // return cell areas
  double * temp_cell_areas =
    xmalloc(num_required_cell_areas * sizeof(*temp_cell_areas));
  yac_alltoallv_p2p(
    requested_cell_areas, recvcounts, rdispls,
    temp_cell_areas, sendcounts, sdispls+1,
    sizeof(*requested_cell_areas), MPI_DOUBLE, comm, routine, __LINE__);
  free(requested_cell_areas);

  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);

  // unpack cell areas
  for (size_t i = 0; i < num_required_cell_areas; ++i)
    cell_areas[reorder_idx[i]] = temp_cell_areas[i];

  free(temp_cell_areas);
  free(reorder_idx);
}

/**
 * Gets cell areas for scaling either from file or compute them
 * (depending on configuration)
 * @param[in] basic_grid_data
 * @param[in] type             used for error message
 *                             (either "source" or "target");
 * @param[in] cell_area_config determines whether areas are written from file
 *                             or computed and contains additional
 *                             configuration specific to respective type
 * @param[in] comm             MPI communicator (required for parallel reading
 *                             of areas from file)
 * @param[in] required_points  local ids of cell for which the areas are
 *                             required
 * @param[in] num_required_points
 * @return cell areas
 */
static double * get_cell_areas(
  struct yac_const_basic_grid_data * basic_grid_data,
  char const * type, struct yac_spmap_cell_area_config const * cell_area_config,
  MPI_Comm comm, size_t const * required_points, size_t num_required_points) {

  size_t num_cells = basic_grid_data->count[YAC_LOC_CELL];

  // determine which cell areas are actually required
  int * required_cell_areas = xcalloc(num_cells, sizeof(*required_cell_areas));
  for (size_t i = 0; i < num_required_points; ++i)
    required_cell_areas[required_points[i]] = 1;

  // compute and scale cell areas
  double * cell_areas = xmalloc(num_cells * sizeof(*cell_areas));

  switch (cell_area_config->type) {

    YAC_UNREACHABLE_DEFAULT_F(
      "ERROR(get_cell_areas): unsupported %s cell area origin", type);

    case (YAC_INTERP_SPMAP_CELL_AREA_YAC): {

      double area_scale =
        cell_area_config->yac.sphere_radius *
        cell_area_config->yac.sphere_radius;
      compute_cell_areas(
        basic_grid_data, required_cell_areas, area_scale, type,
        cell_areas);
      break;
    }
    case (YAC_INTERP_SPMAP_CELL_AREA_FILE): {

      read_cell_areas(
        basic_grid_data, cell_area_config->file_config, comm,
        required_cell_areas, type, cell_areas);
    }
  };

  free(required_cell_areas);

  return cell_areas;
}

/**
 * Scale weights based on source/target cell areas depending on used
 * configuration.
 * @remark this is required if for example the source field contains data in
 *         kg/m^2 (a flux) but the target expects it in kg
 */
static void scale_weights(
  struct yac_interp_grid * interp_grid,
  struct yac_spmap_scale_config const * scale_config,
  size_t num_src_points, size_t const * src_points,
  size_t const * num_tgt_per_src, size_t const * tgt_points,
  size_t total_num_weights, double * weights) {

  // if there is no scaling
  if (scale_config->type == YAC_INTERP_SPMAP_NONE) return;

  // get cell areas if they are required
  double const * src_cell_areas =
    ((scale_config->type == YAC_INTERP_SPMAP_SRCAREA) ||
     (scale_config->type == YAC_INTERP_SPMAP_FRACAREA))?
      get_cell_areas(
        yac_interp_grid_get_basic_grid_data_src(interp_grid),
        "source", scale_config->src, yac_interp_grid_get_MPI_Comm(interp_grid),
        src_points, total_num_weights):NULL;
  double const * tgt_cell_areas =
    ((scale_config->type == YAC_INTERP_SPMAP_INVTGTAREA) ||
     (scale_config->type == YAC_INTERP_SPMAP_FRACAREA))?
      get_cell_areas(
        yac_interp_grid_get_basic_grid_data_tgt(interp_grid),
        "target", scale_config->tgt, yac_interp_grid_get_MPI_Comm(interp_grid),
        tgt_points, total_num_weights):NULL;

#define SCALE_WEIGHTS( \
  SRC_CELL_AREA, TGT_CELL_AREA) \
{ \
  for (size_t i = 0, offset = 0; i < num_src_points; ++i) { \
    size_t curr_num_tgt = num_tgt_per_src[i]; \
    for (size_t j = 0; j < curr_num_tgt; ++j, ++offset) { \
      weights[offset] *= SRC_CELL_AREA / TGT_CELL_AREA; \
    } \
  } \
}

  switch (scale_config->type) {
    YAC_UNREACHABLE_DEFAULT(
      "ERROR(scale_weights): invalid scale_type");
    case(YAC_INTERP_SPMAP_SRCAREA):
      SCALE_WEIGHTS(src_cell_areas[src_points[offset]], 1.0)
      break;
    case(YAC_INTERP_SPMAP_INVTGTAREA):
      SCALE_WEIGHTS(1.0, tgt_cell_areas[tgt_points[offset]])
      break;
    case(YAC_INTERP_SPMAP_FRACAREA):
      SCALE_WEIGHTS(
        src_cell_areas[src_points[offset]], tgt_cell_areas[tgt_points[offset]])
      break;
  }

  free((void*)src_cell_areas);
  free((void*)tgt_cell_areas);
}

/**
 * Applies spreading of source data based on initial search results
 * (single target per source)
 * @param[in]     interp_grid
 * @param[in]     spmap_config
 * @param[in]     num_src_points     number of source points
 * @param[in,out] src_points_        local ids of all source points, will be
 *                                   inflated to the total number of links and
 *                                   contain duplication of the source points
 *                                   (each source points is repeated by the
 *                                   number of links associated with it)
 * @param[in,out] tgt_result_points_ local ids of initial target search results
 *                                   (one target per source), will be updated
 *                                   by this routine and then contain multiple
 *                                   targets per source
 * @param[out]    weights_           weights for all links
 * @param[out]    total_num_weights_ total number of links
 */
static void spread_src_data(
  struct yac_interp_grid * interp_grid,
  struct yac_interp_spmap_config const * spmap_config,
  size_t num_src_points, size_t ** src_points_,
  size_t ** tgt_result_points_, double ** weights_,
  size_t * total_num_weights_) {

  // shortcut in case there is only a single "1.0" weight for all source points
  if ((spmap_config->spread_distance <= 0) &&
      (spmap_config->scale_config->type == YAC_INTERP_SPMAP_NONE)) {

    *weights_ = xmalloc(num_src_points * sizeof(**weights_));
    for (size_t i = 0; i < num_src_points; ++i) (*weights_)[i] = 1.0;
    *total_num_weights_ = num_src_points;
    return;
  }

  size_t * src_points = *src_points_;
  size_t * tgt_result_points = *tgt_result_points_;

  size_t * num_tgt_per_src =
    xmalloc(num_src_points * sizeof(*num_tgt_per_src));
  size_t total_num_weights;

  // search for additional target points if spread distance is bigger than 0.0
  if (spmap_config->spread_distance > 0.0) {

    double sin_spread_distance, cos_spread_distance;
    compute_sin_cos(
      spmap_config->spread_distance, &sin_spread_distance, &cos_spread_distance);
    struct sin_cos_angle spread_distance =
      sin_cos_angle_new(sin_spread_distance, cos_spread_distance);
    yac_const_coordinate_pointer tgt_field_coords =
      yac_interp_grid_get_tgt_field_coords(interp_grid);

    struct bounding_circle * search_bnd_circles =
      xmalloc(num_src_points * sizeof(*search_bnd_circles));
    for (size_t i = 0; i < num_src_points; ++i) {
      memcpy(
        search_bnd_circles[i].base_vector,
        tgt_field_coords[tgt_result_points[i]], sizeof(*tgt_field_coords));
      search_bnd_circles[i].inc_angle = spread_distance;
      search_bnd_circles[i].sq_crd = DBL_MAX;
    }
    size_t * spread_tgt_result_points = NULL;

    // do bounding circle search around found tgt points
    // (ensures that all required target points are locally available)
    yac_interp_grid_do_bnd_circle_search_tgt(
      interp_grid, search_bnd_circles, num_src_points,
      &spread_tgt_result_points, num_tgt_per_src);
    free(search_bnd_circles);

    // remove target points which exceed the spread distance and only keep
    // targets that are connected to the original target point or other
    // target that have already been selected
    total_num_weights =
      check_tgt_result_points(
        interp_grid, spread_distance, num_src_points,
        tgt_result_points, num_tgt_per_src, spread_tgt_result_points);
    free((void*)tgt_result_points);
    tgt_result_points =
      xrealloc(
        spread_tgt_result_points,
        total_num_weights * sizeof(*spread_tgt_result_points));

    // adjust src_points (one source per target)
    size_t * new_src_points =
      xmalloc(total_num_weights * sizeof(*new_src_points));
    for (size_t i = 0, offset = 0; i < num_src_points; ++i)
      for (size_t j = 0, curr_src_point = src_points[i];
            j < num_tgt_per_src[i]; ++j, ++offset)
        new_src_points[offset] = curr_src_point;
    free((void*)src_points);
    src_points = new_src_points;

  } else {

    for (size_t i = 0; i < num_src_points; ++i) num_tgt_per_src[i] = 1;
    total_num_weights = num_src_points;
  }

  // compute weights
  double * weights =
    compute_weights(
      interp_grid, spmap_config->weight_type, num_src_points, src_points,
      num_tgt_per_src, total_num_weights, tgt_result_points);

  // scale weights
  scale_weights(
    interp_grid, spmap_config->scale_config, num_src_points, src_points,
    num_tgt_per_src, tgt_result_points, total_num_weights, weights);

  free(num_tgt_per_src);

  // set return values
  *tgt_result_points_ = tgt_result_points;
  *src_points_ = src_points;
  *weights_ = weights;
  *total_num_weights_ = total_num_weights;
}

/**
 * Applies spmap to a set of given source points for a single spmap_config
 */
static void do_search_spmap_(
  struct yac_interp_grid * interp_grid,
  struct yac_interp_spmap_config const * config,
  size_t * src_points, size_t num_src_points, yac_coordinate_pointer src_coords,
  InterpLink ** combined_result, size_t * combined_result_count) {

  // search for matching tgt points
  size_t * tgt_result_points =
    xmalloc(num_src_points * sizeof(*tgt_result_points));
  yac_interp_grid_do_nnn_search_tgt(
    interp_grid, src_coords, num_src_points, 1, tgt_result_points,
    (config->max_search_distance == 0.0)?M_PI:config->max_search_distance);

  // remove source points for which matching target point was found
  {
    size_t new_num_src_points = 0;
    for (size_t i = 0; i < num_src_points; ++i) {
      if (tgt_result_points[i] != SIZE_MAX) {
        if (i != new_num_src_points) {
          src_points[new_num_src_points] = src_points[i];
          tgt_result_points[new_num_src_points] =
            tgt_result_points[i];
        }
        ++new_num_src_points;
      }
    }
    num_src_points = new_num_src_points;
  }

  // spread the data from each source point to multiple target points
  double * weight_data;
  size_t total_num_weights;
  spread_src_data(
    interp_grid, config, num_src_points, &src_points, &tgt_result_points,
    &weight_data, &total_num_weights);

  // relocate source-target-point-pairs to dist owners of the respective
  // target points (these processes have the required target points in their
  // local list of target points that have to be interpolated)
  size_t result_count = total_num_weights;
  int to_tgt_owner = 1;
  yac_interp_grid_relocate_src_tgt_pairs(
    interp_grid, to_tgt_owner,
    0, &src_points, &tgt_result_points, &weight_data, &result_count);
  total_num_weights = result_count;

  // add interpolation links (src_point, tgt_point, weight) to combined list of
  // results
  *combined_result =
    xrealloc(
      *combined_result,
      (*combined_result_count + result_count) * sizeof(**combined_result));

  struct yac_const_basic_grid_data * src_basic_grid_data =
    yac_interp_grid_get_basic_grid_data_src(interp_grid);
  yac_int const * src_global_ids = src_basic_grid_data->ids[YAC_LOC_CELL];
  InterpLink * temp_result = *combined_result + *combined_result_count;
  for (size_t i = 0; i < result_count; ++i) {
    temp_result[i].src_point.local = src_points[i];
    temp_result[i].src_point.global = src_global_ids[src_points[i]];
    temp_result[i].tgt_point = tgt_result_points[i];
    temp_result[i].weight = weight_data[i];
  }
  *combined_result_count += result_count;

  free(tgt_result_points);
  free(weight_data);
  free(src_points);
}

static void src_point_selection_apply(
  struct yac_point_selection const * src_point_selection,
  size_t * src_points, size_t num_src_points,
  yac_coordinate_pointer src_coords,
  size_t ** selected_src_points, size_t * num_selected_src_points,
  yac_coordinate_pointer * selected_src_coords) {

  // get all matching source points
  yac_point_selection_apply(
    src_point_selection, src_coords, src_points, num_src_points,
    num_selected_src_points);

  // get the coordinates and point indices of the selected source points
  *selected_src_points =
    xmalloc(*num_selected_src_points * sizeof(**selected_src_points));
  memcpy(
    *selected_src_points, src_points + num_src_points - *num_selected_src_points,
    *num_selected_src_points * sizeof(**selected_src_points));
  *selected_src_coords = src_coords + num_src_points - *num_selected_src_points;
}

/**
 * Compare routine for objects of type InterpLink to be used by std::qsort
 */
static int compare_interp_link_tgt_point(void const * a, void const * b) {

  InterpLink const * link_a = (InterpLink const *)a;
  InterpLink const * link_b = (InterpLink const *)b;

  int ret =
    (link_a->tgt_point > link_b->tgt_point) -
    (link_a->tgt_point < link_b->tgt_point);
  if (ret) return ret;

  ret =
    (link_a->src_point.global > link_b->src_point.global) -
    (link_a->src_point.global < link_b->src_point.global);
  if (ret) return ret;

  return
    (link_a->weight > link_b->weight) -
    (link_a->weight < link_b->weight);
}

/// @brief basic routine for the computation of the interpolation weights
/// @param[in]     method                 abstract pointer to interpolation
///                                       method data
/// @param[in,out] interp_grid            interpolation grids information
/// @param[in,out] tgt_points             local indices of all target points
///                                       that are to be interpolated
/// @param[in]     tgt_point_count        number of entries in tgt_points
/// @param[in,out] weights                data structure containing the weights
/// @param[in]     interpolation_complete if it is one, no weights are to be
///                                       computed
/// @return number of target points for which an interpolation was generated
/// @remark entries in tgt_points are reordered such that successfully
///         interpolated target points are moved to the beginning of the array
static size_t do_search_spmap (struct interp_method * method,
                               struct yac_interp_grid * interp_grid,
                               size_t * tgt_points, size_t tgt_point_count,
                               struct yac_interp_weights * weights,
                               int * interpolation_complete) {

  if (*interpolation_complete) return 0;

  YAC_ASSERT(
    yac_interp_grid_get_num_src_fields(interp_grid) == 1,
    "ERROR(do_search_spmap): invalid number of source fields")
  YAC_ASSERT(
    yac_interp_grid_get_src_field_location(interp_grid, 0) == YAC_LOC_CELL,
    "ERROR(do_search_spmap): "
    "invalid source field location (has to be YAC_LOC_CELL)")
  YAC_ASSERT(
    yac_interp_grid_get_tgt_field_location(interp_grid) == YAC_LOC_CELL,
    "ERROR(do_search_spmap): "
    "invalid target field location (has to be YAC_LOC_CELL)")

  // get coordinates of all source points
  size_t * src_points;
  size_t num_src_points;
  yac_interp_grid_get_src_points(
    interp_grid, 0, &src_points, &num_src_points);
  yac_coordinate_pointer src_coords = xmalloc(num_src_points * sizeof(*src_coords));
  yac_interp_grid_get_src_coordinates(
    interp_grid, src_points, num_src_points, 0, src_coords);

  struct yac_spmap_overwrite_config const * const * overwrite_configs =
    (struct yac_spmap_overwrite_config const * const *)
      ((struct interp_method_spmap*)method)->overwrite_configs;

  InterpLink * combined_result = NULL;
  size_t combined_result_count = 0;

  // for all alternative configurations
  for (size_t i = 0;
       (overwrite_configs != NULL) && (overwrite_configs[i] != NULL); ++i) {

    // extract source points matching the current criteria
    size_t * selected_src_points;
    size_t num_selected_src_points;
    yac_coordinate_pointer selected_src_coords;
    src_point_selection_apply(
      overwrite_configs[i]->src_point_selection,
      src_points, num_src_points, src_coords,
      &selected_src_points, &num_selected_src_points, &selected_src_coords);
    num_src_points -= num_selected_src_points;

    // apply source to target mapping to selected source points
    do_search_spmap_(
      interp_grid, overwrite_configs[i]->config, selected_src_points,
      num_selected_src_points, selected_src_coords,
      &combined_result, &combined_result_count);
  }

  // apply default configuration to remaining source points
  do_search_spmap_(
    interp_grid, ((struct interp_method_spmap*)method)->default_config,
    src_points, num_src_points, src_coords, &combined_result,
    &combined_result_count);
  free(src_coords);

  // sort source-target-point-pairs by target points
  qsort(
    combined_result, combined_result_count,
    sizeof(*combined_result), compare_interp_link_tgt_point);

  // generate num_src_per_tgt and determine list of unique target points
  size_t * num_src_per_tgt =
    xmalloc(combined_result_count * sizeof(*num_src_per_tgt));
  size_t * unique_result_tgt_points =
    xmalloc(combined_result_count * sizeof(*unique_result_tgt_points));
  size_t num_unique_result_tgt_points = 0;
  for (size_t i = 0; i < combined_result_count;) {
    size_t prev_i = i;
    size_t curr_tgt = combined_result[i].tgt_point;
    while (
      (i < combined_result_count) &&
      (curr_tgt == combined_result[i].tgt_point)) ++i;
    num_src_per_tgt[num_unique_result_tgt_points] = i - prev_i;
    unique_result_tgt_points[num_unique_result_tgt_points] = curr_tgt;
    ++num_unique_result_tgt_points;
  }
  num_src_per_tgt =
    xrealloc(
      num_src_per_tgt, num_unique_result_tgt_points * sizeof(*num_src_per_tgt));
  unique_result_tgt_points =
    xrealloc(
      unique_result_tgt_points,
      num_unique_result_tgt_points * sizeof(*unique_result_tgt_points));

  // match unique_result_tgt_points with target points that are supposed to
  // be interpolated
  qsort(tgt_points, tgt_point_count, sizeof(*tgt_points), compare_size_t);
  int * unused_tgt_point_flag =
    xmalloc(tgt_point_count * sizeof(*unused_tgt_point_flag));
  {
    size_t j = 0;
    for (size_t i = 0; i < num_unique_result_tgt_points; ++i) {
      size_t curr_result_tgt = unique_result_tgt_points[i];
      while ((j < tgt_point_count) && (tgt_points[j] < curr_result_tgt)) {
        unused_tgt_point_flag[j++] = 1;
      }
      YAC_ASSERT(
        (j < tgt_point_count) && (curr_result_tgt == tgt_points[j]),
        "ERROR(do_search_spmap): "
        "required target points already in use or not available")
      unused_tgt_point_flag[j++] = 0;
    }
    for (; j < tgt_point_count; ++j) unused_tgt_point_flag[j] = 1;
  }

  // sort used target points to the beginning of the array
  yac_flag_sort_size_t(
    tgt_points, unused_tgt_point_flag, num_unique_result_tgt_points);
  free(unused_tgt_point_flag);

  size_t * combined_result_src_points =
    xmalloc(combined_result_count * sizeof(*combined_result_src_points));
  double * combined_result_weights =
    xmalloc(combined_result_count * sizeof(*combined_result_weights));
  for (size_t i = 0; i < combined_result_count; ++i) {
    combined_result_src_points[i] = combined_result[i].src_point.local;
    combined_result_weights[i] = combined_result[i].weight;
  }
  free(combined_result);

  struct remote_point * srcs =
    yac_interp_grid_get_src_remote_points(
      interp_grid, 0, combined_result_src_points, combined_result_count);
  struct remote_points tgts = {
    .data =
      yac_interp_grid_get_tgt_remote_points(
        interp_grid, unique_result_tgt_points, num_unique_result_tgt_points),
    .count = num_unique_result_tgt_points};
  free(unique_result_tgt_points);

  // store results
  yac_interp_weights_add_wsum(
    weights, &tgts, num_src_per_tgt, srcs, combined_result_weights);

  // cleanup
  free(combined_result_src_points);
  free(combined_result_weights);
  free(tgts.data);
  free(srcs);
  free(num_src_per_tgt);

  return num_unique_result_tgt_points;
}

static struct yac_spmap_cell_area_config * yac_spmap_cell_area_config_copy(
  struct yac_spmap_cell_area_config const * cell_area_config,
  char const * desc) {

  if (cell_area_config == NULL) cell_area_config = &cell_area_config_default;

  struct yac_spmap_cell_area_config * copy;

  switch (cell_area_config->type) {
    YAC_UNREACHABLE_DEFAULT_F(
      "ERROR(yac_spmap_cell_area_config_copy): "
      "invalid %s cell area provider type", desc);
    case(YAC_INTERP_SPMAP_CELL_AREA_YAC): {
      copy =
        yac_spmap_cell_area_config_yac_new(cell_area_config->yac.sphere_radius);
      break;
    }
    case(YAC_INTERP_SPMAP_CELL_AREA_FILE): {
      copy =
        yac_spmap_cell_area_config_file_new(
          cell_area_config->file_config.filename,
          cell_area_config->file_config.varname,
          cell_area_config->file_config.min_global_id);
      break;
    }
  }

  return copy;
}

static struct yac_spmap_scale_config * yac_spmap_scale_config_copy(
  struct yac_spmap_scale_config const * scale_config) {

  if (scale_config == NULL) scale_config = &spmap_scale_config_default;

  return
    yac_spmap_scale_config_new(
      scale_config->type, scale_config->src, scale_config->tgt);
}

struct yac_interp_spmap_config * yac_interp_spmap_config_copy(
  struct yac_interp_spmap_config const * spmap_config) {

  if (spmap_config == NULL) spmap_config = &spmap_config_default;

  return
    yac_interp_spmap_config_new(
      spmap_config->spread_distance, spmap_config->max_search_distance,
      spmap_config->weight_type, spmap_config->scale_config);
}

struct yac_spmap_overwrite_config * yac_spmap_overwrite_config_copy(
  struct yac_spmap_overwrite_config const * overwrite_config) {

  if (overwrite_config == NULL) overwrite_config = &overwrite_config_default;

  struct yac_spmap_overwrite_config * copy = xmalloc(1 * sizeof(*copy));

  copy->src_point_selection =
    yac_point_selection_copy(overwrite_config->src_point_selection);
  copy->config = yac_interp_spmap_config_copy(overwrite_config->config);

  return copy;
}

struct yac_spmap_overwrite_config ** yac_spmap_overwrite_configs_copy(
  struct yac_spmap_overwrite_config const * const * overwrite_configs) {

  size_t overwrite_config_count = 0;

  for (;(overwrite_configs != NULL) &&
        (overwrite_configs[overwrite_config_count] != NULL);
       ++overwrite_config_count);

  struct yac_spmap_overwrite_config ** copy =
    (overwrite_config_count > 0)?
      xcalloc(overwrite_config_count + 1, sizeof(*copy)):NULL;

  for (size_t i = 0; i < overwrite_config_count; ++i)
    copy[i] = yac_spmap_overwrite_config_copy(overwrite_configs[i]);

  return copy;
}

struct interp_method * yac_interp_method_spmap_new(
  struct yac_interp_spmap_config const * default_config,
  struct yac_spmap_overwrite_config const * const * overwrite_configs) {

  struct interp_method_spmap * method = xmalloc(1 * sizeof(*method));

  method->vtable = &interp_method_spmap_vtable;
  method->default_config = yac_interp_spmap_config_copy(default_config);
  method->overwrite_configs =
    yac_spmap_overwrite_configs_copy(overwrite_configs);

  return (struct interp_method*)method;
}

struct yac_spmap_cell_area_config * yac_spmap_cell_area_config_yac_new(
  double sphere_radius) {

  YAC_ASSERT_F(
    sphere_radius > 0.0,
    "ERROR(yac_spmap_cell_area_config_yac_new): "
    "invalid sphere_radius %lf (has to be >= 0.0)", sphere_radius);

  struct yac_spmap_cell_area_config * cell_area_config =
    xmalloc(1 * sizeof(*cell_area_config));

  cell_area_config->type = YAC_INTERP_SPMAP_CELL_AREA_YAC;
  cell_area_config->yac.sphere_radius = sphere_radius;

  return cell_area_config;
}

struct yac_spmap_cell_area_config * yac_spmap_cell_area_config_file_new(
  char const * filename, char const * varname, yac_int min_global_id) {

  YAC_ASSERT(
    (filename != NULL) && (strlen(filename) > 0),
    "ERROR(yac_spmap_cell_area_config_file_new): invalid filename for areas");
  YAC_ASSERT(
    (varname != NULL) && (strlen(varname) > 0),
    "ERROR(yac_spmap_cell_area_config_file_new): invalid varname for areas");

  struct yac_spmap_cell_area_config * cell_area_config =
    xmalloc(1 * sizeof(*cell_area_config));

  cell_area_config->type = YAC_INTERP_SPMAP_CELL_AREA_FILE;
  cell_area_config->file_config.filename = strdup(filename);
  cell_area_config->file_config.varname = strdup(varname);
  cell_area_config->file_config.min_global_id = min_global_id;

  return cell_area_config;
}

struct yac_spmap_cell_area_config * yac_spmap_cell_area_config_file_new_f2c(
  char const * filename, char const * varname, int min_global_id) {

  YAC_ASSERT_F(
    (min_global_id >= XT_INT_MIN) && (min_global_id <= XT_INT_MAX),
    "ERROR(yac_spmap_cell_area_config_file_new_f2c): "
    "min_global_id %d is outside of the valid range "
    "(%" XT_INT_FMT " <= min_global_id <= %" XT_INT_FMT ")",
    min_global_id, XT_INT_MIN, XT_INT_MAX);

  return
    yac_spmap_cell_area_config_file_new(
      filename, varname, (yac_int)min_global_id);
}

void yac_spmap_cell_area_config_delete(
  struct yac_spmap_cell_area_config * cell_area_config) {

  if ((cell_area_config != &cell_area_config_default) &&
      (cell_area_config != NULL)) {

    switch(cell_area_config->type) {
      YAC_UNREACHABLE_DEFAULT(
        "ERROR(yac_spmap_cell_area_config_delete): "
        "invalid cell area configuration type");
      case(YAC_INTERP_SPMAP_CELL_AREA_YAC):
        // nothing to be done
        break;
      case(YAC_INTERP_SPMAP_CELL_AREA_FILE): {
        free((void*)cell_area_config->file_config.filename);
        free((void*)cell_area_config->file_config.varname);
        break;
      }
    }
    free(cell_area_config);
  }
}

enum yac_interp_spmap_cell_area_provider yac_spmap_cell_area_config_get_type(
  struct yac_spmap_cell_area_config const * cell_area_config) {

  if (cell_area_config == NULL) cell_area_config = &cell_area_config_default;

  return cell_area_config->type;
}

double yac_spmap_cell_area_config_get_sphere_radius(
  struct yac_spmap_cell_area_config const * cell_area_config) {

  if (cell_area_config == NULL) cell_area_config = &cell_area_config_default;

  YAC_ASSERT_F(
    cell_area_config->type == YAC_INTERP_SPMAP_CELL_AREA_YAC,
    "ERROR(yac_spmap_cell_area_config_get_sphere_radius): "
    "invalid cell area configuration type %d "
    "(has to be YAC_INTERP_SPMAP_CELL_AREA_YAC)", (int)cell_area_config->type);

  return cell_area_config->yac.sphere_radius;
}

char const * yac_spmap_cell_area_config_get_filename(
  struct yac_spmap_cell_area_config const * cell_area_config) {

  if (cell_area_config == NULL) cell_area_config = &cell_area_config_default;

  YAC_ASSERT_F(
    cell_area_config->type == YAC_INTERP_SPMAP_CELL_AREA_FILE,
    "ERROR(yac_spmap_cell_area_config_get_filename): "
    "invalid cell area configuration type %d "
    "(has to be YAC_INTERP_SPMAP_CELL_AREA_FILE)", (int)cell_area_config->type);

  return cell_area_config->file_config.filename;
}

char const * yac_spmap_cell_area_config_get_varname(
  struct yac_spmap_cell_area_config const * cell_area_config) {

  if (cell_area_config == NULL) cell_area_config = &cell_area_config_default;

  YAC_ASSERT_F(
    cell_area_config->type == YAC_INTERP_SPMAP_CELL_AREA_FILE,
    "ERROR(yac_spmap_cell_area_config_get_varname): "
    "invalid cell area configuration type %d "
    "(has to be YAC_INTERP_SPMAP_CELL_AREA_FILE)", (int)cell_area_config->type);

  return cell_area_config->file_config.varname;
}

yac_int yac_spmap_cell_area_config_get_min_global_id(
  struct yac_spmap_cell_area_config const * cell_area_config) {

  if (cell_area_config == NULL) cell_area_config = &cell_area_config_default;

  YAC_ASSERT_F(
    cell_area_config->type == YAC_INTERP_SPMAP_CELL_AREA_FILE,
    "ERROR(yac_spmap_cell_area_config_get_min_global_id): "
    "invalid cell area configuration type %d "
    "(has to be YAC_INTERP_SPMAP_CELL_AREA_FILE)", (int)cell_area_config->type);

  return cell_area_config->file_config.min_global_id;
}

struct yac_spmap_scale_config * yac_spmap_scale_config_new(
  enum yac_interp_spmap_scale_type scale_type,
  struct yac_spmap_cell_area_config const * source_cell_area_config,
  struct yac_spmap_cell_area_config const * target_cell_area_config) {

  struct yac_spmap_scale_config * scale_config =
    xmalloc(1 * sizeof(*scale_config));

  YAC_ASSERT_F(
    (scale_type == YAC_INTERP_SPMAP_NONE) ||
    (scale_type == YAC_INTERP_SPMAP_SRCAREA) ||
    (scale_type == YAC_INTERP_SPMAP_INVTGTAREA) ||
    (scale_type == YAC_INTERP_SPMAP_FRACAREA),
    "ERROR(yac_spmap_scale_config_new): "
    "invalid scale configuration type (%d)", (int)scale_type);

  scale_config->type = scale_type;
  scale_config->src =
    yac_spmap_cell_area_config_copy(source_cell_area_config, "source");
  scale_config->tgt =
    yac_spmap_cell_area_config_copy(target_cell_area_config, "target");

  return scale_config;
}

struct yac_spmap_scale_config * yac_spmap_scale_config_new_f2c(
  int scale_type,
  struct yac_spmap_cell_area_config * source_cell_area_config,
  struct yac_spmap_cell_area_config * target_cell_area_config) {

  YAC_ASSERT_F(
    (scale_type == YAC_INTERP_SPMAP_NONE) ||
    (scale_type == YAC_INTERP_SPMAP_SRCAREA) ||
    (scale_type == YAC_INTERP_SPMAP_INVTGTAREA) ||
    (scale_type == YAC_INTERP_SPMAP_FRACAREA),
    "ERROR(yac_spmap_scale_config_new_f2c): "
    "invalid scale configuration type (%d)", scale_type);

  return
    yac_spmap_scale_config_new(
      (enum yac_interp_spmap_scale_type)scale_type,
      source_cell_area_config, target_cell_area_config);
}

void yac_spmap_scale_config_delete(
  struct yac_spmap_scale_config * scale_config) {

  if ((scale_config != &spmap_scale_config_default) &&
      (scale_config != NULL)) {
    yac_spmap_cell_area_config_delete(scale_config->src);
    yac_spmap_cell_area_config_delete(scale_config->tgt);
    free(scale_config);
  }
}

enum yac_interp_spmap_scale_type yac_spmap_scale_config_get_type(
  struct yac_spmap_scale_config const * scale_config) {

  if (scale_config == NULL) scale_config = &spmap_scale_config_default;

  return scale_config->type;
}

struct yac_spmap_cell_area_config const *
  yac_spmap_scale_config_get_src_cell_area_config(
    struct yac_spmap_scale_config const * scale_config) {

  if (scale_config == NULL) scale_config = &spmap_scale_config_default;

  return scale_config->src;
}

struct yac_spmap_cell_area_config const *
  yac_spmap_scale_config_get_tgt_cell_area_config(
    struct yac_spmap_scale_config const * scale_config) {

  if (scale_config == NULL) scale_config = &spmap_scale_config_default;

  return scale_config->tgt;
}

struct yac_interp_spmap_config * yac_interp_spmap_config_new(
  double spread_distance, double max_search_distance,
  enum yac_interp_spmap_weight_type weight_type,
  struct yac_spmap_scale_config const * scale_config) {

  struct yac_interp_spmap_config * spmap_config =
    xmalloc(1 * sizeof(*spmap_config));

  YAC_ASSERT_F(
    (spread_distance >= 0.0) && (spread_distance <= M_PI_2),
    "ERROR(yac_interp_spmap_config_new): invalid spread_distance "
    "(has to be >= 0 and <= PI/2 (%lf)", spread_distance);

  YAC_ASSERT_F(
    (max_search_distance >= 0.0) && (max_search_distance <= M_PI),
    "ERROR(yac_interp_spmap_config_new): invalid max_search_distance "
    "(has to be >= 0 and <= PI (%lf)", max_search_distance);

  YAC_ASSERT_F(
    (weight_type == YAC_INTERP_SPMAP_AVG) ||
    (weight_type == YAC_INTERP_SPMAP_DIST),
    "ERROR(yac_interp_spmap_config_new): invalid weight type (%d)",
    (int)weight_type);

  spmap_config->spread_distance = spread_distance;
  spmap_config->max_search_distance = max_search_distance;
  spmap_config->weight_type = weight_type;
  spmap_config->scale_config = yac_spmap_scale_config_copy(scale_config);

  return spmap_config;
}

struct yac_interp_spmap_config * yac_interp_spmap_config_new_f2c(
  double spread_distance, double max_search_distance,
  int weight_type, struct yac_spmap_scale_config * scale_config) {

  YAC_ASSERT_F(
    (weight_type == YAC_INTERP_SPMAP_AVG) ||
    (weight_type == YAC_INTERP_SPMAP_DIST),
    "ERROR(yac_interp_spmap_config_new_f2c): invalid weight type (%d)",
    weight_type);

  return
    yac_interp_spmap_config_new(
      spread_distance, max_search_distance,
      (enum yac_interp_spmap_weight_type)weight_type, scale_config);
}

void yac_interp_spmap_config_delete(
  struct yac_interp_spmap_config * config) {

  if ((config != &spmap_config_default) && (config != NULL)) {
    yac_spmap_scale_config_delete(config->scale_config);
    free(config);
  }
}

double yac_interp_spmap_config_get_spread_distance(
  struct yac_interp_spmap_config const * spmap_config) {

  if (spmap_config == NULL) spmap_config = &spmap_config_default;

  return spmap_config->spread_distance;
}

double yac_interp_spmap_config_get_max_search_distance(
  struct yac_interp_spmap_config const * spmap_config) {

  if (spmap_config == NULL) spmap_config = &spmap_config_default;

  return spmap_config->max_search_distance;
}

enum yac_interp_spmap_weight_type yac_interp_spmap_config_get_weight_type(
  struct yac_interp_spmap_config const * spmap_config) {

  if (spmap_config == NULL) spmap_config = &spmap_config_default;

  return spmap_config->weight_type;
}

struct yac_spmap_scale_config const * yac_interp_spmap_config_get_scale_config(
  struct yac_interp_spmap_config const * spmap_config) {

  if (spmap_config == NULL) spmap_config = &spmap_config_default;

  return spmap_config->scale_config;
}

struct yac_spmap_overwrite_config * yac_spmap_overwrite_config_new(
  struct yac_point_selection const * src_point_selection,
  struct yac_interp_spmap_config const * config) {

  struct yac_spmap_overwrite_config * overwrite_config =
    xmalloc(1 * sizeof(*overwrite_config));

  overwrite_config->src_point_selection =
    yac_point_selection_copy(src_point_selection);
  overwrite_config->config = yac_interp_spmap_config_copy(config);

  return overwrite_config;
}

void yac_spmap_overwrite_config_delete(
  struct yac_spmap_overwrite_config * overwrite_config) {

  if ((overwrite_config != &overwrite_config_default) &&
      (overwrite_config != NULL)) {
    yac_point_selection_delete(overwrite_config->src_point_selection);
    yac_interp_spmap_config_delete(overwrite_config->config);
    free(overwrite_config);
  }
}

void yac_spmap_overwrite_configs_delete(
  struct yac_spmap_overwrite_config ** overwrite_configs) {

  for (size_t i = 0;
        (overwrite_configs != NULL) && (overwrite_configs[i] != NULL); ++i)
    yac_spmap_overwrite_config_delete(overwrite_configs[i]);
  free(overwrite_configs);
}

struct yac_point_selection const *
  yac_spmap_overwrite_config_get_src_point_selection(
    struct yac_spmap_overwrite_config const * overwrite_config) {

  if (overwrite_config == NULL) overwrite_config = &overwrite_config_default;

  return overwrite_config->src_point_selection;
}

struct yac_interp_spmap_config const *
  yac_spmap_overwrite_config_get_spmap_config(
    struct yac_spmap_overwrite_config const * overwrite_config) {

  if (overwrite_config == NULL) overwrite_config = &overwrite_config_default;

  return overwrite_config->config;
}

static void delete_spmap(struct interp_method * method) {

  struct interp_method_spmap * method_spmap =
    (struct interp_method_spmap*)(method);

  yac_interp_spmap_config_delete(method_spmap->default_config);
  yac_spmap_overwrite_configs_delete(method_spmap->overwrite_configs);

  free(method);
}

static int yac_spmap_cell_area_config_compare(
  struct yac_spmap_cell_area_config const * a,
  struct yac_spmap_cell_area_config const * b) {

  int ret;
  if ((ret = ((int)(a->type) > (int)(b->type)) -
             ((int)(a->type) < (int)(b->type))))
    return ret;

  switch (a->type) {
    default:
    case (YAC_INTERP_SPMAP_CELL_AREA_YAC): {
      return
        ((int)(a->yac.sphere_radius) > (int)(b->yac.sphere_radius)) -
        ((int)(a->yac.sphere_radius) < (int)(b->yac.sphere_radius));
    }
    case (YAC_INTERP_SPMAP_CELL_AREA_FILE): {
      if ((ret = strcmp(a->file_config.filename, b->file_config.filename)))
        return ret;
      if ((ret = strcmp(a->file_config.varname, b->file_config.varname)))
        return ret;
      return
        ((int)(a->file_config.min_global_id) >
         (int)(b->file_config.min_global_id)) -
        ((int)(a->file_config.min_global_id) <
         (int)(b->file_config.min_global_id));
    }
  }
}

int yac_spmap_scale_config_compare(
  struct yac_spmap_scale_config const * a,
  struct yac_spmap_scale_config const * b) {

  int ret;
  if ((ret = ((int)(a->type) > (int)(b->type)) -
             ((int)(a->type) < (int)(b->type))))
    return ret;
  if ((ret = yac_spmap_cell_area_config_compare(a->src, b->src))) return ret;
  if ((ret = yac_spmap_cell_area_config_compare(a->tgt, b->tgt))) return ret;
  return 0;
}

int yac_interp_spmap_config_compare(
  struct yac_interp_spmap_config const * a,
  struct yac_interp_spmap_config const * b) {
  if (fabs(a->spread_distance - b->spread_distance) > yac_angle_tol)
    return (a->spread_distance > b->spread_distance) -
           (a->spread_distance < b->spread_distance);
  if (fabs(a->max_search_distance - b->max_search_distance) > yac_angle_tol)
    return (a->max_search_distance > b->max_search_distance) -
           (a->max_search_distance < b->max_search_distance);
  if (a->weight_type != b->weight_type)
    return
      (a->weight_type > b->weight_type) - (a->weight_type < b->weight_type);
  return yac_spmap_scale_config_compare(a->scale_config, b->scale_config);
}

int yac_spmap_overwrite_config_compare(
  struct yac_spmap_overwrite_config const * a,
  struct yac_spmap_overwrite_config const * b) {

  int ret =
    yac_point_selection_compare(a->src_point_selection, b->src_point_selection);
  if (ret) return ret;

  return yac_interp_spmap_config_compare(a->config, b->config);
}

static size_t yac_interp_stack_config_get_string_pack_size(
  char const * string, MPI_Comm comm) {

  int strlen_pack_size, string_pack_size;
  yac_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &strlen_pack_size), comm);

  YAC_ASSERT(
    string != NULL, "ERROR(yac_interp_stack_config_get_string_pack_size): "
    "string is NULL");

  yac_mpi_call(
    MPI_Pack_size(
      (int)(strlen(string)), MPI_CHAR, comm, &string_pack_size), comm);

  return (size_t)strlen_pack_size + (size_t)string_pack_size;
}

static size_t yac_spmap_cell_area_config_get_pack_size(
  struct yac_spmap_cell_area_config const * cell_area_config, MPI_Comm comm) {

  int int_pack_size, dbl_pack_size;
  yac_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &int_pack_size), comm);
  yac_mpi_call(MPI_Pack_size(1, MPI_DOUBLE, comm, &dbl_pack_size), comm);

  size_t pack_size = (size_t)int_pack_size;  // cell_area_provider

  if (cell_area_config == NULL) cell_area_config = &cell_area_config_default;

  switch (cell_area_config->type) {
    YAC_UNREACHABLE_DEFAULT(
      "ERROR(yac_spmap_cell_area_config_get_pack_size): "
      "invalid cell area provider type");
    case (YAC_INTERP_SPMAP_CELL_AREA_YAC): {
      pack_size += (size_t)dbl_pack_size; // sphere_radius
      break;
    }
    case (YAC_INTERP_SPMAP_CELL_AREA_FILE): {
      pack_size +=
        yac_interp_stack_config_get_string_pack_size(
          cell_area_config->file_config.filename, comm) + // filename
        yac_interp_stack_config_get_string_pack_size(
          cell_area_config->file_config.varname, comm) + // varname
        int_pack_size;                                   // min_global_id
      break;
    }
  }

  return pack_size;
}

size_t yac_spmap_scale_config_get_pack_size(
  struct yac_spmap_scale_config const * scale_config, MPI_Comm comm) {

  int int_pack_size;
  yac_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &int_pack_size), comm);

  return (size_t)int_pack_size + // type
         yac_spmap_cell_area_config_get_pack_size(scale_config->src, comm) +
         yac_spmap_cell_area_config_get_pack_size(scale_config->tgt, comm);
}

size_t yac_interp_spmap_config_get_pack_size(
  struct yac_interp_spmap_config const * spmap_config, MPI_Comm comm) {

  int int_pack_size, dbl_pack_size;
  yac_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &int_pack_size), comm);
  yac_mpi_call(MPI_Pack_size(1, MPI_DOUBLE, comm, &dbl_pack_size), comm);

  return (size_t)dbl_pack_size + // spread_distance
         (size_t)dbl_pack_size + // max_search_distance
         (size_t)int_pack_size + // weight_type
         yac_spmap_scale_config_get_pack_size(spmap_config->scale_config, comm);
}

size_t yac_spmap_overwrite_config_get_pack_size(
  struct yac_spmap_overwrite_config const * overwrite_config, MPI_Comm comm) {

  return
    yac_point_selection_get_pack_size(
      overwrite_config->src_point_selection, comm) +
    yac_interp_spmap_config_get_pack_size(
      overwrite_config->config, comm);
}

size_t yac_spmap_overwrite_configs_get_pack_size(
  struct yac_spmap_overwrite_config const * const * overwrite_configs,
  MPI_Comm comm) {

  size_t over_write_configs_pack_size = 0;
  for (size_t i = 0;
       (overwrite_configs != NULL) && (overwrite_configs[i] != NULL); ++i)
    over_write_configs_pack_size +=
      yac_spmap_overwrite_config_get_pack_size(overwrite_configs[i], comm);

  int size_t_pack_size;
  yac_mpi_call(
    MPI_Pack_size(1, YAC_MPI_SIZE_T, comm, &size_t_pack_size), comm);

  return over_write_configs_pack_size +
         size_t_pack_size; // overwrite_config_count
}

static void yac_interp_stack_config_pack_string(
  char const * string, void * buffer, int buffer_size, int * position,
  MPI_Comm comm) {

  size_t len = (string == NULL)?0:strlen(string);

  YAC_ASSERT(
    len <= INT_MAX, "ERROR(yac_interp_stack_config_pack_string): string too long")

  int len_int = (int)len;

  yac_mpi_call(
    MPI_Pack(
      &len_int, 1, MPI_INT, buffer, buffer_size, position, comm), comm);

  if (len > 0)
    yac_mpi_call(
      MPI_Pack(
        string, len_int, MPI_CHAR, buffer, buffer_size, position, comm),
      comm);
}

static void yac_spmap_cell_area_config_pack(
  struct yac_spmap_cell_area_config const * cell_area_config,
  void * buffer, int buffer_size, int * position, MPI_Comm comm) {

  if (cell_area_config == NULL) cell_area_config = &cell_area_config_default;

  int type = (int)(cell_area_config->type);
  yac_mpi_call(
    MPI_Pack(&type, 1, MPI_INT, buffer, buffer_size, position, comm), comm);

  switch (cell_area_config->type) {
    YAC_UNREACHABLE_DEFAULT(
      "ERROR(yac_spmap_cell_area_config_pack): "
      "invalid cell area provider type");
    case(YAC_INTERP_SPMAP_CELL_AREA_YAC): {
      yac_mpi_call(
        MPI_Pack(
          &(cell_area_config->yac.sphere_radius), 1, MPI_DOUBLE,
          buffer, buffer_size, position, comm), comm);
      break;
    }
    case(YAC_INTERP_SPMAP_CELL_AREA_FILE): {
      yac_interp_stack_config_pack_string(
        cell_area_config->file_config.filename,
        buffer, buffer_size, position, comm);
      yac_interp_stack_config_pack_string(
        cell_area_config->file_config.varname,
        buffer, buffer_size, position, comm);
      YAC_ASSERT(
        (cell_area_config->file_config.min_global_id >= -INT_MAX) &&
        (cell_area_config->file_config.min_global_id <= INT_MAX),
        "ERRROR(yac_spmap_cell_area_config_pack): invalid minimum global id");
      int min_global_id = (int)(cell_area_config->file_config.min_global_id);
      yac_mpi_call(
        MPI_Pack(
          &min_global_id, 1, MPI_INT, buffer, buffer_size, position, comm),
        comm);
      break;
    }
  }
}

void yac_spmap_scale_config_pack(
  struct yac_spmap_scale_config const * scale_config,
  void * buffer, int buffer_size, int * position, MPI_Comm comm) {

  int scale_type = (int)(scale_config->type);
  yac_mpi_call(
    MPI_Pack(
      &scale_type, 1, MPI_INT, buffer, buffer_size, position, comm),
    comm);
  yac_spmap_cell_area_config_pack(
    scale_config->src, buffer, buffer_size, position, comm);
  yac_spmap_cell_area_config_pack(
    scale_config->tgt, buffer, buffer_size, position, comm);
}

void yac_interp_spmap_config_pack(
  struct yac_interp_spmap_config const * spmap_config,
  void * buffer, int buffer_size, int * position, MPI_Comm comm) {

  yac_mpi_call(
    MPI_Pack(
      &(spmap_config->spread_distance), 1, MPI_DOUBLE,
      buffer, buffer_size, position, comm), comm);
  yac_mpi_call(
    MPI_Pack(
      &(spmap_config->max_search_distance), 1, MPI_DOUBLE,
      buffer, buffer_size, position, comm), comm);
  int weight_type = (int)(spmap_config->weight_type);
  yac_mpi_call(
    MPI_Pack(
      &weight_type, 1, MPI_INT, buffer, buffer_size, position, comm),
    comm);
  yac_spmap_scale_config_pack(
    spmap_config->scale_config, buffer, buffer_size, position, comm);
}

static void yac_spmap_overwrite_config_pack(
  struct yac_spmap_overwrite_config const * overwrite_config,
  void * buffer, int buffer_size, int * position, MPI_Comm comm) {

  yac_point_selection_pack(
    overwrite_config->src_point_selection,
    buffer, buffer_size, position, comm);

  yac_interp_spmap_config_pack(
    overwrite_config->config, buffer, buffer_size, position, comm);
}

void yac_spmap_overwrite_configs_pack(
  struct yac_spmap_overwrite_config const * const * overwrite_configs,
  void * buffer, int buffer_size, int * position, MPI_Comm comm) {

  size_t overwrite_config_count = 0;
  for (; (overwrite_configs != NULL) &&
         (overwrite_configs[overwrite_config_count] != NULL);
       ++overwrite_config_count);

  yac_mpi_call(
    MPI_Pack(
      &overwrite_config_count, 1, YAC_MPI_SIZE_T, buffer, buffer_size,
      position, comm), comm);

  for (size_t i = 0; i < overwrite_config_count; ++i)
    yac_spmap_overwrite_config_pack(
      overwrite_configs[i], buffer, buffer_size, position, comm);
}

static char * yac_interp_stack_config_unpack_string(
  void * buffer, int buffer_size, int * position, MPI_Comm comm) {

  int string_len;
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &string_len, 1, MPI_INT, comm), comm);

  YAC_ASSERT(
    string_len >= 0,
    "ERROR(yac_interp_stack_config_unpack_string): invalid string length")

  char * string = NULL;
  if (string_len > 0) {
    string = xmalloc((size_t)(string_len + 1) * sizeof(*string));
    yac_mpi_call(
      MPI_Unpack(
        buffer, buffer_size, position, string, string_len, MPI_CHAR, comm),
      comm);
    string[string_len] = '\0';
  }
  return string;
}

static void yac_spmap_cell_area_config_unpack(
  void * buffer, int buffer_size, int * position,
  struct yac_spmap_cell_area_config ** cell_area_config, MPI_Comm comm) {

  *cell_area_config = xmalloc(1 * sizeof(**cell_area_config));

  int cell_area_provider;
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &cell_area_provider, 1, MPI_INT, comm),
    comm);
  (*cell_area_config)->type =
    (enum yac_interp_spmap_cell_area_provider)cell_area_provider;

  switch ((*cell_area_config)->type) {
    YAC_UNREACHABLE_DEFAULT(
      "ERROR(yac_spmap_cell_area_config_unpack): "
      "invalid cell area provider type");
    case(YAC_INTERP_SPMAP_CELL_AREA_YAC): {
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position,
          &((*cell_area_config)->yac.sphere_radius), 1, MPI_DOUBLE, comm),
        comm);
      break;
    }
    case(YAC_INTERP_SPMAP_CELL_AREA_FILE): {
      (*cell_area_config)->file_config.filename =
        yac_interp_stack_config_unpack_string(
          buffer, buffer_size, position, comm);
      (*cell_area_config)->file_config.varname =
        yac_interp_stack_config_unpack_string(
          buffer, buffer_size, position, comm);
      int min_global_id;
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position, &min_global_id, 1, MPI_INT, comm),
        comm);
      (*cell_area_config)->file_config.min_global_id = (size_t)min_global_id;
    }
  }
}

void yac_spmap_scale_config_unpack(
  void * buffer, int buffer_size, int * position,
  struct yac_spmap_scale_config ** scale_config, MPI_Comm comm) {

  *scale_config = xmalloc(1 * sizeof(**scale_config));
  int scale_type;
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &scale_type, 1, MPI_INT, comm),
    comm);
  (*scale_config)->type =
    (enum yac_interp_spmap_scale_type)scale_type;
  yac_spmap_cell_area_config_unpack(
    buffer, buffer_size, position, &((*scale_config)->src), comm);
  yac_spmap_cell_area_config_unpack(
    buffer, buffer_size, position, &((*scale_config)->tgt), comm);
}

void yac_interp_spmap_config_unpack(
  void * buffer, int buffer_size, int * position,
  struct yac_interp_spmap_config ** spmap_config, MPI_Comm comm) {

  *spmap_config = xmalloc(1 * sizeof(**spmap_config));
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position,
      &((*spmap_config)->spread_distance), 1, MPI_DOUBLE, comm), comm);
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position,
      &((*spmap_config)->max_search_distance), 1, MPI_DOUBLE, comm), comm);
  int weight_type;
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &weight_type, 1, MPI_INT, comm),
    comm);
  (*spmap_config)->weight_type =
    (enum yac_interp_spmap_weight_type)weight_type;
  yac_spmap_scale_config_unpack(
    buffer, buffer_size, position, &((*spmap_config)->scale_config), comm);
}

static struct yac_spmap_overwrite_config * yac_spmap_overwrite_config_unpack(
  void * buffer, int buffer_size, int * position, MPI_Comm comm) {

  struct yac_spmap_overwrite_config * overwrite_config =
    xmalloc(1 * sizeof(*overwrite_config));

  overwrite_config->src_point_selection =
    yac_point_selection_unpack(buffer, buffer_size, position, comm);
  yac_interp_spmap_config_unpack(
    buffer, buffer_size, position, &(overwrite_config->config), comm);

  return overwrite_config;
}

void yac_spmap_overwrite_configs_unpack(
  void * buffer, int buffer_size, int * position,
  struct yac_spmap_overwrite_config *** overwrite_configs, MPI_Comm comm) {

  size_t overwrite_config_count;
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &overwrite_config_count, 1,
      YAC_MPI_SIZE_T, comm), comm);

  *overwrite_configs =
    (overwrite_config_count > 0)?
      xcalloc(overwrite_config_count + 1, sizeof(**overwrite_configs)):NULL;

  for (size_t i = 0; i < overwrite_config_count; ++i)
    (*overwrite_configs)[i] =
      yac_spmap_overwrite_config_unpack(buffer, buffer_size, position, comm);
}
