// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
// Get the definition of the 'restrict' keyword.
#include "config.h"
#endif

#include "interp_method_internal.h"
#include "interp_method_check.h"
#include "interp_method_callback.h"
#include "yac_mpi_common.h"

struct yac_interp_weights * yac_interp_method_do_search(
  struct interp_method ** method, struct yac_interp_grid * interp_grid) {

  size_t num_src_fields = yac_interp_grid_get_num_src_fields(interp_grid);
  MPI_Comm interp_grid_comm = yac_interp_grid_get_MPI_Comm(interp_grid);
  enum yac_location src_field_locations[num_src_fields];
  for (size_t i = 0; i < num_src_fields; ++i)
    src_field_locations[i] =
      yac_interp_grid_get_src_field_location(interp_grid, i);
  struct yac_interp_weights * weights =
    yac_interp_weights_new(
      interp_grid_comm, yac_interp_grid_get_tgt_field_location(interp_grid),
      src_field_locations, num_src_fields);

  if (*method == NULL) return weights;

  size_t temp_count;
  size_t * tgt_points;
  yac_interp_grid_get_tgt_points(interp_grid, &tgt_points, &temp_count);

  size_t final_count = 0;
  int interpolation_complete = 0;
  while (*method != NULL) {
    final_count +=
      (*method)->vtable->do_search(
        *method, interp_grid, tgt_points + final_count,
        temp_count - final_count, weights, &interpolation_complete);
    yac_mpi_call(
      MPI_Allreduce(
        MPI_IN_PLACE, &interpolation_complete, 1, MPI_INT,
        MPI_MAX, interp_grid_comm), interp_grid_comm);
    ++method;
  }

  free(tgt_points);

  return weights;
}

void yac_interp_method_delete(struct interp_method ** method) {

  while (*method != NULL) {
    struct interp_method * curr_method = *method;
    curr_method->vtable->delete(curr_method);
    ++method;
  }
}

void yac_interp_method_cleanup() {

  yac_interp_method_callback_buf_free();
  yac_interp_method_check_buf_free();
}
