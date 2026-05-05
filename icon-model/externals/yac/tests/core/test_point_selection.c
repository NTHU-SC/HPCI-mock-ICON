// Copyright (c) 2025 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <mpi.h>
#include <math.h>
#include <string.h>

#include "tests.h"
#include "test_common.h"
#include "point_selection.h"
#include "geometry.h"

/** \file test_point_selection.c
 *  \test
* Tests the point selection type.
*/

#define YAC_RAD (0.01745329251994329576923690768489) // M_PI / 180

static void utest_check_point_selection(
  struct yac_point_selection * point_selection,
  yac_const_coordinate_pointer point_coords, size_t * point_indices,
  size_t num_points, size_t * ref_selected_points,
  size_t ref_num_selected_points);
static void utest_check_compare(
  struct yac_point_selection * a, struct yac_point_selection * b);

int main(void) {

  MPI_Init(NULL, NULL);

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  YAC_ASSERT(comm_size == 1, "ERROR wrong number of processes (has to be 1)")

  enum {NLON = 36, NLAT = 19};
  yac_coordinate_pointer point_coords =
    malloc(NLON * NLAT * sizeof(*point_coords));
  for (int i = 0, k = 0; i < NLAT; ++i)
    for (int j = 0; j < NLON; ++j, ++k)
      LLtoXYZ_deg((double)(10 * j), (double)(-90 + 10 * i), point_coords[k]);

  { // trivial test with empty point selection

    struct yac_point_selection * point_selection = NULL;

    enum {NUM_POINTS = 100};
    size_t point_indices[NUM_POINTS];
    for (size_t i = 0; i < NUM_POINTS; ++i) point_indices[i] = i;
    enum {REF_NUM_SELECTED_POINTS = 0};
    size_t * ref_selected_points = NULL;

    utest_check_point_selection(
      point_selection, (yac_const_coordinate_pointer)point_coords,
      point_indices, NUM_POINTS, ref_selected_points, REF_NUM_SELECTED_POINTS);

    yac_point_selection_delete(point_selection);
  }

  { // testing bounding circle point selection
    // * covers north pole

    double center_lon = 0.0;
    double center_lat = M_PI_2;
    double inc_angle = 11.0 * YAC_RAD;

    struct yac_point_selection * point_selection =
      yac_point_selection_bnd_circle_new(center_lon, center_lat, inc_angle);

    enum {NUM_POINTS = 3 * NLAT};
    size_t point_indices[NUM_POINTS];
    for (size_t i = 0, k = 0; i < 3; ++i)
      for (size_t j = 0; j < NLAT; ++j, ++k)
        point_indices[k] = i * 10 + j * NLON;
    enum {REF_NUM_SELECTED_POINTS = 6};
    size_t ref_selected_points[REF_NUM_SELECTED_POINTS] =
      {(NLAT - 2) * NLON + 0 * 10, (NLAT - 1) * NLON + 0 * 10,
       (NLAT - 2) * NLON + 1 * 10, (NLAT - 1) * NLON + 1 * 10,
       (NLAT - 2) * NLON + 2 * 10, (NLAT - 1) * NLON + 2 * 10};

    utest_check_point_selection(
      point_selection, (yac_const_coordinate_pointer)point_coords,
      point_indices, NUM_POINTS, ref_selected_points, REF_NUM_SELECTED_POINTS);

    yac_point_selection_delete(point_selection);
  }

  { // check yac_point_selection_compare routine
    struct yac_point_selection * point_selections[] = {
        NULL,
        yac_point_selection_bnd_circle_new(0.0, 0.0, 0.0),
        yac_point_selection_bnd_circle_new(0.0, 0.1, 0.0),
        yac_point_selection_bnd_circle_new(0.0, 0.0, 0.1),
        yac_point_selection_bnd_circle_new(0.1, 0.1, 0.0),
        yac_point_selection_bnd_circle_new(0.1, 0.0, 0.1),
        yac_point_selection_bnd_circle_new(0.0, 0.1, 0.1),
        yac_point_selection_bnd_circle_new(0.1, 0.1, 0.1)};
    enum {
      NUM_POINT_SELECTIONS =
        sizeof(point_selections)/sizeof(point_selections[0])
    };

    for (size_t i = 0; i < NUM_POINT_SELECTIONS; ++i) {
      for (size_t j = i + 1; j < NUM_POINT_SELECTIONS; ++j)
        utest_check_compare(point_selections[i], point_selections[j]);
      yac_point_selection_delete(point_selections[i]);
    }
  }

  free(point_coords);

  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static inline int compare_size_t(const void * a, const void * b) {

  size_t const * a_ = a, * b_ = b;

  return (*a_ > *b_) - (*b_ > *a_);
}

static void utest_check_point_selection_(
  struct yac_point_selection * point_selection,
  yac_const_coordinate_pointer point_coords, size_t * point_indices,
  size_t num_points, size_t * ref_sorted_points,
  size_t * ref_sorted_selected_points, size_t ref_num_selected_points) {

  size_t * point_indices_copy =
    malloc(num_points * sizeof(*point_indices_copy));
  memcpy(
    point_indices_copy, point_indices, num_points * sizeof(*point_indices));
  yac_coordinate_pointer selected_point_coords =
    malloc(num_points * sizeof(*selected_point_coords));
  for (size_t i = 0; i < num_points; ++i)
    memcpy(
      selected_point_coords[i], point_coords[point_indices[i]],
      sizeof(*selected_point_coords));

  size_t num_selected_points;

  yac_point_selection_apply(
    point_selection, selected_point_coords, point_indices_copy, num_points,
    &num_selected_points);

  if (num_selected_points == ref_num_selected_points) {

    // check selected points that were moved to the end of the
    // original list of points
    qsort(
      point_indices_copy + num_points - num_selected_points,
      num_selected_points, sizeof(*point_indices_copy), compare_size_t);
    for (size_t i = 0, j = num_points - num_selected_points;
         i < num_selected_points; ++i, ++j)
      if (point_indices_copy[j] != ref_sorted_selected_points[i])
        PUT_ERR("error in selected point within original point list");

    // check original list of points
    qsort(
      point_indices_copy, num_points, sizeof(*point_indices_copy),
      compare_size_t);
    for (size_t i = 0; i < num_points; ++i)
      if (point_indices_copy[i] != ref_sorted_points[i])
        PUT_ERR("error in original point list");

  } else {
    PUT_ERR("wrong number of selected points");
  }

  free(selected_point_coords);
  free(point_indices_copy);
}

static void utest_check_point_selection(
  struct yac_point_selection * point_selection,
  yac_const_coordinate_pointer point_coords, size_t * point_indices,
  size_t num_points, size_t * ref_selected_points,
  size_t ref_num_selected_points) {

  size_t * sorted_point_indices =
    malloc(num_points * sizeof(*sorted_point_indices));
  memcpy(
    sorted_point_indices, point_indices, num_points * sizeof(*point_indices));
  qsort(
    sorted_point_indices, num_points, sizeof(*sorted_point_indices),
    compare_size_t);
  qsort(
    ref_selected_points, ref_num_selected_points, sizeof(*ref_selected_points),
    compare_size_t);

  { // check original point selection
    utest_check_point_selection_(
      point_selection, point_coords, point_indices, num_points,
      sorted_point_indices, ref_selected_points, ref_num_selected_points);
  }

  { // check point selection generated by yac_point_selection_copy
    struct yac_point_selection * point_selection_copy =
      yac_point_selection_copy(point_selection);
    utest_check_point_selection_(
      point_selection_copy, point_coords, point_indices, num_points,
      sorted_point_indices, ref_selected_points, ref_num_selected_points);
    yac_point_selection_delete(point_selection_copy);
  }

  { // check point selection that was packed and unpacked
    size_t pack_size =
      yac_point_selection_get_pack_size(point_selection, MPI_COMM_WORLD);
    void * pack_buffer = malloc(pack_size);
    int position = 0;
    yac_point_selection_pack(
      point_selection, pack_buffer, (int)pack_size, &position, MPI_COMM_WORLD);
    position = 0;
    struct yac_point_selection * point_selection_copy =
      yac_point_selection_unpack(
        pack_buffer, (int)pack_size, &position, MPI_COMM_WORLD);
    free(pack_buffer);
    utest_check_point_selection_(
      point_selection_copy, point_coords, point_indices, num_points,
      sorted_point_indices, ref_selected_points, ref_num_selected_points);
    yac_point_selection_delete(point_selection_copy);
  }

  free(sorted_point_indices);
}

static void utest_check_compare(
  struct yac_point_selection * a, struct yac_point_selection * b) {

  if (yac_point_selection_compare(a, a))
    PUT_ERR("error in yac_point_selection_compare (a != a)")
  if (yac_point_selection_compare(b, b))
    PUT_ERR("error in yac_point_selection_compare (b != b)")

  int cmp_a = yac_point_selection_compare(a, b);
  int cmp_b = yac_point_selection_compare(b, a);

  if (cmp_a != -cmp_b)
    PUT_ERR("error in yac_point_selection_compare (cmp_a != -cmp_b)")
  if (!cmp_a) PUT_ERR("error in yac_point_selection_compare (cmp_a == 0)")
  if (!cmp_b) PUT_ERR("error in yac_point_selection_compare (cmp_b == 0)")
}
