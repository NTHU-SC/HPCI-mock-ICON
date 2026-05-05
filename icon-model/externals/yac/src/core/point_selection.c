// Copyright (c) 2025 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include "string.h"

#include "point_selection.h"
#include "utils_core.h"
#include "yac_mpi_internal.h"
#include "geometry.h"

struct yac_point_selection {
  enum yac_point_selection_type type;
  union {
    struct {
      double center_lon;
      double center_lat;
      double inc_angle;
    } bnd_circle;
  } data;
};

struct yac_point_selection * yac_point_selection_bnd_circle_new(
  double center_lon, double center_lat, double inc_angle) {

  YAC_ASSERT_F(
    (center_lat >= - M_PI_2) && (center_lat <= M_PI_2),
    "ERROR(yac_point_selection_bnd_circle_new): "
    "invalid latitude for bounding circle center (%lf) "
    "(has to be in the range of [-M_PI_2;M_PI_2])", center_lat);
  YAC_ASSERT_F(
    (inc_angle >= 0.0) && (inc_angle < M_PI),
    "ERROR(yac_point_selection_bnd_circle_new): "
    "invalid angle for bounding circle (%lf) "
    "(has to be in the range of [0;M_PI[)", inc_angle);

  struct yac_point_selection * point_select =
    xmalloc(1 * sizeof(*point_select));

  point_select->type = YAC_POINT_SELECTION_TYPE_BND_CIRCLE;
  point_select->data.bnd_circle.center_lon = center_lon;
  point_select->data.bnd_circle.center_lat = center_lat;
  point_select->data.bnd_circle.inc_angle = inc_angle;

  return point_select;
}

struct yac_point_selection * yac_point_selection_copy(
  struct yac_point_selection const * point_select) {

  enum yac_point_selection_type type =
    yac_point_selection_get_type(point_select);

  struct yac_point_selection * point_select_copy;
  switch (type) {
    YAC_UNREACHABLE_DEFAULT(
      "ERROR(yac_point_selection_copy): invalid point selection type");
    case(YAC_POINT_SELECTION_TYPE_EMPTY):
      point_select_copy = NULL;
      break;
    case(YAC_POINT_SELECTION_TYPE_BND_CIRCLE):
      point_select_copy =
        yac_point_selection_bnd_circle_new(
          point_select->data.bnd_circle.center_lon,
          point_select->data.bnd_circle.center_lat,
          point_select->data.bnd_circle.inc_angle);
      break;
  }

  return point_select_copy;
}

void yac_point_selection_delete(struct yac_point_selection * point_select) {
  free(point_select);
}

size_t yac_point_selection_get_pack_size(
  struct yac_point_selection const * point_select, MPI_Comm comm) {

  int int_pack_size;
  yac_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &int_pack_size), comm);

  enum yac_point_selection_type type =
    yac_point_selection_get_type(point_select);

  size_t pack_size = int_pack_size; // type

  switch (type) {
    YAC_UNREACHABLE_DEFAULT(
      "ERROR(yac_point_selection_get_pack_size): invalid point selection type");
    case(YAC_POINT_SELECTION_TYPE_EMPTY):
      break;
    case(YAC_POINT_SELECTION_TYPE_BND_CIRCLE): {
      int dble_pack_size;
      yac_mpi_call(
        MPI_Pack_size(1, MPI_DOUBLE, comm, &dble_pack_size), comm);
      pack_size += dble_pack_size + // center_lon
                   dble_pack_size + // center_lat
                   dble_pack_size;  // inc_angle
    };
  }

  return pack_size;
}

void yac_point_selection_pack(
  struct yac_point_selection const * point_select,
  void * buffer, int buffer_size, int * position, MPI_Comm comm) {

  int int_pack_size;
  yac_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &int_pack_size), comm);

  int type = (int)yac_point_selection_get_type(point_select);

  yac_mpi_call(
    MPI_Pack(&type, 1, MPI_INT, buffer, buffer_size, position, comm), comm);

  switch (type) {
    YAC_UNREACHABLE_DEFAULT(
      "ERROR(yac_point_selection_pack): invalid point selection type");
    case(YAC_POINT_SELECTION_TYPE_EMPTY):
      break;
    case(YAC_POINT_SELECTION_TYPE_BND_CIRCLE): {
      yac_mpi_call(
        MPI_Pack(
          &(point_select->data.bnd_circle.center_lon), 1, MPI_DOUBLE,
          buffer, buffer_size, position, comm), comm);
      yac_mpi_call(
        MPI_Pack(
          &(point_select->data.bnd_circle.center_lat), 1, MPI_DOUBLE,
          buffer, buffer_size, position, comm), comm);
      yac_mpi_call(
        MPI_Pack(
          &(point_select->data.bnd_circle.inc_angle), 1, MPI_DOUBLE,
          buffer, buffer_size, position, comm), comm);
    }
  }
}

struct yac_point_selection * yac_point_selection_unpack(
  void const * buffer, int buffer_size, int * position, MPI_Comm comm) {

  int type;
  yac_mpi_call(
    MPI_Unpack(
      buffer, buffer_size, position, &type, 1, MPI_INT, comm), comm);

  struct yac_point_selection * point_select;

  switch (type) {
    YAC_UNREACHABLE_DEFAULT(
      "ERROR(yac_point_selection_unpack): invalid point selection type");
    case(YAC_POINT_SELECTION_TYPE_EMPTY):
      point_select = NULL;
      break;
    case(YAC_POINT_SELECTION_TYPE_BND_CIRCLE): {
      double center_lon, center_lat, inc_angle;
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position, &center_lon, 1,
          MPI_DOUBLE, comm), comm);
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position, &center_lat, 1,
          MPI_DOUBLE, comm), comm);
      yac_mpi_call(
        MPI_Unpack(
          buffer, buffer_size, position, &inc_angle, 1,
          MPI_DOUBLE, comm), comm);
      point_select =
        yac_point_selection_bnd_circle_new(center_lon, center_lat, inc_angle);
      break;
    }
  }

  return point_select;
}

int yac_point_selection_compare(
  struct yac_point_selection const * a, struct yac_point_selection const * b) {

  enum yac_point_selection_type type_a = yac_point_selection_get_type(a);
  enum yac_point_selection_type type_b = yac_point_selection_get_type(b);

  int ret = (type_a > type_b) - (type_a < type_b);

  if (!ret) {

    switch (type_a) {
      YAC_UNREACHABLE_DEFAULT(
        "ERROR(yac_point_selection_compare): invalid point selection type");
      case(YAC_POINT_SELECTION_TYPE_EMPTY):
        ret = 0;
        break;
      case(YAC_POINT_SELECTION_TYPE_BND_CIRCLE): {
        ret =
          (a->data.bnd_circle.center_lon > b->data.bnd_circle.center_lon) -
          (a->data.bnd_circle.center_lon < b->data.bnd_circle.center_lon);
        if (ret) break;
        ret =
          (a->data.bnd_circle.center_lat > b->data.bnd_circle.center_lat) -
          (a->data.bnd_circle.center_lat < b->data.bnd_circle.center_lat);
        if (ret) break;
        ret =
          (a->data.bnd_circle.inc_angle > b->data.bnd_circle.inc_angle) -
          (a->data.bnd_circle.inc_angle < b->data.bnd_circle.inc_angle);
        break;
      }
    }
  }
  return ret;
}

// Sorts the provided arrays based on a flag-array (containing the
// values "0" and "!= 0"). After the sort, all array elements whose associated
// flag value is "0" are the front of the array.
//
// This sort is:
//   * not stable
//   * has a time complexity of O(n)
static void flag_sort(
  size_t * array_size_t, yac_coordinate_pointer array_coord,
  int * flag, size_t false_count/*, size_t true_count*/) {

  // The number of "true" elements in the 0...false_count-1 range of the
  // array is identical to the number of "false" elements in the
  // false_count...false_count+true_count-1. We just have to find matching
  // pairs and swap them.
  for (size_t i = 0, j = false_count; i < false_count; ++i) {
    // if there is a wrongfully placed "true" element
    if (flag[i]) {
      // find a wrongfully place "false" element
      for (; flag[j]; ++j);
      // swap elements
      double temp_coord[3];
      memcpy(temp_coord, array_coord[i], sizeof(temp_coord));
      memcpy(array_coord[i], array_coord[j], sizeof(temp_coord));
      memcpy(array_coord[j], temp_coord, sizeof(temp_coord));
      size_t temp_size_t = array_size_t[i];
      array_size_t[i] = array_size_t[j];
      array_size_t[j] = temp_size_t;
      // set to next element in "true" list
      ++j;
    }
  }
}

void yac_point_selection_apply(
  struct yac_point_selection const * point_select,
  yac_coordinate_pointer point_coords, size_t * point_indices,
  size_t num_points, size_t * num_selected_points) {

  enum yac_point_selection_type type =
    yac_point_selection_get_type(point_select);

  switch (type) {
    YAC_UNREACHABLE_DEFAULT(
      "ERROR(yac_point_selection_apply): invalid point selection type");
    case(YAC_POINT_SELECTION_TYPE_EMPTY):
      *num_selected_points = 0;
      break;
    case(YAC_POINT_SELECTION_TYPE_BND_CIRCLE): {

      // generate bounding bounding circle
      struct bounding_circle bnd_circle = {.sq_crd = DBL_MAX};
      LLtoXYZ(
        point_select->data.bnd_circle.center_lon,
        point_select->data.bnd_circle.center_lat,
        bnd_circle.base_vector);
      double sin_inc_angle, cos_inc_angle;
      compute_sin_cos(
        point_select->data.bnd_circle.inc_angle,
        &sin_inc_angle, &cos_inc_angle);
      bnd_circle.inc_angle = sin_cos_angle_new(sin_inc_angle, cos_inc_angle);

      int * match_flag = xmalloc(num_points * sizeof(*match_flag));

      // determine matching points and count them
      size_t match_count = 0;
      for (size_t i = 0; i < num_points; ++i) {

        int point_in_bounding_circle =
          yac_point_in_bounding_circle_vec(
            (double*)(point_coords[i]), &bnd_circle);
        match_flag[i] = point_in_bounding_circle;
        if (point_in_bounding_circle) ++match_count;
      }

      // sort selected points to the end of the list
      flag_sort(
        point_indices, point_coords, match_flag, num_points - match_count);
      free(match_flag);

      *num_selected_points = match_count;
      break;
    }
  }
}

enum yac_point_selection_type yac_point_selection_get_type(
  struct yac_point_selection const * point_select) {

  return
    (point_select != NULL)?point_select->type:YAC_POINT_SELECTION_TYPE_EMPTY;
}

void yac_point_selection_bnd_circle_get_config(
  struct yac_point_selection const * point_selection,
  double * center_lon, double * center_lat, double * inc_angle) {

  YAC_ASSERT(
    yac_point_selection_get_type(point_selection) ==
    YAC_POINT_SELECTION_TYPE_BND_CIRCLE,
    "ERROR(yac_point_selection_bnd_circle_get_config): "
    "invalid point selection type");

  *center_lon = point_selection->data.bnd_circle.center_lon;
  *center_lat = point_selection->data.bnd_circle.center_lat;
  *inc_angle = point_selection->data.bnd_circle.inc_angle;
}
