// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include "basic_grid_data.h"
#include "grid_reg2d_common.h"
#include "geometry.h"
#include "utils_common.h"

void yac_rotate_coordinates(
  yac_coordinate_pointer coordinates, size_t num_coordinates,
  double north_pole[3]) {

  // compute angle between original and new north pol
  struct sin_cos_angle rot_angle =
    get_vector_angle_2(north_pole, (double[]){0.0, 0.0, 1.0});

  YAC_ASSERT_F(
    compare_angles(rot_angle, SIN_COS_TOL) > 0,
    "ERROR(yac_rotate_coordinates): "
    "new north pole is the original north pole "
    "(new north pole (%.3lf; %.3lf; %.3lf); angle (sin: %e, cos: %e)",
    north_pole[0], north_pole[1], north_pole[2], rot_angle.sin, rot_angle.cos);

  YAC_ASSERT_F(
    compare_angles(rot_angle, SIN_COS_M_PI_2) <= 0,
    "ERROR(yac_rotate_coordinates): "
    "new north pole is on the southern hemisphere "
    "(new north pole (%.3lf; %.3lf; %.3lf); angle (sin: %e, cos: %e)",
    north_pole[0], north_pole[1], north_pole[2], rot_angle.sin, rot_angle.cos);

  // compute rotatation axis
  // (has an angle of 90 degree between original and new north pole)
  double rot_axis[3];
  crossproduct_kahan((double[]){0.0, 0.0, 1.0}, north_pole, rot_axis);
  normalise_vector(rot_axis);

  for (size_t i = 0; i < num_coordinates; ++i) {
    double temp[3];
    rotate_vector2(rot_axis, rot_angle, coordinates[i], temp);
    coordinates[i][0] = temp[0];
    coordinates[i][1] = temp[1];
    coordinates[i][2] = temp[2];
  }
}

static struct yac_basic_grid_data yac_generate_basic_grid_data_reg_2d_rot_(
  size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices,
  double north_pol_lon, double north_pol_lat,
  void (*LLtoXYZ_ptr)(double, double, double[])) {

  YAC_ASSERT(
    !cyclic[1],
    "ERROR(yac_generate_basic_grid_data_reg_2d_rot): "
    "cyclic[1] != 0 not yet supported")

  size_t num_cells_2d[2] =
    {nbr_vertices[0] - (cyclic[0]?0:1), nbr_vertices[1] - (cyclic[1]?0:1)};
  size_t num_vertices_2d[2] = {num_cells_2d[0] + (cyclic[0]?0:1), num_cells_2d[1] + 1};
  size_t num_vertices = num_vertices_2d[0] * num_vertices_2d[1];
  size_t num_edges =
    (num_cells_2d[0] + (cyclic[0]?0:1)) * num_cells_2d[1] +
    num_cells_2d[0] * (num_cells_2d[1] + 1);

  yac_coordinate_pointer vertex_coordinates =
    xmalloc(num_vertices * sizeof(*vertex_coordinates));
  for (size_t i = 0, k = 0; i < num_vertices_2d[1]; ++i)
    for (size_t j = 0; j < num_vertices_2d[0]; ++j, ++k)
      LLtoXYZ_ptr(lon_vertices[j], lat_vertices[i], vertex_coordinates[k]);

  // rotate all coordinates according to provided north pole
  double north_pole[3];
  LLtoXYZ_ptr(north_pol_lon, north_pol_lat, north_pole);
  yac_rotate_coordinates(vertex_coordinates, num_vertices, north_pole);

  enum yac_edge_type * edge_type = xmalloc(num_edges * sizeof(*edge_type));
  for (size_t i = 0; i < num_edges; ++i) edge_type[i] = YAC_GREAT_CIRCLE_EDGE;

  struct yac_basic_grid_data grid =
    yac_generate_basic_grid_data_reg2d_common(nbr_vertices, cyclic);
  grid.vertex_coordinates = vertex_coordinates;
  grid.edge_type = edge_type;
  return grid;
}

struct yac_basic_grid_data yac_generate_basic_grid_data_reg_2d_rot(
  size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices,
  double north_pol_lon, double north_pol_lat) {

  return
    yac_generate_basic_grid_data_reg_2d_rot_(
      nbr_vertices, cyclic, lon_vertices, lat_vertices,
      north_pol_lon, north_pol_lat, LLtoXYZ);
}

struct yac_basic_grid_data yac_generate_basic_grid_data_reg_2d_rot_deg(
  size_t nbr_vertices[2], int cyclic[2],
  double *lon_vertices, double *lat_vertices,
  double north_pol_lon, double north_pol_lat) {

  return
    yac_generate_basic_grid_data_reg_2d_rot_(
      nbr_vertices, cyclic, lon_vertices, lat_vertices,
      north_pol_lon, north_pol_lat, LLtoXYZ_deg);
}
