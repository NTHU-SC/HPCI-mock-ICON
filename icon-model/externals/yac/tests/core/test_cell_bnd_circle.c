// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tests.h"
#include "geometry.h"
#include "grid_cell.h"
#include "test_common.h"

/** \file test_cell_bnd_circle.c
 *  \test
 * These are some examples on how to use \ref yac_get_cell_bounding_circle.
 */

double const tol = 1.0e-10;

static void utest_check_latlon_cell(double * coordinates_x, double * coordinates_y);
static void utest_check_gc_triangle(double * coordinates_x, double * coordinates_y);
static void utest_check_gc_quad(double * coordinates_x, double * coordinates_y);
static void utest_check_gc(double * coordinates_x, double * coordinates_y, size_t num_corners);
// tests whether all corners of cell are within the bounding circle
static void utest_test_circle(struct yac_grid_cell cell, struct bounding_circle circle);
// tests whether a point is within the bounding circle
static unsigned utest_point_in_circle(double point[3], struct bounding_circle circle);

int main (void) {

   { // test without corners
      struct yac_grid_cell cell = {
         .coordinates_xyz = NULL,
         .edge_type = NULL,
         .num_corners = 0,
         .array_size = 0
      };

      struct bounding_circle bnd_circle;
      yac_get_cell_bounding_circle(cell, &bnd_circle);
   }

   { // test regular cell

      double coordinates_x[] = {-1.0, 1.0, 1.0, -1.0};
      double coordinates_y[] = {-1.0, -1.0, 1.0, 1.0};
      utest_check_latlon_cell(coordinates_x, coordinates_y);
   }

   { // test regular cell

      double coordinates_x[] = {40.0, 45.0, 45.0, 40.0};
      double coordinates_y[] = {20.0, 20.0, 25.0, 25.0};
      utest_check_latlon_cell(coordinates_x, coordinates_y);
   }

   { // test regular cell

      double coordinates_x[] = { 175.0, -175.0, -175.0, 175.0};
      double coordinates_y[] = {-5.0, -5.0, 5.0, 5.0};
      utest_check_latlon_cell(coordinates_x, coordinates_y);
   }

   { // test regular cell


      double coordinates_x[] = {30.0, 40.0, 40.0, 30.0};
      double coordinates_y[] = {80.0, 80.0, 85.0, 85.0};
      utest_check_latlon_cell(coordinates_x, coordinates_y);
   }

   { // test regular cell

      double coordinates_x[] = {30.0, 40.0, 40.0, 30.0};
      double coordinates_y[] = {80.0, 80.0, 90.0, 90.0};
      utest_check_latlon_cell(coordinates_x, coordinates_y);
   }

   { // test regular cell


      double coordinates_x[] = {0.0, 1.0, 1.0, 0.0};
      double coordinates_y[] = {89.0, 89.0, 90.0, 90.0};
      utest_check_latlon_cell(coordinates_x, coordinates_y);
   }

   { // test regular cell


      double coordinates_x[] = {0.0, 1.0, 0.0, 0.0};
      double coordinates_y[] = {89.0, 89.0, 90.0, 90.0};
      utest_check_latlon_cell(coordinates_x, coordinates_y);
   }

   { // test triangle

      double coordinates_x[] = {-5.0, 5.0, 0.0};
      double coordinates_y[] = {-5.0, 5.0, -5.0};
      utest_check_gc_triangle(coordinates_x, coordinates_y);
   }

   { // test triangle

      double coordinates_x[] = {0.0, 120.0, -120.0};
      double coordinates_y[] = {85.0, 85.0, 85.0};
      utest_check_gc_triangle(coordinates_x, coordinates_y);
   }

   { // test triangle

      double coordinates_x[] = {0.0, 120.0, -120.0};
      double coordinates_y[] = {-85.0, -85.0, -85.0};
      utest_check_gc_triangle(coordinates_x, coordinates_y);
   }

   { // test triangle

      double coordinates_x[] = {-5.0, 5.0, 1.0};
      double coordinates_y[] = {0.0, 0.0, 1.0};
      utest_check_gc_triangle(coordinates_x, coordinates_y);
   }

   { // test triangle

      double coordinates_x[] = {0.0, 170.0, 260.0};
      double coordinates_y[] = {85.0, 85.0, 89.0};
      utest_check_gc_triangle(coordinates_x, coordinates_y);
   }

   { // test great circle quad

      double coordinates_x[] = {-1.0, 1.0, 1.0, -1.0};
      double coordinates_y[] = {-1.0, -1.0, 1.0, 1.0};
      utest_check_gc_quad(coordinates_x, coordinates_y);
   }

   { // test great circle quad

      double coordinates_x[] = {0.0, 90.0, 180.0, 270.0};
      double coordinates_y[] = {85.0, 85.0, 85.0, 85.0};
      utest_check_gc_quad(coordinates_x, coordinates_y);
   }

   { // test great circle quad

      double coordinates_x[] = {0.0, 90.0, 180.0, 270.0};
      double coordinates_y[] = {-85.0, -85.0, -85.0, -85.0};
      utest_check_gc_quad(coordinates_x, coordinates_y);
   }

   { // test great circle quad


      double coordinates_x[] = {0.0, 10.0, 0.0, -10.0};
      double coordinates_y[] = {-10.0, 0.0, 10.0, 0.0};
      utest_check_gc_quad(coordinates_x, coordinates_y);
   }

   { // test great circle pentagon


      double coordinates_x[] = {-1.0, -0.75, 0.75, 1.0, 0.0};
      double coordinates_y[] = {0.25, -1.0, -1.0, 0.25, 1.0};
      utest_check_gc(coordinates_x, coordinates_y, 5);
   }

   { // test great circle pentagon


      double coordinates_x[] = {0.0, 72.0, 144.0, 216.0, 288.0};
      double coordinates_y[] = {88.0, 88.0, 88.0, 88.0, 88.0};
      utest_check_gc(coordinates_x, coordinates_y, 5);
   }

   { // test great circle hexagon


      double coordinates_x[] = {-1.0, -0.75, 0.75, 1.0, 0.75, -0.75};
      double coordinates_y[] = {0.0, -1.0, -1.0, 0.0, 1.0, 1.0};
      utest_check_gc(coordinates_x, coordinates_y, 6);
   }

   { // test great circle hexagon


      double coordinates_x[] = {0.0, 60.0, 120.0, 180.0, 240.0, 300.0};
      double coordinates_y[] = {88.0, 88.0, 88.0, 88.0, 88.0, 88.0};
      utest_check_gc(coordinates_x, coordinates_y, 6);
   }

   { // test great circle polygon with 7 corners


      double coordinates_x[] = {0.0, 52.0, 104.0, 156.0, 208.0, 260.0, 312.0};
      double coordinates_y[] = {88.0, 88.0, 88.0, 88.0, 88.0, 88.0};
      utest_check_gc(coordinates_x, coordinates_y, 7);
   }

   return TEST_EXIT_CODE;
}

static void utest_check_gc_triangle(double * coordinates_x, double * coordinates_y) {

   enum yac_edge_type edges[] = {
     YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE};

   for (int order = -1; order <= 1; order += 2) {

      for (int start = 0; start < 3; ++start) {

         double temp_coordinates_x[3];
         double temp_coordinates_y[3];
         double coords[3][3];

         for (int i = 0; i < 3; ++i) {
            temp_coordinates_x[i] = coordinates_x[(3+i*order+start)%3];
            temp_coordinates_y[i] = coordinates_y[(3+i*order+start)%3];
            LLtoXYZ(temp_coordinates_x[i]*YAC_RAD,
                    temp_coordinates_y[i]*YAC_RAD, coords[i]);
         }

         struct yac_grid_cell cell =
           generate_cell_deg(temp_coordinates_x, temp_coordinates_y, edges, 3);

         struct bounding_circle bnd_circle;

         yac_get_cell_bounding_circle(cell, &bnd_circle);

         utest_test_circle(cell, bnd_circle);

         yac_free_grid_cell(&cell);
      }
   }
}

static void utest_check_gc_quad(double * coordinates_x, double * coordinates_y) {

   enum yac_edge_type edges[] = {YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
                                 YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE};

   for (unsigned i = 0; i < 4; ++i) {
     coordinates_x[i] *= YAC_RAD;
     coordinates_y[i] *= YAC_RAD;
   }

   for (int order = -1; order <= 1; order += 2) {

      for (int start = 0; start < 4; ++start) {

         double temp_coordinates_x[4];
         double temp_coordinates_y[4];

         for (int i = 0; i < 4; ++i) {
            temp_coordinates_x[i] = coordinates_x[(4+i*order+start)%4];
            temp_coordinates_y[i] = coordinates_y[(4+i*order+start)%4];
         }

         struct yac_grid_cell cell =
           generate_cell_deg(temp_coordinates_x, temp_coordinates_y, edges, 4);

         struct bounding_circle bnd_circle;

         yac_get_cell_bounding_circle(cell, &bnd_circle);

         utest_test_circle(cell, bnd_circle);

         yac_free_grid_cell(&cell);
      }
   }
}

static void utest_check_gc(
  double * coordinates_x, double * coordinates_y, size_t num_corners) {

   enum yac_edge_type edges[num_corners];

   for (size_t i = 0; i < num_corners; ++i) {
     coordinates_x[i] *= YAC_RAD;
     coordinates_y[i] *= YAC_RAD;
     edges[i] = YAC_GREAT_CIRCLE_EDGE;
   }

   for (int order = -1; order <= 1; order += 2) {

      for (size_t start = 0; start < num_corners; ++start) {

         double temp_coordinates_x[num_corners];
         double temp_coordinates_y[num_corners];

         for (size_t i = 0; i < num_corners; ++i) {
            temp_coordinates_x[i] =
              coordinates_x[(num_corners+i*order+start)%num_corners];
            temp_coordinates_y[i] =
              coordinates_y[(num_corners+i*order+start)%num_corners];
         }

         struct yac_grid_cell cell =
           generate_cell_deg(
             temp_coordinates_x, temp_coordinates_y, edges, num_corners);

         struct bounding_circle bnd_circle;

         yac_get_cell_bounding_circle(cell, &bnd_circle);

         utest_test_circle(cell, bnd_circle);

         yac_free_grid_cell(&cell);
      }
   }
}

static void utest_check_latlon_cell(double * coordinates_x, double * coordinates_y) {

   for (int order = -1; order <= 1; order += 2) {

      for (int start = 0; start < 4; ++start) {

         double temp_coordinates_x[4];
         double temp_coordinates_y[4];
         enum yac_edge_type edges[4];
         double coords[4][3];

         for (int i = 0; i < 4; ++i) {
            temp_coordinates_x[i] = coordinates_x[(4+i*order+start)%4];
            temp_coordinates_y[i] = coordinates_y[(4+i*order+start)%4];
            LLtoXYZ(temp_coordinates_x[i]*YAC_RAD,
                    temp_coordinates_y[i]*YAC_RAD, coords[i]);
         }


         for (int i = 0; i < 4; ++i) {
            int temp =
             fabs(temp_coordinates_y[i] - temp_coordinates_y[(i+3)%4]) > 0.0;
            edges[i] = (temp)?(YAC_LAT_CIRCLE_EDGE):(YAC_LON_CIRCLE_EDGE);
         }

         struct yac_grid_cell cell =
           generate_cell_deg(temp_coordinates_x, temp_coordinates_y, edges, 4);

         struct bounding_circle bnd_circle;

         yac_get_cell_bounding_circle(cell, &bnd_circle);

         utest_test_circle(cell, bnd_circle);

         yac_free_grid_cell(&cell);
      }
   }
}

/** normalises a vector while keeping z (lat) constant */
static inline void normalise_vector_lat(double v[]) {

  double temp = v[0] * v[0] + v[1] * v[1];

  if (fabs(temp) < 1e-9) {

    v[0] = 0.0;
    v[1] = 0.0;
    v[2] = copysign(1.0, v[2]);

  } else {

    double norm;

    if (fabs(v[2]) < 1e-9) {
      norm = 1.0;
      v[2] = 0.0;
    } else {
      norm = 1.0 - v[2] * v[2];
    }

    norm = sqrt(norm/temp);

    v[0] *= norm;
    v[1] *= norm;
  }
}

static void utest_test_circle(struct yac_grid_cell cell, struct bounding_circle circle) {

  for (size_t i = 0; i < cell.num_corners; ++i) {

    if (!utest_point_in_circle(cell.coordinates_xyz[i], circle))
      PUT_ERR("point is not in bounding circle\n");
    if (!yac_point_in_bounding_circle_vec(cell.coordinates_xyz[i], &circle))
      PUT_ERR("point is not in bounding circle\n");

    double edge_middle_point[3];
    size_t j = (i + 1) % cell.num_corners;
    edge_middle_point[0] =
      cell.coordinates_xyz[i][0] + cell.coordinates_xyz[j][0];
    edge_middle_point[1] =
      cell.coordinates_xyz[i][1] + cell.coordinates_xyz[j][1];
    edge_middle_point[2] =
      cell.coordinates_xyz[i][2] + cell.coordinates_xyz[j][2];

    switch (cell.edge_type[i]) {
      case(YAC_GREAT_CIRCLE_EDGE):
      case(YAC_LON_CIRCLE_EDGE):
        edge_middle_point[2] =
          cell.coordinates_xyz[i][2] + cell.coordinates_xyz[j][2];
        normalise_vector(edge_middle_point);
        break;
      case(YAC_LAT_CIRCLE_EDGE):
        edge_middle_point[2] = cell.coordinates_xyz[i][2];
        normalise_vector_lat(edge_middle_point);
        break;
      default:
        PUT_ERR("invalid edge type");
        continue;
    }

    if (!utest_point_in_circle(edge_middle_point, circle))
      PUT_ERR("point is not in bounding circle\n");
    if (!yac_point_in_bounding_circle_vec(edge_middle_point, &circle))
      PUT_ERR("point is not in bounding circle\n");
  }
}

static unsigned utest_point_in_circle(double point[3], struct bounding_circle circle) {

   // double const tol = 1.0e-12;
   // double const pi = 3.14159265358979323846;

   struct sin_cos_angle inc_angle =
      sum_angles_no_check(
        *(struct sin_cos_angle*)&(circle.inc_angle), SIN_COS_TOL);

   // if (circle.inc_angle + tol >= pi) return 1 == 1;
   if (compare_angles(inc_angle, SIN_COS_M_PI) >= 0) return 1 == 1;

   return
      compare_angles(
         get_vector_angle_2(circle.base_vector, point), inc_angle) <= 0;
}

