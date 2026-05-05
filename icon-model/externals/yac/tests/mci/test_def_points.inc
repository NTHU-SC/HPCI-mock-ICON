! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

!> \file test_def_points.inc
!! \test
!! This shows how to use yac_fdef_points.

#include "test_macros.inc"

#ifndef TEST_PRECISION
#error "TEST_PRECISION is not defined"
#endif

  use utest
  use yac

  implicit none

  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND(6, 37)   !< single precision
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15, 307) !< double precision
  INTEGER, PARAMETER :: wp = TEST_PRECISION              !< selected working precision

  REAL(wp), PARAMETER :: YAC_RAD = 0.017453292519943295769_wp ! M_PI / 180

  integer :: comp_id
  integer :: grid_id
  integer :: point_id

  integer, parameter :: nbr_vertices(2) = (/3,3/)
  integer, parameter :: nbr_cells(2) = (/2,2/)
  integer, parameter :: nbr_edges = 12
  integer, parameter :: cyclic(2) = (/0,0/)

  real(wp), parameter :: &
    x_vertices(nbr_vertices(1)) = real((/-1.0, 0.0, 1.0/), wp) * YAC_RAD
  real(wp), parameter :: &
    y_vertices(nbr_vertices(2)) = real((/-1.0, 0.0, 1.0/), wp) * YAC_RAD

  real(wp), parameter :: &
    x_cells(nbr_cells(1)) = real((/-0.5, 0.5/), wp) * YAC_RAD
  real(wp), parameter :: &
    y_cells(nbr_cells(2)) = real((/-0.5, 0.5/), wp) * YAC_RAD
  real(wp), parameter :: &
    x_cells_2d(nbr_cells(1), nbr_cells(2)) = &
      reshape(real((/-0.5, 0.5, -0.5, 0.5/), wp) * YAC_RAD, nbr_cells)
  real(wp), parameter :: &
    y_cells_2d(nbr_cells(1), nbr_cells(2)) = &
      reshape(real((/-0.5, -0.5, 0.5, 0.5/), wp) * YAC_RAD, nbr_cells)

  real(wp), parameter ::  &
    x_edges(nbr_edges) =  &
      real(               &
        (/-0.5, 0.5,      &
          -1.0, 0.0, 1.0, &
          -0.5, 0.5,      &
          -1.0, 0.0, 1.0, &
          -0.5, 0.5/), wp) * YAC_RAD
  real(wp), parameter ::    &
    y_edges(nbr_edges) =    &
      real(                 &
        (/-1.0, 1.0,        &
          -0.5, -0.5, -0.5, &
          0.0, 0.0,         &
          0.5, 0.5, 0.5,    &
          1.0, 1.0/), wp) * YAC_RAD

  call yac_finit ( )

  ! set up dummy component
  call yac_fdef_comp ( "dummy-comp", comp_id )

  call test ( comp_id /= -99 )

  ! regular 2d grid

  !  7-------8-------9
  !  |       |       |
  !  |   3   |   4   |
  !  |       |       |
  !  4-------5-------6
  !  |       |       |
  !  |   1   |   2   |
  !  |       |       |
  !  1-------2-------3

  grid_id = -99

  call yac_fdef_grid ( 'reg2d_grid', &
                       nbr_vertices, &
                       cyclic,       &
                       x_vertices,   &
                       y_vertices,   &
                       grid_id )

  call test(grid_id /= -99)
  call test(yac_fget_grid_size(yac_location_cell, grid_id) == product(nbr_cells))
  call test(yac_fget_grid_size(yac_location_corner, grid_id) == product(nbr_vertices))
  call test(yac_fget_grid_size(yac_location_edge, grid_id) == nbr_edges)

  call start_test("def_points")

  point_id = -99

  ! regular 2d points

  call yac_fdef_points ( grid_id,             &
                         nbr_vertices,        &
                         yac_location_corner, &
                         x_vertices,          &
                         y_vertices,          &
                         point_id )

  call test(point_id /= -99)
  call test(yac_fget_points_size(point_id) == product(nbr_vertices))

  ! uniform unstructed 2d points

  call yac_fdef_points ( grid_id,           &
                         nbr_cells,         &
                         yac_location_cell, &
                         x_cells_2d,        &
                         y_cells_2d,        &
                         point_id )

  call test(point_id /= -99)
  call test(yac_fget_points_size(point_id) == product(nbr_cells))

  ! unstructed points

  call yac_fdef_points ( grid_id,           &
                         nbr_edges,         &
                         yac_location_edge, &
                         x_edges,           &
                         y_edges,           &
                         point_id )

  call test(point_id /= -99)
  call test(yac_fget_points_size(point_id) == nbr_edges)

  ! regular 2d rotatet points (cell points slightly rotated to the north)

  call yac_fdef_points ( grid_id,            &
                         nbr_cells,          &
                         yac_location_cell,  &
                         x_cells,            &
                         y_cells,            &
                         180.0_wp * YAC_RAD, &
                         89.7_wp * YAC_RAD,  &
                         point_id )

  call test(point_id /= -99)
  call test(yac_fget_points_size(point_id) == product(nbr_cells))

  call yac_ffinalize ( )

  call stop_test

  call exit_tests
