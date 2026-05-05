! Copyright (c) 2024 The YAC Authors
!
! SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
! Get the definition of the 'YAC_MPI_FINT_FC_KIND' macro.
#include "config.h"
#endif

module yac_core
  use, intrinsic :: iso_c_binding, only : c_int, c_long, &
                                        & c_long_long, c_short, c_char, &
                                        & c_size_t, c_double, c_null_char, &
                                        & c_ptr, c_null_ptr

  implicit none

  public

  !------------------------------------------------
  ! Constants for Fortran-C interoperability
  !------------------------------------------------

  integer, parameter :: YAC_MPI_FINT_KIND = YAC_MPI_FINT_FC_KIND
  real(kind=c_double), parameter :: &
    YAC_FRAC_MASK_NO_VALUE = 133713371337.0_c_double

  !---------------
  ! Location enums
  !---------------

  enum, bind(c)
   enumerator :: YAC_LOC_CELL   = 0
   enumerator :: YAC_LOC_CORNER = 1
   enumerator :: YAC_LOC_EDGE   = 2
  end enum

  !----------------------------------
  ! Handling of existing weight files
  !----------------------------------

  enum, bind(c)
   enumerator :: YAC_WEIGHT_FILE_ERROR     = 0
   enumerator :: YAC_WEIGHT_FILE_KEEP      = 1
   enumerator :: YAC_WEIGHT_FILE_OVERWRITE = 2
  end enum

  !---------------------------
  ! Interpolation method enums
  !---------------------------

  enum, bind(c)
     enumerator :: YAC_INTERP_AVG_ARITHMETIC = 0
     enumerator :: YAC_INTERP_AVG_DIST       = 1
     enumerator :: YAC_INTERP_AVG_BARY       = 2
  end enum

  enum, bind(c)
     enumerator :: YAC_INTERP_NCC_AVG  = 0
     enumerator :: YAC_INTERP_NCC_DIST = 1
  end enum

  enum, bind(c)
     enumerator :: YAC_INTERP_NNN_AVG    = 0
     enumerator :: YAC_INTERP_NNN_DIST   = 1
     enumerator :: YAC_INTERP_NNN_GAUSS  = 2
     enumerator :: YAC_INTERP_NNN_RBF    = 3
     enumerator :: YAC_INTERP_NNN_ZERO   = 4
  end enum

  enum, bind(c)
     enumerator :: YAC_INTERP_CONSERV_DESTAREA = 0
     enumerator :: YAC_INTERP_CONSERV_FRACAREA = 1
  end enum

  enum, bind(c)
     enumerator :: YAC_INTERP_SPMAP_AVG  = 0
     enumerator :: YAC_INTERP_SPMAP_DIST = 1
  end enum

  enum, bind(c)
     enumerator :: YAC_INTERP_SPMAP_NONE       = 0
     enumerator :: YAC_INTERP_SPMAP_SRCAREA    = 1
     enumerator :: YAC_INTERP_SPMAP_INVTGTAREA = 2
     enumerator :: YAC_INTERP_SPMAP_FRACAREA   = 3
  end enum

  enum, bind(c)
    enumerator :: YAC_INTERP_FILE_MISSING_ERROR = 0 !<  abort on missing file
    enumerator :: YAC_INTERP_FILE_MISSING_CONT  = 1 !< continue on missing file
  end enum

  enum, bind(c)
    enumerator :: YAC_INTERP_FILE_SUCCESS_STOP = 0 !< prevents following interpolation
                                                   !! method from computating further
                                                   !! weights
    enumerator :: YAC_INTERP_FILE_SUCCESS_CONT = 1 !< continue weight computation with
                                                   !! following interpolation methods
  end enum

  !------------------------------------
  ! Interpolation method default values
  !------------------------------------

  ! average
  integer(kind=c_int), parameter :: YAC_INTERP_AVG_WEIGHT_TYPE_DEFAULT_F = 0_c_int
  integer(kind=c_int), parameter :: YAC_INTERP_AVG_PARTIAL_COVERAGE_DEFAULT_F = 0_c_int

  ! nearest corner cells
  integer(kind=c_int), parameter :: YAC_INTERP_NCC_WEIGHT_TYPE_DEFAULT_F = 0_c_int
  integer(kind=c_int), parameter :: YAC_INTERP_NCC_PARTIAL_COVERAGE_DEFAULT_F = 0_c_int

  ! conserv
  integer(kind=c_int), parameter :: YAC_INTERP_CONSERV_ORDER_DEFAULT_F = 1_c_int
  integer(kind=c_int), parameter :: YAC_INTERP_CONSERV_ENFORCED_CONSERV_DEFAULT_F = 0_c_int
  integer(kind=c_int), parameter :: YAC_INTERP_CONSERV_PARTIAL_COVERAGE_DEFAULT_F = 0_c_int
  integer(kind=c_int), parameter :: YAC_INTERP_CONSERV_NORMALISATION_DEFAULT_F = 0_c_int

  ! creep
  integer(kind=c_int), parameter :: YAC_INTERP_CREEP_DISTANCE_DEFAULT_F = -1_c_int

  ! fixed
  real(kind=c_double), parameter :: YAC_INTERP_FIXED_VALUE_DEFAULT_F = HUGE(1.0_c_double)

  ! n-nearest-neighbour
  integer(kind=c_int), parameter    :: YAC_INTERP_NNN_WEIGHTED_DEFAULT_F = 0_c_int
  integer(kind=c_size_t), parameter :: YAC_INTERP_NNN_N_DEFAULT_F = 1_c_size_t
  real(kind=c_double), parameter    :: YAC_INTERP_NNN_MAX_SEARCH_DISTANCE_DEFAULT_F = 0.0_c_double
  real(kind=c_double), parameter    :: YAC_INTERP_NNN_GAUSS_SCALE_DEFAULT_F = 0.1_c_double

  ! rbf
  integer(kind=c_size_t), parameter :: YAC_INTERP_RBF_N_DEFAULT_F = 9_c_size_t
  real(kind=c_double), parameter    :: YAC_INTERP_RBF_MAX_SEARCH_DISTANCE_DEFAULT_F = 0.0_c_double
  real(kind=c_double), parameter    :: YAC_INTERP_RBF_SCALE_DEFAULT_F = 1.487973e+01_c_double
  integer(kind=c_int), parameter    :: YAC_INTERP_RBF_KERNEL_DEFAULT_F = 0_c_int

  ! source-point mapping
  real(kind=c_double), parameter :: YAC_INTERP_SPMAP_SPREAD_DISTANCE_DEFAULT_F = 0.0_c_double
  real(kind=c_double), parameter :: YAC_INTERP_SPMAP_MAX_SEARCH_DISTANCE_DEFAULT_F = 0.0_c_double
  integer(kind=c_int), parameter :: YAC_INTERP_SPMAP_WEIGHTED_DEFAULT_F = 0_c_int
  integer(kind=c_int), parameter :: YAC_INTERP_SPMAP_SCALE_TYPE_DEFAULT_F = 0_c_int
  real(kind=c_double), parameter :: YAC_INTERP_SPMAP_SPHERE_RADIUS_DEFAULT_F = 1.0_c_double
  character(kind=c_char,len=*), parameter :: YAC_INTERP_SPMAP_FILENAME_DEFAULT_F = "" // c_null_char
  character(kind=c_char,len=*), parameter :: YAC_INTERP_SPMAP_VARNAME_DEFAULT_F = "" // c_null_char
  integer(kind=c_int), parameter :: YAC_INTERP_SPMAP_MIN_GLOBAL_ID_DEFAULT_F = 0_c_int
  type(c_ptr), parameter :: YAC_INTERP_SPMAP_CONFIG_DEFAULT_F = c_null_ptr
  type(c_ptr), parameter :: YAC_INTERP_SPMAP_OVERWRITE_DEFAULT_F = c_null_ptr

  ! check
  character(kind=c_char,len=*), parameter :: YAC_INTERP_CHECK_CONSTRUCTOR_KEY_DEFAULT_F = "" // c_null_char
  character(kind=c_char,len=*), parameter :: YAC_INTERP_CHECK_DO_SEARCH_KEY_DEFAULT_F = "" // c_null_char

  ! user file
  integer(kind=c_int), parameter :: YAC_INTERP_FILE_ON_MISSING_FILE_DEFAULT_F = 0_c_int
  integer(kind=c_int), parameter :: YAC_INTERP_FILE_ON_SUCCESS_DEFAULT_F = 1_c_int

  !---------------------------
  ! Interpolation weight enums
  !---------------------------

  enum, bind(c)
    enumerator :: YAC_MAPPING_ON_SRC = 0
    enumerator :: YAC_MAPPING_ON_TGT = 1
  end enum

  !------------------------------------------------
  ! Initialisation and finalization of MPI and yaxt
  !------------------------------------------------

  interface

    ! initialise MPI (if not already initialised)
    subroutine yac_mpi_init_c () &
      bind ( c, name='yac_mpi_init' )
    end subroutine yac_mpi_init_c

    ! initialise yaxt
    subroutine yac_yaxt_init_c (comm) &
      bind ( c, name='yac_yaxt_init_f2c' )

      import :: YAC_MPI_FINT_KIND

      integer(kind=YAC_MPI_FINT_KIND), value :: comm
    end subroutine yac_yaxt_init_c

    ! check whether MPI is initialised
    function yac_mpi_is_initialised_c () &
      bind ( c, name='yac_mpi_is_initialised' )

      use, intrinsic :: iso_c_binding

      integer(kind=c_int) :: yac_mpi_is_initialised_c
    end function yac_mpi_is_initialised_c

    ! free internal buffers and finalize yaxt
    subroutine yac_mpi_cleanup_c () &
      bind ( c, name="yac_mpi_cleanup" )
    end subroutine yac_mpi_cleanup_c

    ! finalize MPI (if initialised by YAC)
    subroutine yac_mpi_finalize_c () &
      bind ( c, name='yac_mpi_finalize' )
    end subroutine yac_mpi_finalize_c

  end interface

  !-----------------------
  ! Handling of basic grid
  !-----------------------

  interface

    function yac_basic_grid_empty_new_c(name) &
      bind ( C, name='yac_basic_grid_empty_new' )

      use, intrinsic :: iso_c_binding

      character(kind=c_char) :: name(*)

      type(c_ptr) :: yac_basic_grid_empty_new_c
    end function yac_basic_grid_empty_new_c

    subroutine yac_basic_grid_delete_c(grid) &
      bind ( C, name='yac_basic_grid_delete' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: grid
    end subroutine yac_basic_grid_delete_c

    function yac_basic_grid_get_data_size_c(grid, location) &
      bind ( C, name='yac_basic_grid_get_data_size_f2c')

      use, intrinsic :: iso_c_binding

      type(c_ptr), value         :: grid
      integer(kind=c_int), value :: location

      integer(kind=c_size_t) :: yac_basic_grid_get_data_size_c
    end function yac_basic_grid_get_data_size_c

    function yac_basic_grid_add_coordinates_c( &
      grid, location, coordinates, count) &
      bind ( C, name='yac_basic_grid_add_coordinates_f2c')

      use, intrinsic :: iso_c_binding

      type(c_ptr), value            :: grid
      integer(kind=c_int), value    :: location
      real(kind=c_double)           :: coordinates(*)
      integer(kind=c_size_t), value :: count

      integer(kind=c_size_t) :: yac_basic_grid_add_coordinates_c
    end function yac_basic_grid_add_coordinates_c

    function yac_basic_grid_add_coordinates_nocpy_c( &
      grid, location, coordinates) &
      bind ( C, name='yac_basic_grid_add_coordinates_nocpy_f2c')

      use, intrinsic :: iso_c_binding

      type(c_ptr), value            :: grid
      integer(kind=c_int), value    :: location
      real(kind=c_double)           :: coordinates(*)

      integer(kind=c_size_t) :: yac_basic_grid_add_coordinates_nocpy_c
    end function yac_basic_grid_add_coordinates_nocpy_c

    function yac_basic_grid_add_mask_c( &
      grid, location, mask, count, name) &
      bind ( C, name='yac_basic_grid_add_mask_f2c')

      use, intrinsic :: iso_c_binding

      type(c_ptr), value            :: grid
      integer(kind=c_int), value    :: location
      integer(kind=c_int)           :: mask(*)
      integer(kind=c_size_t), value :: count
      character(kind=c_char)        :: name(*)

      integer(kind=c_size_t) :: yac_basic_grid_add_mask_c
    end function yac_basic_grid_add_mask_c

    function yac_basic_grid_add_mask_nocpy_c( &
      grid, location, mask, name) &
      bind ( C, name='yac_basic_grid_add_mask_nocpy_f2c')

      use, intrinsic :: iso_c_binding

      type(c_ptr), value            :: grid
      integer(kind=c_int), value    :: location
      integer(kind=c_int)           :: mask(*)
      character(kind=c_char)        :: name(*)

      integer(kind=c_size_t) :: yac_basic_grid_add_mask_nocpy_c
    end function yac_basic_grid_add_mask_nocpy_c

    function yac_basic_grid_reg_2d_new_c( &
      name, nbr_vertices, cyclic, lon_vertices, lat_vertices) &
      bind ( C, name='yac_basic_grid_reg_2d_new')

      use, intrinsic :: iso_c_binding

      character(kind=c_char) :: name(*)
      integer(kind=c_size_t) :: nbr_vertices(2)
      integer(kind=c_int)    :: cyclic(2)
      real(kind=c_double)    :: lon_vertices(*)
      real(kind=c_double)    :: lat_vertices(*)

      ! struct yac_basic_grid *
      type(c_ptr) :: yac_basic_grid_reg_2d_new_c
    end function yac_basic_grid_reg_2d_new_c

    function yac_basic_grid_reg_2d_deg_new_c( &
      name, nbr_vertices, cyclic, lon_vertices, lat_vertices) &
      bind ( C, name='yac_basic_grid_reg_2d_deg_new')

      use, intrinsic :: iso_c_binding

      character(kind=c_char) :: name(*)
      integer(kind=c_size_t) :: nbr_vertices(2)
      integer(kind=c_int)    :: cyclic(2)
      real(kind=c_double)    :: lon_vertices(*)
      real(kind=c_double)    :: lat_vertices(*)

      ! struct yac_basic_grid *
      type(c_ptr) :: yac_basic_grid_reg_2d_deg_new_c
    end function yac_basic_grid_reg_2d_deg_new_c

    function yac_basic_grid_curve_2d_new_c( &
      name, nbr_vertices, cyclic, lon_vertices, lat_vertices) &
      bind ( C, name='yac_basic_grid_curve_2d_new')

      use, intrinsic :: iso_c_binding

      character(kind=c_char) :: name(*)
      integer(kind=c_size_t) :: nbr_vertices(2)
      integer(kind=c_int)    :: cyclic(2)
      real(kind=c_double)    :: lon_vertices(*)
      real(kind=c_double)    :: lat_vertices(*)

      ! struct yac_basic_grid *
      type(c_ptr) :: yac_basic_grid_curve_2d_new_c
    end function yac_basic_grid_curve_2d_new_c

    function yac_basic_grid_curve_2d_deg_new_c( &
      name, nbr_vertices, cyclic, lon_vertices, lat_vertices) &
      bind ( C, name='yac_basic_grid_curve_2d_deg_new')

      use, intrinsic :: iso_c_binding

      character(kind=c_char) :: name(*)
      integer(kind=c_size_t) :: nbr_vertices(2)
      integer(kind=c_int)    :: cyclic(2)
      real(kind=c_double)    :: lon_vertices(*)
      real(kind=c_double)    :: lat_vertices(*)

      ! struct yac_basic_grid *
      type(c_ptr) :: yac_basic_grid_curve_2d_deg_new_c
    end function yac_basic_grid_curve_2d_deg_new_c

    function yac_basic_grid_unstruct_new_c( &
      name, nbr_vertices, nbr_cells, num_vertices_per_cell, &
      x_vertices, y_vertices, cell_to_vertex) &
      bind ( C, name='yac_basic_grid_unstruct_new')

      use, intrinsic :: iso_c_binding

      character(kind=c_char)        :: name(*)
      integer(kind=c_size_t), value :: nbr_vertices
      integer(kind=c_size_t), value :: nbr_cells
      integer(kind=c_int)           :: num_vertices_per_cell(*)
      real(kind=c_double)           :: x_vertices(*)
      real(kind=c_double)           :: y_vertices(*)
      integer(kind=c_int)           :: cell_to_vertex(*)

      ! struct yac_basic_grid *
      type(c_ptr) :: yac_basic_grid_unstruct_new_c
    end function yac_basic_grid_unstruct_new_c

    function yac_basic_grid_unstruct_deg_new_c( &
      name, nbr_vertices, nbr_cells, num_vertices_per_cell, &
      x_vertices, y_vertices, cell_to_vertex) &
      bind ( C, name='yac_basic_grid_unstruct_deg_new')

      use, intrinsic :: iso_c_binding

      character(kind=c_char)        :: name(*)
      integer(kind=c_size_t), value :: nbr_vertices
      integer(kind=c_size_t), value :: nbr_cells
      integer(kind=c_int)           :: num_vertices_per_cell(*)
      real(kind=c_double)           :: x_vertices(*)
      real(kind=c_double)           :: y_vertices(*)
      integer(kind=c_int)           :: cell_to_vertex(*)

      ! struct yac_basic_grid *
      type(c_ptr) :: yac_basic_grid_unstruct_deg_new_c
    end function yac_basic_grid_unstruct_deg_new_c

    function yac_basic_grid_unstruct_ll_new_c( &
      name, nbr_vertices, nbr_cells, num_vertices_per_cell, &
      x_vertices, y_vertices, cell_to_vertex) &
      bind ( C, name='yac_basic_grid_unstruct_ll_new')

      use, intrinsic :: iso_c_binding

      character(kind=c_char)        :: name(*)
      integer(kind=c_size_t), value :: nbr_vertices
      integer(kind=c_size_t), value :: nbr_cells
      integer(kind=c_int)           :: num_vertices_per_cell(*)
      real(kind=c_double)           :: x_vertices(*)
      real(kind=c_double)           :: y_vertices(*)
      integer(kind=c_int)           :: cell_to_vertex(*)

      ! struct yac_basic_grid *
      type(c_ptr) :: yac_basic_grid_unstruct_ll_new_c
    end function yac_basic_grid_unstruct_ll_new_c

    function yac_basic_grid_unstruct_ll_deg_new_c( &
      name, nbr_vertices, nbr_cells, num_vertices_per_cell, &
      x_vertices, y_vertices, cell_to_vertex) &
      bind ( C, name='yac_basic_grid_unstruct_ll_deg_new')

      use, intrinsic :: iso_c_binding

      character(kind=c_char)        :: name(*)
      integer(kind=c_size_t), value :: nbr_vertices
      integer(kind=c_size_t), value :: nbr_cells
      integer(kind=c_int)           :: num_vertices_per_cell(*)
      real(kind=c_double)           :: x_vertices(*)
      real(kind=c_double)           :: y_vertices(*)
      integer(kind=c_int)           :: cell_to_vertex(*)

      ! struct yac_basic_grid *
      type(c_ptr) :: yac_basic_grid_unstruct_ll_deg_new_c
    end function yac_basic_grid_unstruct_ll_deg_new_c

    function yac_basic_grid_unstruct_edge_new_c( &
      name, nbr_vertices, nbr_cells, nbr_edges, num_edges_per_cell, &
      x_vertices, y_vertices, cell_to_edge, edge_to_vertex) &
      bind ( C, name='yac_basic_grid_unstruct_edge_new')

      use, intrinsic :: iso_c_binding

      character(kind=c_char)        :: name(*)
      integer(kind=c_size_t), value :: nbr_vertices
      integer(kind=c_size_t), value :: nbr_cells
      integer(kind=c_size_t), value :: nbr_edges
      integer(kind=c_int)           :: num_edges_per_cell(*)
      real(kind=c_double)           :: x_vertices(*)
      real(kind=c_double)           :: y_vertices(*)
      integer(kind=c_int)           :: cell_to_edge(*)
      integer(kind=c_int)           :: edge_to_vertex(*)

      ! struct yac_basic_grid *
      type(c_ptr) :: yac_basic_grid_unstruct_edge_new_c
    end function yac_basic_grid_unstruct_edge_new_c

    function yac_basic_grid_unstruct_edge_deg_new_c( &
      name, nbr_vertices, nbr_cells, nbr_edges, num_edges_per_cell, &
      x_vertices, y_vertices, cell_to_edge, edge_to_vertex) &
      bind ( C, name='yac_basic_grid_unstruct_edge_deg_new')

      use, intrinsic :: iso_c_binding

      character(kind=c_char)        :: name(*)
      integer(kind=c_size_t), value :: nbr_vertices
      integer(kind=c_size_t), value :: nbr_cells
      integer(kind=c_size_t), value :: nbr_edges
      integer(kind=c_int)           :: num_edges_per_cell(*)
      real(kind=c_double)           :: x_vertices(*)
      real(kind=c_double)           :: y_vertices(*)
      integer(kind=c_int)           :: cell_to_edge(*)
      integer(kind=c_int)           :: edge_to_vertex(*)

      ! struct yac_basic_grid *
      type(c_ptr) :: yac_basic_grid_unstruct_edge_deg_new_c
    end function yac_basic_grid_unstruct_edge_deg_new_c

    function yac_basic_grid_unstruct_edge_ll_new_c( &
      name, nbr_vertices, nbr_cells, nbr_edges, num_edges_per_cell, &
      x_vertices, y_vertices, cell_to_edge, edge_to_vertex) &
      bind ( C, name='yac_basic_grid_unstruct_edge_ll_new')

      use, intrinsic :: iso_c_binding

      character(kind=c_char)        :: name(*)
      integer(kind=c_size_t), value :: nbr_vertices
      integer(kind=c_size_t), value :: nbr_cells
      integer(kind=c_size_t), value :: nbr_edges
      integer(kind=c_int)           :: num_edges_per_cell(*)
      real(kind=c_double)           :: x_vertices(*)
      real(kind=c_double)           :: y_vertices(*)
      integer(kind=c_int)           :: cell_to_edge(*)
      integer(kind=c_int)           :: edge_to_vertex(*)

      ! struct yac_basic_grid *
      type(c_ptr) :: yac_basic_grid_unstruct_edge_ll_new_c
    end function yac_basic_grid_unstruct_edge_ll_new_c

    function yac_basic_grid_unstruct_edge_ll_deg_new_c( &
      name, nbr_vertices, nbr_cells, nbr_edges, num_edges_per_cell, &
      x_vertices, y_vertices, cell_to_edge, edge_to_vertex) &
      bind ( C, name='yac_basic_grid_unstruct_edge_ll_deg_new')

      use, intrinsic :: iso_c_binding

      character(kind=c_char)        :: name(*)
      integer(kind=c_size_t), value :: nbr_vertices
      integer(kind=c_size_t), value :: nbr_cells
      integer(kind=c_size_t), value :: nbr_edges
      integer(kind=c_int)           :: num_edges_per_cell(*)
      real(kind=c_double)           :: x_vertices(*)
      real(kind=c_double)           :: y_vertices(*)
      integer(kind=c_int)           :: cell_to_edge(*)
      integer(kind=c_int)           :: edge_to_vertex(*)

      ! struct yac_basic_grid *
      type(c_ptr) :: yac_basic_grid_unstruct_edge_ll_deg_new_c
    end function yac_basic_grid_unstruct_edge_ll_deg_new_c

    function yac_basic_grid_cloud_new_c( &
      name, nbr_points, x_points, y_points) &
      bind ( C, name='yac_basic_grid_cloud_new')

      use, intrinsic :: iso_c_binding

      character(kind=c_char)        :: name(*)
      integer(kind=c_size_t), value :: nbr_points
      real(kind=c_double)           :: x_points(*)
      real(kind=c_double)           :: y_points(*)

      ! struct yac_basic_grid *
      type(c_ptr) :: yac_basic_grid_cloud_new_c
    end function yac_basic_grid_cloud_new_c

    function yac_basic_grid_cloud_deg_new_c( &
      name, nbr_points, x_points, y_points) &
      bind ( C, name='yac_basic_grid_cloud_deg_new')

      use, intrinsic :: iso_c_binding

      character(kind=c_char)        :: name(*)
      integer(kind=c_size_t), value :: nbr_points
      real(kind=c_double)           :: x_points(*)
      real(kind=c_double)           :: y_points(*)

      ! struct yac_basic_grid *
      type(c_ptr) :: yac_basic_grid_cloud_deg_new_c
    end function yac_basic_grid_cloud_deg_new_c

    function yac_basic_grid_reg_2d_rot_new_c( &
      name, nbr_vertices, cyclic, lon_vertices, lat_vertices, &
      north_pol_lon, north_pol_lat) &
      bind ( C, name='yac_basic_grid_reg_2d_rot_new')

      use, intrinsic :: iso_c_binding

      character(kind=c_char)     :: name(*)
      integer(kind=c_size_t)     :: nbr_vertices(2)
      integer(kind=c_int)        :: cyclic(2)
      real(kind=c_double)        :: lon_vertices(*)
      real(kind=c_double)        :: lat_vertices(*)
      real(kind=c_double), value :: north_pol_lon
      real(kind=c_double), value :: north_pol_lat

      type(c_ptr) :: yac_basic_grid_reg_2d_rot_new_c
    end function yac_basic_grid_reg_2d_rot_new_c

    function yac_basic_grid_reg_2d_rot_deg_new_c( &
      name, nbr_vertices, cyclic, lon_vertices, lat_vertices, &
      north_pol_lon, north_pol_lat) &
      bind ( C, name='yac_basic_grid_reg_2d_rot_deg_new')

      use, intrinsic :: iso_c_binding

      character(kind=c_char)     :: name(*)
      integer(kind=c_size_t)     :: nbr_vertices(2)
      integer(kind=c_int)        :: cyclic(2)
      real(kind=c_double)        :: lon_vertices(*)
      real(kind=c_double)        :: lat_vertices(*)
      real(kind=c_double), value :: north_pol_lon
      real(kind=c_double), value :: north_pol_lat

      type(c_ptr) :: yac_basic_grid_reg_2d_rot_deg_new_c
    end function yac_basic_grid_reg_2d_rot_deg_new_c

    subroutine yac_basic_grid_compute_cell_areas_c(grid, cell_areas) &
      bind ( C, name='yac_basic_grid_compute_cell_areas')

      use, intrinsic :: iso_c_binding

      type(c_ptr), value  :: grid
      real(kind=c_double) :: cell_areas(*)
    end subroutine yac_basic_grid_compute_cell_areas_c

  end interface

  !------------------------------
  ! Handling of distributed grids
  !------------------------------

  interface

    function yac_dist_grid_pair_new_c(grid_a, grid_b, comm) &
      bind ( C, name='yac_dist_grid_pair_new_f2c' )

      use, intrinsic :: iso_c_binding
      import :: YAC_MPI_FINT_KIND

      type(c_ptr), value                     :: grid_a
      type(c_ptr), value                     :: grid_b
      integer(kind=YAC_MPI_FINT_KIND), value :: comm

      ! struct yac_dist_grid_pair *
      type(c_ptr) :: yac_dist_grid_pair_new_c
    end function yac_dist_grid_pair_new_c

    subroutine yac_dist_grid_pair_delete_c(grid_pair) &
      bind ( C, name='yac_dist_grid_pair_delete' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: grid_pair
    end subroutine yac_dist_grid_pair_delete_c

  end interface

  !-------------------------------
  ! Handling of interpolation grid
  !-------------------------------

  interface

    function yac_interp_grid_new_c( dist_grid_pair,            &
                                    src_grid_name,             &
                                    tgt_grid_name,             &
                                    num_src_fields,            &
                                    src_field_locations,       &
                                    src_field_coordinate_idxs, &
                                    src_field_masks_idxs,      &
                                    tgt_field_location,        &
                                    tgt_field_coordinate_idx,  &
                                    tgt_field_masks_idx)        &
      bind ( C, name='yac_interp_grid_new_f2c' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value            :: dist_grid_pair
      character(kind=c_char)        :: src_grid_name(*)
      character(kind=c_char)        :: tgt_grid_name(*)
      integer(kind=c_size_t), value :: num_src_fields
      integer(kind=c_int)           :: src_field_locations(*)
      integer(kind=c_size_t)        :: src_field_coordinate_idxs(*)
      integer(kind=c_size_t)        :: src_field_masks_idxs(*)
      integer(kind=c_int), value    :: tgt_field_location
      integer(kind=c_size_t), value :: tgt_field_coordinate_idx
      integer(kind=c_size_t), value :: tgt_field_masks_idx

      ! struct yac_interp_grid *
      type(c_ptr) :: yac_interp_grid_new_c
    end function yac_interp_grid_new_c

    subroutine yac_interp_grid_delete_c(interp_grid) &
      bind ( C, name='yac_interp_grid_delete' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: interp_grid
    end subroutine yac_interp_grid_delete_c

  end interface

  !---------------------------------------------------------------------
  ! Handling of datastructures required for extended spmap configuration
  !---------------------------------------------------------------------

  interface

    function yac_point_selection_bnd_circle_new_c( &
      center_lon, center_lat, inc_angle)           &
      bind ( C, name='yac_point_selection_bnd_circle_new' )

      use, intrinsic :: iso_c_binding

      real(kind=c_double), value :: center_lon
      real(kind=c_double), value :: center_lat
      real(kind=c_double), value :: inc_angle

      ! struct yac_point_selection *
      type(c_ptr) :: yac_point_selection_bnd_circle_new_c
    end function yac_point_selection_bnd_circle_new_c

    subroutine yac_point_selection_delete_c(point_select) &
      bind ( C, name='yac_point_selection_delete' )

      use, intrinsic :: iso_c_binding

      ! struct yac_point_selection *
      type(c_ptr), value :: point_select
    end subroutine yac_point_selection_delete_c

    function yac_interp_spmap_config_new_c( &
      spread_distance, max_search_distance, &
      weight_type, scale_config)            &
      bind ( C, name='yac_interp_spmap_config_new_f2c' )

      use, intrinsic :: iso_c_binding

      real(kind=c_double), value :: spread_distance
      real(kind=c_double), value :: max_search_distance
      integer(kind=c_int), value :: weight_type
      ! struct yac_spmap_scale_config *
      type(c_ptr), value         :: scale_config

      ! struct yac_interp_spmap_config *
      type(c_ptr) :: yac_interp_spmap_config_new_c
    end function yac_interp_spmap_config_new_c

    subroutine yac_interp_spmap_config_delete_c(spmap_config) &
      bind ( C, name='yac_interp_spmap_config_delete' )

      use, intrinsic :: iso_c_binding

      ! struct yac_interp_spmap_config *
      type(c_ptr), value :: spmap_config
    end subroutine yac_interp_spmap_config_delete_c

    function yac_spmap_scale_config_new_c(                          &
      scale_type, source_cell_area_config, target_cell_area_config) &
      bind ( C, name='yac_spmap_scale_config_new_f2c' )

      use, intrinsic :: iso_c_binding

      integer(kind=c_int), value :: scale_type
      ! struct yac_spmap_cell_area_config *
      type(c_ptr), value         :: source_cell_area_config
      ! struct yac_spmap_cell_area_config *
      type(c_ptr), value         :: target_cell_area_config

      ! struct yac_spmap_scale_config *
      type(c_ptr) :: yac_spmap_scale_config_new_c
    end function yac_spmap_scale_config_new_c

    subroutine yac_spmap_scale_config_delete_c(scale_config) &
      bind ( C, name='yac_spmap_scale_config_delete' )

      use, intrinsic :: iso_c_binding

      ! struct yac_spmap_scale_config *
      type(c_ptr), value :: scale_config
    end subroutine yac_spmap_scale_config_delete_c

    function yac_spmap_cell_area_config_yac_new_c(sphere_radius) &
      bind ( C, name='yac_spmap_cell_area_config_yac_new' )

      use, intrinsic :: iso_c_binding

      real(kind=c_double), value :: sphere_radius

      ! struct yac_spmap_cell_area_config *
      type(c_ptr) :: yac_spmap_cell_area_config_yac_new_c
    end function yac_spmap_cell_area_config_yac_new_c

    function yac_spmap_cell_area_config_file_new_c( &
      filename, varname, min_global_id) &
      bind ( C, name='yac_spmap_cell_area_config_file_new_f2c' )

      use, intrinsic :: iso_c_binding

      character(kind=c_char)     :: filename(*)
      character(kind=c_char)     :: varname(*)
      integer(kind=c_int), value :: min_global_id

      ! struct yac_spmap_cell_area_config *
      type(c_ptr) :: yac_spmap_cell_area_config_file_new_c
    end function yac_spmap_cell_area_config_file_new_c

    subroutine yac_spmap_cell_area_config_delete_c(cell_area_config) &
      bind ( C, name='yac_spmap_cell_area_config_delete' )

      use, intrinsic :: iso_c_binding

      ! struct yac_spmap_cell_area_config *
      type(c_ptr), value :: cell_area_config
    end subroutine yac_spmap_cell_area_config_delete_c

    function yac_spmap_overwrite_config_new_c(src_point_selection, config) &
      bind ( C, name='yac_spmap_overwrite_config_new' )

      use, intrinsic :: iso_c_binding

      ! struct yac_point_selection *
      type(c_ptr), value :: src_point_selection
      ! struct yac_interp_spmap_config *
      type(c_ptr), value :: config

      ! struct yac_spmap_overwrite_config *
      type(c_ptr) :: yac_spmap_overwrite_config_new_c
    end function yac_spmap_overwrite_config_new_c

    subroutine yac_spmap_overwrite_config_delete_c(overwrite_config) &
      bind ( C, name='yac_spmap_overwrite_config_delete' )

      use, intrinsic :: iso_c_binding

      ! struct yac_spmap_overwrite_config *
      type(c_ptr), value :: overwrite_config
    end subroutine yac_spmap_overwrite_config_delete_c

    subroutine yac_spmap_overwrite_configs_delete_c(overwrite_configs) &
      bind ( C, name='yac_spmap_overwrite_configs_delete' )

      use, intrinsic :: iso_c_binding

      ! struct yac_spmap_overwrite_config **
      type(c_ptr), value :: overwrite_configs
    end subroutine yac_spmap_overwrite_configs_delete_c

  end interface

  !--------------------------------
  ! Handling of interpolation stack
  !--------------------------------

  interface

    function yac_interp_stack_config_new_c() &
      bind ( C, name='yac_interp_stack_config_new' )

      use, intrinsic :: iso_c_binding

      type(c_ptr) :: yac_interp_stack_config_new_c
    end function yac_interp_stack_config_new_c

    function yac_interp_stack_config_copy_c(interp_stack_config) &
      bind ( C, name='yac_interp_stack_config_copy' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: interp_stack_config

      type(c_ptr) :: yac_interp_stack_config_copy_c
    end function yac_interp_stack_config_copy_c

    subroutine yac_interp_stack_config_delete_c(interp_stack_config) &
      bind ( C, name='yac_interp_stack_config_delete' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: interp_stack_config
    end subroutine yac_interp_stack_config_delete_c

    subroutine yac_interp_stack_config_add_average_c( &
      interp_stack_config, reduction_type, partial_coverage) &
      bind ( C, name='yac_interp_stack_config_add_average_f2c' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value         :: interp_stack_config
      integer(kind=c_int), value :: reduction_type
      integer(kind=c_int), value :: partial_coverage
    end subroutine yac_interp_stack_config_add_average_c

    subroutine yac_interp_stack_config_add_ncc_c( &
      interp_stack_config, weight_type, partial_coverage) &
      bind ( C, name='yac_interp_stack_config_add_ncc_f2c' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value         :: interp_stack_config
      integer(kind=c_int), value :: weight_type
      integer(kind=c_int), value :: partial_coverage
    end subroutine yac_interp_stack_config_add_ncc_c

    subroutine yac_interp_stack_config_add_nnn_c( &
      interp_stack_config, type, n, max_search_distance, scale) &
      bind ( C, name='yac_interp_stack_config_add_nnn_f2c' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value            :: interp_stack_config
      integer(kind=c_int), value    :: type
      integer(kind=c_size_t), value :: n
      real(kind=c_double), value    :: max_search_distance
      real(kind=c_double), value    :: scale
    end subroutine yac_interp_stack_config_add_nnn_c

    subroutine yac_interp_stack_config_add_rbf_c( &
      interp_stack_config, n, max_search_distance, scale) &
      bind ( C, name='yac_interp_stack_config_add_rbf' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value            :: interp_stack_config
      integer(kind=c_size_t), value :: n
      real(kind=c_double), value    :: max_search_distance
      real(kind=c_double), value    :: scale
    end subroutine yac_interp_stack_config_add_rbf_c

    subroutine yac_interp_stack_config_add_conservative_c(            &
      interp_stack_config, order, enforced_conserv, partial_coverage, &
      normalisation) &
      bind ( C, name='yac_interp_stack_config_add_conservative_f2c' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value         :: interp_stack_config
      integer(kind=c_int), value :: order
      integer(kind=c_int), value :: enforced_conserv
      integer(kind=c_int), value :: partial_coverage
      integer(kind=c_int), value :: normalisation
    end subroutine yac_interp_stack_config_add_conservative_c

    subroutine yac_interp_stack_config_add_spmap_ext_c( &
      interp_stack_config, default_config, overwrite_configs) &
      bind ( C, name='yac_interp_stack_config_add_spmap_ext' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: interp_stack_config
      type(c_ptr), value :: default_config
      type(c_ptr), value :: overwrite_configs
    end subroutine yac_interp_stack_config_add_spmap_ext_c

    subroutine yac_interp_stack_config_add_spmap_c( &
      interp_stack_config, spread_distance, max_search_distance, &
      weight_type, scale_type, &
      src_sphere_radius, src_filename, src_varname, src_min_global_id, &
      tgt_sphere_radius, tgt_filename, tgt_varname, tgt_min_global_id) &
      bind ( C, name='yac_interp_stack_config_add_spmap_f2c' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value         :: interp_stack_config
      real(kind=c_double), value :: spread_distance
      real(kind=c_double), value :: max_search_distance
      integer(kind=c_int), value :: weight_type
      integer(kind=c_int), value :: scale_type
      real(kind=c_double), value :: src_sphere_radius
      character(kind=c_char)     :: src_filename(*)
      character(kind=c_char)     :: src_varname(*)
      integer(kind=c_int), value :: src_min_global_id
      real(kind=c_double), value :: tgt_sphere_radius
      character(kind=c_char)     :: tgt_filename(*)
      character(kind=c_char)     :: tgt_varname(*)
      integer(kind=c_int), value :: tgt_min_global_id
    end subroutine yac_interp_stack_config_add_spmap_c

    subroutine yac_interp_stack_config_add_hcsbb_c(interp_stack_config) &
      bind ( C, name='yac_interp_stack_config_add_hcsbb' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: interp_stack_config
    end subroutine yac_interp_stack_config_add_hcsbb_c

    subroutine yac_interp_stack_config_add_user_file_c( &
      interp_stack_config, filename, on_missing_file, on_success) &
      bind ( C, name='yac_interp_stack_config_add_user_file_f2c' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value         :: interp_stack_config
      character(kind=c_char)     :: filename(*)
      integer(kind=c_int), value :: on_missing_file
      integer(kind=c_int), value :: on_success
    end subroutine yac_interp_stack_config_add_user_file_c

    subroutine yac_interp_stack_config_add_fixed_c( &
      interp_stack_config, value) &
      bind ( C, name='yac_interp_stack_config_add_fixed' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value         :: interp_stack_config
      real(kind=c_double), value :: value
    end subroutine yac_interp_stack_config_add_fixed_c

    subroutine yac_interp_stack_config_add_creep_c( &
      interp_stack_config, n) &
      bind ( C, name='yac_interp_stack_config_add_creep' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value         :: interp_stack_config
      integer(kind=c_int), value :: n
    end subroutine yac_interp_stack_config_add_creep_c

    subroutine yac_interp_stack_config_add_check_c( &
      interp_stack_config, constructor_key, do_search_key) &
      bind ( C, name='yac_interp_stack_config_add_check' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value     :: interp_stack_config
      character(kind=c_char) :: constructor_key(*)
      character(kind=c_char) :: do_search_key(*)
    end subroutine yac_interp_stack_config_add_check_c

    function yac_interp_stack_config_compare_c(a, b) &
      bind ( C, name='yac_interp_stack_config_compare' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: a
      type(c_ptr), value :: b

      integer(kind=c_int) :: yac_interp_stack_config_compare_c
    end function yac_interp_stack_config_compare_c

    function yac_interp_stack_config_generate_c(interp_stack_config) &
      bind ( C, name='yac_interp_stack_config_generate' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: interp_stack_config

      ! struct interp_method **
      type(c_ptr) :: yac_interp_stack_config_generate_c
    end function yac_interp_stack_config_generate_c

  end interface

  !----------------------------------
  ! Handling of interpolation methods
  !----------------------------------

  interface

    function yac_interp_method_do_search_c(interp_method_stack, interp_grid) &
      bind ( C, name='yac_interp_method_do_search' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: interp_method_stack
      type(c_ptr), value :: interp_grid

      ! struct yac_interp_weights *
      type(c_ptr) :: yac_interp_method_do_search_c
    end function yac_interp_method_do_search_c

    subroutine yac_interp_method_delete_c(interp_method_stack) &
      bind ( C, name='yac_interp_method_delete' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: interp_method_stack
    end subroutine yac_interp_method_delete_c

    subroutine yac_interp_method_cleanup_c() &
      bind ( C, name='yac_interp_method_cleanup' )

    end subroutine yac_interp_method_cleanup_c

  end interface

  !----------------------------------
  ! Handling of interpolation weights
  !----------------------------------

  interface

    subroutine yac_interp_weights_delete_c(interp_weights) &
      bind ( C, name='yac_interp_weights_delete' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: interp_weights
    end subroutine yac_interp_weights_delete_c

    subroutine yac_interp_weights_write_to_file_c(interp_weights, &
                                                  filename,       &
                                                  src_grid_name,  &
                                                  tgt_grid_name,  &
                                                  src_grid_size,  &
                                                  tgt_grid_size,  &
                                                  on_existing)    &
      bind ( C, name='yac_interp_weights_write_to_file' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value            :: interp_weights
      character(kind=c_char)        :: filename(*)
      character(kind=c_char)        :: src_grid_name(*)
      character(kind=c_char)        :: tgt_grid_name(*)
      integer(kind=c_size_t), value :: src_grid_size
      integer(kind=c_size_t), value :: tgt_grid_size
      integer(kind=c_int), value    :: on_existing
    end subroutine yac_interp_weights_write_to_file_c

    function yac_interp_weights_get_interpolation_c( &
      interp_weights, reorder, collection_size, frac_mask_fallback_value, &
      scaling_factor, scaling_summand, yaxt_exchanger_name, &
      is_source, is_target) &
      bind ( C, name='yac_interp_weights_get_interpolation_f2c' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value            :: interp_weights
      integer(kind=c_int), value    :: reorder
      integer(kind=c_size_t), value :: collection_size
      real(kind=c_double), value    :: frac_mask_fallback_value
      real(kind=c_double), value    :: scaling_factor
      real(kind=c_double), value    :: scaling_summand
      character(kind=c_char)        :: yaxt_exchanger_name
      integer(kind=c_int), value    :: is_source
      integer(kind=c_int), value    :: is_target

      ! struct yac_interpolation *
      type(c_ptr) :: yac_interp_weights_get_interpolation_c
    end function yac_interp_weights_get_interpolation_c

  end interface

  !---------------------------
  ! Handling of interpolations
  !---------------------------

  interface

    subroutine yac_interpolation_delete_c(interp) &
      bind ( C, name='yac_interpolation_delete' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: interp
    end subroutine yac_interpolation_delete_c

    subroutine yac_interpolation_execute_c( &
      interp, src_fields, tgt_field) &
      bind ( C, name='yac_interpolation_execute' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: interp
      type(c_ptr), value :: src_fields
      type(c_ptr), value :: tgt_field
    end subroutine yac_interpolation_execute_c

    subroutine yac_interpolation_execute_frac_c( &
      interp, src_fields, src_frac_masks, tgt_field) &
      bind ( C, name='yac_interpolation_execute_frac' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: interp
      type(c_ptr), value :: src_fields
      type(c_ptr), value :: src_frac_masks
      type(c_ptr), value :: tgt_field
    end subroutine yac_interpolation_execute_frac_c

    subroutine yac_interpolation_execute_put_c(interp, src_fields) &
      bind ( C, name='yac_interpolation_execute_put' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: interp
      type(c_ptr), value :: src_fields
    end subroutine yac_interpolation_execute_put_c

    subroutine yac_interpolation_execute_put_frac_c( &
      interp, src_fields, src_frac_masks) &
      bind ( C, name='yac_interpolation_execute_put_frac' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: interp
      type(c_ptr), value :: src_fields
      type(c_ptr), value :: src_frac_masks
    end subroutine yac_interpolation_execute_put_frac_c

    subroutine yac_interpolation_execute_get_c(interp, tgt_field) &
      bind ( C, name='yac_interpolation_execute_get' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: interp
      type(c_ptr), value :: tgt_field
    end subroutine yac_interpolation_execute_get_c

    subroutine yac_interpolation_execute_get_async_c(interp, tgt_field) &
      bind ( C, name='yac_interpolation_execute_get_async' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: interp
      type(c_ptr), value :: tgt_field
    end subroutine yac_interpolation_execute_get_async_c

    function yac_interpolation_execute_put_test_c(interp) &
      bind ( C, name='yac_interpolation_execute_put_test' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: interp

      integer(kind=c_int) :: yac_interpolation_execute_put_test_c
    end function yac_interpolation_execute_put_test_c

    function yac_interpolation_execute_get_test_c(interp) &
      bind ( C, name='yac_interpolation_execute_get_test' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: interp

      integer(kind=c_int) :: yac_interpolation_execute_get_test_c
    end function yac_interpolation_execute_get_test_c

    subroutine yac_interpolation_execute_wait_c(interp) &
      bind ( C, name='yac_interpolation_execute_wait' )

      use, intrinsic :: iso_c_binding

      type(c_ptr), value :: interp
    end subroutine yac_interpolation_execute_wait_c

    function yac_interpolation_get_const_frac_mask_no_value_c() &
      bind ( C, name='yac_interpolation_get_const_frac_mask_no_value_c2f' )

      use, intrinsic :: iso_c_binding

      real(kind=c_double) :: &
        yac_interpolation_get_const_frac_mask_no_value_c
    end function yac_interpolation_get_const_frac_mask_no_value_c

    function yac_interpolation_get_const_frac_mask_undef_c() &
      bind ( C, name='yac_interpolation_get_const_frac_mask_undef_c2f' )

      use, intrinsic :: iso_c_binding

      real(kind=c_double) :: &
        yac_interpolation_get_const_frac_mask_undef_c
    end function yac_interpolation_get_const_frac_mask_undef_c

  end interface

  !--------------------
  ! Abort functionality
  !--------------------

  interface

     subroutine yac_abort_message_c ( text, file, line ) &
           bind ( c, name='yac_abort_message' )

       use, intrinsic :: iso_c_binding, only : c_int, c_char

       character ( kind=c_char ), dimension(*) :: text       !< [IN] error message
       character ( kind=c_char ), dimension(*) :: file       !< [IN] error message
       integer ( kind=c_int ), value           :: line       !< [IN] line number

     end subroutine yac_abort_message_c

  end interface

  !---------------------------------
  ! Setting of callback routines for
  ! check interpolation method
  !---------------------------------

  interface

     subroutine yac_interp_method_check_add_constructor_callback_c ( &
       constructor_callback, user_data, key ) &
       bind ( c, name='yac_interp_method_check_add_constructor_callback' )

       use, intrinsic :: iso_c_binding, only : c_funptr, c_ptr, c_char

       type(c_funptr), value                   :: constructor_callback
       type(c_ptr), value                      :: user_data
       character ( kind=c_char ), dimension(*) :: key

     end subroutine yac_interp_method_check_add_constructor_callback_c

     subroutine yac_interp_method_check_add_do_search_callback_c ( &
       do_search_callback, user_data, key ) &
       bind ( c, name='yac_interp_method_check_add_do_search_callback' )

       use, intrinsic :: iso_c_binding, only : c_funptr, c_ptr, c_char

       type(c_funptr), value                   :: do_search_callback
       type(c_ptr), value                      :: user_data
       character ( kind=c_char ), dimension(*) :: key

     end subroutine yac_interp_method_check_add_do_search_callback_c

  end interface

end module yac_core
