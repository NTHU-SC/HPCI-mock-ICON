! --------------------------------------------------------------------
!> Example plugin for the ICON Community Interface (ComIn)
!  with basic (not MPI-parallel) callbacks and accessing variables and
!  descriptive data structures.
!
!  Note that in order to demonstrate ComIn's language interoperability,
!  a similary plugin has been implemented in C, see the subdirectory
!  "simple_c".
!
!
!  @authors 01/2023 :: ICON Community INTERFACE  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  Please see the file LICENSE in the root of the source tree for this code.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.
! --------------------------------------------------------------------
MODULE simple_fortran_plugin
  USE mpi
  USE iso_c_binding,           ONLY : C_INT, c_ptr, c_f_pointer, C_ASSOCIATED
  USE comin_plugin_interface,  ONLY : comin_callback_register,                                          &
    &                                 comin_var_get, comin_parallel_get_host_mpi_comm,                  &
    &                                 t_comin_var_descriptor, t_comin_var_handle,                       &
    &                                 comin_var_request_add,                                            &
    &                                 comin_descrdata_get_domain, t_comin_descrdata_domain,             &
    &                                 comin_descrdata_get_global, t_comin_descrdata_global,             &
    &                                 comin_descrdata_get_simulation_interval,                          &
    &                                 t_comin_descrdata_simulation_interval,                            &
    &                                 t_comin_setup_version_info, comin_setup_get_version,              &
    &                                 EP_SECONDARY_CONSTRUCTOR, EP_DESTRUCTOR,                          &
    &                                 EP_ATM_PHYSICS_BEFORE, EP_ATM_WRITE_OUTPUT_BEFORE,                &
    &                                 COMIN_FLAG_READ, COMIN_FLAG_WRITE, COMIN_ZAXIS_2D,                &
    &                                 comin_parallel_get_host_mpi_rank, comin_current_get_domain_id,    &
    &                                 COMIN_DOMAIN_OUTSIDE_LOOP, comin_callback_get_ep_name,            &
    &                                 t_comin_plugin_info, comin_current_get_plugin_info,               &
    &                                 comin_plugin_finish, comin_metadata_set,                          &
    &                                 comin_metadata_get,                                               &
    &                                 comin_descrdata_get_timesteplength, COMIN_HGRID_UNSTRUCTURED_EDGE,&
    &                                 comin_error_check, comin_var_get_descr_list_head,                 &
    &                                 comin_var_get_descr_list_next, comin_var_get_descr_list_var_desc, &
    &                                 comin_error_check, comin_print_info,                              &
    &                                 COMIN_DIM_SEMANTICS_NPROMA, COMIN_DIM_SEMANTICS_LEVEL,            &
    &                                 COMIN_DIM_SEMANTICS_BLOCK, COMIN_DIM_SEMANTICS_UNUSED,            &
    &                                 COMIN_VAR_DATATYPE_DOUBLE

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: pluginname = "simple_fortran_plugin"

  !> working precision (will be compared to ComIn's and ICON's)
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)
  TYPE(t_comin_setup_version_info) :: version

  TYPE(t_comin_var_handle)           :: pres, vn, simple_fortran_var, simple_fortran_tracer
  INTEGER                         :: rank

  !> access descriptive data structures
  TYPE(t_comin_descrdata_domain),     POINTER   :: p_patch
  TYPE(t_comin_descrdata_global),     POINTER   :: p_global
  TYPE(t_comin_descrdata_simulation_interval), POINTER   :: p_simulation_interval

  TYPE(t_comin_var_handle),  ALLOCATABLE :: qv(:)

  CHARACTER(LEN=200)            :: text

CONTAINS

  ! --------------------------------------------------------------------
  ! ComIn primary constructor.
  ! --------------------------------------------------------------------
  SUBROUTINE comin_main()  BIND(C)
    !
    TYPE(t_comin_plugin_info)     :: this_plugin
    TYPE(t_comin_var_descriptor)  :: simple_fortran_d, simple_fortran_tracer_d
    REAL(wp)                      :: dtime_1, dtime_2

    rank = comin_parallel_get_host_mpi_rank()
    CALL comin_print_info("setup")

    version = comin_setup_get_version()
    IF (version%version_no_major > 1)  THEN
      CALL comin_plugin_finish("comin_main (simple_fortran_plugin)", "incompatible version!")
    END IF

    !> check plugin id
    CALL comin_current_get_plugin_info(this_plugin)
    WRITE (text,'(a,i4)') "     plugin  id: ", this_plugin%id
    CALL comin_print_info(text)

    !> add requests for additional ICON variables

    ! request host model to add variable simple_fortran_var
    simple_fortran_d = t_comin_var_descriptor(id = 1, name = "simple_fortran_var")
    CALL comin_var_request_add(simple_fortran_d, .FALSE.)
    CALL comin_metadata_set(simple_fortran_d, "tracer", .FALSE.)
    CALL comin_metadata_set(simple_fortran_d, "restart", .FALSE.)

    ! request host model to add tracer simple_fortran_tracer
    simple_fortran_tracer_d = t_comin_var_descriptor( id = -1, name = "simple_fortran_tracer" )
    CALL comin_var_request_add(simple_fortran_tracer_d, .FALSE.)
    CALL comin_metadata_set(simple_fortran_tracer_d, "tracer", .TRUE.)
    CALL comin_metadata_set(simple_fortran_tracer_d, "restart", .FALSE.)

    ! register callbacks
    CALL comin_callback_register(EP_SECONDARY_CONSTRUCTOR,   simple_fortran_constructor)
    CALL comin_callback_register(EP_ATM_WRITE_OUTPUT_BEFORE, simple_fortran_diagfct)
    CALL comin_callback_register(EP_DESTRUCTOR,              simple_fortran_destructor)

    ! get descriptive data structures
    p_patch    => comin_descrdata_get_domain(1)
    p_global   => comin_descrdata_get_global()
    p_simulation_interval    => comin_descrdata_get_simulation_interval()

    dtime_1 = comin_descrdata_get_timesteplength(1)
    dtime_2 = comin_descrdata_get_timesteplength(2)
    WRITE(text,"(A,F5.0,F5.0)") "     timesteplength from comin_descrdata_get_timesteplength", dtime_1, dtime_2
    CALL comin_print_info(text)
  END SUBROUTINE comin_main

  ! --------------------------------------------------------------------
  ! ComIn secondary constructor.
  ! --------------------------------------------------------------------
  SUBROUTINE simple_fortran_constructor()  BIND(C)
    TYPE(t_comin_var_descriptor)          :: var_desc
    INTEGER                               :: jg, hgrid_id, datatype
    LOGICAL                               :: tracer, multi_timelevel
    REAL(WP), POINTER                     :: tracer_slice(:,:,:)
    TYPE(c_ptr)                           :: it
    TYPE(t_comin_var_descriptor)          :: descriptor

    CALL comin_print_info("third party callback: secondary constructor.")
    CALL comin_print_info("third party callback: iterate over variable list:")

    it = comin_var_get_descr_list_head()
    DO WHILE (C_ASSOCIATED(it))
      CALL comin_var_get_descr_list_var_desc(it, descriptor)
      WRITE (text,*) "Variable found: ", TRIM(descriptor%name), &
        "(", descriptor%id, ")"
      CALL comin_print_info(text)
      it = comin_var_get_descr_list_next(it)
    END DO

    WRITE (text,*) "     ", pluginname, " - register some variables in some context"
    CALL comin_print_info(text)
    var_desc = t_comin_var_descriptor('pres', 1)
    CALL comin_var_get([EP_ATM_WRITE_OUTPUT_BEFORE], &
      &                 var_desc, COMIN_FLAG_READ, pres)
    IF (.NOT. pres%valid()) &
      &  CALL comin_plugin_finish("simple_fortran_constructor", "Internal error!")

    IF (ANY(pres%dim_semantics() /= [COMIN_DIM_SEMANTICS_NPROMA, &
         &                           COMIN_DIM_SEMANTICS_LEVEL,  &
         &                           COMIN_DIM_SEMANTICS_BLOCK,  &
         &                           COMIN_DIM_SEMANTICS_UNUSED,  &
         &                           COMIN_DIM_SEMANTICS_UNUSED ])) &
      &  CALL comin_plugin_finish("simple_fortran_constructor", "Dimension check failed!")

    var_desc = t_comin_var_descriptor('vn', 1)
    CALL comin_var_get([EP_ATM_WRITE_OUTPUT_BEFORE], &
      &                 var_desc, COMIN_FLAG_READ, vn)
    IF (.NOT. vn%valid()) &
      &  CALL comin_plugin_finish("simple_fortran_constructor vn", "Internal error!")

    CALL comin_metadata_get(var_desc, "hgrid_id", hgrid_id)
    IF (hgrid_id/=COMIN_HGRID_UNSTRUCTURED_EDGE) &
      &  CALL comin_plugin_finish("comin_var_get_metadata_hgrid", "Internal error!")
    CALL comin_metadata_get(var_desc, "multi_timelevel", multi_timelevel)

    ALLOCATE(qv(p_global%n_dom))
    DO jg =1, p_global%n_dom
      CALL comin_var_get([EP_ATM_WRITE_OUTPUT_BEFORE], &
           &             t_comin_var_descriptor(name='qv', id=jg), &
           &             COMIN_FLAG_READ, qv(jg))
      CALL comin_metadata_get(t_comin_var_descriptor(name='qv', id=jg), &
           &             "datatype", datatype)
      IF (datatype/=COMIN_VAR_DATATYPE_DOUBLE) &
           &  CALL comin_plugin_finish("simple_fortran_constructor", "Internal error!")
      IF (.NOT. qv(jg)%valid()) &
        &  CALL comin_plugin_finish("simple_fortran_constructor", "Internal error!")

      IF (ANY(qv(jg)%dim_semantics() /= [ COMIN_DIM_SEMANTICS_NPROMA, &
           &                              COMIN_DIM_SEMANTICS_LEVEL,  &
           &                              COMIN_DIM_SEMANTICS_BLOCK,  &
           &                              COMIN_DIM_SEMANTICS_UNUSED,  &
           &                              COMIN_DIM_SEMANTICS_UNUSED ] )) &
      &  CALL comin_plugin_finish("simple_fortran_constructor", "Dimension check failed!")
    END DO

    CALL comin_var_get([EP_ATM_WRITE_OUTPUT_BEFORE], &
         &                t_comin_var_descriptor(name="simple_fortran_var", id=1), &
         &                COMIN_FLAG_WRITE, simple_fortran_var)
    IF (.NOT. simple_fortran_var%valid()) &
      &  CALL comin_plugin_finish("simple_fortran_constructor", "Internal error!")

    CALL comin_metadata_get(t_comin_var_descriptor(name="simple_fortran_var", id=1), &
      &                     "tracer", tracer)

    CALL comin_var_get([EP_ATM_WRITE_OUTPUT_BEFORE], &
         &             t_comin_var_descriptor(name="simple_fortran_tracer", id=1), &
         &             COMIN_FLAG_WRITE, simple_fortran_tracer)
    IF (.NOT. simple_fortran_tracer%valid()) &
      &  CALL comin_plugin_finish("simple_fortran_constructor", "Internal error!")

    ! access the "simple_fortran_tracer" via a 3D array pointer:
    CALL simple_fortran_tracer%to_3d(tracer_slice)
    tracer_slice(1,9,1) = 1.0
  END SUBROUTINE simple_fortran_constructor

  ! --------------------------------------------------------------------
  ! ComIn callback function.
  ! --------------------------------------------------------------------
  SUBROUTINE simple_fortran_diagfct()  BIND(C)
    INTEGER                         :: ierr, domain_id, jg, comm, root
    REAL(wp)                        :: local_max, global_max
    REAL(wp), POINTER, DIMENSION(:,:,:,:,:) :: simple_fortran_ptr, &
         &                                     pres_ptr, qv_ptr

    CALL comin_print_info("third party callback: before output.")

    domain_id = comin_current_get_domain_id()

    IF (domain_id == COMIN_DOMAIN_OUTSIDE_LOOP) THEN
      CALL comin_print_info("currently not in domain loop")
    ELSE
      WRITE(text,'(a,i0)')  "      currently on domain ", domain_id
      CALL comin_print_info(text)
    END IF

    CALL simple_fortran_var%get_ptr(simple_fortran_ptr)
    CALL pres%get_ptr(pres_ptr)
    simple_fortran_ptr = pres_ptr + 42

    DO jg =1, p_global%n_dom
      comm      = comin_parallel_get_host_mpi_comm()
      CALL qv(jg)%get_ptr(qv_ptr)
      local_max = MAXVAL(qv_ptr)
      root      = 0
      CALL MPI_REDUCE(local_max, global_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, root, comm, ierr)
      WRITE(text,"(A,I8,A,F8.6)") "domain ", jg, ": global max = ", global_max
      CALL comin_print_info(text)
    END DO
  END SUBROUTINE simple_fortran_diagfct

  ! --------------------------------------------------------------------
  ! ComIn callback function.
  ! --------------------------------------------------------------------
  SUBROUTINE simple_fortran_destructor() BIND(C)
    CALL comin_print_info("third party callback: destructor.")
  END SUBROUTINE simple_fortran_destructor

END MODULE simple_fortran_plugin
