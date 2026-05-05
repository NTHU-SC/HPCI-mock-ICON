! --------------------------------------------------------------------
!> Example plugin for the ICON Community Interface (ComIn)
!  using the YAXT library to gather one ICON variable (temperature) on
!  task zero to calculate the global average near surface temperature.
!
!  @authors 08/2021 :: ICON Community Interface  <icon@dwd.de>
!
!  Note that in order to demonstrate ComIn's language interoperability,
!  a similar plugin has been implemented in C, see the subdirectory
!  "yaxt_c".
! --------------------------------------------------------------------
MODULE yaxt_fortran_plugin

  USE iso_c_binding,           ONLY : C_INT
  USE comin_plugin_interface,  ONLY : comin_callback_register    &
                                      , comin_var_get              &
                                      , t_comin_var_descriptor     &
                                      , t_comin_var_handle            &
                                      , comin_descrdata_get_domain &
                                      , t_comin_descrdata_domain   &
                                      , comin_descrdata_get_global &
                                      , t_comin_descrdata_global   &
                                      , t_comin_setup_version_info     &
                                      , comin_setup_get_version        &
                                      , EP_SECONDARY_CONSTRUCTOR       &
                                      , EP_DESTRUCTOR                  &
                                      , EP_ATM_WRITE_OUTPUT_BEFORE     &
                                      , COMIN_FLAG_READ                &
                                      , comin_parallel_get_host_mpi_rank &
                                      , comin_current_get_domain_id      &
                                      , comin_parallel_get_plugin_mpi_comm &
                                      , COMIN_DOMAIN_OUTSIDE_LOOP        &
                                      , t_comin_plugin_info              &
                                      , comin_current_get_plugin_info    &
                                      , comin_plugin_finish              &
                                      , comin_metadata_get               &
                                      , comin_error_check, comin_print_info    &
                                      , comin_print_debug
  USE yaxt,                     ONLY: xt_redist, xi => xt_int_kind

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: pluginname = "yaxt_fortran_plugin"

  !> working precision (will be compared to ComIn's and ICON's)
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12,307)
  TYPE(t_comin_setup_version_info) :: version

  TYPE(t_comin_var_handle)  :: temp ! temperature
  INTEGER                         :: rank
  CHARACTER(LEN=:), ALLOCATABLE   :: units

  !> access descriptive data structures
  TYPE(t_comin_descrdata_domain),     POINTER   :: p_patch
  TYPE(t_comin_descrdata_global),     POINTER   :: p_global

  !> yaxt related variables
  TYPE(xt_redist) :: yaxt_redist

  CHARACTER(LEN=120)            :: text

CONTAINS

  ! --------------------------------------------------------------------
  ! ComIn primary constructor.
  ! --------------------------------------------------------------------
  SUBROUTINE comin_main()  BIND(C)

    USE yaxt, ONLY: xt_initialize, xt_idxlist, xt_idxlist_delete &
                    , xt_xmap, xt_xmap_delete, xt_idxstripes_new, xt_idxvec_new &
                    , xt_stripe, xt_idxempty_new, xt_xmap_all2all_new &
                    , xt_redist_p2p_new, xt_initialized
    USE mpi,  ONLY: MPI_DOUBLE_PRECISION

    CHARACTER(LEN=*), PARAMETER   :: substr = 'comin_main (yaxt_fortran_plugin)'
    TYPE(t_comin_plugin_info)     :: this_plugin
    INTEGER                       :: p_all_comm ! communicator of all ICON tasks
    TYPE(xt_idxlist)              :: src_idxlist, tgt_idxlist
    TYPE(xt_xmap)                 :: xmap
    INTEGER(kind=xi), DIMENSION(:), ALLOCATABLE :: idxvec, idxmap
    INTEGER                       :: k

    !> get the rank of the current process and say hello to the world
    rank = comin_parallel_get_host_mpi_rank()
    CALL comin_print_info("setup")

    !> check, if the ComIn library version is compatible
    version = comin_setup_get_version()
    IF (version%version_no_major > 1)  THEN
      CALL comin_plugin_finish(substr, "incompatible ComIn library version!")
    END IF

    !> check plugin id
    CALL comin_current_get_plugin_info(this_plugin)
    WRITE (text,'(a,a,a,i4)') "     plugin " &
      , TRIM(this_plugin%name), " has id: ", this_plugin%id
    CALL comin_print_info(text)

    !> add requests for additional ICON variables
    !  not applicable for this example

    !> register callbacks
    CALL comin_callback_register(EP_SECONDARY_CONSTRUCTOR &
                                 , yaxt_fortran_constructor)
    CALL comin_callback_register(EP_ATM_WRITE_OUTPUT_BEFORE &
                                 , yaxt_fortran_gather)
    CALL comin_callback_register(EP_DESTRUCTOR &
                                 , yaxt_fortran_destructor)

    !> get descriptive data structures
    p_patch    => comin_descrdata_get_domain(1)
    p_global   => comin_descrdata_get_global()

    !> setup yaxt
    p_all_comm = comin_parallel_get_plugin_mpi_comm()
    IF (.NOT. xt_initialized()) THEN
      CALL comin_print_info("Initialize yaxt...")
      CALL xt_initialize(p_all_comm)
    ENDIF

    !> construct yaxt variables ...
    ! ... get halo info using decomp_domain :-
    ! ... 0=core, 1=shared edge with owned, 2=shared vertex with owned, <0: undefined
    idxmap = RESHAPE(p_patch%cells%decomp_domain &
                     , (/ SIZE(p_patch%cells%decomp_domain) /))
    ! ... get local ids of all core cells
    idxvec = INT(PACK( [(k,k=1,p_patch%cells%ncells)], idxmap == 0 ), xi)
    ! ... convert local ids to global ids
    idxvec = p_patch%cells%glb_index(idxvec)
    ! ... generate idxlist for all core cells
    src_idxlist = xt_idxvec_new(idxvec)

    IF (rank == 0) THEN
      tgt_idxlist = xt_idxstripes_new( &
                    (/ xt_stripe(1, 1, p_patch%cells%ncells_global) /))
    ELSE
      ! ... empty on all other pe
      tgt_idxlist = xt_idxempty_new()
    ENDIF

    ! ... create exchange map ...
    xmap = xt_xmap_all2all_new(src_idxlist, tgt_idxlist, p_all_comm)

    ! ... create redistribution instance for DP ...
    yaxt_redist = xt_redist_p2p_new(xmap, MPI_DOUBLE_PRECISION)

    ! ... clean up
    CALL xt_xmap_delete(xmap)
    CALL xt_idxlist_delete(src_idxlist)
    CALL xt_idxlist_delete(tgt_idxlist)
    DEALLOCATE(idxvec)
    DEALLOCATE(idxmap)

  END SUBROUTINE comin_main

  ! --------------------------------------------------------------------
  ! ComIn secondary constructor.
  ! --------------------------------------------------------------------
  SUBROUTINE yaxt_fortran_constructor()  BIND(C)

    USE yaxt, ONLY: xt_idxlist, xt_idxlist_delete &
                    , xt_xmap, xt_xmap_delete &
                    , xt_redist

    CHARACTER(LEN=*), PARAMETER :: substr = 'yaxt_fortran_constructor (yaxt_fortran_plugin)'
    TYPE(t_comin_var_descriptor)         :: var_desc

    CALL comin_print_info("secondary constructor")

    CALL comin_print_info("request temperature")
    var_desc = t_comin_var_descriptor('temp', 1)
    CALL comin_var_get([EP_ATM_WRITE_OUTPUT_BEFORE], &
                       var_desc, COMIN_FLAG_READ, temp)

    CALL comin_metadata_get(var_desc, 'units', units)

  END SUBROUTINE yaxt_fortran_constructor

  ! --------------------------------------------------------------------
  ! ComIn callback function.
  ! --------------------------------------------------------------------
  SUBROUTINE yaxt_fortran_gather()  BIND(C)

    USE yaxt,          ONLY: xt_redist_s_exchange1
    USE iso_c_binding, ONLY: C_LOC

    CHARACTER(LEN=*), PARAMETER :: substr = 'yaxt_fortran_gather (yaxt_fortran_plugin)'
    TYPE(t_comin_plugin_info)                :: this_plugin
    INTEGER                                  :: domain_id
    REAL(kind=wp), DIMENSION(:,:),   POINTER :: src
    REAL(kind=wp), DIMENSION(:,:,:), POINTER :: src3d
    REAL(kind=wp), DIMENSION(:),     POINTER :: tgt, area

    CALL comin_print_info("callback before output")

    !> check plugin id
    CALL comin_current_get_plugin_info(this_plugin)

    domain_id = comin_current_get_domain_id()
    IF (domain_id == COMIN_DOMAIN_OUTSIDE_LOOP) THEN
      CALL comin_print_debug("currently not in domain loop")
    ELSE
      WRITE(text,'(a,a,i0)') "currently on domain ", domain_id
      CALL comin_print_debug(text)
    END IF

    !> reset pointers
    NULLIFY(src, src3d, tgt, area)

    CALL temp%to_3d(src3d)

    !> extract near surface temperature
    ALLOCATE(src(p_global%nproma,p_patch%cells%nblks))
    src(:,:) = src3d(:,p_patch%nlev,:)

    !> allocate local space to gather global information
    ALLOCATE(tgt(p_patch%cells%ncells_global))
    ALLOCATE(area(p_patch%cells%ncells_global))

    !> gather information in rank zero
    CALL xt_redist_s_exchange1(yaxt_redist, C_LOC(src), C_LOC(tgt))
    CALL xt_redist_s_exchange1(yaxt_redist, C_LOC(p_patch%cells%area) &
                               , C_LOC(area) )

    IF (rank == 0) &
      WRITE(0,*) substr, ': global average temperature is ' &
      , SUM(tgt*area)/SUM(area), TRIM(units)

    !> clean up memory
    DEALLOCATE(src)
    DEALLOCATE(tgt)
    DEALLOCATE(area)

  END SUBROUTINE yaxt_fortran_gather

  ! --------------------------------------------------------------------
  ! ComIn callback function.
  ! --------------------------------------------------------------------
  SUBROUTINE yaxt_fortran_destructor() BIND(C)

    USE yaxt, ONLY: xt_finalize, xt_redist_delete

    CALL comin_print_info("destructor")

    !> free yaxt related memory
    CALL xt_redist_delete(yaxt_redist)

    !> finalize yaxt
    CALL xt_finalize()

  END SUBROUTINE yaxt_fortran_destructor

END MODULE yaxt_fortran_plugin
