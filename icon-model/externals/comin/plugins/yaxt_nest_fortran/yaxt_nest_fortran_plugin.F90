! --------------------------------------------------------------------
!> Example plugin for the ICON Community Interface (ComIn)
!  using the YAXT library to set up child/parent communication and
!  get the owner PE of cells...
!
!  @authors 07/2025 :: ICON Community Interface <comin@icon-model.org>
!
! --------------------------------------------------------------------
MODULE yaxt_nest_fortran_plugin

  USE comin_plugin_interface, ONLY: comin_callback_register &
                                    , comin_var_get &
                                    , t_comin_var_descriptor &
                                    , t_comin_var_handle &
                                    , comin_descrdata_get_domain &
                                    , t_comin_descrdata_domain &
                                    , comin_descrdata_get_global &
                                    , t_comin_descrdata_global &
                                    , t_comin_setup_version_info &
                                    , comin_setup_get_version &
                                    , EP_SECONDARY_CONSTRUCTOR &
                                    , EP_DESTRUCTOR &
                                    , EP_ATM_INTEGRATE_END &
                                    , COMIN_FLAG_READ &
                                    , COMIN_FLAG_WRITE &
                                    , COMIN_ZAXIS_2D &
                                    , comin_parallel_get_host_mpi_rank &
                                    , comin_current_get_domain_id &
                                    , comin_parallel_get_plugin_mpi_comm &
                                    , COMIN_DOMAIN_OUTSIDE_LOOP &
                                    , t_comin_plugin_info &
                                    , comin_current_get_plugin_info &
                                    , comin_plugin_finish &
                                    , comin_metadata_get &
                                    , comin_metadata_set &
                                    , comin_var_request_add &
                                    , comin_error_check, comin_print_info &
                                    , comin_print_debug
  USE yaxt, ONLY: xt_redist

  IMPLICIT NONE

  CHARACTER(LEN=*), PARAMETER :: pluginname = "yaxt_nest_fortran_plugin"

  !> working precision (will be compared to ComIn's and ICON's)
  INTEGER, PARAMETER :: wp = SELECTED_REAL_KIND(12, 307)
  TYPE(t_comin_setup_version_info) :: version

  TYPE(t_comin_var_handle), ALLOCATABLE  :: temp(:) ! temperature
  TYPE(t_comin_var_handle), ALLOCATABLE  :: temp_diag(:) ! temperature
  INTEGER                         :: rank
  INTEGER                         :: n_dom, nproma

  !> access descriptive data structures
  TYPE(t_comin_descrdata_global), POINTER   :: p_global
  TYPE t_patch
    TYPE(t_comin_descrdata_domain), POINTER   :: ptr
  END TYPE t_patch
  TYPE(t_patch), ALLOCATABLE :: p_patch(:)

  !> yaxt related variables
  TYPE(xt_redist), ALLOCATABLE :: yaxt_redist(:)

  CHARACTER(LEN=120)            :: text

CONTAINS

  ! --------------------------------------------------------------------
  ! ComIn primary constructor.
  ! --------------------------------------------------------------------
  SUBROUTINE comin_main() BIND(C)

    USE yaxt, ONLY: xt_initialize, xt_idxlist, xt_idxlist_delete &
                    , xt_xmap, xt_xmap_delete, xt_idxstripes_new, xt_idxvec_new &
                    , xt_stripe, xt_idxempty_new, xt_xmap_dist_dir_new &
                    , xt_redist_p2p_new, xt_initialized, xi => xt_int_kind
    USE mpi, ONLY: MPI_DOUBLE_PRECISION

    CHARACTER(LEN=*), PARAMETER   :: substr = 'comin_main (yaxt_nest_fortran_plugin)'
    TYPE(t_comin_plugin_info)     :: this_plugin
    TYPE(t_comin_var_descriptor)  :: temp_d
    INTEGER                       :: p_all_comm ! communicator of all ICON tasks
    TYPE(xt_idxlist)              :: src_idxlist, tgt_idxlist
    TYPE(xt_xmap)                 :: xmap
    INTEGER                       :: ii, ic, cid
    INTEGER, DIMENSION(:), ALLOCATABLE :: src_indices, parent_indices, tgt_indices

    !> get the rank of the current process and say hello to the world
    rank = comin_parallel_get_host_mpi_rank()
    CALL comin_print_info("setup")

    !> check, if the ComIn library version is compatible
    version = comin_setup_get_version()
    IF (version%version_no_major > 1) THEN
      CALL comin_plugin_finish(substr, "incompatible ComIn library version!")
    END IF

    !> check plugin id
    CALL comin_current_get_plugin_info(this_plugin)
    WRITE (text, '(a,a,a,i4)') "     plugin " &
      , TRIM(this_plugin%name), " has id: ", this_plugin%id
    CALL comin_print_info(text)

    !> add requests for additional ICON variables
    !  not applicable for this example

    !> register callbacks
    CALL comin_callback_register(EP_SECONDARY_CONSTRUCTOR &
                                 , yaxt_nest_fortran_constructor)
    CALL comin_callback_register(EP_ATM_INTEGRATE_END &
                                 , yaxt_nest_fortran_exchange)
    CALL comin_callback_register(EP_DESTRUCTOR &
                                 , yaxt_nest_fortran_destructor)

    !> get descriptive data structures
    p_global => comin_descrdata_get_global()
    nproma = p_global%nproma
    n_dom = p_global%n_dom
    IF (n_dom == 1) THEN
      CALL comin_plugin_finish(substr, "only applicable for nested domain setups")
    END IF
    ALLOCATE (p_patch(n_dom))
    DO ii = 1, n_dom
      p_patch(ii)%ptr => comin_descrdata_get_domain(ii)
    END DO

    !> setup yaxt
    p_all_comm = comin_parallel_get_plugin_mpi_comm()
    IF (.NOT. xt_initialized()) THEN
      CALL comin_print_info("Initialize yaxt...")
      CALL xt_initialize(p_all_comm)
    END IF

    !> construct yaxt variables ...
    ALLOCATE (yaxt_redist(n_dom))
    DO ii = 1, n_dom
      DO ic = 1, p_patch(ii)%ptr%n_childdom
        cid = p_patch(ii)%ptr%child_id(ic)
        ALLOCATE (src_indices(COUNT(p_patch(ii)%ptr%cells%decomp_domain == 0)))
        ALLOCATE (parent_indices(SIZE(p_patch(cid)%ptr%cells%parent_glb_idx)))
        ALLOCATE (tgt_indices(COUNT(p_patch(cid)%ptr%cells%decomp_domain >= 0)))
        ! assign values to source indices: if reshaped array of child_id for this
        ! parent equals the current child id (cid): use global index, else: -1
        src_indices = MERGE(p_patch(ii)%ptr%cells%glb_index, -1_xi              &
          &              , RESHAPE(p_patch(ii)%ptr%cells%child_id             &
          &                , (/SIZE(p_patch(ii)%ptr%cells%child_id)/)) == cid)
        ! compute 1D global parent index for each cell
        parent_indices = idx_1d(RESHAPE(p_patch(cid)%ptr%cells%parent_glb_idx   &
          &                , (/SIZE(p_patch(cid)%ptr%cells%parent_glb_idx)/))   &
          &              , RESHAPE(p_patch(cid)%ptr%cells%parent_glb_blk        &
          &                , (/SIZE(p_patch(cid)%ptr%cells%parent_glb_blk)/)))
        ! restrict target indices to the cells owned by the task
        tgt_indices = pack(parent_indices                                       &
          &              , RESHAPE(p_patch(cid)%ptr%cells%decomp_domain         &
          &                , (/SIZE(p_patch(cid)%ptr%cells%decomp_domain)/)) >= 0)

        src_idxlist = xt_idxvec_new(src_indices, SIZE(src_indices, 1))
        tgt_idxlist = xt_idxvec_new(tgt_indices, SIZE(tgt_indices, 1))

        ! ... create exchange map ...
        xmap = xt_xmap_dist_dir_new(src_idxlist, tgt_idxlist, p_all_comm)

        ! ... create redistribution instance for DP ...
        yaxt_redist(cid) = xt_redist_p2p_new(xmap, MPI_DOUBLE_PRECISION)

        ! ... clean up
        CALL xt_xmap_delete(xmap)
        CALL xt_idxlist_delete(src_idxlist)
        CALL xt_idxlist_delete(tgt_idxlist)
        DEALLOCATE (src_indices, parent_indices, tgt_indices)
      END DO

      temp_d = t_comin_var_descriptor(id=ii, name="temp_yaxt_diag")
      CALL comin_var_request_add(temp_d, .TRUE.)
      CALL comin_metadata_set(temp_d, "zaxis_id", COMIN_ZAXIS_2D)
      CALL comin_metadata_set(temp_d, "tracer", .FALSE.)
      CALL comin_metadata_set(temp_d, "restart", .FALSE.)
      CALL comin_metadata_set(temp_d, "units", "K")
    END DO

  END SUBROUTINE comin_main

  ! --------------------------------------------------------------------
  ! ComIn secondary constructor.
  ! --------------------------------------------------------------------
  SUBROUTINE yaxt_nest_fortran_constructor() BIND(C)

    USE yaxt, ONLY: xt_idxlist, xt_idxlist_delete &
                    , xt_xmap, xt_xmap_delete &
                    , xt_redist

    INTEGER :: ii, pid

    CALL comin_print_info("secondary constructor")

    ALLOCATE (temp(n_dom), temp_diag(n_dom))
    DO ii = 1, n_dom
      CALL comin_print_info("request temperature")
      CALL comin_var_get([EP_ATM_INTEGRATE_END], &
        &                t_comin_var_descriptor(name='temp', id=ii), COMIN_FLAG_READ, temp(ii))
      CALL comin_print_info("request diagnostic temperature")
      pid = p_patch(ii)%ptr%parent_id
      IF (pid <= 0) CYCLE
      CALL comin_var_get([EP_ATM_INTEGRATE_END], &
        &                t_comin_var_descriptor(name='temp_yaxt_diag', id=ii), COMIN_FLAG_WRITE, temp_diag(ii))
    END DO
  END SUBROUTINE yaxt_nest_fortran_constructor

  ! --------------------------------------------------------------------
  ! ComIn callback function.
  ! --------------------------------------------------------------------
  SUBROUTINE yaxt_nest_fortran_exchange() BIND(C)

    USE yaxt, ONLY: xt_redist_s_exchange1
    USE iso_c_binding, ONLY: C_LOC

    INTEGER :: domain_id, pid
    REAL(kind=wp), DIMENSION(:,:),     POINTER :: src
    REAL(kind=wp), DIMENSION(:, :, :), POINTER :: temp3d, tgt

    CALL comin_print_info("callback before output")

    domain_id = comin_current_get_domain_id()
    IF (domain_id == COMIN_DOMAIN_OUTSIDE_LOOP) THEN
      CALL comin_print_debug("currently not in domain loop")
      RETURN
    ELSE
      WRITE (text, '(a,a,i0)') "currently on domain ", domain_id
      CALL comin_print_debug(text)
    END IF

    pid = p_patch(domain_id)%ptr%parent_id
    IF (pid <= 0) RETURN

    CALL temp_diag(domain_id)%to_3d(tgt)
    tgt = 0._wp

    CALL temp(pid)%to_3d(temp3d)
    ALLOCATE (src(p_global%nproma,p_patch(pid)%ptr%cells%nblks))
    ! select the near surface temperature as source variable field
    src(:,:) = temp3d(:,p_patch(pid)%ptr%nlev,:)

    !> exchange parent to child
    CALL xt_redist_s_exchange1(yaxt_redist(domain_id), C_LOC(src), C_LOC(tgt))
    ! find the min/max temperatures, filter out 0 values
    WRITE(text, '(a,2(1x,f7.3))') ': min/max temperature src' &
      , MINVAL(PACK(src, src > 0._wp)), MAXVAL(PACK(src, src > 0._wp))
    CALL comin_print_info(text)
    WRITE(text, '(a,2(1x,f7.3))') ': min/max temperature tgt' &
      , MINVAL(PACK(tgt, tgt > 0._wp)), MAXVAL(PACK(tgt, tgt > 0._wp))
    CALL comin_print_info(text)
    DEALLOCATE(src)

  END SUBROUTINE yaxt_nest_fortran_exchange

  ! --------------------------------------------------------------------
  ! ComIn callback function.
  ! --------------------------------------------------------------------
  SUBROUTINE yaxt_nest_fortran_destructor() BIND(C)

    USE yaxt, ONLY: xt_finalize, xt_redist_delete

    INTEGER :: ii, pid

    CALL comin_print_info("destructor")

    !> free yaxt related memory
    DO ii = 1, n_dom
      pid = p_patch(ii)%ptr%parent_id
      IF (pid <= 0) CYCLE
      CALL xt_redist_delete(yaxt_redist(ii))
    END DO
    DEALLOCATE (temp, temp_diag, p_patch)

    !> finalize yaxt
    CALL xt_finalize()

  END SUBROUTINE yaxt_nest_fortran_destructor

  ELEMENTAL INTEGER FUNCTION idx_1d(jl, jb)
    INTEGER, INTENT(IN) :: jl, jb
    IF (jb <= 0) THEN
      idx_1d = 0 ! This covers the special case nproma==1,jb=0,jl=1
      ! All other cases are invalid and get also a 0 returned
    ELSE
      idx_1d = SIGN((jb - 1)*nproma + ABS(jl), jl)
    END IF
  END FUNCTION idx_1d

END MODULE yaxt_nest_fortran_plugin
