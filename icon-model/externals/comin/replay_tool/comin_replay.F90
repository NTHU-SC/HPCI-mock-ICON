!> @file comin_replay.F90
!! @brief Imitates a host model by reading data from a netcdf file
!! (created with the `comin_run_recorder_plugin`) and plays the data back to comin.
!
!  @authors 01/2024 :: ICON Community Interface  <comin@icon-model.org>
!
!  SPDX-License-Identifier: BSD-3-Clause
!
!  See LICENSES for license information.
!  Where software is supplied by third parties, it is indicated in the
!  headers of the routines.

MODULE mo_comin_replay
  USE mpi,                     ONLY: MPI_INIT, MPI_COMM_WORLD, MPI_COMM_RANK, MPI_COMM_SIZE,        &
       &                             MPI_FINALIZE, MPI_SUCCESS, MPI_ABORT, MPI_DOUBLE_PRECISION,    &
       &                             MPI_FLOAT, MPI_INTEGER
  USE comin_ftnlist_ifc,       ONLY: comin_ftnlist_new, comin_ftnlist_push_back, comin_ftnlist_iterator_begin, &
       &                             comin_ftnlist_iterator_next, comin_ftnlist_iterator_value,                &
       &                             comin_ftnlist_iterator_delete, comin_ftnlist_is_end, comin_ftnlist_delete
  USE comin_host_interface,    ONLY: t_comin_var_handle, t_comin_setup_version_info,                   &
       &                             t_comin_descrdata_global, t_comin_descrdata_domain,            &
       &                             t_comin_descrdata_simulation_interval,                         &
       &                             t_comin_plugin_description, comin_setup_init,                  &
       &                             comin_setup_errhandler, comin_setup_set_verbosity_level,       &
       &                             comin_setup_check, comin_setup_get_version,                    &
       &                             comin_descrdata_set_global, comin_descrdata_set_domain,        &
       &                             comin_descrdata_set_timesteplength,                            &
       &                             comin_descrdata_set_simulation_interval,                       &
       &                             comin_parallel_mpi_handshake, comin_plugin_primaryconstructor, &
       &                             comin_current_set_datetime, comin_callback_context_call,       &
       &                             t_comin_var_descriptor, comin_request_get_list,                &
       &                             comin_var_list_append, comin_metadata_set, COMIN_ZAXIS_2D,     &
       &                             COMIN_ZAXIS_3D, COMIN_ZAXIS_3D_HALF,mpi_handshake, EP_FINISH,  &
       &                             COMIN_DOMAIN_OUTSIDE_LOOP, EP_ATM_YAC_DEFCOMP_AFTER,           &
       &                             EP_ATM_YAC_SYNCDEF_AFTER, EP_ATM_YAC_ENDDEF_AFTER,             &
       &                             EP_ATM_TIMELOOP_BEFORE, t_comin_var_metadata_iterator,         &
       &                             comin_metadata_get_iterator, COMIN_METADATA_TYPEID_INTEGER,    &
       &                             COMIN_METADATA_TYPEID_REAL, COMIN_METADATA_TYPEID_CHARACTER,   &
       &                             COMIN_METADATA_TYPEID_LOGICAL, comin_metadata_get_or,          &
       &                             t_comin_request_item,                                          &
       &                             comin_descrdata_set_fct_glb2loc_cell,                          &
       &                             comin_var_set_sync_halo, comin_var_list_finalize,              &
       &                             comin_var_request_list_finalize, comin_var_descr_list_finalize,&
       &                             comin_descrdata_finalize, comin_setup_finalize,                &
       &                             COMIN_HGRID_UNSTRUCTURED_CELL, COMIN_HGRID_UNSTRUCTURED_EDGE,  &
       &                             COMIN_HGRID_UNSTRUCTURED_VERTEX,                               &
       &                             COMIN_VAR_DATATYPE_DOUBLE, COMIN_VAR_DATATYPE_FLOAT,           &
       &                             COMIN_VAR_DATATYPE_INT,                                        &
       &                             COMIN_DIM_SEMANTICS_NPROMA, COMIN_DIM_SEMANTICS_LEVEL,         &
       &                             COMIN_DIM_SEMANTICS_BLOCK, COMIN_DIM_SEMANTICS_UNUSED,         &
       &                             COMIN_DIM_SEMANTICS_CONTAINER, comin_metadata_get
  USE netcdf,                  ONLY: nf90_open, nf90_get_att, nf90_inq_ncid, nf90_get_var,          &
       &                             nf90_close, NF90_GLOBAL, NF90_NOWRITE
  USE netcdf_utils,            ONLY: nf90, nf90_utils_def_var, nf90_utils_get_shape
  USE comin_descrdata_load,    ONLY: comin_descrdata_load_domain, comin_descrdata_load_global
  USE iso_c_binding,           ONLY: C_DOUBLE, C_NULL_PTR, c_ptr, c_f_pointer, C_LOC, C_FLOAT, C_INT
  USE utils,                   ONLY: int2string
#ifdef ENABLE_YAC
  USE yac,                     ONLY: yac_finit_comm, yac_fdef_calendar, YAC_PROLEPTIC_GREGORIAN, &
       &                             yac_fdef_datetime, yac_fdef_comp, yac_fsync_def, yac_fenddef
#endif

  USE global_index_lookup,     ONLY: t_glb2loc_index_lookup, setup_glb2loc, glb2loc_lookup

#ifdef ENABLE_YAXT
  USE yaxt,                    ONLY: xt_redist, xi => xt_int_kind
#endif

  IMPLICIT NONE

  PUBLIC:: main

  TYPE :: replay_var
    REAL(C_DOUBLE), POINTER :: data_ptr_dp(:,:,:,:,:)
    REAL(C_FLOAT), POINTER :: data_ptr_sp(:,:,:,:,:)
    INTEGER(C_INT), POINTER :: data_ptr_i(:,:,:,:,:)
  END TYPE replay_var

  INTEGER, PARAMETER :: wp = C_DOUBLE
  INTEGER :: ierr
  INTEGER :: ncid, grp_ncid, varid
  INTEGER :: host_comm, comin_comm, host_rank, host_size
  INTEGER :: file_host_rank, file_host_size
  INTEGER :: i, nplugins
  CHARACTER(len=256) :: namelist_filename, filename
  TYPE(t_comin_setup_version_info) :: comin_version
  INTEGER :: file_comin_version(3)
  INTEGER :: output_unit
  CHARACTER(:), ALLOCATABLE :: output_filename

  TYPE(t_comin_descrdata_global)                        :: comin_global
  TYPE(t_comin_descrdata_domain), ALLOCATABLE, SAVE     :: comin_domain(:)
  TYPE(t_comin_descrdata_simulation_interval)           :: comin_simulation_interval

  TYPE(replay_var), TARGET, ALLOCATABLE :: vars(:)

  TYPE(t_comin_plugin_description) :: plugin_list(16) !< list of dynamic libs (max: 16)
  NAMELIST /comin_nml/ plugin_list

  CHARACTER(len=256) :: replay_data_path
  INTEGER :: msg_level = 13
  CHARACTER(len=256) :: output_dir = ""
  NAMELIST /replay_tool_nml/ replay_data_path, msg_level, output_dir

  INTEGER :: shap(2)
  INTEGER, ALLOCATABLE :: current_ep(:)
  INTEGER, ALLOCATABLE :: current_domain_id(:)
  CHARACTER(LEN=32), ALLOCATABLE :: current_datetime(:)
  REAL(C_DOUBLE) :: dt
#ifdef ENABLE_YAC
  INTEGER :: yac_instance_id, yac_comp_id
#endif

  TYPE(t_glb2loc_index_lookup), ALLOCATABLE, SAVE :: glb2loc_lookup_cell(:)

#ifdef ENABLE_YAXT
  TYPE(xt_redist),  DIMENSION(:), ALLOCATABLE, SAVE    :: yaxt_redist_dp
  TYPE(xt_redist),  DIMENSION(:), ALLOCATABLE, SAVE    :: yaxt_redist_sp
  TYPE(xt_redist),  DIMENSION(:), ALLOCATABLE, SAVE    :: yaxt_redist_i
  TYPE :: t_idxvec
    INTEGER(kind=xi), DIMENSION(:), ALLOCATABLE        :: idxvec
  END TYPE t_idxvec
  TYPE(t_idxvec), DIMENSION(:), ALLOCATABLE, SAVE      :: idxvec_domain
#endif

CONTAINS
  SUBROUTINE main()

    INTEGER :: jg

    CALL start_mpi()

    IF ( COMMAND_ARGUMENT_COUNT() == 1 ) THEN
      CALL GET_COMMAND_ARGUMENT(1, namelist_filename)
    ELSE
      CALL finish("replay", "Exactly one argument is necessary. Please provide an input namelist.")
    ENDIF

    OPEN(1, file=namelist_filename)
    READ(1, nml=replay_tool_nml)
    READ(1, nml=comin_nml)
    CLOSE(1)

    IF (LEN_TRIM(output_dir) > 0) THEN
      output_filename = TRIM(output_dir) // "/" // int2string(host_rank)
      OPEN(newunit=output_unit, file=output_filename, status="replace", action="write")
      CALL comin_setup_init(.TRUE., output_unit)
    ELSE
      output_unit = 0
      CALL comin_setup_init(host_rank == 0)
    END IF

    filename = TRIM(replay_data_path) // int2string(host_rank) // ".nc"
    IF(host_rank == 0) WRITE (output_unit,*) "reading ", TRIM(filename)
    CALL nf90(nf90_open(TRIM(filename), NF90_NOWRITE, ncid))

    CALL comin_setup_errhandler(finish)
    CALL comin_setup_set_verbosity_level(iverbosity=msg_level)
    CALL comin_setup_check("replay", wp)

    ! check version of file and library
    comin_version = comin_setup_get_version()
    CALL nf90(nf90_get_att(ncid, NF90_GLOBAL, "comin_version", file_comin_version))
    IF(host_rank == 0) WRITE (output_unit,*) "using comin version          ", comin_version
    IF(host_rank == 0) WRITE (output_unit,*) "file was written with version ", file_comin_version
    IF(file_comin_version(1) /= comin_version%version_no_major .OR. &
       file_comin_version(2) /= comin_version%version_no_minor) &
      CALL finish("replay", "Incompatible versions")

    ! check host comm
    CALL nf90(nf90_get_att(ncid, NF90_GLOBAL, "host_comm_rank", file_host_rank))
    IF(host_rank /= file_host_rank) CALL finish("comin_replay", "Wrong host rank number in file")
    CALL nf90(nf90_get_att(ncid, NF90_GLOBAL, "host_comm_size", file_host_size))
    IF(host_size /= file_host_size) CALL finish("comin_replay", "Wrong host comm size in file")

    CALL comin_descrdata_load_global(ncid, comin_global)
#ifdef ENABLE_YAC
    IF (comin_global%yac_instance_id /= -1) THEN
      comin_global%yac_instance_id = yac_instance_id
    ENDIF
#endif

    ! This will be replaced if the replay tool also supports devices
    comin_global%device_name      = "Test device"
    comin_global%device_vendor    = "Dummy vendor"
    comin_global%device_driver    = "Dummy driver"

    CALL comin_descrdata_set_global(comin_global)
    ALLOCATE(comin_domain(comin_global%n_dom))
    ALLOCATE(glb2loc_lookup_cell(comin_global%n_dom))
    DO jg = 1,comin_global%n_dom
      CALL nf90(nf90_inq_ncid(ncid, "domain_"//int2string(jg), grp_ncid))
      CALL comin_descrdata_load_domain(grp_ncid, comin_domain(jg))
      CALL nf90(nf90_get_att(grp_ncid, NF90_GLOBAL, "timestep_length", dt))
      CALL comin_descrdata_set_timesteplength(jg, dt)
      CALL setup_glb2loc(comin_domain(jg)%cells%glb_index, glb2loc_lookup_cell(jg))
    END DO
    CALL comin_descrdata_set_fct_glb2loc_cell(glb2loc_lookup_wrapper_cell)

    CALL comin_descrdata_set_domain(comin_domain)

    CALL nf90(nf90_get_att(ncid, NF90_GLOBAL, "exp_start", comin_simulation_interval%exp_start))
    CALL nf90(nf90_get_att(ncid, NF90_GLOBAL, "exp_stop", comin_simulation_interval%exp_stop))
    CALL nf90(nf90_get_att(ncid, NF90_GLOBAL, "run_start", comin_simulation_interval%run_start))
    CALL nf90(nf90_get_att(ncid, NF90_GLOBAL, "run_stop", comin_simulation_interval%run_stop))
    CALL comin_descrdata_set_simulation_interval(comin_simulation_interval)
#ifdef ENABLE_YAC
    CALL yac_fdef_datetime(yac_instance_id, &
         & TRIM(comin_simulation_interval%run_start), &
         & TRIM(comin_simulation_interval%run_stop))
#endif

    DO i=1,SIZE(plugin_list)
      ! plugin and primary_constructor are optional. Only if both are not specified the plugin is considered as empty
      IF (len_TRIM(plugin_list(i)%plugin_library) == 0 .AND. &
           &   TRIM(plugin_list(i)%primary_constructor) == "comin_main") EXIT
    END DO

    nplugins = i-1
    CALL comin_parallel_mpi_handshake(comin_comm, plugin_list(1:nplugins)%comm, "atmo") ! currently the minimial example only emulates the atmo
#ifdef ENABLE_YAXT
    ALLOCATE(yaxt_redist_dp(comin_global%n_dom))
    ALLOCATE(yaxt_redist_sp(comin_global%n_dom))
    ALLOCATE(yaxt_redist_i(comin_global%n_dom))
    ALLOCATE(idxvec_domain(comin_global%n_dom))
    CALL set_yaxt_redist()
    CALL comin_var_set_sync_halo(replay_halo_sync_variable)
#endif
    CALL comin_plugin_primaryconstructor(plugin_list(1:nplugins))

    CALL add_variables(comin_global, comin_domain)
    ! run callback loop
    shap(1:1) = nf90_utils_get_shape(ncid, "current_ep", 1, varid)
    ALLOCATE(current_ep(shap(1)))
    CALL nf90(nf90_get_var(ncid, varid, current_ep))
    shap(1:1) = nf90_utils_get_shape(ncid, "current_domain_id", 1, varid)
    ALLOCATE(current_domain_id(shap(1)))
    CALL nf90(nf90_get_var(ncid, varid, current_domain_id))

    shap = nf90_utils_get_shape(ncid, "current_datetime", 2, varid)
    ALLOCATE(current_datetime(shap(2)))
    CALL nf90(nf90_get_var(ncid, varid, current_datetime))

    DO i=1,SIZE(current_ep)
      !    WRITE(output_unit,*) "time is " // TRIM(current_datetime(i)) // &
      !         & " calling EP " // int2string(current_ep(i)) // &
      !         & " on domain " // int2string(current_domain_id(i))
      CALL comin_current_set_datetime(trim_null(current_datetime(i)))
      CALL yac_routines(current_ep(i))
      CALL comin_callback_context_call(current_ep(i), current_domain_id(i), .FALSE.)
    END DO

    CALL comin_var_list_finalize()
    CALL comin_var_request_list_finalize()
    CALL comin_var_descr_list_finalize()
    CALL comin_descrdata_finalize()
    CALL comin_setup_finalize()

    CALL nf90(nf90_close(ncid))

    CALL MPI_FINALIZE(ierr); CALL handle_mpi_errcode(ierr)

  CONTAINS

    SUBROUTINE start_mpi()
      INTEGER :: ierr
#ifdef ENABLE_YAC
      CHARACTER(LEN=256), PARAMETER :: group_names(3) = [&
                                       "replay", &
                                       "comin ", &
                                       "yac   "]
#else
      CHARACTER(LEN=256), PARAMETER :: group_names(2) = [&
                                       "replay", &
                                       "comin "]
#endif
      INTEGER :: group_comms(3)

      CALL MPI_INIT (ierr); CALL handle_mpi_errcode(ierr)
      CALL mpi_handshake(MPI_COMM_WORLD, group_names, group_comms)

#ifdef ENABLE_YAC
      CALL yac_fdef_calendar(YAC_PROLEPTIC_GREGORIAN)
      CALL yac_finit_comm(group_comms(3), yac_instance_id)
#endif

      host_comm = group_comms(1)
      comin_comm = group_comms(2)
      CALL MPI_COMM_RANK(host_comm, host_rank, ierr); CALL handle_mpi_errcode(ierr)

      CALL MPI_COMM_SIZE(host_comm, host_size, ierr); CALL handle_mpi_errcode(ierr)
      IF (host_rank == 0) THEN
        WRITE (output_unit,*) "running with ", host_size, " MPI tasks"
      END IF
    END SUBROUTINE start_mpi

    !> Utility function.
    SUBROUTINE handle_mpi_errcode(errcode)
      INTEGER, INTENT(IN) :: errcode
      IF (errcode .NE. MPI_SUCCESS) THEN
        CALL finish("replay", "Error in MPI program. Terminating.")
      END IF
    END SUBROUTINE handle_mpi_errcode

    SUBROUTINE finish(routine, text)
      CHARACTER(LEN=*), INTENT(IN) :: routine
      CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: text
      INTEGER :: ierr, iexit

      CALL comin_callback_context_call(EP_FINISH, COMIN_DOMAIN_OUTSIDE_LOOP, .FALSE.)

      iexit = -1
      IF (host_rank == 0)  WRITE (output_unit,*) routine, ": ", text, iexit
      CALL MPI_ABORT(MPI_COMM_WORLD, iexit, ierr)
    END SUBROUTINE finish

    SUBROUTINE add_variables(global, domain)
      TYPE(t_comin_descrdata_global), INTENT(IN) :: global
      TYPE(t_comin_descrdata_domain), INTENT(IN) :: domain(:)
      TYPE(t_comin_var_descriptor) :: descr
      TYPE(t_comin_var_metadata_iterator) :: metadata_it
      INTEGER :: no_vars = 0, no_tracers = 0, dim_semantics(5), dimshape(5), dom_idx, ncontained
      INTEGER :: tracer_counter(global%n_dom), jg
      LOGICAL :: tracer, logbuff
      INTEGER :: zaxis_id, intbuff, hgrid_id, nblks, datatype
      REAL(C_DOUBLE) :: realbuff
      CHARACTER(LEN=:),ALLOCATABLE :: current_metadata_key, charbuff
      TYPE(c_ptr) :: list, it, cptr
      TYPE(t_comin_request_item), POINTER :: item
      TYPE(t_comin_var_handle) :: dummy

      list = comin_request_get_list()
      CALL comin_ftnlist_iterator_begin(list, it)
      DO WHILE (.NOT. comin_ftnlist_is_end(list,it))
        CALL comin_ftnlist_iterator_value(it, cptr)
        CALL c_f_POINTER(cptr, item)
        CALL comin_metadata_get_or(item%metadata, "tracer", tracer, .FALSE.)

        IF (tracer) THEN
          no_tracers = no_tracers+1
        END IF
        no_vars = no_vars+1
        CALL comin_ftnlist_iterator_next(it)
      END DO
      CALL comin_ftnlist_iterator_delete(it)

      ALLOCATE(vars(no_vars + global%n_dom))  ! plus 1 for the "tracer" container variable (on each domain)

      no_vars = 1
      ! add the "tracer" variables
      DO jg=1,global%n_dom
        ! We only support tracers in double precision
        ALLOCATE(vars(jg)%data_ptr_dp(global%nproma, domain(jg)%nlev, &
             & domain(jg)%cells%nblks, 1, no_tracers))
        descr = t_comin_var_descriptor("container",jg)
        CALL comin_var_list_append(descr, &
             & C_LOC(vars(no_vars)%data_ptr_dp), &
             & device_ptr=C_NULL_PTR, &
             & array_shape = SHAPE(vars(no_vars)%data_ptr_dp), &
             & type_id = COMIN_VAR_DATATYPE_DOUBLE, &
             & dim_semantics = [COMIN_DIM_SEMANTICS_NPROMA, &
             & COMIN_DIM_SEMANTICS_LEVEL,  &
             & COMIN_DIM_SEMANTICS_BLOCK,  &
             & COMIN_DIM_SEMANTICS_UNUSED, &
             & COMIN_DIM_SEMANTICS_CONTAINER ], &
             & lcontainer = .TRUE.,  &
             & ncontained=no_tracers, &
             & var_handle=dummy)
        no_vars = no_vars + 1
      END DO

      tracer_counter = 0

      CALL comin_ftnlist_iterator_begin(list, it)
      do while (.not. comin_ftnlist_is_end(list,it))
        CALL comin_ftnlist_iterator_value(it, cptr)
        CALL c_f_POINTER(cptr, item)
        dom_idx = item%descriptor%id
        CALL comin_metadata_get_or(item%metadata, "zaxis_id", zaxis_id, COMIN_ZAXIS_3D)
        CALL comin_metadata_get_or(item%metadata, "hgrid_id", hgrid_id, COMIN_HGRID_UNSTRUCTURED_CELL)
        SELECT CASE(hgrid_id)
        CASE (COMIN_HGRID_UNSTRUCTURED_CELL)
          nblks = domain(dom_idx)%cells%nblks
        CASE (COMIN_HGRID_UNSTRUCTURED_EDGE)
          nblks = domain(dom_idx)%edges%nblks
        CASE (COMIN_HGRID_UNSTRUCTURED_VERTEX)
          nblks = domain(dom_idx)%verts%nblks
        CASE DEFAULT
          CALL finish("comin_replay", "Unknown HGRID")
        END SELECT
        IF (zaxis_id == COMIN_ZAXIS_3D) THEN
          dimshape=[global%nproma, domain(dom_idx)%nlev, nblks, 1, 1]
          dim_semantics = [COMIN_DIM_SEMANTICS_NPROMA, &
               &           COMIN_DIM_SEMANTICS_LEVEL, &
               &           COMIN_DIM_SEMANTICS_BLOCK, &
               &           COMIN_DIM_SEMANTICS_UNUSED, &
               &           COMIN_DIM_SEMANTICS_UNUSED ]
        ELSE IF(zaxis_id == COMIN_ZAXIS_3D_HALF) THEN
          dimshape=[global%nproma, domain(dom_idx)%nlev + 1, nblks, 1, 1]
          dim_semantics = [COMIN_DIM_SEMANTICS_NPROMA, &
               &           COMIN_DIM_SEMANTICS_LEVEL, &
               &           COMIN_DIM_SEMANTICS_BLOCK, &
               &           COMIN_DIM_SEMANTICS_UNUSED, &
               &           COMIN_DIM_SEMANTICS_UNUSED ]
        ELSE IF(zaxis_id == COMIN_ZAXIS_2D) THEN
          dimshape=[global%nproma, nblks, 1, 1, 1]
          dim_semantics = [COMIN_DIM_SEMANTICS_NPROMA, &
               &           COMIN_DIM_SEMANTICS_BLOCK, &
               &           COMIN_DIM_SEMANTICS_UNUSED, &
               &           COMIN_DIM_SEMANTICS_UNUSED, &
               &           COMIN_DIM_SEMANTICS_UNUSED ]
        ELSE
          CALL finish("comin_replay", "Unknown ZAXIS")
        END IF

        CALL comin_metadata_get_or(item%metadata, "tracer", tracer, .FALSE.)
        CALL comin_metadata_get_or(item%metadata, "datatype", datatype, COMIN_VAR_DATATYPE_DOUBLE)

        ! ALLOCATE(vars(no_vars)%comin_var_ptr)
        IF (tracer) THEN
          IF( datatype /= COMIN_VAR_DATATYPE_DOUBLE) &
            CALL finish("comin_replay", "tracers are only supported in double precision")
          tracer_counter(dom_idx) = tracer_counter(dom_idx) + 1
          ncontained = tracer_counter(dom_idx)
          vars(no_vars)%data_ptr_dp => vars(dom_idx)%data_ptr_dp(:,:,:,:,ncontained:ncontained)
          cptr = C_LOC(vars(no_vars)%data_ptr_dp)
        ELSE
          SELECT CASE(datatype)
          CASE (COMIN_VAR_DATATYPE_DOUBLE)
            ALLOCATE(vars(no_vars)%data_ptr_dp(dimshape(1), dimshape(2), &
                 & dimshape(3), dimshape(4), dimshape(5)))
            cptr = C_LOC(vars(no_vars)%data_ptr_dp)
          CASE (COMIN_VAR_DATATYPE_FLOAT)
            ALLOCATE(vars(no_vars)%data_ptr_sp(dimshape(1), dimshape(2), &
                 & dimshape(3), dimshape(4), dimshape(5)))
            cptr = C_LOC(vars(no_vars)%data_ptr_sp)
          CASE (COMIN_VAR_DATATYPE_INT)
            ALLOCATE(vars(no_vars)%data_ptr_i(dimshape(1), dimshape(2), &
                 & dimshape(3), dimshape(4), dimshape(5)))
            cptr = C_LOC(vars(no_vars)%data_ptr_i)
          CASE DEFAULT
            CALL finish("comin_replay", "Unknown datatype")
          END SELECT
          ncontained = 0
        END IF

        CALL comin_var_list_append(item%descriptor, &
             & cptr, &
             & device_ptr=C_NULL_PTR,  &
             & array_shape = dimshape, &
             & type_id = datatype,     &
             & dim_semantics = dim_semantics, &
             & lcontainer = .FALSE.,  &
             & ncontained= ncontained, &
             & var_handle=dummy &
             )

        ! set metadata defaults (may be overridden below)
        ! see Metadat section of icon_comin_doc.md for default values
        CALL comin_metadata_set(item%descriptor, "zaxis_id", zaxis_id)
        CALL comin_metadata_set(item%descriptor, "hgrid_id", hgrid_id)
        CALL comin_metadata_set(item%descriptor, "restart", .FALSE.)
        CALL comin_metadata_set(item%descriptor, "multi_timelevel", .FALSE.)
        CALL comin_metadata_set(item%descriptor, "tracer", tracer)
        CALL comin_metadata_set(item%descriptor, "tracer_turb", .FALSE.)
        CALL comin_metadata_set(item%descriptor, "tracer_conv", .FALSE.)
        CALL comin_metadata_set(item%descriptor, "tracer_hlimit", 4)
        CALL comin_metadata_set(item%descriptor, "tracer_vlimit", 1)
        CALL comin_metadata_set(item%descriptor, "tracer_hadv", 2)
        CALL comin_metadata_set(item%descriptor, "tracer_vadv", 3)
        CALL comin_metadata_set(item%descriptor, "units", "")
        CALL comin_metadata_set(item%descriptor, "standard_name", "")
        CALL comin_metadata_set(item%descriptor, "long_name", "")
        CALL comin_metadata_set(item%descriptor, "short_name", "")
        CALL comin_metadata_set(item%descriptor, "datatype", datatype)

        ! Iterate through request list metadata, forward metadata to var list
        CALL item%metadata%get_iterator(metadata_it)
        DO WHILE(.NOT. metadata_it%is_end())
          current_metadata_key = metadata_it%key()
          SELECT CASE(item%metadata%query(current_metadata_key))
          CASE (COMIN_METADATA_TYPEID_INTEGER)
            CALL item%metadata%get(TRIM(current_metadata_key), intbuff)
            CALL comin_metadata_set(item%descriptor, current_metadata_key, intbuff)
          CASE (COMIN_METADATA_TYPEID_REAL)
            CALL item%metadata%get(current_metadata_key, realbuff)
            CALL comin_metadata_set(item%descriptor, current_metadata_key, realbuff)
          CASE (COMIN_METADATA_TYPEID_CHARACTER)
            CALL item%metadata%get(current_metadata_key, charbuff)
            CALL comin_metadata_set(item%descriptor, current_metadata_key, charbuff)
          CASE (COMIN_METADATA_TYPEID_LOGICAL)
            CALL item%metadata%get(current_metadata_key, logbuff)
            CALL comin_metadata_set(item%descriptor, current_metadata_key, logbuff)
          END SELECT
          CALL metadata_it%next()
        ENDDO
        CALL metadata_it%delete()

        no_vars = no_vars+1
        CALL comin_ftnlist_iterator_next(it)
      END DO
      CALL comin_ftnlist_iterator_delete(it)
    END SUBROUTINE add_variables

    FUNCTION trim_null( string)
      CHARACTER(len=*) :: string
      CHARACTER(len=LEN(string)) :: trim_null

      INTEGER :: pos

      pos = INDEX( string, ACHAR(0) )
      IF ( pos > 0 ) THEN
        trim_null = string(1:pos-1)
      ELSE
        trim_null = string
      ENDIF
    END FUNCTION trim_null

    SUBROUTINE yac_routines(current_ep)
      INTEGER, INTENT(IN) :: current_ep
#ifdef ENABLE_YAC
      IF (current_ep == EP_ATM_YAC_DEFCOMP_AFTER) THEN
        CALL yac_fdef_comp(yac_instance_id, "comin_replay", yac_comp_id)
      ENDIF
      IF (current_ep == EP_ATM_YAC_SYNCDEF_AFTER) THEN
        CALL yac_fsync_def(yac_instance_id)
      ENDIF
      IF (current_ep == EP_ATM_YAC_ENDDEF_AFTER) THEN
        CALL yac_fenddef(yac_instance_id)
      ENDIF
#endif
    END SUBROUTINE yac_routines

#ifdef ENABLE_YAXT
    SUBROUTINE set_yaxt_redist()
      USE yaxt, ONLY: xt_idxlist, xt_xmap, xt_redist, &
           &     xt_initialized, xt_initialize,       &
           &     xt_xmap_all2all_new, xt_idxvec_new,  &
           &     xt_redist_p2p_new, xt_xmap_delete,   &
           &     xt_idxlist_delete
      INTEGER(kind=xi), DIMENSION(:), ALLOCATABLE :: idxmap
      TYPE(xt_idxlist)                            :: src_idxlist, tgt_idxlist
      TYPE(xt_xmap)                               :: xmap
      INTEGER                                     :: k, jg
      IF (.NOT. xt_initialized()) THEN
        !initialize yaxt
        CALL xt_initialize(host_comm)
      END IF
      !construct yaxt variables
      DO jg = 1, comin_global%n_dom
        idxmap = RESHAPE(comin_domain(jg)%cells%decomp_domain &
                         , (/ SIZE(comin_domain(jg)%cells%decomp_domain) /))
        idxvec_domain(jg)%idxvec = INT(PACK( [(k,k=1,comin_domain(jg)%cells%ncells)], idxmap == 0 ), xi)
        src_idxlist = xt_idxvec_new(comin_domain(jg)%cells%glb_index(idxvec_domain(jg)%idxvec))
        tgt_idxlist = xt_idxvec_new(comin_domain(jg)%cells%glb_index)
        !create exchange map
        xmap = xt_xmap_all2all_new(src_idxlist, tgt_idxlist, host_comm)
        !create redistribution instance for DP
        yaxt_redist_dp(jg) = xt_redist_p2p_new(xmap, MPI_DOUBLE_PRECISION)
        yaxt_redist_sp(jg) = xt_redist_p2p_new(xmap, MPI_FLOAT)
        yaxt_redist_i(jg) = xt_redist_p2p_new(xmap, MPI_INTEGER)
        !clean up
        CALL xt_xmap_delete(xmap)
        CALL xt_idxlist_delete(src_idxlist)
        CALL xt_idxlist_delete(tgt_idxlist)
      END DO
      DEALLOCATE(idxmap)

    END SUBROUTINE set_yaxt_redist

    SUBROUTINE replay_halo_sync_variable_dp(var_ptr, halo_sync_mode)
      USE yaxt, ONLY: xt_redist_s_exchange1
      TYPE(t_comin_var_handle), INTENT(IN) :: var_ptr
      INTEGER, INTENT(IN) :: halo_sync_mode
      REAL(C_DOUBLE), DIMENSION(:,:),   POINTER   :: tgt
      REAL(C_DOUBLE), DIMENSION(:),   POINTER     :: buffer1d
      REAL(C_DOUBLE), DIMENSION(:,:,:), POINTER   :: src3d
      INTEGER                                    :: jk,jg
      TYPE(t_comin_var_descriptor) :: descr
      descr = var_ptr%descriptor()
      jg = descr%id
      ALLOCATE(tgt(comin_global%nproma,comin_domain(jg)%cells%nblks))
      ALLOCATE(buffer1d(comin_global%nproma*comin_domain(jg)%cells%nblks))
      CALL var_ptr%to_3d(src3d)
      !do the halo exchange for var pointer
      IF(halo_sync_mode==COMIN_ZAXIS_2D) THEN
        buffer1d(:) = RESHAPE(src3d, [comin_global%nproma*comin_domain(jg)%cells%nblks])
        buffer1d(:SIZE(idxvec_domain(jg)%idxvec)) = buffer1d(idxvec_domain(jg)%idxvec)
        CALL xt_redist_s_exchange1(yaxt_redist_dp(jg), C_LOC(buffer1d), C_LOC(tgt))
        src3d(:,:,1) = tgt(:,:)
      ELSE IF (halo_sync_mode==COMIN_ZAXIS_3D) THEN
        DO jk = 1, comin_domain(jg)%nlev
          buffer1d(:) = RESHAPE(src3d(:,jk,:), [comin_global%nproma*comin_domain(jg)%cells%nblks])
          buffer1d(:SIZE(idxvec_domain(jg)%idxvec)) = buffer1d(idxvec_domain(jg)%idxvec)
          CALL xt_redist_s_exchange1(yaxt_redist_dp(jg), C_LOC(buffer1d), C_LOC(tgt))
          src3d(:,jk,:) = tgt(:,:)
        END DO
      END IF
      DEALLOCATE(buffer1d)

    END SUBROUTINE replay_halo_sync_variable_dp

    SUBROUTINE replay_halo_sync_variable_sp(var_ptr, halo_sync_mode)
      USE yaxt, ONLY: xt_redist_s_exchange1
      TYPE(t_comin_var_handle), INTENT(IN) :: var_ptr
      INTEGER, INTENT(IN) :: halo_sync_mode
      REAL(C_FLOAT), DIMENSION(:,:),   POINTER   :: tgt
      REAL(C_FLOAT), DIMENSION(:),   POINTER     :: buffer1d
      REAL(C_FLOAT), DIMENSION(:,:,:), POINTER   :: src3d
      INTEGER                                    :: jk,jg
      TYPE(t_comin_var_descriptor) :: descr
      descr = var_ptr%descriptor()
      jg = descr%id
      ALLOCATE(tgt(comin_global%nproma,comin_domain(jg)%cells%nblks))
      ALLOCATE(buffer1d(comin_global%nproma*comin_domain(jg)%cells%nblks))
      CALL var_ptr%to_3d(src3d)
      !do the halo exchange for var pointer
      IF(halo_sync_mode==COMIN_ZAXIS_2D) THEN
        buffer1d(:) = RESHAPE(src3d, [comin_global%nproma*comin_domain(jg)%cells%nblks])
        buffer1d(:SIZE(idxvec_domain(jg)%idxvec)) = buffer1d(idxvec_domain(jg)%idxvec)
        CALL xt_redist_s_exchange1(yaxt_redist_sp(jg), C_LOC(buffer1d), C_LOC(tgt))
        src3d(:,:,1) = tgt(:,:)
      ELSE IF (halo_sync_mode==COMIN_ZAXIS_3D) THEN
        DO jk = 1, comin_domain(jg)%nlev
          buffer1d(:) = RESHAPE(src3d(:,jk,:), [comin_global%nproma*comin_domain(jg)%cells%nblks])
          buffer1d(:SIZE(idxvec_domain(jg)%idxvec)) = buffer1d(idxvec_domain(jg)%idxvec)
          CALL xt_redist_s_exchange1(yaxt_redist_sp(jg), C_LOC(buffer1d), C_LOC(tgt))
          src3d(:,jk,:) = tgt(:,:)
        END DO
      END IF
      DEALLOCATE(buffer1d)

    END SUBROUTINE replay_halo_sync_variable_sp

    SUBROUTINE replay_halo_sync_variable_i(var_ptr, halo_sync_mode)
      USE yaxt, ONLY: xt_redist_s_exchange1
      TYPE(t_comin_var_handle), INTENT(IN) :: var_ptr
      INTEGER, INTENT(IN) :: halo_sync_mode
      INTEGER(C_INT), DIMENSION(:,:),   POINTER   :: tgt
      INTEGER(C_INT), DIMENSION(:),   POINTER     :: buffer1d
      INTEGER(C_INT), DIMENSION(:,:,:), POINTER   :: src3d
      INTEGER                                    :: jk,jg
      TYPE(t_comin_var_descriptor) :: descr
      descr = var_ptr%descriptor()
      jg = descr%id
      ALLOCATE(tgt(comin_global%nproma,comin_domain(jg)%cells%nblks))
      ALLOCATE(buffer1d(comin_global%nproma*comin_domain(jg)%cells%nblks))
      CALL var_ptr%to_3d(src3d)
      !do the halo exchange for var pointer
      IF(halo_sync_mode==COMIN_ZAXIS_2D) THEN
        buffer1d(:) = RESHAPE(src3d, [comin_global%nproma*comin_domain(jg)%cells%nblks])
        buffer1d(:SIZE(idxvec_domain(jg)%idxvec)) = buffer1d(idxvec_domain(jg)%idxvec)
        CALL xt_redist_s_exchange1(yaxt_redist_i(jg), C_LOC(buffer1d), C_LOC(tgt))
        src3d(:,:,1) = tgt(:,:)
      ELSE IF (halo_sync_mode==COMIN_ZAXIS_3D) THEN
        DO jk = 1, comin_domain(jg)%nlev
          buffer1d(:) = RESHAPE(src3d(:,jk,:), [comin_global%nproma*comin_domain(jg)%cells%nblks])
          buffer1d(:SIZE(idxvec_domain(jg)%idxvec)) = buffer1d(idxvec_domain(jg)%idxvec)
          CALL xt_redist_s_exchange1(yaxt_redist_i(jg), C_LOC(buffer1d), C_LOC(tgt))
          src3d(:,jk,:) = tgt(:,:)
        END DO
      END IF
      DEALLOCATE(buffer1d)

    END SUBROUTINE replay_halo_sync_variable_i

    SUBROUTINE replay_halo_sync_variable(var_ptr, halo_sync_mode)
      TYPE(t_comin_var_handle), INTENT(IN) :: var_ptr
      INTEGER, INTENT(IN) :: halo_sync_mode
      INTEGER :: datatype
      CALL comin_metadata_get(var_ptr%descriptor(), "datatype", datatype)
      SELECT CASE(datatype)
      CASE (COMIN_VAR_DATATYPE_DOUBLE)
        CALL replay_halo_sync_variable_dp(var_ptr, halo_sync_mode)
      CASE (COMIN_VAR_DATATYPE_FLOAT)
        CALL replay_halo_sync_variable_sp(var_ptr, halo_sync_mode)
      CASE (COMIN_VAR_DATATYPE_INT)
        CALL replay_halo_sync_variable_i(var_ptr, halo_sync_mode)
      CASE DEFAULT
        CALL finish("comin_replay", "Unknown datatype!")
      END SELECT
    END SUBROUTINE replay_halo_sync_variable
#endif

    INTEGER FUNCTION glb2loc_lookup_wrapper_cell(jg, glb) RESULT(loc)
      INTEGER, INTENT(IN) :: jg
      INTEGER, INTENT(IN) :: glb
      loc = glb2loc_lookup(glb2loc_lookup_cell(jg), glb)
    END FUNCTION glb2loc_lookup_wrapper_cell

  END SUBROUTINE main

END MODULE mo_comin_replay

PROGRAM comin_replay
  USE mo_comin_replay, ONLY: main

  CALL main
END PROGRAM comin_replay
