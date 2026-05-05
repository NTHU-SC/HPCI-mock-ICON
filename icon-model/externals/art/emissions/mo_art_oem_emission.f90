!
! Source Module for OEM: Online Anthropogenic Emissions
!--------------------------------------------------------------------------------------
!
! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2025, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE mo_art_oem_emission

  !------------------------------------------------------------------------------------
  ! Description:
  !   This module contains subroutines for computatoin of the emissions by 
  !   scaling the gridded emissions with temporal and vertical profiles and
  !   adds them to the tracers.
  
  !! Modifications: 
  !! 2024: Arash Hamzehloo, Empa
  !! - OEM was substantially refactored & ported to GPUs.
  !!
  !====================================================================================
    ! ICON 
    USE mo_parallel_config,        ONLY: nproma
    USE mo_kind,                   ONLY: wp, i8
    USE mo_model_domain,           ONLY: t_patch
    USE mo_nonhydro_state,         ONLY: p_nh_state
  
    USE mo_time_config,            ONLY: time_config
    USE mtime,                     ONLY: MAX_DATETIME_STR_LEN,      &
                                     &   datetime,                  &
                                     &   timedelta,                 &
                                     &   newTimedelta,              &
                                     &   getdayofyearfromdatetime,  &
                                     &   OPERATOR(-), OPERATOR(+),  &
                                     &   OPERATOR(>=)

    USE mo_util_mtime,             ONLY: getElapsedSimTimeInSeconds
    USE mo_loopindices,            ONLY: get_indices_c
    USE mo_impl_constants_grf,     ONLY: grf_bdywidth_c
    USE mo_limarea_config,         ONLY: latbc_config
    USE mo_async_latbc_types,      ONLY: t_latbc_data
    USE mo_exception,              ONLY: message
    USE mo_util_string,            ONLY: tolower

    USE mo_nwp_phy_state,          ONLY: prm_diag

  
    ! ART
    USE mo_art_atmo_data,          ONLY: t_art_atmo
    USE mo_art_data,               ONLY: p_art_data
    USE mo_art_wrapper_routines,   ONLY: art_get_indices_c
  
    ! OEM
    USE mo_art_oem_types,          ONLY: p_art_oem_data
    USE mo_oem_config,             ONLY: restart_init_time

    ! VPRM
    USE mo_oem_config,             ONLY: vprm_par,          & 
                                     &   vprm_lambda,       & 
                                     &   vprm_alpha,        & 
                                     &   vprm_beta,         & 
                                     &   vprm_tmin,         & 
                                     &   vprm_tmax,         & 
                                     &   vprm_topt,         & 
                                     &   vprm_tlow,         & 
                                     &   lcut_area,         & 
                                     &   lon_cut_start,     &
                                     &   lon_cut_end,       &
                                     &   lat_cut_start,     & 
                                     &   lat_cut_end
  
  !---------------------------------------------------------------------------------
  
    IMPLICIT NONE
  
    PRIVATE
    PUBLIC ::                       &
      &  art_oem_compute_emissions, &
      &  art_oem_extract_time_information
  
    ! Constant variable
    INTEGER,  PARAMETER :: tp_param_hourofday = 24
    INTEGER,  PARAMETER :: tp_param_dayofweek = 7
    INTEGER,  PARAMETER :: tp_param_monthofyear = 12
    INTEGER(KIND=2), PARAMETER :: tp_param_hour = 8784
  
  
    CHARACTER(LEN=3), DIMENSION(7) :: day_of_week = &
      & (/ 'MON', 'TUE', 'WED', 'THU', 'FRI', 'SAT', 'SUN' /)
  
    INTEGER :: global_iteration = 0, hour_iteration = 0, hour_sec = 0, global_sec = 0

    INTEGER, ALLOCATABLE :: hod(:), dow(:), moy(:), hoy(:), buff(:), is_m(:), ie_m(:), is_b_m(:), ie_b_m(:)

    INTEGER, ALLOCATABLE :: grd_index(:,:), vp_cat_idx(:,:), tp_cat_idx(:,:), tp_country_idx(:,:)
 
    REAL(KIND=wp), ALLOCATABLE :: tsf_now(:,:,:,:), tsf_next(:,:,:,:), total_weight(:,:,:,:,:)

    TYPE(datetime) :: next_latbc_time, next_full_hour 


    INTEGER :: ncat_max = 200

    CHARACTER(LEN=16), ALLOCATABLE :: buffer_name(:), ensemble_name(:)

  !====================================================================================
  ! Module procedures
  !====================================================================================
   
  CONTAINS
  
    SUBROUTINE art_oem_compute_emissions(p_tracer_now,p_patch,dtime,mtime_current,ierror,yerrmsg)
  
  !-----------------------------------------------------------------------------------
  ! Description: This subroutine uses the temporal and vertical profiles to compute
  ! the gridded emissions and add them to the OEM-tracers.
  !-----------------------------------------------------------------------------------
  
    IMPLICIT NONE
  
    REAL(wp), INTENT(INOUT) :: &
     &  p_tracer_now(:,:,:,:)     !< tracer mixing ratio
  
    TYPE(t_patch), INTENT(IN) :: &
     &  p_patch                  !< patch on which computation is performed
  
    REAL(wp), INTENT(IN) :: &
     &  dtime                     !< time step    
  
    TYPE(datetime), POINTER ::  &
     &  mtime_current             !< current datetime
     
    INTEGER, INTENT(INOUT)             :: ierror
    CHARACTER(LEN= *), INTENT(INOUT)   :: yerrmsg
  
    !---------------------------------------------------------------------------------
    ! Local variables
    INTEGER ::                 &
      &  nc, nt, is, ie, jg,   &
      &  nlev, k, trcr_idx, jb, jc,          &
      &  nblks_c, i_startblk, i_endblk, idx, &
      &  ii, jj, kk, ll, &
      &  min, &
      &  nr, nens, ens_int_idx, table_nr, ens_count, nt_old, &
      &  i_startblk_bdry, i_endblk_bdry, name_idx, nc_ensemble
  
    REAL(KIND=wp) ::                         &
      &  n1, n2, &
      &  temp_scaling_fact, lambda, temp_tracer
  
    CHARACTER(LEN=3) :: numstring
       
    TYPE(t_art_atmo), POINTER :: &
      &  art_atmo
  
    TYPE(datetime) :: datetime_next
    TYPE(timedelta), POINTER :: mtime_td

    TYPE(t_latbc_data), TARGET :: latbc

    INTEGER, PARAMETER :: n_vprm_cats = 16 ! number of VPRM categories

    CHARACTER(*), PARAMETER :: routine = "art_oem_compute_emissions"

  ! End of header
  !==================================================================================== 
    ierror = 0
    yerrmsg = '   '
    jg   = p_patch%id
    nblks_c = p_art_data(jg)%atmo%nblks
    nlev = p_art_data(jg)%atmo%nlev

  !------------------------------------------------------------------------------------
  
    IF (p_art_oem_data%configure%emis_tracer>0) THEN

      ! Read start- and end-block for this PE:
      i_startblk = p_art_data(jg)%atmo%i_startblk
      i_endblk   = p_art_data(jg)%atmo%i_endblk

      i_startblk_bdry = p_patch%cells%start_blk(1,1)
      i_endblk_bdry   = p_patch%cells%end_blk(grf_bdywidth_c,1)

      global_iteration = global_iteration + 1
      global_sec = global_iteration*dtime


      IF (global_iteration==1) THEN

        next_full_hour = mtime_current
        
        ALLOCATE(hod(2), dow(2), moy(2), hoy(2))

      ENDIF 

      IF (mtime_current>=next_full_hour) THEN

        ! Extract the different time information for this timestep
        CALL art_oem_extract_time_information(time_config%tc_current_date,hod(1),dow(1),moy(1),hoy(1))
        ! Extract the time information for one hour later
        mtime_td => newTimedelta("PT01H")
        datetime_next = time_config%tc_current_date + mtime_td
        CALL art_oem_extract_time_information(datetime_next,hod(2),dow(2),moy(2),hoy(2))

        hour_iteration = 0

      END IF

      hour_iteration = hour_iteration + 1
      hour_sec = hour_iteration*dtime


!------------------------------------------------------------------------------------
! Section 1:
!     
! The following section is executed only once at the beginning of the 
! simulation. Here, indices of the grid emission, vertical profile, temporal 
! profile, and countries are updated. 
!
! After that, the total "weight" due to the "gridded" emission and vertical 
! scaling are evaluated. 
!------------------------------------------------------------------------------------

      IF (global_iteration == 1) THEN

        next_latbc_time = mtime_current
  
        ALLOCATE(grd_index(p_art_oem_data%configure%emis_tracer, ncat_max))
        ALLOCATE(vp_cat_idx(p_art_oem_data%configure%emis_tracer, ncat_max))
        ALLOCATE(tp_cat_idx(p_art_oem_data%configure%emis_tracer, ncat_max))
        ALLOCATE(buffer_name(SIZE(p_art_oem_data%data_fields%p_latbc_data%buffer%name_tracer(:))))
        ALLOCATE(ensemble_name(SIZE(p_art_oem_data%ensemble%ens_name(:))))
        ALLOCATE(buff(SIZE(p_art_oem_data%data_fields%p_latbc_data%buffer%idx_tracer(:))))


        DO nt = 1, p_art_oem_data%configure%emis_tracer 
          DO nc = 1, ncat_max

            IF (p_art_oem_data%configure%ycatl_l(nt,nc) /= "") THEN

              grd_index(nt,nc) = 0

              DO ii = 1, SIZE(p_art_oem_data%configure%gridded_emissions_idx)
                IF(p_art_oem_data%configure%ycatl_l(nt,nc) == p_art_oem_data%configure%gridded_emissions_idx(ii)) THEN
                  grd_index(nt,nc) = ii
                  EXIT
                END IF
              END DO

              vp_cat_idx(nt,nc) = 0

              DO jj = 1, SIZE(p_art_oem_data%configure%vp_category)
                IF(p_art_oem_data%configure%yvpl_l(nt,nc) == p_art_oem_data%configure%vp_category(jj)) THEN
                  vp_cat_idx(nt,nc) = jj
                  EXIT
                END IF
              END DO

              tp_cat_idx(nt,nc) = 0

              DO kk = 1, SIZE(p_art_oem_data%configure%tp_category)
                IF(p_art_oem_data%configure%ytpl_l(nt,nc) == p_art_oem_data%configure%tp_category(kk)) THEN
                  tp_cat_idx(nt,nc) = kk
                  EXIT
                END IF
              END DO
            END IF
          END DO
        END DO

        DO jb = i_startblk, i_endblk

          CALL art_get_indices_c(jg, jb, is, ie)

          IF (jb == i_startblk) THEN

            ALLOCATE(tp_country_idx(nproma, i_endblk))
            ALLOCATE(total_weight(nproma, nlev, p_art_oem_data%configure%emis_tracer,ncat_max, i_endblk))
            ALLOCATE(tsf_now(nproma, p_art_oem_data%configure%emis_tracer,ncat_max, i_endblk))
            ALLOCATE(tsf_next(nproma, p_art_oem_data%configure%emis_tracer,ncat_max, i_endblk))  
            ALLOCATE(is_m(i_endblk))
            ALLOCATE(ie_m(i_endblk))
            ALLOCATE(is_b_m(i_endblk_bdry))
            ALLOCATE(ie_b_m(i_endblk_bdry))

          ENDIF

          is_m(jb) = is
          ie_m(jb) = ie

          DO nt = 1, p_art_oem_data%configure%emis_tracer 
            DO nc = 1, ncat_max

              IF (p_art_oem_data%configure%ycatl_l(nt,nc) /= "") THEN

                DO k=1, nlev 
                  DO jc = is, ie

                    total_weight(jc,k,nt,nc,jb) = dtime * p_art_oem_data%data_fields%gridded_emissions(jc,jb,grd_index(nt,nc)) &
                    & * p_art_oem_data%data_fields%vert_scaling_fact(jc,k,jb,vp_cat_idx(nt,nc))

                    tp_country_idx(jc,jb) = 0
                    DO ll = 1, p_art_oem_data%configure%tp_ncountry
                      IF(p_art_oem_data%data_fields%country_ids(jc,jb) == p_art_oem_data%data_fields%tp_countryid(ll)) THEN
                        tp_country_idx(jc,jb) = ll
                      END IF
                    END DO
                    
                  END DO

                END DO

              END IF

            END DO
            p_art_oem_data%configure%emis_name(nt) = tolower(p_art_oem_data%configure%emis_name(nt))
          END DO

        END DO

        DO nens = 1,p_art_oem_data%configure%ens_tracer

          p_art_oem_data%ensemble%bg_ens_name(nens) = tolower(p_art_oem_data%ensemble%bg_ens_name(nens))
          p_art_oem_data%ensemble%ens_name(nens) = tolower(p_art_oem_data%ensemble%ens_name(nens))

          write(numstring,'(I3.3)') nens
          ensemble_name(nens) = TRIM(p_art_oem_data%ensemble%ens_name(nens)) // '-' // numstring

        ENDDO

        DO idx=1,SIZE(p_art_oem_data%data_fields%p_latbc_data%buffer%name_tracer(:))

          buffer_name(idx) = TRIM(p_art_oem_data%data_fields%p_latbc_data%buffer%name_tracer(idx))

        ENDDO 


        DO jb = i_startblk_bdry, i_endblk_bdry

        
          IF (jb == i_startblk_bdry) THEN
            is = MAX(1,p_patch%cells%start_index(1))
            ie   = nproma
            IF (jb == i_endblk_bdry) ie = p_patch%cells%end_index(grf_bdywidth_c)
          ELSE IF (jb == i_endblk_bdry) THEN
            is = 1
            ie   = p_patch%cells%end_index(grf_bdywidth_c)
          ELSE
            is = 1
            ie = nproma
          ENDIF


          is_b_m(jb) = is
          ie_b_m(jb) = ie

        ENDDO

        buff(:) = p_art_oem_data%data_fields%p_latbc_data%buffer%idx_tracer(:)

        !$ACC ENTER DATA COPYIN(p_nh_state, p_art_oem_data, p_patch, p_art_data(jg), &  
        !$ACC           & p_art_oem_data%configure%ycatl_l, &
        !$ACC           & p_art_oem_data%configure%emis_tracer, &
        !$ACC           & p_art_oem_data%ensemble%ens_name, p_art_oem_data%configure%emis_name, p_art_oem_data%data_fields%reg_map, &
        !$ACC           & p_art_oem_data%ensemble%ens_table, p_art_oem_data%data_fields%lambda_mat, &
        !$ACC           & total_weight, ncat_max, restart_init_time, &
        !$ACC           & p_art_oem_data%configure%restart_tracer, p_art_oem_data%configure%l_restarted, &
        !$ACC           & p_art_oem_data%data_fields%chem_init_3D, p_art_oem_data%configure%restart_idx, &
        !$ACC           & p_art_oem_data%data_fields%boundary_lambdas, &
        !$ACC           & p_art_oem_data%ensemble%bg_ens_name, &
        !$ACC           & p_art_oem_data%data_fields%lambda_categories_ids, &
        !$ACC           & p_art_oem_data%configure%ens_tracer, p_art_oem_data%configure%emis_diag_out, &
        !$ACC           & p_art_oem_data%data_fields%tp_hourofday, p_art_oem_data%data_fields%tp_dayofweek, &
        !$ACC           & p_art_oem_data%data_fields%tp_monthofyear, p_art_oem_data%data_fields%tp_hourofyear, &
        !$ACC           & p_art_oem_data%configure%itype_tscale_l, tp_cat_idx, tp_country_idx, tsf_now, tsf_next, p_tracer_now, is_m, ie_m)


        ! VPRM-related data
#if defined( _OPENACC )
        IF (p_art_oem_data%configure%vprm_tracer>0) THEN
          !$ACC ENTER DATA COPYIN(prm_diag, p_art_oem_data%configure%vprm_tracer, p_art_oem_data%configure%vprm_idx,  &
          !$ACC           & p_art_oem_data%data_fields%newflux_vprm , p_art_data(jg)%atmo%t_2m, p_art_oem_data%configure%vprm_flux_type, &
          !$ACC           & p_art_data(jg)%atmo%lon, p_art_data(jg)%atmo%lat, p_art_oem_data%data_fields%vprm_lu_class_fraction, &
          !$ACC           & p_art_oem_data%configure%vprm_name, p_art_oem_data%ensemble%vprm_bg_ens, &
          !$ACC           & p_art_oem_data%data_fields%i_vprm_lc_shrub, p_art_oem_data%data_fields%i_vprm_lc_grass, &
          !$ACC           & p_art_oem_data%data_fields%lswi_max, p_art_oem_data%data_fields%lswi_min, p_art_oem_data%data_fields%lswi, & 
          !$ACC           & p_art_oem_data%data_fields%i_vprm_lc_evergreen, p_art_oem_data%data_fields%i_vprm_lc_savanna, &
          !$ACC           & p_art_oem_data%data_fields%evi_min, p_art_oem_data%data_fields%evi_max, p_art_oem_data%data_fields%evi, &
          !$ACC           & p_art_oem_data%ensemble%vprm_bg_ens_name, p_art_data(jg)%atmo%nlev, &
          !$ACC           & vprm_par, vprm_lambda, vprm_alpha, vprm_beta, vprm_tmin, vprm_tmax, vprm_topt, vprm_tlow, lcut_area, &
          !$ACC           & lon_cut_start, lon_cut_end, lat_cut_start, lat_cut_end)
        ENDIF
#endif

        CALL message(routine, "OEM's section 1 executed")
      ENDIF ! global_iteration == 1

!------------------------------------------------------------------------------------
! End of Section 1
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
! Section 2:
!     
! The following section is executed at the update time of the time-dependant BCs. 
! Here, tracer tendencies at the lateral boundaries are updated.
!------------------------------------------------------------------------------------

        IF (mtime_current>=next_latbc_time) THEN
          !$ACC UPDATE HOST(p_nh_state%diag%grf_tend_tracer) 

          DO jb = i_startblk_bdry, i_endblk_bdry

            is = is_b_m(jb)
            ie = ie_b_m(jb) 
  
            DO idx=1,SIZE(p_art_oem_data%data_fields%p_latbc_data%buffer%name_tracer(:))
  
              IF ( ANY( p_art_oem_data%ensemble%bg_ens_name==buffer_name(idx) ) ) THEN

                DO nens = 1,p_art_oem_data%configure%ens_tracer
  
                  IF (p_art_oem_data%ensemble%bg_ens_name(nens) == buffer_name(idx)) THEN
  

                    DO name_idx=1,SIZE(p_art_oem_data%data_fields%p_latbc_data%buffer%name_tracer(:))
                      IF ( ensemble_name(nens)  == buffer_name(name_idx) ) EXIT
                    END DO
  

                    DO k = 1,nlev
                      DO jc = is, ie
    
                        p_nh_state(1)%diag%grf_tend_tracer(jc,k,jb,buff(name_idx)) = p_nh_state(1)%diag%grf_tend_tracer(jc,k,jb,buff(idx))*p_art_oem_data%data_fields%boundary_lambdas(jc,jb,nens)
    
                      ENDDO
                    ENDDO
                  ENDIF
                ENDDO
              ENDIF
            ENDDO
          ENDDO
   
          next_latbc_time = next_latbc_time + latbc_config%intv(1)%dtime_latbc_mtime    ! no variable latbc interval supported
          !$ACC UPDATE DEVICE(p_nh_state%diag%grf_tend_tracer)

          CALL message(routine, "OEM's section 2 executed. The time-dependant BCs updated")
  
        END IF

!------------------------------------------------------------------------------------
! End of Section 2
!------------------------------------------------------------------------------------
      !$ACC DATA COPYIN(hod, dow, moy, hoy)
!------------------------------------------------------------------------------------
! Section 3:
!     
! The following section is executed at one hour (physical time) intervals. 
! Here, the temporal scaling factors are evaluated.
!------------------------------------------------------------------------------------
 
      IF (mtime_current>=next_full_hour) THEN
        DO jb = i_startblk, i_endblk
          is = is_m(jb)
          ie = ie_m(jb) 

          DO nc = 1, ncat_max
            DO nt = 1, p_art_oem_data%configure%emis_tracer
              IF (p_art_oem_data%configure%ycatl_l(nt,nc) /= "") THEN

                !$ACC PARALLEL DEFAULT(PRESENT)
                !$ACC LOOP GANG VECTOR
                DO jc = is, ie

                  SELECT CASE (p_art_oem_data%configure%itype_tscale_l(nt))
  
                  CASE(0)

                    tsf_now(jc,nt,nc,jb) =  1._wp
                    tsf_next(jc,nt,nc,jb) =  1._wp

                  CASE(1)

                    tsf_now(jc,nt,nc,jb) = p_art_oem_data%data_fields%tp_hourofday(hod(1), tp_cat_idx(nt,nc), tp_country_idx(jc,jb)) * &
                      &                       p_art_oem_data%data_fields%tp_dayofweek(dow(1), tp_cat_idx(nt,nc), tp_country_idx(jc,jb)) * &
                      &                       p_art_oem_data%data_fields%tp_monthofyear(moy(1), tp_cat_idx(nt,nc), tp_country_idx(jc,jb))

                    tsf_next(jc,nt,nc,jb) = p_art_oem_data%data_fields%tp_hourofday(hod(2), tp_cat_idx(nt,nc), tp_country_idx(jc,jb)) * &
                      &                        p_art_oem_data%data_fields%tp_dayofweek(dow(2), tp_cat_idx(nt,nc), tp_country_idx(jc,jb)) * &
                      &                        p_art_oem_data%data_fields%tp_monthofyear(moy(2), tp_cat_idx(nt,nc), tp_country_idx(jc,jb))

                  CASE(2)

                    tsf_now(jc,nt,nc,jb) = p_art_oem_data%data_fields%tp_hourofyear(hoy(1), tp_cat_idx(nt,nc), tp_country_idx(jc,jb))
                    tsf_next(jc,nt,nc,jb) = p_art_oem_data%data_fields%tp_hourofyear(hoy(2), tp_cat_idx(nt,nc), tp_country_idx(jc,jb))

                  END SELECT

                END DO ! jc = is, ie
                !$ACC END PARALLEL
              END IF ! p_art_oem_data%configure%ycatl_l(nt,nc) /= ""
            END DO ! nt = 1, p_art_oem_data%configure%emis_tracer
          END DO ! nc = 1, ncat_max
        END DO ! jb = i_startblk, i_endblk 

        next_full_hour = next_full_hour + mtime_td

        CALL message(routine, "OEM's section 3 executed. The temporal scaling factors updated")

      END IF ! mtime_current>=next_full_hour
         
!------------------------------------------------------------------------------------
! End of Section 3
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
! Section 4:
!
! This section is exectuted at every time-step and includes the main loops over 
! the tracers. The final temporal scaling facatore is calculated using a linear 
! interpolation to the values obtained in section 2.
! 
! Please note that in the code provided, the locations of the 'gang' and 'vector',
! (and 'worker', if needed) directives are case-sensitive. The GPU porting approach 
! outlined below is optimised for scenarios involving a limited number of tracers  
! and categories, with a high number of ensemble members (typically exceeding 190).      
!------------------------------------------------------------------------------------
      nt_old = 0
      ens_count = 0
      n2 = REAL(hour_sec, wp)/3600._wp
      n1 = 1._wp - n2
      
      DO jb = i_startblk, i_endblk    
        is = is_m(jb)
        ie = ie_m(jb)     
        DO  nt = 1, p_art_oem_data%configure%emis_tracer 

          trcr_idx = p_art_oem_data%configure%emis_idx(nt)

          IF ( ANY( p_art_oem_data%configure%emis_diag_out==trcr_idx ) ) THEN
            p_tracer_now(:,:,jb,trcr_idx) = 0._wp
            
            !$ACC UPDATE DEVICE(p_tracer_now(:,:,jb,trcr_idx)) 
          END IF

          DO nc = 1, ncat_max

            IF (p_art_oem_data%configure%ycatl_l(nt,nc) /= "") THEN
              IF (p_art_oem_data%configure%ens_tracer>0) THEN
                ! Evaluate the ensemble tracers 
                !$ACC PARALLEL DEFAULT(PRESENT)
                IF ( ANY( p_art_oem_data%ensemble%ens_name==p_art_oem_data%configure%emis_name(nt) ) ) THEN
                  nc_ensemble = p_art_oem_data%data_fields%lambda_categories_ids(nc+n_vprm_cats)
                  IF (nt_old/=nt) THEN
                    ens_count = ens_count+1
                    nt_old = nt
                  ENDIF
  
                  !$ACC LOOP GANG COLLAPSE(2) 
                  DO nens = 1,SIZE(p_art_oem_data%data_fields%lambda_mat, dim=4)  
                    DO table_nr=1,p_art_oem_data%configure%ens_tracer
                      IF (p_art_oem_data%ensemble%ens_table(1,table_nr)==nens .AND. p_art_oem_data%ensemble%ens_name(table_nr)==p_art_oem_data%configure%emis_name(nt)) THEN

                        !$ACC LOOP VECTOR COLLAPSE(2) PRIVATE(nr, lambda, temp_scaling_fact, temp_tracer, ens_int_idx)
                        DO k = 1,nlev
                          DO jc = is, ie

                            ens_int_idx = p_art_oem_data%ensemble%ens_table(2,table_nr)

                            nr = p_art_oem_data%data_fields%reg_map(jc,jb)
                            lambda = p_art_oem_data%data_fields%lambda_mat(ens_count,nc_ensemble,nr,nens)
                            temp_scaling_fact = n1*tsf_now(jc,nt,nc,jb)+n2*tsf_next(jc,nt,nc,jb)

                            temp_tracer = total_weight(jc,k,nt,nc,jb) & 
                              &                               * temp_scaling_fact &
                              &                               * lambda &  
                              &                              / p_nh_state(jg)%diag%airmass_now(jc,k,jb)

                            p_tracer_now(jc,k,jb,ens_int_idx) = p_tracer_now(jc,k,jb,ens_int_idx) + temp_tracer

                          ENDDO ! jc = is, ie 
                        ENDDO ! k = 1,nlev
                      ENDIF ! ens_table
                    ENDDO ! table_nr=1,p_art_oem_data%configure%ens_tracer
                  ENDDO ! nens = 1,SIZE(lambda_mat, dim=4)
                ENDIF ! lens(nt)==.TRUE. 
                !$ACC END PARALLEL
              ENDIF 
              
              !$ACC PARALLEL DEFAULT(PRESENT)
              IF ( ANY( p_art_oem_data%configure%emis_diag_out==trcr_idx ) ) THEN  
                !$ACC LOOP GANG VECTOR COLLAPSE(2)   
                DO k = 1,nlev
                  DO jc = is, ie

                    p_tracer_now(jc,k,jb,trcr_idx) = p_tracer_now(jc,k,jb,trcr_idx) + total_weight(jc,k,nt,nc,jb) &
                      &                               * (n1*tsf_now(jc,nt,nc,jb)+n2*tsf_next(jc,nt,nc,jb)) &
                      &                              / dtime
                  ENDDO
                ENDDO

              ELSE
                ! Evaluate the emission tracers
                !$ACC LOOP GANG VECTOR COLLAPSE(2) 
                DO k = 1,nlev
                  DO jc = is, ie
                    
                    p_tracer_now(jc,k,jb,trcr_idx) = p_tracer_now(jc,k,jb,trcr_idx) + total_weight(jc,k,nt,nc,jb) &
                      &                               * (n1*tsf_now(jc,nt,nc,jb)+n2*tsf_next(jc,nt,nc,jb)) &  
                      &                              / p_nh_state(jg)%diag%airmass_now(jc,k,jb)
                  ENDDO
                ENDDO

              ENDIF
              !$ACC END PARALLEL

            ENDIF ! p_art_oem_data%configure%ycatl_l(nt,nc) /= ""

          ENDDO ! nc 
        ENDDO ! nt
      ENDDO ! jb = i_startblk, i_endblk

!------------------------------------------------------------------------------------
! End of Section 4
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
! Section 5:
!     
! The following section is executed just once at a certain time. 
! Here, concentrations are replaced by values obtained from a previous simulation.
!------------------------------------------------------------------------------------
       !$ACC PARALLEL DEFAULT(PRESENT)
       IF (p_art_oem_data%configure%restart_tracer>0) THEN
         IF (.NOT. p_art_oem_data%configure%l_restarted) THEN
           IF (global_sec>=restart_init_time) THEN
             !$ACC LOOP SEQ
             DO jb = i_startblk, i_endblk
               is = is_m(jb) 
               ie = ie_m(jb)
               !$ACC LOOP SEQ
               DO nt=1,p_art_oem_data%configure%restart_tracer
                 trcr_idx = p_art_oem_data%configure%restart_idx(nt)
                 !$ACC LOOP GANG VECTOR COLLAPSE(2)
                 DO k = 1,nlev
                   DO jc = is, ie
                     p_tracer_now(jc,k,jb,trcr_idx) = p_art_oem_data%data_fields%chem_init_3D(jc,k,jb,nt)
                   ENDDO 
                 ENDDO 
               ENDDO
             ENDDO
             p_art_oem_data%configure%l_restarted = .TRUE.
           ENDIF ! global_sec>=restart_init_time
         ENDIF ! .NOT. l_restarted
       ENDIF ! oem_config%restart_tracer>0
       !$ACC END PARALLEL

!------------------------------------------------------------------------------------
! End of Section 5
!------------------------------------------------------------------------------------

      !$ACC END DATA
      CALL message(routine, "OEM executed")
    ENDIF ! p_art_oem_data%configure%emis_tracer>0

  !-----------------------------------------------------------------------------------
  ! End of the Subroutine
  !------------------------------------------------------------------------------------
  
    END SUBROUTINE art_oem_compute_emissions
  
  !====================================================================================
  
    SUBROUTINE art_oem_extract_time_information(date, hour_of_day, day_of_week, month_of_year, &
      &                                         hour_of_year)
  
      IMPLICIT NONE
  
      ! Parameters
      TYPE(datetime), INTENT(IN) :: date
      INTEGER, INTENT(INOUT) :: hour_of_day
      INTEGER, INTENT(INOUT) :: day_of_week
      INTEGER, INTENT(INOUT) :: month_of_year
      INTEGER, INTENT(INOUT) :: hour_of_year
  
      ! Local variables
      INTEGER :: errno, y, d, k, j, moy, doy

      ! Extract data from the datetime object
      y = date%date%year
      month_of_year = date%date%month
      d = date%date%day
      hour_of_day = date%time%hour + 1

      ! We use Zeller's Congruence formula for the day_of_week calculation:
      ! Adjust January and February to be the 13th and 14th months of the previous year
      moy = month_of_year
      IF (moy < 3) THEN
        moy = moy + 12
        y = y - 1
      END IF
      k = MOD(y, 100)
      j = y / 100
      day_of_week = MOD(d + INT((13 * (moy + 1)) / 5) + k + INT(k / 4) + INT(j / 4) + 5 * j, 7)
      ! Adjust to make Monday = 1, ..., Sunday = 7
      day_of_week = MOD(day_of_week + 5, 7) + 1

      ! For the hour_of_year we compute its offset relative to 01 Jan T00 of this year:
      doy = getdayofyearfromdatetime(date)
      hour_of_year = 24 * (doy - 1) + hour_of_day

    END SUBROUTINE art_oem_extract_time_information

  !====================================================================================
  
  END MODULE mo_art_oem_emission
