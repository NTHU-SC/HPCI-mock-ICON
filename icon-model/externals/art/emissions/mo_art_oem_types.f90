!
! mo_art_oem_types
! This module provides datastructures for the
! online emisison module OEM
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

MODULE mo_art_oem_types

! ICON
  USE mo_kind,                          ONLY: wp
  USE mo_async_latbc_types,             ONLY: t_latbc_data
  
  IMPLICIT NONE
  
  PRIVATE

  PUBLIC :: p_art_oem_data, &
            t_art_oem, &
            t_art_oem_data, &
            t_art_oem_config, &
            t_art_oem_ensemble

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  TYPE t_art_oem_data

    INTEGER, DIMENSION(:,:), ALLOCATABLE :: &
      & country_ids                 ! EMEP country code for each domain (nproma,nblks)

    REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
      & gridded_emissions,     & ! "2D" emissions fields (nproma,nblks,ncat)
      & tp_dayofweek,          & ! day-of-week scaling factor (ndays,ncat,ncountries)
      & tp_monthofyear,        & ! seasonal scaling factor
      & tp_hourofday,          & ! diurnal scaling factor
      & tp_hourofyear,         & ! hourly scaling factor
      & boundary_lambdas         ! lambdas for background scaling (nproma,nblks,ntracer,nens)

    REAL(KIND=wp), DIMENSION(:,:,:,:), ALLOCATABLE :: &
      & lambda_mat,            & ! lambdas for ensembles
      & chem_init_3D             ! restart-fields (nproma,nlev,nnblks,n_restart_tracer)

    REAL(KIND=wp), DIMENSION(:,:,:,:), ALLOCATABLE :: &
      & vert_scaling_fact        ! vertical scale factors on model grid

    REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
      & lswi,                  & ! LSWI field (jc,jb,nday)
      & evi,                   & ! EVI field (jc,jb,nday)
      & vprm_lu_class_fraction   ! VPRM land-use class fractions (jc,jb,num_vprm_lu_classes)

    REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
      & lswi_min,              & ! Minimum of LSWI field (jc,jb)
      & evi_min,               & ! Minimum of EVI field (jc,jb)
      & lswi_max,              & ! Minimum of LSWI field (jc,jb)
      & evi_max                  ! Minumum of EVI field (jc,jb)

    INTEGER :: &
      & i_vprm_lc_evergreen,   & !< VPRM land-use class 'evergreen'
      & i_vprm_lc_deciduous,   & !< VPRM land-use class 'deciduous'
      & i_vprm_lc_mixed,       & !< VPRM land-use class 'mixed forest'
      & i_vprm_lc_shrub,       & !< VPRM land-use class 'shrubland'
      & i_vprm_lc_savanna,     & !< VPRM land-use class 'savanna'
      & i_vprm_lc_crop,        & !< VPRM land-use class 'cropland'
      & i_vprm_lc_grass,       & !< VPRM land-use class 'grassland'
      & i_vprm_lc_urban          !< VPRM land-use class 'urban area'

    INTEGER, DIMENSION(:), ALLOCATABLE :: &
      & lambda_categories_ids,            & ! Lambda indeces for different emission categories 
      & lambda_temp_ids                     ! Lambda indeces for different hour of the day  

    REAL(KIND=wp), DIMENSION(:), ALLOCATABLE :: &
      & newflux_vprm                 ! An array of bioflux tendencies
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: &
      & reg_map                  ! region masks for grid cells for ensembles (nreg,nrpoma,nblks)

    INTEGER, DIMENSION(:), ALLOCATABLE :: &
      & tp_countryid             ! EMEP country code

      TYPE(t_latbc_data), POINTER :: p_latbc_data

  END TYPE t_art_oem_data
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  TYPE t_art_oem_config

    INTEGER :: tp_ncountry, & ! number of countries in dataset
      &        emis_tracer, & ! number of emisison tracer
      &        emis_diag,   & ! check for outputting emissions diag field
      &        emis_diag_index,   & ! number of tracers for which we want to output emissions diag field
      &        bg_tracer,   & ! number of background tracer
      &        ens_tracer,  & ! number of ensemble tracer
      &        restart_tracer, & ! number of tracer with field from previous sim.
      &        vprm_tracer    ! number of VPRM tracer

    CHARACTER(LEN=20), DIMENSION(:), ALLOCATABLE   :: &
      & gridded_emissions_idx,  & ! name of the annual mean emissions fields
      & vp_category,            & ! category names of vertical profiles
      & tp_category,            & ! category names of temporal profiles
      & bg_name,                & ! name of the background tracers
      & emis_name,              & ! name of the emission tracer
      & restart_name,           & ! name of the restart tracers
      & vprm_name,              & ! name of the VPRM tracer
      & vprm_flux_type            ! name of the VPRM flux type ('resp' or 'gpp')

    CHARACTER(LEN=256), DIMENSION(:), ALLOCATABLE   :: &
      & latbc_buffer_names


    CHARACTER(LEN=20), DIMENSION(:,:), ALLOCATABLE :: &
      & ycatl_l, yvpl_l, ytpl_l   ! categories, vertical- and temporal profiles (ntracer,ncat)

    INTEGER, DIMENSION(:), ALLOCATABLE             :: &
      & emis_idx,               & ! indices of OEM-tracers with emissions within the ICON container
      & emis_diag_out,          & ! indices of OEM-tracers with emissions within the ICON container for which we want to output emissions
      & bg_idx,                 & ! indices of OEM-tracers as background within the ICON container
      & restart_idx,            & ! indices of OEM-tracers which are restarted within the ICON container
      & tend_tracer_idx,        & ! indices of tracers in latbc tendency data structure
      & vprm_idx,               & ! indices of OEM-tracers with VPRM within the ICON container
      & itype_tscale_l            ! type of temporal scaling

    REAL(KIND=wp), DIMENSION(:), ALLOCATABLE       :: &
      & ra_lifetime               ! lifetimes of radioacive decay

    LOGICAL, DIMENSION(:), ALLOCATABLE :: radioact ! indices of OEM emission tracers with radiactive decay

    LOGICAL :: l_first_tend_scaled, l_boundary_scaling, l_restarted

    REAL(KIND=wp) :: offset_val   ! offset value for initialization


  END TYPE t_art_oem_config

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  TYPE t_art_oem_ensemble

    CHARACTER(LEN=20), DIMENSION(:), ALLOCATABLE  :: &
      & ens_name,               & ! name of the reference tracer, where the ensemble belongs to
      & bg_ens_name               ! name of the reference tracer, where the ensemble belongs to

    CHARACTER(LEN=20), DIMENSION(:,:), ALLOCATABLE  :: &
      & vprm_bg_ens_name          ! name of the 2 reference tracers, where the VPRM BG ensemble belongs to (RA, then GPP, order matters!!!)

    LOGICAL, DIMENSION(:), ALLOCATABLE  :: &
      & vprm_bg_ens               ! check for vprm_bg ensemble tracers

    INTEGER, DIMENSION(:,:), ALLOCATABLE          :: &
      & ens_table                 ! indices of OEM-tracers with emissions within the ICON container


  END TYPE t_art_oem_ensemble

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
  TYPE t_art_oem

    TYPE(t_art_oem_data)               :: data_fields
    TYPE(t_art_oem_config)             :: configure
    TYPE(t_art_oem_ensemble)           :: ensemble

  END TYPE t_art_oem

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

  TYPE(t_art_oem),TARGET  :: &
    &  p_art_oem_data
 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!


END MODULE mo_art_oem_types