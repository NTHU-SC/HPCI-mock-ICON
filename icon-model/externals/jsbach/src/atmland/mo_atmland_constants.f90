!> Atmosphere-land constants, conversion factors, plus init routine
!>
!> ICON-Land
!>
!> ---------------------------------------
!> Copyright (C) 2013-2026, MPI-M, MPI-BGC
!>
!> Contact: icon-model.org
!> Authors: AUTHORS.md
!> See LICENSES/ for license information
!> SPDX-License-Identifier: BSD-3-Clause
!> ---------------------------------------
!>
MODULE mo_atmland_constants
#ifndef __NO_JSBACH__

  USE mo_kind,            ONLY: wp

  IMPLICIT NONE

  ! solar constant
  REAL(wp),SAVE :: frac_vis_swtot_srf          !< fraction of solar shortwave radiation in the VIS band (here: 320-700nm for albedo)
  REAL(wp),SAVE :: frac_par_swvis_srf          !< fraction of visible solar shortwave radiation in the PAR band (400-700nm)

  ! wind speed
  REAL(wp),SAVE :: min_wind                    !< minimum wind speed allowed to prevent total decoupling of the surface (m/s)
  REAL(wp),SAVE :: max_wind                    !< maximum wind speed allowed to prevent excess cooling of the surface (m/s)
  REAL(wp),SAVE :: def_wind                    !< wind (m/s)

  ! conversion factors
  REAL(wp),SAVE :: standard_press_srf          !< standard pressure in Pa
  REAL(wp),SAVE :: eps_vpd                     !< eps_vpd used for unit conversion Pa -> g g-1

  ! default values of model forcing for model testing:
  !  @TODO some of the 'default values of model forcing' might not be used (and can be removed ?)
  REAL(wp),SAVE :: def_t_air                   !< temperature in deg K
  REAL(wp),SAVE :: def_vpd                     !< vpd in Pa
  REAL(wp),SAVE :: def_swdown                  !< PAR in W m-2
  REAL(wp),SAVE :: def_fdiffuse                !< diffuse fraction
  REAL(wp),SAVE :: def_cos_angle               !< cosine of zenith angle
  REAL(wp),SAVE :: def_co2_mixing_ratio        !< atmospheric CO2 in ppm
  REAL(wp),SAVE :: def_co2_mixing_ratio_C13    !< atmospheric 13CO2 in ppm
  REAL(wp),SAVE :: def_co2_mixing_ratio_C14    !< atmospheric 14CO2 in 10^-12 ppm
  REAL(wp),SAVE :: def_co2_deltaC13            !< molar mixing ratio of 13C / 12C of atmospheric CO2
  REAL(wp),SAVE :: def_co2_deltaC14            !< molar mixing ratio of 14C / 12C of atmospheric CO2

  PUBLIC

  CHARACTER(len=*), PARAMETER, PRIVATE :: modname = 'mo_atmland_constants'

  !$ACC DECLARE CREATE(frac_vis_swtot_srf, frac_par_swvis_srf)
  !$ACC DECLARE CREATE(min_wind, max_wind, def_wind, standard_press_srf, eps_vpd)
  !$ACC DECLARE CREATE(def_t_air, def_vpd, def_swdown, def_fdiffuse, def_cos_angle)
  !$ACC DECLARE CREATE(def_co2_mixing_ratio, def_co2_deltaC13, def_co2_mixing_ratio_C13)
  !$ACC DECLARE CREATE(def_co2_deltaC14, def_co2_mixing_ratio_C14)

CONTAINS

  !-----------------------------------------------------------------------------------------------------
  !> initialise atmland constants
  !! called by jsbach_setup_models() !
  !!
  !! initialise Earth's orbit; solar angle calculation, solar constant, diffuse fraction calculation
  !! diurnal temperature calculation, conversion factors, default values
  !-----------------------------------------------------------------------------------------------------
  SUBROUTINE init_atmland_constants

    USE mo_isotope_util,          ONLY: calc_mixing_ratio_C13C12, calc_mixing_ratio_C14C

    IMPLICIT NONE
    ! ---------------------------
    ! 0.1 InOut
    !
    ! ---------------------------
    ! 0.2 Local
    CHARACTER(len=*), PARAMETER :: routine = TRIM(modname)//':init_atmland_constants'

    !> initialise solar constant
    !!
    frac_vis_swtot_srf        = 0.54_wp    ! UV-A+PAR / UV-A+PAR+NIR; Monteith & Unsworth 1995
    frac_par_swvis_srf        = 0.9_wp     ! PAR/PAR+UV-A

    !> wind speed
    !!
    min_wind                  = 0.3_wp
    max_wind                  = 5.0_wp
    def_wind                  = 1.0_wp

    !> conversion factors
    !!
    standard_press_srf        = 101300.0_wp
    eps_vpd                   = 0.622_wp

    !> default values
    !!
    def_t_air                 = 298.15_wp
    def_vpd                   = 1000._wp
    def_swdown                = 652.173913_wp
    def_fdiffuse              = 0.2_wp
    def_cos_angle             = 1.0_wp
    def_co2_mixing_ratio      = 380.0_wp
    def_co2_deltaC13          = -6.6_wp     ! CMIP7 PI control 1850 reference value, cf. Graven, H. (2025)
                                            ! doi: 10.25981/ESGF.input4MIPs.CMIP7/2583902.
    def_co2_mixing_ratio_C13  = def_co2_mixing_ratio / (1._wp + 1._wp/calc_mixing_ratio_C13C12(def_co2_deltaC13))
    def_co2_deltaC14          = -3.2_wp     ! CMIP7 PI control 1850 reference value for NH, cf. Graven, H. (2025)
                                            ! doi: 10.25981/ESGF.input4MIPs.CMIP7/2583902.
    def_co2_mixing_ratio_C14  = calc_mixing_ratio_C14C(def_co2_deltaC13, def_co2_deltaC14) * def_co2_mixing_ratio

    !$ACC UPDATE DEVICE(frac_vis_swtot_srf, frac_par_swvis_srf) ASYNC(1)
    !$ACC UPDATE DEVICE(min_wind, max_wind, def_wind, standard_press_srf, eps_vpd) ASYNC(1)
    !$ACC UPDATE DEVICE(def_t_air, def_vpd, def_swdown, def_fdiffuse, def_cos_angle) ASYNC(1)
    !$ACC UPDATE DEVICE(def_co2_mixing_ratio, def_co2_deltaC13, def_co2_mixing_ratio_C13) ASYNC(1)
    !$ACC UPDATE DEVICE(def_co2_deltaC14, def_co2_mixing_ratio_C14) ASYNC(1)

  END SUBROUTINE init_atmland_constants

#endif
END MODULE mo_atmland_constants
