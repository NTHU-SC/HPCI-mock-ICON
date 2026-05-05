! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

MODULE mo_wave_io_config
  USE mo_name_list_output_config, ONLY:is_variable_in_output_dom

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_wave_var_in_output
  PUBLIC :: init_wave_var_in_output
  PUBLIC :: wave_var_in_output

  ! Derived type to collect logical variables indicating if optional diagnostics are requested for output
  TYPE t_wave_var_in_output
     !
     ! diagnostics for Stokes drift vertical profile
     LOGICAL :: last_idx_depth = .FALSE. !< Flag. TRUE if the storage is required
     LOGICAL :: kbar           = .FALSE. ! --//--
     LOGICAL :: T_stokes       = .FALSE. ! --//--
     LOGICAL :: u3d_stokes     = .FALSE. ! --//--
     LOGICAL :: v3d_stokes     = .FALSE. ! --//--

  END type t_wave_var_in_output

  TYPE(t_wave_var_in_output), ALLOCATABLE :: wave_var_in_output(:)

CONTAINS

  !! Precomputation of derived type collecting logical variables indicating whether
  !! optional diagnostics are requested in the output namelists
  !!
  !! Replaces repeated calculations of the same that used to be scattered around various places in the model code
  !!
  SUBROUTINE init_wave_var_in_output(n_dom)

    INTEGER, INTENT(in)  :: n_dom  ! number of model domains

    INTEGER :: jg

    ALLOCATE(wave_var_in_output(n_dom))

    ! diagnostics for Stokes drift vertical profile
    DO jg=1,n_dom
      wave_var_in_output(jg)%last_idx_depth = is_variable_in_output_dom(var_name="last_idx_depth", jg=jg)
      wave_var_in_output(jg)%kbar           = is_variable_in_output_dom(var_name="kbar"          , jg=jg)
      wave_var_in_output(jg)%T_stokes       = is_variable_in_output_dom(var_name="T_stokes"      , jg=jg)
      wave_var_in_output(jg)%u3d_stokes     = is_variable_in_output_dom(var_name="u3d_stokes"    , jg=jg)
      wave_var_in_output(jg)%v3d_stokes     = is_variable_in_output_dom(var_name="v3d_stokes"    , jg=jg)
    END DO

  END SUBROUTINE init_wave_var_in_output

END MODULE mo_wave_io_config
