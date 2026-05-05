!>
!! @file xt_xmap_dist_dir_bucket_gen_f.f90
!! @brief Fortran interface to yaxt bucket generator declarations
!!
!! @copyright Copyright  (C)  2024 Jörg Behrens <behrens@dkrz.de>
!!                                 Moritz Hanke <hanke@dkrz.de>
!!                                 Thomas Jahns <jahns@dkrz.de>
!!
!! @author Jörg Behrens <behrens@dkrz.de>
!!         Moritz Hanke <hanke@dkrz.de>
!!         Thomas Jahns <jahns@dkrz.de>
!!

!
! Keywords:
! Maintainer: Jörg Behrens <behrens@dkrz.de>
!             Moritz Hanke <hanke@dkrz.de>
!             Thomas Jahns <jahns@dkrz.de>
! URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
!
! Redistribution and use in source and binary forms, with or without
! modification, are  permitted provided that the following conditions are
! met:
!
! Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! Neither the name of the DKRZ GmbH nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
#include "fc_feature_defs.inc"
MODULE xt_xmap_dist_dir_bucket_gen
  USE iso_c_binding, ONLY: c_ptr, c_null_ptr, c_associated
  IMPLICIT NONE
  PRIVATE

  TYPE, PUBLIC, BIND(c) :: xt_xmdd_bucket_gen
#ifndef __G95__
    PRIVATE
#endif
    TYPE(c_ptr) :: cptr = c_null_ptr
  END TYPE xt_xmdd_bucket_gen

  INTERFACE
    ! this function must not be implemented in Fortran because
    ! PGI 11.x chokes on that
    FUNCTION xt_xmdd_bucket_gen_f2c(gen) &
         BIND(c, name='xt_xmdd_bucket_gen_f2c') RESULT(p)
      IMPORT :: c_ptr, xt_xmdd_bucket_gen
      IMPLICIT NONE
      TYPE(xt_xmdd_bucket_gen), INTENT(in) :: gen
      TYPE(c_ptr) :: p
    END FUNCTION xt_xmdd_bucket_gen_f2c
  END INTERFACE

  PUBLIC :: xt_xmdd_bucket_gen_new, xt_xmdd_bucket_gen_delete
  PUBLIC :: xt_xmdd_bucket_gen_c2f, xt_xmdd_bucket_gen_f2c

  INTERFACE xt_is_null
    MODULE PROCEDURE xt_xmdd_bucket_gen_is_null
  END INTERFACE xt_is_null
  PUBLIC :: xt_is_null

CONTAINS

  FUNCTION xt_xmdd_bucket_gen_new() RESULT(gen)
    TYPE(xt_xmdd_bucket_gen) :: gen
    INTERFACE
      FUNCTION xt_xmdd_bucket_gen_new_c() RESULT(gen) &
           BIND(c, name='xt_xmdd_bucket_gen_new')
        IMPORT :: c_ptr
        IMPLICIT NONE
        TYPE(c_ptr) :: gen
      END FUNCTION xt_xmdd_bucket_gen_new_c
    END INTERFACE
    gen%cptr = xt_xmdd_bucket_gen_new_c()
  END FUNCTION xt_xmdd_bucket_gen_new

  SUBROUTINE xt_xmdd_bucket_gen_delete(gen)
    TYPE(xt_xmdd_bucket_gen), INTENT(in) :: gen
    INTERFACE
      SUBROUTINE xt_xmdd_bucket_gen_delete_c(gen) &
           BIND(c, name='xt_xmdd_bucket_gen_delete')
        IMPORT :: c_ptr
        IMPLICIT NONE
        TYPE(c_ptr), VALUE, INTENT(in) :: gen
      END SUBROUTINE xt_xmdd_bucket_gen_delete_c
    END INTERFACE
    CALL xt_xmdd_bucket_gen_delete_c(gen%cptr)
  END SUBROUTINE xt_xmdd_bucket_gen_delete

  FUNCTION xt_xmdd_bucket_gen_c2f(gen) RESULT(p)
    TYPE(c_ptr), INTENT(in) :: gen
    TYPE(xt_xmdd_bucket_gen) :: p
    p%cptr = gen
  END FUNCTION xt_xmdd_bucket_gen_c2f

  FUNCTION xt_xmdd_bucket_gen_is_null(bucket_gen) RESULT(p)
    TYPE(xt_xmdd_bucket_gen), INTENT(in) :: bucket_gen
    LOGICAL :: p
    p = .NOT. C_ASSOCIATED(bucket_gen%cptr)
  END FUNCTION xt_xmdd_bucket_gen_is_null

END MODULE xt_xmap_dist_dir_bucket_gen
!
! Local Variables:
! f90-continuation-indent: 5
! coding: utf-8
! indent-tabs-mode: nil
! show-trailing-whitespace: t
! require-trailing-newline: t
! End:
!
