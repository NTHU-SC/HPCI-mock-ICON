dnl acx_fc_bind_c_integer_align.m4 --- determine alignment of BIND(C) variable
dnl
dnl Copyright  (C)  2023  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords:
dnl Author: Thomas Jahns <jahns@dkrz.de>
dnl Maintainer: Thomas Jahns <jahns@dkrz.de>
dnl URL: https://www.dkrz.de/redmine/projects/scales-ppm
dnl
dnl Redistribution and use in source and binary forms, with or without
dnl modification, are  permitted provided that the following conditions are
dnl met:
dnl
dnl Redistributions of source code must retain the above copyright notice,
dnl this list of conditions and the following disclaimer.
dnl
dnl Redistributions in binary form must reproduce the above copyright
dnl notice, this list of conditions and the following disclaimer in the
dnl documentation and/or other materials provided with the distribution.
dnl
dnl Neither the name of the DKRZ GmbH nor the names of its contributors
dnl may be used to endorse or promote products derived from this software
dnl without specific prior written permission.
dnl
dnl THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
dnl IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
dnl TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
dnl PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
dnl OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
dnl EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
dnl PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
dnl PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
dnl LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
dnl NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
dnl SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
dnl
dnl
dnl ACX_FC_BIND_C_ALIGN([TYPE-DECL],[VAR-TO-ASSIGN],
dnl                     [ACTION-IF-FOUND],[ACTION-IF-NOT-FOUND],
dnl                     [OPTIONAL-HEADER],[OPTIONAL-DEFAULT-VALUE])
dnl
dnl This macro tests the alignment of variables or rather the common
dnl block created by their BIND(C) declaration. This is meant to set a
dnl corresponding alignment in C sources. If VAR-TO-ASSIGN already has
dnl a valid value, uses that, otherwise assigns it from the result of
dnl a compiler check.
dnl
dnl This currently does not work if only MS dumpbin is available, but
dnl there also does not seem to be a relevant platform that only has
dnl that, so implementing that has to wait for a later day.
dnl
m4_define([_ACX_IF_IFVAL],[m4_ifval([$1],[AS_VAR_SET_IF([$1],[$2],[$3])],[$3])])
AC_DEFUN([ACX_FC_BIND_C_ALIGN],
  [AC_REQUIRE([AC_PROG_FC])dnl
   AS_VAR_PUSHDEF([acx_cv_fc_bind_c_align],[acx_cv_fc_bind_c_align_$1])dnl
   AC_MSG_CHECKING([alignment chosen by $FC for $1, PUBLIC, BIND(c)])
   _ACX_IF_IFVAL([$2],[AS_VAR_COPY([acx_cv_fc_bind_c_align],[$2])],
      [AC_CACHE_VAL([acx_cv_fc_bind_c_align],
        [AC_LANG_PUSH([Fortran])dnl
         AC_COMPILE_IFELSE([AC_LANG_SOURCE(
           [      MODULE CONFTEST_ALIGN
      $5
      IMPLICIT NONE
      $1, PUBLIC, BIND(c, name='acx_align_test_variable') :: acx_align_test_variable
      END MODULE CONFTEST_ALIGN])],
           [acx_temp=`"$ac_pwd"/libtool --mode=clean --silent \
              ls conftest.lo 2>/dev/null | grep '\.'"$OBJEXT"'$' | xargs $NM \
              | sed -n 'h
/ C _\{0,1\}acx_align_test_variable$/s/^0*\(@<:@0-9a-f@:>@*\) C .*/0x\1/p
g
/ C _\{0,1\}acx_align_test_variable$/q'`
            AS_CASE([$acx_temp],
              [0x],[# in ipo configurations, alignment will be reported as 0
],
              [0x*],
              [AS_VAR_SET([acx_cv_fc_bind_c_align],[$acx_temp])])])
         AC_LANG_POP([Fortran])])])
   AS_VAR_SET_IF([acx_cv_fc_bind_c_align],
     [m4_ifval([$3],[$3
])dnl
      AS_VAR_COPY([acx_temp],[acx_cv_fc_bind_c_align])],
     [m4_ifval([$4],[$4
])dnl
      m4_ifval([$6],
        [AS_VAR_SET([acx_cv_fc_bind_c_align],[$6])
         AS_VAR_COPY([acx_temp],[acx_cv_fc_bind_c_align])
         acx_temp="using default of $acx_temp"],
        [acx_temp='not found!'])])
   AC_MSG_RESULT([$acx_temp])
   m4_ifval([$2],[AS_VAR_COPY([$2],[acx_cv_fc_bind_c_align])
])dnl
   AS_VAR_POPDEF([acx_cv_fc_bind_c_align])])

dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
