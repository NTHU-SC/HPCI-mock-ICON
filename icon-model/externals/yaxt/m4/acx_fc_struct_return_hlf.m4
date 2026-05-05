dnl acx_fc_struct_return_hlf.m4 --- check if Fortran compiler
dnl                                 correctly handles higher-level
dnl                                 functions whose argument returns
dnl                                 a struct.
dnl
dnl
dnl Copyright  (C)  2024  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Keywords:
dnl Author: Thomas Jahns <jahns@dkrz.de>
dnl Maintainer: Thomas Jahns <jahns@dkrz.de>
dnl URL: https://swprojects.dkrz.de/redmine/projects/scales-ppm
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
dnl Commentary:
dnl
dnl
dnl
dnl Code:
dnl
# ACX_FC_STRUCT_RETURN_HLF([ACTION-IF-WORKING],[ACTION-IF-NOT-WORKING])
# -----------------
# Checks whether the compiler supports procedures that accept function
# arguments which in turn return a struct. This is unfortunately buggy
# in some older versions of the NAG compiler.
#
# This macro depends, of course, on the Fortran compiler producing
# module files. See comment to AC_FC_MOD_PATH_FLAG.
#
#
AC_DEFUN([ACX_FC_STRUCT_RETURN_HLF],
  [AC_CACHE_CHECK([for struct returning function arguments],
     [acx_cv_fc_struct_return_hlf],
     [AC_LANG_PUSH([Fortran])
      AC_COMPILE_IFELSE(
[      MODULE conftest_module
      TYPE apple
      END TYPE
      INTERFACE
        SUBROUTINE b(c, d)
        IMPORT
          INTERFACE
            FUNCTION c() RESULT(banana)
              IMPORT apple
              TYPE(apple) :: banana
            END FUNCTION c
          END INTERFACE
          INTEGER, VALUE :: d
        END SUBROUTINE b
      END INTERFACE
      END MODULE conftest_module
],
        [AC_COMPILE_IFELSE([      PROGRAM conftest
      USE conftest_module
      END PROGRAM conftest],
        [acx_cv_fc_struct_return_hlf=yes],
        [acx_cv_fc_struct_return_hlf=no])],
        [acx_cv_fc_struct_return_hlf=no])
      rm -f conftest_module* CONFTEST_MODULE*
dnl Some Fortran compilers create module files not in the current working directory but
dnl in the directory with the object file, therefore we try to delete everything:
   AS_IF([expr "$ac_compile" : '.*/libtool --mode=compile' >/dev/null],
     [AS_IF([test -n "$objdir"],
        [rm -f "$objdir"/conftest_module* "$objdir"/CONFTEST_MODULE*])])
   AC_LANG_POP([Fortran])])
dnl
   AS_IF([test x"$acx_cv_fc_struct_return_hlf" = xyes],
     [$1],
     [m4_default([$2],
       [AC_MSG_WARN([Could not build Fortran program using module procedure with struct-returning function argument!])])])
dnl
])dnl ACX_FC_STRUCT_RETURN_HLF
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://swprojects.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
