dnl acx_mpi_type_create_args.m4 --- check whether MPI functions to
dnl                                 create derived datatypes use
dnl                                 const qualified array arguments
dnl
dnl Copyright  (C)  2023  Thomas Jahns <jahns@dkrz.de>
dnl
dnl Keywords: configure configure.ac autoconf MPI mpirun mpiexec
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
dnl ACX_MPI_TYPE_CREATE_CONST_ARRAY_ARG([ACTION-IF-ACCEPTS-CONST-ARRAY-ARGS],
dnl   [ACTION-IF-NOT-ACCEPTS-CONST-ARRAY-ARGS])
dnl
AC_DEFUN([ACX_MPI_TYPE_CREATE_CONST_ARRAY_ARGS],
  [AC_CACHE_CHECK([whether MPI_Type_create_hindexed declares array arguments const],
     [acx_cv_mpi_type_create_hindexed_const_array_args],
     [AC_LANG_PUSH([C])
      AC_LANG_CONFTEST([AC_LANG_SOURCE([@%:@include <stdlib.h>
@%:@include <mpi.h>

@%:@define xmpi(ret)           \\
  do {                      \\
    if ((ret) != MPI_SUCCESS) \\
      exit(EXIT_FAILURE);   \\
  } while (0)

extern int
MPI_Type_create_hindexed(int count,
                         const int array_of_blocklengths@<:@@:>@,
                         const MPI_Aint array_of_displacements@<:@@:>@,
                         MPI_Datatype oldtype, MPI_Datatype *newtype);

int main(int argc, char **argv)
{
  xmpi(MPI_Init(&argc, &argv));
  static const int foo@<:@@:>@ = { 1 };
  static const MPI_Aint bar@<:@@:>@ = { 0 };
  MPI_Datatype dt;
  xmpi(MPI_Type_create_hindexed(1, foo, bar, MPI_INT, &dt));
  xmpi(MPI_Type_free(&dt));
  xmpi(MPI_Finalize());
  return EXIT_SUCCESS;
}
])])
      AC_COMPILE_IFELSE(,
        [# compilation worked without error,
         # inspect if removing const errors out or creates extra warnings next
         acx_temp=`cat conftest.err | wc -l`
         sed 's/const //' conftest.c >conftest.er1 ; mv conftest.er1 conftest.c
         AC_COMPILE_IFELSE(,
           [AS_IF([test "$acx_temp" -lt `cat conftest.err | wc -l`],
              [acx_cv_mpi_type_create_hindexed_const_array_args=yes],
              [acx_cv_mpi_type_create_hindexed_const_array_args=no])],
           [acx_cv_mpi_type_create_hindexed_const_array_args=yes])],
        [acx_cv_mpi_type_create_hindexed_const_array_args=no])
      AC_LANG_POP([C])])
    AS_IF([test x"$acx_cv_mpi_type_create_hindexed_const_array_args" = xyes],[$1],
     [m4_default([$2],
        [AC_MSG_FAILURE([MPI_Type_create_hindexed does not accept const qualified array arguments])])])
  ])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
