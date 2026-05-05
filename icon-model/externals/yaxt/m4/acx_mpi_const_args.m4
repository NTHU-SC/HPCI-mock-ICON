dnl acx_mpi_const_args.m4 --- check whether some MPI functions take
dnl                           const-qualified arguments
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
dnl _ACX_MPI_ARGCHECK_TEST([prototype-to-check])
m4_define([_ACX_MPI_ARGCHECK_TEST],
  [AC_LANG_CONFTEST([AC_LANG_SOURCE([@%:@include <stdlib.h>
@%:@include <mpi.h>

@%:@define xmpi(ret)           \\
  do {                      \\
    if (ret != MPI_SUCCESS) \\
      exit(EXIT_FAILURE);   \\
  } while (0)

$1

int main(int argc, char **argv)
{
  xmpi(MPI_Init(&argc, &argv));
  static const int foo = 1;
  int rank, baz;
  xmpi(MPI_Comm_rank(MPI_COMM_WORLD, &rank));
  if (rank == 0)
    xmpi(MPI_Send(&foo, 1, MPI_INT, 1, 1, MPI_COMM_WORLD));
  else if (rank == 1)
    xmpi(MPI_Recv(&baz, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE));
  xmpi(MPI_Finalize());
  return EXIT_SUCCESS;
}
])])])
dnl _ACX_CHECK_MPI_CONST_COMPILE([CACHE_VAR_TO_SET])
m4_define([_ACX_CHECK_MPI_CONST_COMPILE],
  [AC_COMPILE_IFELSE(,
     [# compilation worked without error,
      # inspect if removing const errors out or creates extra warnings next
      acx_temp=`cat conftest.err | wc -l`
      sed 's/const //' conftest.c >conftest.er1 ; mv conftest.er1 conftest.c
      AC_COMPILE_IFELSE(,
        [AS_IF([test "$acx_temp" -lt `cat conftest.err | wc -l`],
           [$1=yes],
           [$1=no])],
        [$1=yes])],
     [$1=no])])
dnl
dnl _ACX_CHECK_MPI_CONST_ARGS([FUNCNAME],[TEST_DECL],[CACHE_VAR_TO_SET],
dnl                           [CHECK_DESCRIPTION],
dnl                           [ACTION-IF-ACCEPTS-CONST],
dnl                           [ACTION-IF-NOT-ACCEPTS-CONST])
m4_define([_ACX_CHECK_MPI_CONST_ARGS],
  [AC_CACHE_CHECK([whether $1 accepts $4],
     [$3],
     [AC_LANG_PUSH([C])
      _ACX_MPI_ARGCHECK_TEST([$2])
      _ACX_CHECK_MPI_CONST_COMPILE([$3])
      AC_LANG_POP([C])])
    AS_IF([test x"$$3" = xyes],[$5],
     [m4_default([$6],
        [AC_MSG_FAILURE([$1 does not accept $4])])])])
dnl
dnl ACX_MPI_SEND_CONST_VOID_P_BUF_ARG([ACTION-IF-ACCEPTS-CONST-VOID-P],
dnl   [ACTION-IF-NOT-ACCEPTS-CONST-VOID-P])
dnl
AC_DEFUN([ACX_MPI_SEND_CONST_VOID_P_BUF_ARG],
  [_ACX_CHECK_MPI_CONST_ARGS([MPI_Send],
    [extern int MPI_Send(const void *buf, int count, MPI_Datatype datatype,
                    int dest, int tag, MPI_Comm comm);],
    [acx_cv_mpi_send_takes_const_void],
    [const void * as first argument], [$1], [$2])])
dnl
dnl
dnl ACX_MPI_GET_ADDRESS_CONST_VOID_P_LOCATION_ARG(
dnl   [ACTION-IF-ACCEPTS-CONST-VOID-P],
dnl   [ACTION-IF-NOT-ACCEPTS-CONST-VOID-P])
dnl
AC_DEFUN([ACX_MPI_GET_ADDRESS_CONST_VOID_P_LOCATION_ARG],
  [_ACX_CHECK_MPI_CONST_ARGS([MPI_Get_address],
     [extern int MPI_Get_address(const void *location, MPI_Aint *address);],
     [acx_cv_mpi_get_address_takes_const_void],
     [const void * as first argument], [$1], [$2])])
dnl
dnl
dnl ACX_MPI_TYPE_CREATE_STRUCT_CONST_ARRAY_ARGS(
dnl   [ACTION-IF-ACCEPTS-CONST-ARRAYS],
dnl   [ACTION-IF-NOT-ACCEPTS-CONST-ARRAYS])
dnl
AC_DEFUN([ACX_MPI_TYPE_CREATE_STRUCT_CONST_ARRAY_ARGS],
  [_ACX_CHECK_MPI_CONST_ARGS([MPI_Type_create_struct],
     [extern int MPI_Type_create_struct(int count, const int array_of_blocklengths@<:@@:>@,
            const MPI_Aint array_of_displacements@<:@@:>@, const MPI_Datatype array_of_types@<:@@:>@,
            MPI_Datatype *newtype);],
     [acx_cv_mpi_type_create_struct_takes_const_arrays],
     [const array arguments], [$1], [$2])])
dnl
dnl ACX_MPI_TYPE_CREATE_HINDEXED_BLOCK_CONST_DISP(
dnl   [ACTION-IF-ACCEPTS-CONST-DISPLACEMENTS],
dnl   [ACTION-IF-NOT-ACCEPTS-CONST-DISPLACEMENTS])
dnl
AC_DEFUN([ACX_MPI_TYPE_CREATE_HINDEXED_BLOCK_CONST_DISP],
  [_ACX_CHECK_MPI_CONST_ARGS([MPI_Type_create_hindexed_block],
     [extern int MPI_Type_create_hindexed_block(int count, int blocklength,
            const MPI_Aint array_of_displacements@<:@@:>@, MPI_Datatype oldtype,
            MPI_Datatype *newtype);],
     [acx_cv_mpi_type_create_hindexed_block_takes_const_disp],
     [const-qualified displacements], [$1], [$2])])
dnl
dnl ACX_MPI_TYPE_CREATE_INDEXED_BLOCK_CONST_DISP(
dnl   [ACTION-IF-ACCEPTS-CONST-DISPLACEMENTS],
dnl   [ACTION-IF-NOT-ACCEPTS-CONST-DISPLACEMENTS])
dnl
AC_DEFUN([ACX_MPI_TYPE_CREATE_INDEXED_BLOCK_CONST_DISP],
  [_ACX_CHECK_MPI_CONST_ARGS([MPI_Type_create_indexed_block],
     [extern int MPI_Type_create_indexed_block(int count, int blocklength,
            const int array_of_displacements@<:@@:>@, MPI_Datatype oldtype,
            MPI_Datatype *newtype);],
     [acx_cv_mpi_type_create_indexed_block_takes_const_disp],
     [const-qualified displacements], [$1], [$2])])
dnl
dnl ACX_MPI_TYPE_INDEXED_CONST_ARRAY_ARGS(
dnl   [ACTION-IF-ACCEPTS-CONST-ARRAY-ARGS],
dnl   [ACTION-IF-NOT-ACCEPTS-CONST-ARRAY-ARGS])
dnl
AC_DEFUN([ACX_MPI_TYPE_INDEXED_CONST_ARRAY_ARGS],
  [_ACX_CHECK_MPI_CONST_ARGS([MPI_Type_indexed],
     [extern int MPI_Type_indexed(int count, const int array_of_blocklengths@<:@@:>@,
            const int array_of_displacements@<:@@:>@, MPI_Datatype oldtype,
            MPI_Datatype *newtype);],
     [acx_cv_mpi_type_indexed_takes_const_arrays],
     [const array arguments], [$1], [$2])])
dnl
dnl ACX_MPI_TYPE_CREATE_HINDEXED_CONST_ARRAY_ARGS(
dnl   [ACTION-IF-ACCEPTS-CONST-ARRAY-ARGS],
dnl   [ACTION-IF-NOT-ACCEPTS-CONST-ARRAY-ARGS])
dnl
AC_DEFUN([ACX_MPI_TYPE_CREATE_HINDEXED_CONST_ARRAY_ARGS],
  [_ACX_CHECK_MPI_CONST_ARGS([MPI_Type_create_hindexed],
     [extern int MPI_Type_create_hindexed(int count, const int array_of_blocklengths@<:@@:>@,
            const MPI_Aint array_of_displacements@<:@@:>@, MPI_Datatype oldtype,
            MPI_Datatype *newtype);],
     [acx_cv_mpi_type_create_hindexed_takes_const_arrays],
     [const array arguments], [$1], [$2])])
dnl
dnl ACX_MPI_TYPE_CREATE_SUBARRAY_CONST_ARRAY_ARGS(
dnl   [ACTION-IF-ACCEPTS-CONST-ARRAY-ARGS],
dnl   [ACTION-IF-NOT-ACCEPTS-CONST-ARRAY-ARGS])
dnl
AC_DEFUN([ACX_MPI_TYPE_CREATE_SUBARRAY_CONST_ARRAY_ARGS],
  [_ACX_CHECK_MPI_CONST_ARGS([MPI_Type_create_subarray],
     [extern int MPI_Type_create_subarray(int ndims, const int array_of_sizes@<:@@:>@,
            const int array_of_subsizes@<:@@:>@, const int array_of_starts@<:@@:>@,
            int order, MPI_Datatype oldtype, MPI_Datatype *newtype);],
     [acx_cv_mpi_type_create_subarray_takes_const_arrays],
     [const array arguments], [$1], [$2])])
dnl
dnl Local Variables:
dnl mode: autoconf
dnl license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
dnl license-default: "bsd"
dnl End:
