/**
 * @file xt_mpi_ddt_wrap.h
 *
 * @copyright Copyright  (C)  2023 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef XT_MPI_DDT_WRAP_H
#define XT_MPI_DDT_WRAP_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef XT_MPI_TYPE_CREATE_INDEXED_BLOCK_CONST_DISP
#  define Xt_Type_create_indexed_block(count, blocklength, disp,        \
                                       oldtype, newtype, comm)          \
  xt_mpi_call(MPI_Type_create_indexed_block(count, blocklength, disp,   \
                                            oldtype, newtype), comm)
#else
#  define Xt_Type_create_indexed_block(count, blocklength, disp,        \
                                       oldtype, newtype, comm)          \
  xt_mpi_call(MPI_Type_create_indexed_block(count, blocklength,         \
                                            (int *)(intptr_t)disp,      \
                                            oldtype, newtype), comm)
#endif

#if MPI_VERSION < 3
static inline int
XtMPI_Type_create_hindexed_block(int count, int blocklength,
                                 const MPI_Aint array_of_displacements[],
                                 MPI_Datatype oldtype, MPI_Datatype *newtype)
{
  size_t count_ = count > 0 ? (size_t)count : 0;
  enum { bl_auto_max = 32 };
  int blocklengths_auto[bl_auto_max];
  int *restrict blocklengths
    = count_ > bl_auto_max
    ? xmalloc(count_ * sizeof (*blocklengths))
    : blocklengths_auto;
  for (size_t i = 0; i < count_; ++i)
    blocklengths[i] = blocklength;
  int rc = MPI_Type_create_hindexed(count, blocklengths,
                                    CAST_MPI_SEND_BUF(array_of_displacements),
                                    oldtype, newtype);
  if (count_ > bl_auto_max)
    free(blocklengths);
  return rc;
}

#  define MPI_Type_create_hindexed_block XtMPI_Type_create_hindexed_block
#  define XT_MPI_TYPE_CREATE_HINDEXED_BLOCK_CONST_DISP
#  define Xt_Type_create_hindexed_block(count, blocklength,             \
                                        array_of_displacements,         \
                                        oldtype, newtype, comm)         \
  xt_mpi_call(XtMPI_Type_create_hindexed_block(count, blocklength,      \
                                               array_of_displacements,  \
                                               oldtype, newtype), comm)
enum { MPI_COMBINER_HINDEXED_BLOCK=128 };
#elif defined XT_MPI_TYPE_CREATE_HINDEXED_BLOCK_CONST_DISP
#  define Xt_Type_create_hindexed_block(count, blocklength,             \
                                        array_of_displacements,         \
                                        oldtype, newtype, comm)         \
  xt_mpi_call(MPI_Type_create_hindexed_block(count, blocklength,        \
                                             array_of_displacements,    \
                                             oldtype, newtype), comm)
#else
#  define Xt_Type_create_hindexed_block(count, blocklength,             \
                                        array_of_displacements,         \
                                        oldtype, newtype, comm)         \
  xt_mpi_call(MPI_Type_create_hindexed_block(                           \
                count, blocklength,                                     \
                (MPI_Aint *)(intptr_t)array_of_displacements,           \
                oldtype, newtype), comm)
#endif

#ifdef XT_MPI_TYPE_INDEXED_CONST_ARRAY_ARGS
#  define Xt_Type_indexed(count, blocklength, disp, oldtype, newtype,   \
                          comm)                                         \
  xt_mpi_call(MPI_Type_indexed(count, blocklength,                      \
                               disp, oldtype, newtype), comm)
#else
#  define Xt_Type_indexed(count, blocklength, disp, oldtype, newtype,   \
                          comm)                                         \
  xt_mpi_call(MPI_Type_indexed(count, (int *)(intptr_t)blocklength,     \
                               (int *)(intptr_t)disp, oldtype,          \
                               newtype), comm)
#endif

#ifdef XT_MPI_TYPE_CREATE_HINDEXED_CONST_ARRAY_ARGS
#  define Xt_Type_create_hindexed(count, blocklength, disp, oldtype,    \
                                  newtype, comm)                        \
  xt_mpi_call(MPI_Type_create_hindexed(count, blocklength,              \
                                       disp, oldtype, newtype), comm)
#else
#  define Xt_Type_create_hindexed(count, blocklength, disp, oldtype,    \
                                  newtype, comm)                        \
  xt_mpi_call(MPI_Type_create_hindexed(count,                           \
                                       (int *)(intptr_t)blocklength,    \
                                       (MPI_Aint *)(intptr_t)disp,      \
                                       oldtype, newtype), comm)
#endif

#ifdef XT_MPI_TYPE_CREATE_STRUCT_CONST_ARRAY_ARGS
#  define Xt_Type_create_struct(count, array_of_blocklengths,           \
                                array_of_displacements,                 \
                                array_of_types,                         \
                                newtype, comm)                          \
  xt_mpi_call(MPI_Type_create_struct(                                   \
                count, array_of_blocklengths,                           \
                array_of_displacements,                                 \
                array_of_types,                                         \
                newtype), comm)
#else
#  define Xt_Type_create_struct(count, array_of_blocklengths,           \
                                array_of_displacements,                 \
                                array_of_types,                         \
                                newtype, comm)                          \
  xt_mpi_call(MPI_Type_create_struct(                                   \
                count, (int *)(intptr_t)array_of_blocklengths,          \
                (MPI_Aint*)(intptr_t)array_of_displacements,            \
                (MPI_Datatype *)(intptr_t)array_of_types,               \
                newtype), comm)
#endif

#ifdef XT_MPI_TYPE_INDEXED_CONST_ARRAY_ARGS
#  define Xt_Type_create_subarray(ndims, array_of_sizes,                \
                                  array_of_subsizes, array_of_starts,   \
                                  order, oldtype, newtype, comm)        \
  xt_mpi_call(MPI_Type_create_subarray(ndims, array_of_sizes,           \
                                       array_of_subsizes,               \
                                       array_of_starts,                 \
                                       order, oldtype, newtype), comm)
#else
#  define Xt_Type_create_subarray(ndims, array_of_sizes,                \
                                  array_of_subsizes, array_of_starts,   \
                                  order, oldtype, newtype, comm)        \
  xt_mpi_call(MPI_Type_create_subarray(                                 \
    ndims, (int *)(intptr_t)array_of_sizes,                             \
    (int *)(intptr_t)array_of_subsizes,                                 \
    (int *)(intptr_t)array_of_starts,                                   \
    order, oldtype, newtype), comm)
#endif

#endif // XT_MPI_DDT_WRAP_H

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
