/**
 * @file xt_mpi_ddt_gen.c
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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <stdlib.h>

#include <mpi.h>
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "xt_mpi_ddt_cache.h"

//! COMPACT_DT enables the analysis of displacements in order to give a
//! more compact description to the datatype generators of MPI. For strong
//! enough MPI implementations this not required. Then you can undefine
//! COMPACT_DT and save some processing time within yaxt without losing communication
//! performance.
#define COMPACT_DT

#ifndef COMPACT_DT
static MPI_Datatype copy_mpi_datatype(MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_dup(old_type, &datatype), comm);

  return datatype;
}

static MPI_Datatype
gen_mpi_datatype_simple(int displacement, MPI_Datatype old_type, MPI_Comm comm)
{
  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_create_indexed_block(1, 1, &displacement, old_type,
                                                &datatype), comm);
  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static MPI_Datatype
gen_mpi_datatype_contiguous(int displacement, int blocklength,
                            MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  if (displacement == 0)
    xt_mpi_call(MPI_Type_contiguous(blocklength, old_type, &datatype),
                    comm);
  else
    xt_mpi_call(MPI_Type_create_indexed_block(1, blocklength,
                                                  &displacement, old_type,
                                                  &datatype), comm);

  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;

}

static MPI_Datatype
gen_mpi_datatype_vector(int count, int blocklength, int stride,
                        int offset, MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_vector(count, blocklength, stride, old_type,
                              &datatype), comm);
  if (offset != 0) {

    MPI_Datatype datatype_;
    int hindexed_blocklength = 1;
    MPI_Aint old_type_size, old_type_lb;

    xt_mpi_call(MPI_Type_get_extent(old_type, &old_type_lb,
                                    &old_type_size), comm);

    MPI_Aint displacement = offset * old_type_size;

    xt_mpi_call(MPI_Type_create_hindexed(1, &hindexed_blocklength,
                                         &displacement, datatype, &datatype_),
                comm);
    xt_mpi_call(MPI_Type_free(&datatype), comm);
    datatype = datatype_;
  }
  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static MPI_Datatype
gen_mpi_datatype_indexed_block(int const * displacements, int blocklength,
                               int count, MPI_Datatype old_type, MPI_Comm comm)
{
  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_create_indexed_block(count, blocklength,
                                                (void *)displacements,
                                                old_type, &datatype), comm);
  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static MPI_Datatype
gen_mpi_datatype_indexed(const int *displacements, const int *blocklengths,
                         int count, MPI_Datatype old_type, MPI_Comm comm) {

  MPI_Datatype datatype;

  xt_mpi_call(MPI_Type_indexed(count, (int*)blocklengths, (void*)displacements,
                                   old_type, &datatype), comm);
  xt_mpi_call(MPI_Type_commit(&datatype), comm);

  return datatype;
}

static inline int
check_for_vector_type(const int *displacements, const int *blocklengths,
                      int count) {

  int blocklength = blocklengths[0];

  for (int i = 1; i < count; ++i)
    if (blocklengths[i] != blocklength)
      return 0;

  int stride = displacements[1] - displacements[0];

  for (int i = 1; i + 1 < count; ++i)
    if (displacements[i+1] - displacements[i] != stride)
      return 0;

  return 1;
}

static inline int check_for_indexed_block_type(const int *blocklengths,
                                               int count) {

  int blocklength = blocklengths[0];

  for (int i = 1; i < count; ++i)
    if (blocklengths[i] != blocklength)
      return 0;

  return 1;
}
#endif

#ifdef COMPACT_DT
static MPI_Datatype
xt_mpi_generate_compact_datatype_block(const int *disp, const int *blocklengths,
                                       int count, MPI_Datatype old_type,
                                       MPI_Comm comm);

static MPI_Datatype
xt_mpi_generate_compact_datatype(int const *disp, int disp_len,
                                 MPI_Datatype old_type, MPI_Comm comm);
#endif

MPI_Datatype
xt_mpi_generate_datatype_block(const int *displacements,
                               const int *blocklengths,
                               int count, MPI_Datatype old_type,
                               MPI_Comm comm) {

#ifdef COMPACT_DT
  return xt_mpi_generate_compact_datatype_block(displacements, blocklengths,
                                                count, old_type, comm);
#else
  MPI_Datatype datatype;

  if (count == 0)
    datatype = MPI_DATATYPE_NULL;
  else if (count == 1 && blocklengths[0] == 1 && displacements[0] == 0)
    datatype = copy_mpi_datatype(old_type, comm);
  else if (count == 1 && blocklengths[0] == 1)
    datatype = gen_mpi_datatype_simple(displacements[0], old_type, comm);
  else if (count == 1)
    datatype = gen_mpi_datatype_contiguous(displacements[0], blocklengths[0],
                                           old_type, comm);
  else if (check_for_vector_type(displacements, blocklengths, count))
    datatype = gen_mpi_datatype_vector(count, blocklengths[0],
                                       displacements[1] - displacements[0],
                                       displacements[0], old_type, comm);
  else if (check_for_indexed_block_type(blocklengths, count))
    datatype = gen_mpi_datatype_indexed_block(displacements, blocklengths[0],
                                              count, old_type, comm);
  else
    datatype = gen_mpi_datatype_indexed(displacements, blocklengths, count,
                                        old_type, comm);

  return datatype;
#endif
}

MPI_Datatype
xt_mpi_generate_datatype(int const * displacements, int count,
                         MPI_Datatype old_type, MPI_Comm comm)
{
  if (count <= 0)
    return MPI_DATATYPE_NULL;

#ifdef COMPACT_DT
  return xt_mpi_generate_compact_datatype(displacements, count, old_type, comm);
#else
  int * blocklengths = xmalloc((size_t)count * sizeof(*blocklengths));
  int new_count = 0;
  {
    int i = 0;
    do {
      int j = 1;
      while (i + j < count && displacements[i] + j == displacements[i + j])
        ++j;
      blocklengths[new_count++] = j;
      i += j;
    } while (i < count);
  }

  int * tmp_displ = NULL;
  const int *displ;

  if (new_count != count) {

    tmp_displ = xmalloc((size_t)new_count * sizeof(*tmp_displ));

    int offset = 0;

    for (int i = 0; i < new_count; ++i) {

      tmp_displ[i] = displacements[offset];
      offset += blocklengths[i];
    }

    displ = tmp_displ;
  } else
    displ = displacements;

  MPI_Datatype datatype;

  datatype = xt_mpi_generate_datatype_block(displ, blocklengths, new_count,
                                            old_type, comm);

  free(blocklengths);

  free(tmp_displ);

  return datatype;
#endif
}

#define XT_MPI_STRP_PRS_PREFIX
#define XT_MPI_STRP_PRS_UNITSTRIDE 1
#define XT_MPI_STRP_PRS_AOFS_TYPE int
#define XT_MPI_STRP_PRS_DISP_ADJUST(val) ((val) * params->old_type_extent)
#define XT_MPI_STRP_PRS_BLOCK_VEC_CREATE Xt_mpi_ddt_cache_acquire_vector
#define XT_MPI_STRP_PRS_INDEXED_BLOCK_CREATE \
  Xt_mpi_ddt_cache_acquire_indexed_block
#define XT_MPI_STRP_PRS_INDEXED_CREATE Xt_mpi_ddt_cache_acquire_indexed
#include "xt_mpi_stripe_parse_func.h"
#undef XT_MPI_STRP_PRS_PREFIX
#undef XT_MPI_STRP_PRS_UNITSTRIDE
#undef XT_MPI_STRP_PRS_AOFS_TYPE
#undef XT_MPI_STRP_PRS_DISP_ADJUST
#undef XT_MPI_STRP_PRS_BLOCK_VEC_CREATE
#undef XT_MPI_STRP_PRS_INDEXED_BLOCK_CREATE
#undef XT_MPI_STRP_PRS_INDEXED_CREATE

#if MPI_VERSION < 3
static inline int
XtMPI_Type_create_hindexed_block(int count, int blocklength,
                                 const MPI_Aint array_of_displacements[],
                                 MPI_Datatype oldtype, MPI_Datatype *newtype)
{
  size_t count_ = count > 0 ? (size_t)count : 0;
  MPI_Datatype *restrict oldtypes = xmalloc(count_ * sizeof (*oldtypes)
                                   + count_ * sizeof (int));
  int *restrict blocklengths = (int *)(oldtypes + count_);
  for (size_t i = 0; i < count_; ++i) {
    blocklengths[i] = blocklength;
    oldtypes[i] = oldtype;
  }
  int rc = MPI_Type_create_struct(count, blocklengths,
                                  CAST_MPI_SEND_BUF(array_of_displacements),
                                  oldtypes, newtype);
  free(oldtypes);
  return rc;
}

#define MPI_Type_create_hindexed_block XtMPI_Type_create_hindexed_block
#endif

#define XT_MPI_STRP_PRS_PREFIX a
#define XT_MPI_STRP_PRS_UNITSTRIDE params->old_type_extent
#define XT_MPI_STRP_PRS_AOFS_TYPE MPI_Aint
#define XT_MPI_STRP_PRS_DISP_ADJUST(val) (val)
#define XT_MPI_STRP_PRS_BLOCK_VEC_CREATE Xt_mpi_ddt_cache_acquire_hvector
#define XT_MPI_STRP_PRS_INDEXED_BLOCK_CREATE \
  Xt_mpi_ddt_cache_acquire_hindexed_block
#define XT_MPI_STRP_PRS_INDEXED_CREATE Xt_mpi_ddt_cache_acquire_hindexed
#include "xt_mpi_stripe_parse_func.h"


static MPI_Datatype
xt_mpi_generate_compact_datatype_block(const int *disp, const int *blocklengths,
                                       int count, MPI_Datatype old_type,
                                       MPI_Comm comm)
{
  struct Xt_mpi_strp_prs_params params;
  xt_init_mpi_strp_prs_params(&params, old_type, comm);
  MPI_Datatype dt = xt_mpi_ddt_block_gen(count, disp, blocklengths, &params);
  xt_destroy_mpi_strp_prs_params(&params);
  return dt;
}

MPI_Datatype
xt_mpi_ddt_block_gen(int count, const int *disp, const int *blocklengths,
                     struct Xt_mpi_strp_prs_params *params)
{
  size_t count_ = (size_t)0;
  for (int i=0; i<count; ++i)
    count_ += (size_t)(blocklengths[i] > 0);
  if (count_ < 1) return MPI_DATATYPE_NULL;
  struct Xt_offset_ext *restrict v = xmalloc(sizeof(*v) * count_);
  size_t j=0;
  for (size_t i=0; i<(size_t)count; ++i) {
    v[j].start = disp[i];
    v[j].stride = 1;
    int bl = blocklengths[i];
    v[j].size = bl;
    j += (size_t)(bl > 0);
  }
  MPI_Datatype dt = parse_stripe(v, count_, params);
  free(v);
  return dt;
}

MPI_Datatype
xt_mpi_generate_compact_datatype(const int *disp, int disp_len,
                                 MPI_Datatype old_type, MPI_Comm comm)
{
  if (disp_len < 1) return MPI_DATATYPE_NULL;

  size_t vlen = xt_disp2ext_count((size_t)disp_len, disp);
  struct Xt_offset_ext *v = xmalloc(sizeof(*v) * vlen);
  xt_disp2ext((size_t)disp_len, disp, v);
  struct Xt_mpi_strp_prs_params params;
  xt_init_mpi_strp_prs_params(&params, old_type, comm);
  MPI_Datatype dt = parse_stripe(v, vlen, &params);
  xt_destroy_mpi_strp_prs_params(&params);
  free(v);
  return dt;
}



/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
