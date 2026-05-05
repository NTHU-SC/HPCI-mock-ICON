/**
 * @file xt_idxvec.c
 *
 * @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
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
#include <config.h>
#endif

#include <assert.h>
#include <limits.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "xt/xt_core.h"
#include "xt/xt_idxlist.h"
#include "xt_idxlist_internal.h"
#include "xt/xt_idxempty.h"
#include "xt/xt_idxvec.h"
#include "xt_idxvec_internal.h"
#include "xt/xt_idxstripes.h"
#include "xt_idxstripes_internal.h"
#include "xt/xt_mpi.h"
#include "xt_idxlist_unpack.h"
#include "core/ppm_xfuncs.h"
#include "core/core.h"
#include "xt_stripe_util.h"
#include "xt/xt_sort.h"
#include "xt_config_internal.h"
#include "instr.h"

#define MAX(a,b) (((a)>=(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))

static void
idxvec_delete(Xt_idxlist data);

static size_t
idxvec_get_pack_size(Xt_idxlist data, MPI_Comm comm);

static void
idxvec_pack(Xt_idxlist data, void *buffer, int buffer_size,
            int *position, MPI_Comm comm);

static Xt_idxlist
idxvec_copy(Xt_idxlist idxlist);

static Xt_idxlist
idxvec_sorted_copy(Xt_idxlist idxlist, Xt_config config);

static void
idxvec_get_indices(Xt_idxlist idxlist, Xt_int *indices);

static Xt_int const*
idxvec_get_indices_const(Xt_idxlist idxlist);

static int
idxvec_get_num_index_stripes(Xt_idxlist idxlist);

static void
idxvec_get_index_stripes(Xt_idxlist idxlist, struct Xt_stripe *stripes,
                         size_t num_stripes_alloc);

static int
idxvec_get_index_at_position(Xt_idxlist idxlist, int position, Xt_int * index);

static int
idxvec_get_indices_at_positions(Xt_idxlist idxlist, const int *positions,
                                int num, Xt_int *index, Xt_int undef_idx);

static int
idxvec_get_position_of_index(Xt_idxlist idxlist, Xt_int index, int * position);

static int
idxvec_get_position_of_index_off(Xt_idxlist idxlist, Xt_int index,
                                 int * position, int offset);

static size_t
idxvec_get_positions_of_indices(Xt_idxlist idxlist, const Xt_int *indices,
                                size_t num_indices, int *positions,
                                int single_match_only);

static Xt_int
idxvec_get_min_index(Xt_idxlist idxlist);

static Xt_int
idxvec_get_max_index(Xt_idxlist idxlist);

static int
idxvec_get_sorting(Xt_idxlist idxlist);

static const struct xt_idxlist_vtable idxvec_vtable = {
  .delete                      = idxvec_delete,
  .get_pack_size               = idxvec_get_pack_size,
  .pack                        = idxvec_pack,
  .copy                        = idxvec_copy,
  .sorted_copy			= idxvec_sorted_copy,
  .get_indices                 = idxvec_get_indices,
  .get_indices_const           = idxvec_get_indices_const,
  .get_num_index_stripes	= idxvec_get_num_index_stripes,
  .get_index_stripes           = idxvec_get_index_stripes,
  .get_index_at_position       = idxvec_get_index_at_position,
  .get_indices_at_positions    = idxvec_get_indices_at_positions,
  .get_position_of_index       = idxvec_get_position_of_index,
  .get_positions_of_indices    = idxvec_get_positions_of_indices,
  .get_position_of_index_off   = idxvec_get_position_of_index_off,
  .get_positions_of_indices_off = NULL,
  .get_min_index               = idxvec_get_min_index,
  .get_max_index               = idxvec_get_max_index,
  .get_sorting                 = idxvec_get_sorting,
  .get_bounding_box            = NULL,
  .idxlist_pack_code           = VECTOR,
};

static const char filename[] = "xt_idxvec.c";


typedef struct Xt_idxvec_ *Xt_idxvec;

// index vector data structure
struct Xt_idxvec_ {

  struct Xt_idxlist_ parent;
  unsigned flags;
  Xt_int max, min;
  const Xt_int *vector;

  // internal array used to optimise access to vector data
  const Xt_int *sorted_vector;        // sorted version of vector
  int    *sorted_vec_positions; // original positions of the
                                // indices in sorted_vector
  /*
    we have the following relations:
    sorted_vector[i-1] <= sorted_vector[i],
    vector[sorted_vec_positions[i]] = sorted_vector[i]
    iff sorted_vector != vector
   */
};

static struct Xt_vec_alloc
idxvec_alloc_no_init(int num_indices)
{
  size_t vector_size = (size_t)num_indices * sizeof (Xt_int),
    header_size = ((sizeof (struct Xt_idxvec_) + sizeof (Xt_int) - 1)
                   /sizeof (Xt_int)) * sizeof (Xt_int);
  struct Xt_idxvec_ *idxvec_obj = xmalloc(header_size + vector_size);
  Xt_int *vector_assign
    = (Xt_int *)(void *)((unsigned char *)idxvec_obj + header_size);
  idxvec_obj->vector = vector_assign;
  return (struct Xt_vec_alloc) { idxvec_obj, vector_assign };
}

static struct Xt_vec_alloc
idxvec_alloc(int num_indices);

struct Xt_vec_alloc
xt_idxvec_alloc(int num_indices)
{
  return idxvec_alloc(num_indices);
}

static struct Xt_vec_alloc
idxvec_alloc(int num_indices)
{
  struct Xt_vec_alloc vec_alloc = idxvec_alloc_no_init(num_indices);
  Xt_idxlist_init(&vec_alloc.idxvec->parent, &idxvec_vtable, num_indices);
  return vec_alloc;
}

Xt_idxlist xt_idxvec_new(const Xt_int *idxvec, int num_indices) {
  INSTR_DEF(t_idxvec_new,"xt_idxvec_new")
  // ensure that yaxt is initialized
  assert(xt_initialized());

  if (num_indices > 0)
    ;
  else if (num_indices == 0)
    return xt_idxempty_new();
  else
    die("number of indices passed to xt_idxvec_new must not be negative!");

  INSTR_START(t_idxvec_new);
  struct Xt_vec_alloc vec_alloc = idxvec_alloc(num_indices);
  Xt_int *restrict vector = vec_alloc.vector;
  vector[0] = idxvec[0];
  int ntrans_up=0, ntrans_dn=0;
  Xt_int maxv = idxvec[0], minv = idxvec[0];
  for (size_t i = 1; i < (size_t)num_indices; ++i) {
    Xt_int v = idxvec[i];
    vector[i] = v;
    ntrans_up += idxvec[i-1] < v;
    ntrans_dn += idxvec[i-1] > v;
    maxv = MAX(v, maxv);
    minv = MIN(v, minv);
  }
  vec_alloc.idxvec->max = maxv;
  vec_alloc.idxvec->min = minv;
  if (ntrans_dn == 0)
    vec_alloc.idxvec->sorted_vector = vector;
  else
    vec_alloc.idxvec->sorted_vector = NULL;
  vec_alloc.idxvec->sorted_vec_positions = NULL;
  vec_alloc.idxvec->flags = XT_SORT_FLAGS(ntrans_up, ntrans_dn);
  INSTR_STOP(t_idxvec_new);
  return (void *)vec_alloc.idxvec;
}

struct flags_min_max
{
  unsigned flags;
  Xt_int max, min;
};

static struct flags_min_max
get_sort_flags(size_t num_indices, const Xt_int vector[num_indices])
{
  int ntrans_up=0, ntrans_dn=0;
  Xt_int maxv = vector[0], minv = vector[0];
  for (size_t i = 1; i < num_indices; ++i) {
    Xt_int v = vector[i];
    ntrans_up += vector[i-1] < v;
    ntrans_dn += vector[i-1] > v;
    maxv = MAX(v, maxv);
    minv = MIN(v, minv);
  }
  unsigned flags = XT_SORT_FLAGS(ntrans_up, ntrans_dn);
  /* icc 12 and 13 miscompile the above when it's a leaf function */
#if defined __ICC && defined __OPTIMIZE__      \
  && ( __INTEL_COMPILER_BUILD_DATE == 20110811 \
       || __INTEL_COMPILER == 1300 )
  if (num_indices == -1)
    fprintf(stderr, "%s: ntrans_up=%d, ntrans_dn=%d, flags=%u\n",
            filename, ntrans_up, ntrans_dn, flags);
#else
  (void)filename;
#endif
  return (struct flags_min_max){ .flags = flags, .min = minv, .max = maxv };
}

Xt_idxlist
xt_idxvec_congeal(struct Xt_vec_alloc vec_alloc)
{
  assert(vec_alloc.idxvec->vector == vec_alloc.vector);
  Xt_idxlist idxlist;
  Xt_idxvec idxvec = vec_alloc.idxvec;
  size_t num_indices = (size_t)idxvec->parent.num_indices;
  if (num_indices) {
    struct flags_min_max fmm = get_sort_flags(num_indices, vec_alloc.vector);
    unsigned flags = fmm.flags;
    idxvec->min = fmm.min;
    idxvec->max = fmm.max;
    idxvec->flags = flags;
    idxvec->sorted_vector = ((flags & sort_mask) == sort_asc
                             || (flags & sort_mask) == sort_idt)
      ? vec_alloc.vector : NULL;
    idxvec->sorted_vec_positions = NULL;
    idxlist = &vec_alloc.idxvec->parent;
  } else {
    free(idxvec);
    idxlist = xt_idxempty_new();
  }
  return idxlist;
}

Xt_idxlist xt_idxvec_prealloc_new(const Xt_int *idxvec, int num_indices)
{
  if (num_indices > 0)
    ;
  else if (num_indices == 0)
    return xt_idxempty_new();
  else
    die("number of indices passed to xt_idxvec_new must not be negative!");
  struct Xt_idxvec_ *restrict idxvec_obj = xmalloc(sizeof (*idxvec_obj));
  Xt_idxlist_init(&idxvec_obj->parent, &idxvec_vtable, num_indices);
  idxvec_obj->vector = idxvec;
  return xt_idxvec_congeal((struct Xt_vec_alloc){
      .idxvec = idxvec_obj, .vector = (void *)(intptr_t)idxvec });
}

static size_t
decode_stripe(struct Xt_stripe stripe, Xt_int *sorted_vector,
              int * sorted_vec_pos, int pos_offset) {

  Xt_int stride = stripe.stride, start;
  int sign;
  if (stride >= 0) {
    sign = 1;
    start = stripe.start;
  } else {
    sign = -1;
    start = (Xt_int)(stripe.start + (Xt_int)(stripe.nstrides-1) * stride);
    stride = (Xt_int)-stride;
    pos_offset += stripe.nstrides-1;
  }
  for (int i = 0; i < stripe.nstrides; ++i) {
    sorted_vector[i] = (Xt_int)(start + i * stride);
    sorted_vec_pos[i] = pos_offset + i * sign;
  }

  return (size_t)stripe.nstrides;
}

static void
generate_sorted_vector_from_stripes(const struct Xt_stripe stripes[],
                                    int num_stripes_,
                                    Xt_idxvec idxvec,
                                    Xt_config config)
{
  assert(num_stripes_ > 0);
  size_t num_stripes = (size_t)num_stripes_,
    num_indices = (size_t)xt_idxlist_get_num_indices(&idxvec->parent);
  Xt_int *restrict sorted_vector_assign
    = xmalloc(num_indices * sizeof(*(idxvec->sorted_vector)));
  idxvec->sorted_vector = sorted_vector_assign;
  int *sorted_vec_positions
    = idxvec->sorted_vec_positions
    = xmalloc(num_indices * sizeof(*sorted_vec_positions));

  /* stripe_minmax[0][i] is the minimal index in stripe i at first,
   * later of sorted stripe i, stripe_minmax[1][i] is the
   * corresponding maximal index */
  Xt_int (*restrict stripe_minmax)[num_stripes]
    = xmalloc(2 * sizeof(*stripe_minmax));
  int *restrict sorted_stripe_min_pos
    = xmalloc(num_stripes * 3 * sizeof(*sorted_stripe_min_pos));

  Xt_int min = XT_INT_MAX;
  for(size_t i = 0; i < num_stripes; ++i) {
    Xt_int ofs = (Xt_int)(stripes[i].stride * (stripes[i].nstrides - 1)),
      mask = Xt_isign_mask(ofs);
    Xt_int stripe_min = (Xt_int)(stripes[i].start + (ofs & mask));
    stripe_minmax[0][i] = stripe_min;
    min = MIN(min, stripe_min);
  }
  idxvec->min = min;

  xt_assign_id_map_int(num_stripes, sorted_stripe_min_pos, 0);
  config->sort_funcs->sort_xt_int_permutation(stripe_minmax[0], num_stripes,
                                              sorted_stripe_min_pos);

  int *restrict sorted_pos_prefix_sum = sorted_stripe_min_pos + num_stripes,
    *restrict orig_pos_prefix_sum
    = xmalloc(num_stripes * sizeof(*orig_pos_prefix_sum));

  int accum = 0;
  for (size_t i = 0; i < num_stripes; ++i) {
    orig_pos_prefix_sum[i] = accum;
    accum += stripes[i].nstrides;
  }

  Xt_int max = XT_INT_MIN;
  for (size_t i = 0; i < num_stripes; ++i) {
    int sorted_pos = sorted_stripe_min_pos[i];
    sorted_pos_prefix_sum[i] = orig_pos_prefix_sum[sorted_pos];
    Xt_int ofs = (Xt_int)(stripes[sorted_pos].stride
                          * (stripes[sorted_pos].nstrides - 1)),
      mask = Xt_isign_mask(ofs);
    Xt_int stripe_max = (Xt_int)(stripes[sorted_pos].start + (ofs & ~mask));
    stripe_minmax[1][i] = stripe_max;
    max = MAX(max, stripe_max);
  }
  idxvec->max = max;

  free(orig_pos_prefix_sum);

  /* i'th stripe overlaps with overlap_count[i] following stripes, or
   * is part of a stretch of this many overlapping stripes, if
   * overlap_count[i] is > 0, in case overlap_count[i] <= 0, this many
   * non-overlapping stripes follow after negation + 1 */
  int *restrict overlap_count
    = sorted_stripe_min_pos + 2 * num_stripes;
  for (size_t i = 0; i < num_stripes - 1; ++i) {
    bool do_overlap = stripe_minmax[1][i] >= stripe_minmax[0][i + 1];
    size_t j = i + 1;
    if (do_overlap) {
      /* range_max_idx is the maximal index encountered in a rage of
       * overlapping stripes, only stop when a stripe starting at
       * index larger than this is encountered */
      Xt_int range_max_idx = MAX(stripe_minmax[1][i], stripe_minmax[1][i+1]);
      while (j + 1 < num_stripes
             && stripe_minmax[0][j + 1] <= range_max_idx) {
        range_max_idx = MAX(range_max_idx, stripe_minmax[1][j+1]);
        ++j;
      }
      overlap_count[i] = (int)(j - i);
      i = j;
    } else {
      while (j + 1 < num_stripes
             && stripe_minmax[0][j + 1] > stripe_minmax[1][j])
        ++j;
      overlap_count[i] = -(int)(j - i - 1);
      i = j - 1;
    }
  }
  overlap_count[num_stripes - 1] = 0;

  size_t offset = 0;

  size_t i = 0;
  void (*sort_xt_int_permutation)(Xt_int a[], size_t n, int permutation[])
    = config->sort_funcs->sort_xt_int_permutation;
  while (i < num_stripes) {

    bool do_overlap = overlap_count[i] > 0;
    size_t num_selection = (size_t)(abs(overlap_count[i])) + 1;
    size_t sel_size = 0;

    for (size_t j = 0; j < num_selection; ++j)
      sel_size += decode_stripe(stripes[sorted_stripe_min_pos[i+j]],
                                sorted_vector_assign + offset + sel_size,
                                sorted_vec_positions + offset + sel_size,
                                sorted_pos_prefix_sum[i+j]);

    if (do_overlap)
      sort_xt_int_permutation(sorted_vector_assign + offset,
                              sel_size, sorted_vec_positions + offset);

    offset += sel_size;
    i += num_selection;
  }

  free(sorted_stripe_min_pos);
  free(stripe_minmax);
}

Xt_idxlist
xt_idxvec_from_stripes_new(const struct Xt_stripe stripes[],
                           int num_stripes) {
  // ensure that yaxt is initialized
  assert(xt_initialized());
  size_t num_stripes_ = (size_t)(num_stripes > 0 ? num_stripes : 0);
  struct Xt_stripe_summary summa = xt_summarize_stripes(num_stripes_, stripes);
  long long num_indices = summa.num_indices;
  assert((sizeof (long long) > sizeof (int)) & (num_indices <= INT_MAX)
         & (num_indices >= 0));

  Xt_idxlist idxlist;
  if (num_indices > 0) {
    struct Xt_vec_alloc vec_alloc = idxvec_alloc((int)num_indices);
    Xt_int *restrict indices = vec_alloc.vector;
    size_t k = (size_t)-1;
    for (int i = 0; i < num_stripes; ++i)
      for (int j = 0; j < stripes[i].nstrides; ++j)
        indices[++k] = (Xt_int)(stripes[i].start + j * stripes[i].stride);
    vec_alloc.idxvec->flags = summa.flags;
    generate_sorted_vector_from_stripes(stripes, num_stripes, vec_alloc.idxvec,
                                        &xt_default_config);
    idxlist = (Xt_idxlist)vec_alloc.idxvec;
  } else
    idxlist = xt_idxempty_new();
  return idxlist;
}

static void idxvec_delete(Xt_idxlist obj) {

  Xt_idxvec vec_obj = (Xt_idxvec)obj;
  if (vec_obj->sorted_vector != vec_obj->vector)
    free((void *)(intptr_t)vec_obj->sorted_vector);
  free(vec_obj->sorted_vec_positions);
  free(obj);
}

enum {
  pack_header_size = 2,
  unpack_header_size = pack_header_size-1,
};

static size_t idxvec_get_pack_size(Xt_idxlist obj, MPI_Comm comm) {

  int size_xt_idx, size_header;

  xt_mpi_call(MPI_Pack_size(pack_header_size, MPI_INT, comm, &size_header), comm);
  xt_mpi_call(MPI_Pack_size(obj->num_indices, Xt_int_dt, comm,
                            &size_xt_idx), comm);

  return (size_t)size_xt_idx + (size_t)size_header;
}

void idxvec_pack(Xt_idxlist obj, void *buffer, int buffer_size,
                 int *position, MPI_Comm comm) {

  assert(obj && xt_idxlist_get_num_indices(obj) > 0);
  Xt_idxvec idxvec = (Xt_idxvec)obj;
  int header[pack_header_size] = { VECTOR, xt_idxlist_get_num_indices(obj) };
  xt_mpi_call(MPI_Pack(header, pack_header_size, MPI_INT, buffer,
                       buffer_size, position, comm), comm);
  xt_mpi_call(MPI_Pack(CAST_MPI_SEND_BUF(idxvec->vector),
                       xt_idxlist_get_num_indices(obj),
                       Xt_int_dt, buffer,
                       buffer_size, position, comm), comm);
}

Xt_idxlist xt_idxvec_unpack(void *buffer, int buffer_size, int *position,
                            MPI_Comm comm) {

  int num_indices;
  xt_mpi_call(MPI_Unpack(buffer, buffer_size, position,
                         &num_indices, 1, MPI_INT, comm), comm);
  assert(num_indices > 0);
  Xt_idxlist idxvec = NULL;
  struct Xt_vec_alloc vec_alloc = idxvec_alloc(num_indices);
  xt_mpi_call(MPI_Unpack(buffer, buffer_size, position,
                         vec_alloc.vector, num_indices,
                         Xt_int_dt, comm), comm);
  struct flags_min_max fmm
    = get_sort_flags((size_t)num_indices, vec_alloc.vector);
  unsigned flags = fmm.flags;
  vec_alloc.idxvec->sorted_vector
    = ((flags & sort_mask) == sort_asc || (flags & sort_mask) == sort_idt)
    ? vec_alloc.vector : NULL;
  vec_alloc.idxvec->sorted_vec_positions = NULL;
  vec_alloc.idxvec->flags = flags;
  vec_alloc.idxvec->max = fmm.max;
  vec_alloc.idxvec->min = fmm.min;
  idxvec = (Xt_idxlist)vec_alloc.idxvec;
  return idxvec;
}

static const Xt_int *
get_sorted_vector(Xt_idxvec idxvec, Xt_config config);

PPM_DSO_INTERNAL const Xt_int *
xt_idxvec_get_sorted_vector(Xt_idxlist idxvec, Xt_config config)
{
  return get_sorted_vector((Xt_idxvec)idxvec, config);
}

static const Xt_int *
get_sorted_vector(Xt_idxvec idxvec, Xt_config config)
{
  if (idxvec->sorted_vector)
    return idxvec->sorted_vector;

  size_t num_indices
    = (size_t)xt_idxlist_get_num_indices(&idxvec->parent);

  const Xt_int *restrict vector = idxvec->vector;
  if ((idxvec->flags & sort_mask) == sort_asc
      || (idxvec->flags & sort_mask) == sort_idt) {
    /* we are done if vector is already sorted */
    idxvec->sorted_vec_positions = NULL;
    return idxvec->sorted_vector = vector;
  }
  if (XT_CONFIG_GET_FORCE_NOSORT(config))
    return NULL;

  size_t svec_size = num_indices * sizeof(Xt_int);
  Xt_int *sorted_vector = xmalloc(svec_size);
  idxvec->sorted_vec_positions = xmalloc(num_indices *
                                         sizeof(*(idxvec->sorted_vec_positions)));

  memcpy(sorted_vector, vector, svec_size);
/* todo: accelerate for case when vector is sorted but descending */
  xt_assign_id_map_int(num_indices, idxvec->sorted_vec_positions, 0);
  config->sort_funcs->sort_xt_int_permutation(
    sorted_vector, num_indices, idxvec->sorted_vec_positions);

  return idxvec->sorted_vector = sorted_vector;
}

Xt_idxlist
xt_idxvec_get_intersection(Xt_idxlist idxlist_src, Xt_idxlist idxlist_dst,
                           Xt_config config)
{

  // both lists are index vectors:

  Xt_idxvec idxvec_src = (Xt_idxvec)idxlist_src,
    idxvec_dst = (Xt_idxvec)idxlist_dst;


  size_t num_indices_inter = 0,
    num_indices_src = (size_t)xt_idxlist_get_num_indices(idxlist_src),
    num_indices_dst = (size_t)xt_idxlist_get_num_indices(idxlist_dst);

  struct Xt_vec_alloc vec_alloc = idxvec_alloc_no_init((int)num_indices_dst);
  Xt_idxvec inter_vector = vec_alloc.idxvec;
  Xt_int *vector_assign = vec_alloc.vector;

  struct Xt_config_ sort_config = *config;
  XT_CONFIG_UNSET_FORCE_NOSORT(&sort_config);
  // get sorted indices of source and destination
  const Xt_int *restrict sorted_src_vector
    = get_sorted_vector(idxvec_src, &sort_config),
               *restrict sorted_dst_vector
    = get_sorted_vector(idxvec_dst, &sort_config);

  // compute the intersection
  for (size_t i = 0, j = 0; i < num_indices_dst; ++i) {

    while (j < num_indices_src
           && sorted_src_vector[j] < sorted_dst_vector[i]) ++j;
    if (j < num_indices_src) {
      vector_assign[num_indices_inter] = sorted_dst_vector[i];
      num_indices_inter += sorted_src_vector[j] == sorted_dst_vector[i];
    } else
      break;
  }

  if (num_indices_inter) {
    size_t vector_size = num_indices_inter * sizeof (idxvec_dst->vector[0]),
      header_size = ((sizeof (struct Xt_idxvec_) + sizeof (Xt_int) - 1)
                     /sizeof (Xt_int)) * sizeof (Xt_int);
    inter_vector = xrealloc(inter_vector, header_size + vector_size);
    inter_vector->vector
      = (Xt_int *)(void *)((unsigned char *)inter_vector + header_size);
    Xt_idxlist_init(&inter_vector->parent, &idxvec_vtable, (int)num_indices_inter);
    inter_vector->sorted_vector = inter_vector->vector;
    inter_vector->sorted_vec_positions = NULL;
    inter_vector->flags = sort_asc;
    inter_vector->min = inter_vector->vector[0];
    inter_vector->max = inter_vector->vector[num_indices_inter-1];
  } else {
    free(inter_vector);
    inter_vector = (Xt_idxvec)xt_idxempty_new();
  }
  return (Xt_idxlist)inter_vector;
}

PPM_DSO_INTERNAL Xt_idxlist
xt_idxvec_get_idxstripes(Xt_idxlist idxlist)
{
  assert(idxlist->vtable == &idxvec_vtable);
  Xt_idxvec idxvec = (Xt_idxvec)idxlist;
  size_t num_stripes = xt_indices_count_stripes(
    (size_t)idxvec->parent.num_indices,
    idxvec->vector);
  struct Xt_stripes_alloc stripes_alloc
    = xt_idxstripes_alloc(num_stripes);
  xt_convert_indices_to_stripes_buf(
    (size_t)idxvec->parent.num_indices,
    idxvec->vector,
    num_stripes,
    stripes_alloc.stripes);
  return xt_idxstripes_congeal(stripes_alloc);
}


PPM_DSO_INTERNAL Xt_idxlist
xt_idxvec_get_idxstripes_intersection(Xt_idxlist idxlist_src,
                                      Xt_idxlist idxlist_dst,
                                      Xt_config config)
{
  /* destination is index stripes, source index vector */
  assert(idxlist_dst->vtable->idxlist_pack_code == STRIPES
         && idxlist_src->vtable == &idxvec_vtable);
  Xt_idxvec vec_src = (Xt_idxvec)idxlist_src;
  if (XT_CONFIG_GET_FORCE_NOSORT(config)) {
    Xt_idxlist idxstripes_src = xt_idxvec_get_idxstripes(idxlist_src);
    Xt_idxlist intersection
      = xt_idxstripes_get_intersection(idxstripes_src, idxlist_dst, config);
    xt_idxlist_delete(idxstripes_src);
    return intersection;
  }
  const struct Xt_stripe *stripes_dst
    = xt_idxstripes_get_index_stripes_const(idxlist_dst);
  const Xt_int *src_sorted = get_sorted_vector(vec_src, config);
  Xt_int src_min = idxvec_get_min_index(idxlist_src),
    src_max = idxvec_get_max_index(idxlist_src);
  size_t nsrc = (size_t)idxlist_src->num_indices,
    ndst = (size_t)idxlist_dst->num_indices;
  /** @todo:
   * 1. better estimate of matched indices,
   * 2. use stride sign to choose search loop
   */
  struct Xt_vec_alloc vec_alloc = idxvec_alloc_no_init((int)ndst);
  Xt_int *found_indices = vec_alloc.vector;
  size_t ndst_stripes = xt_idxstripes_get_num_index_stripes(idxlist_dst),
    num_indices_inter = 0, last_match_pos = 0;
  for (size_t i = 0; i < ndst_stripes; ++i) {
    size_t nstrides = (size_t)stripes_dst[i].nstrides;
    Xt_int start = stripes_dst[i].start,
      stride = stripes_dst[i].stride,
      end = (Xt_int)(start + (int)(nstrides-1) * stride);
    if (stride < 0) {
      Xt_int temp = start;
      start = end;
      end = temp;
      stride = (Xt_int)-stride;
    }
    if (start > src_max || end < src_min)
      continue;
    size_t adj_first = start >= src_min
      ? (size_t)0 : (size_t)((src_min - start)/stride);
    if (end > src_max) {
      size_t adj_nstrides = (size_t)((end - src_max)/stride);
      nstrides -= adj_nstrides;
    }
    for (size_t k = adj_first; k < nstrides; ++k) {
      Xt_int dst_idx = (Xt_int)(start + stride * (Xt_int)k);
      size_t lb = 0, ub = nsrc,
        guess = last_match_pos;
      do {
        Xt_int src_val = src_sorted[guess];
        if (src_val < dst_idx)
          lb = guess+1;
        else if (src_val > dst_idx)
          ub = guess;
        else {                  /* src_val == dst_idx */
          found_indices[num_indices_inter++] = dst_idx;
          last_match_pos = guess;
          break;
        }
        guess = lb + (ub-lb)/2;
      } while (lb < ub);
    }
  }
  if (num_indices_inter) {
    if (num_indices_inter != ndst) {
      size_t vector_size = num_indices_inter * sizeof (found_indices[0]),
        header_size = ((sizeof (struct Xt_idxvec_) + sizeof (Xt_int) - 1)
                       /sizeof (Xt_int)) * sizeof (Xt_int);
      vec_alloc.idxvec = xrealloc(vec_alloc.idxvec, header_size + vector_size);
      vec_alloc.idxvec->vector = vec_alloc.vector
        = (Xt_int *)(void *)((unsigned char *)vec_alloc.idxvec + header_size);
    }
  } else {
    free(vec_alloc.idxvec);
    return xt_idxempty_new();
  }
  config->sort_funcs->sort_xt_int(vec_alloc.vector, num_indices_inter);
  Xt_idxlist_init(&vec_alloc.idxvec->parent, &idxvec_vtable, (int)num_indices_inter);
  vec_alloc.idxvec->sorted_vector = vec_alloc.idxvec->vector;
  vec_alloc.idxvec->sorted_vec_positions = NULL;
  vec_alloc.idxvec->flags = sort_asc;
  return (Xt_idxlist)vec_alloc.idxvec;
}


static Xt_idxlist
idxvec_copy(Xt_idxlist idxlist) {

   return xt_idxvec_new(((Xt_idxvec)idxlist)->vector,
                        xt_idxlist_get_num_indices(idxlist));
}

static Xt_idxlist
idxvec_sorted_copy(Xt_idxlist idxlist, Xt_config config) {

  Xt_idxvec idxvec = (Xt_idxvec)idxlist;
  struct Xt_config_ sort_config = *config;
  XT_CONFIG_UNSET_FORCE_NOSORT(&sort_config);
  const Xt_int *sorted_vector = get_sorted_vector(idxvec, &sort_config);
  size_t num_indices = (size_t)xt_idxlist_get_num_indices(idxlist);
  if ((int)num_indices > config->idxv_cnv_size) {
    size_t num_stripes_alloc
      = xt_indices_count_stripes(num_indices, sorted_vector);
    struct Xt_stripes_alloc stripes_alloc
      = xt_idxstripes_alloc(num_stripes_alloc);
#ifndef NDEBUG
    size_t num_stripes =
#endif
      xt_convert_indices_to_stripes_buf(
        num_indices, sorted_vector, num_stripes_alloc, stripes_alloc.stripes);
    assert(num_stripes == num_stripes_alloc);
    return xt_idxstripes_congeal(stripes_alloc);
  } else {
    /* fixme: use congeal */
    return xt_idxvec_new(sorted_vector, (int)num_indices);
  }
}

static void
idxvec_get_indices(Xt_idxlist idxlist, Xt_int *indices) {

   memcpy(indices, ((Xt_idxvec)idxlist)->vector,
          (size_t)xt_idxlist_get_num_indices(idxlist) * sizeof(*indices));
}

static Xt_int const*
idxvec_get_indices_const(Xt_idxlist idxlist) {
  Xt_idxvec idxvec = (Xt_idxvec)idxlist;

  return idxvec->vector;
}


static int
idxvec_get_num_index_stripes(Xt_idxlist idxlist)
{
  size_t num_stripes = xt_indices_count_stripes(
    (size_t)xt_idxlist_get_num_indices(idxlist),
    ((Xt_idxvec)idxlist)->vector);
  assert(num_stripes <= INT_MAX);
  return (int)num_stripes;
}

static void
idxvec_get_index_stripes(Xt_idxlist idxlist, struct Xt_stripe *stripes,
                         size_t num_stripes_alloc) {

  xt_convert_indices_to_stripes_buf(
    (size_t)xt_idxlist_get_num_indices(idxlist),
    ((Xt_idxvec)idxlist)->vector, num_stripes_alloc, stripes);
}

static int
idxvec_get_index_at_position(Xt_idxlist idxlist, int position, Xt_int * index) {

  if (position < 0 || position >= xt_idxlist_get_num_indices(idxlist))
    return 1;

  *index = ((Xt_idxvec)idxlist)->vector[position];

  return 0;
}

static int
idxvec_get_indices_at_positions(Xt_idxlist idxlist,
                                const int *restrict positions,
                                int num_pos_, Xt_int *index,
                                Xt_int undef_idx) {

  Xt_idxvec idxvec = (Xt_idxvec)idxlist;
  size_t num_indices = (size_t)idxvec->parent.num_indices;
  const Xt_int *restrict v = idxvec->vector;

  int undef_count = 0;
  size_t num_pos = num_pos_ >= 0 ? (size_t)num_pos_ : (size_t)0;
  for (size_t ip = 0; ip < num_pos; ip++) {
    int p = positions[ip];
    if (p >= 0 && (size_t)p < num_indices) {
      index[ip] = v[p];
    } else {
      index[ip] = undef_idx;
      undef_count++;
    }
  }

  return undef_count;
}

/**
 * \todo check datatype of variables lb, ub and middle
 */
static int
idxvec_get_position_of_index_off(Xt_idxlist idxlist, Xt_int index,
                                 int * position, int offset) {

  Xt_idxvec idxvec_obj = (Xt_idxvec)idxlist;

  *position = -1;

  size_t num_indices = (size_t)idxvec_obj->parent.num_indices;
  if ((offset < 0) || ((size_t)offset >= num_indices))
    return 1;

  struct Xt_config_ sort_config = xt_default_config;
  XT_CONFIG_UNSET_FORCE_NOSORT(&sort_config);

  const Xt_int *sorted_vector
    = get_sorted_vector(idxvec_obj, &sort_config);

  if ((index < sorted_vector[0]) ||
      (index > sorted_vector[num_indices-1]))
    return 1;

  // bisection to find one matching position:
  size_t lb = 0;
  size_t ub = num_indices - 1;

  while (sorted_vector[lb] < index) {

    size_t middle = (ub + lb + 1)/2;

    if (sorted_vector[middle] <= index)
      lb = middle;
    else if (ub == middle)
      return 1;
    else
      ub = middle;
  }

  // find left most match:
  while (lb > 0 && sorted_vector[lb-1] == index) --lb;

  // go forward until offset condition is satisfied:
  const int *sorted_vec_positions = idxvec_obj->sorted_vec_positions;
  if (sorted_vec_positions) {
    while (lb < num_indices - 1            // boundary condition
           && sorted_vec_positions[lb] < offset // ignore postions left of offset
           && sorted_vector[lb] == index) ++lb;  // check if index is valid
  } else {
    while (lb < num_indices - 1            // boundary condition
           && (int)lb < offset // ignore postions left of offset
           && sorted_vector[lb] == index) ++lb;  // check if index is valid
  }
  // check if position is invalid:
  if (lb >= num_indices || sorted_vector[lb] != index)
    return 1; // failure

  // result:
  *position = sorted_vec_positions ? sorted_vec_positions[lb] : (int)lb;
  return 0;
}

static int
idxvec_get_position_of_index(Xt_idxlist idxlist, Xt_int index, int * position) {

  return idxvec_get_position_of_index_off(idxlist, index, position, 0);
}

static bool idx_vec_is_sorted(Xt_int const *idx, size_t n) {

  if (n>=2)
    for (size_t i = 1; i < n; i++)
      if (idx[i] < idx[i-1]) return false;

  return true;
}

static size_t
idxvec_get_positions_of_indices(Xt_idxlist body_idxlist,
                                const Xt_int *selection_idx,
                                size_t num_selection, int *positions,
                                int single_match_only) {

  bool selection_is_ordered = idx_vec_is_sorted(selection_idx, num_selection);
  /// \todo try linear scan of sorted data instead (requires performance test first)

  Xt_int const *sorted_selection;
  int *sorted_selection_pos = NULL;
  Xt_int *tmp_idx = NULL;

  if (selection_is_ordered) {
    sorted_selection = selection_idx;
  } else {
    size_t idx_memsize = num_selection * sizeof(*sorted_selection),
      pos_memsize = num_selection * sizeof(*sorted_selection_pos),
#if defined _CRAYC && _RELEASE_MAJOR < 9
      pos_ofs_roundup = _MAXVL_8,
#else
      pos_ofs_roundup = sizeof(int),
#endif
      /* round pos_memsize up to next multiple of sizeof (int) */
      pos_ofs = ((idx_memsize + pos_ofs_roundup - 1)
                 & ((size_t)-(ssize_t)(pos_ofs_roundup))),
      /* compute size of merged allocation */
      alloc_size = pos_ofs + pos_memsize;

    tmp_idx = xmalloc(alloc_size);
    memcpy(tmp_idx, selection_idx, idx_memsize);

    sorted_selection_pos
      = (void *)((unsigned char *)tmp_idx + pos_ofs);
    xt_assign_id_map_int(num_selection, sorted_selection_pos, 0);
    xt_default_config.sort_funcs->sort_xt_int_permutation(
      tmp_idx, num_selection, sorted_selection_pos);
    sorted_selection = tmp_idx;
  }

  /* motivation for usage of single_match_only:
   *  on the target side we want single_match_only,
   *  on the source side we don't
   */
  Xt_idxvec body_idxvec = (Xt_idxvec)body_idxlist;
  struct Xt_config_ sort_config = xt_default_config;
  XT_CONFIG_UNSET_FORCE_NOSORT(&sort_config);
  const Xt_int *sorted_body
    = get_sorted_vector(body_idxvec, &sort_config);
  const int *sorted_body_pos = body_idxvec->sorted_vec_positions;
  size_t search_end = (size_t)body_idxvec->parent.num_indices - 1;
  size_t num_unmatched = 0;

  // after the match we will move on one step in order to avoid matching the same position again
  size_t post_match_step = single_match_only != 0;

  size_t i=0;
  for (size_t search_start = 0, ub_guess_ofs = 1;
       i < num_selection && search_start<=search_end;
       ++i) {
    size_t selection_pos = selection_is_ordered ? i : (size_t)sorted_selection_pos[i];
    Xt_int isel = sorted_selection[i];
    // bisection to find one matching position:
    size_t ub = MIN(search_start + ub_guess_ofs, search_end);
    size_t lb = search_start;
     /* guess too low? */
    if (sorted_body[ub] < isel) {
      lb = MIN(ub + 1, search_end);
      ub = search_end;
    }
    /* dividing (ub-lb) by 2 gives 0 iff (ub-lb) < 2 but uses less
     * instructions than comparing to literal 1 */
    while ((ub-lb)/16) {
      size_t middle = (ub + lb + 1) / 2;
      /* todo: make branch free with mask/inv mask by predicate */
      if (sorted_body[middle] <= isel)
        lb = middle;
      else
        ub = middle;
    }
    /* use linear scan for last part of search */
    while (sorted_body[lb] < isel && lb < ub)
      ++lb;
    size_t match_pos;
    // search is now narrowed to two positions, select one of them:
    if (isel == sorted_body[lb]) {
      match_pos = lb;
    } else {
      num_unmatched++;
      positions[selection_pos] = -1;
      continue;
    }

    // find left most match >= search_start (bisection can lead to any match >= search_start)
    while (match_pos > search_start && sorted_body[match_pos-1] == isel)
      --match_pos;

    // result:
    // update positions and prepare next search:
    positions[selection_pos]
      = sorted_body_pos ? sorted_body_pos[match_pos] : (int)match_pos;
    ub_guess_ofs = match_pos - search_start;
    search_start = match_pos + post_match_step;
  }
  if (i < num_selection) {
    num_unmatched += num_selection - i;
    if (selection_is_ordered)
      do {
        positions[i] = -1;
      } while (++i < num_selection);
    else
      do {
        positions[sorted_selection_pos[i]] = -1;
      } while (++i < num_selection);
  }
  if (tmp_idx) free(tmp_idx);

  return num_unmatched;
}

static Xt_int
idxvec_get_min_index(Xt_idxlist idxlist) {

  Xt_idxvec idxvec_obj = (Xt_idxvec)idxlist;
  return idxvec_obj->min;
}

static Xt_int
idxvec_get_max_index(Xt_idxlist idxlist) {

  Xt_idxvec idxvec_obj = (Xt_idxvec)idxlist;

  return idxvec_obj->max;
}

static int
idxvec_get_sorting(Xt_idxlist idxlist)
{
  unsigned flags = ((Xt_idxvec)idxlist)->flags;
  return (int)(flags & sort_mask)-((flags & sort_mask) < 3);
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
