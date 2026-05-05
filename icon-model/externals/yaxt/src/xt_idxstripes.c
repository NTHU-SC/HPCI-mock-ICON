/**
 * @file xt_idxstripes.c
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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "xt/xt_core.h"
#include "xt_arithmetic_util.h"
#include "xt_arithmetic_long.h"
#include "xt/xt_idxlist.h"
#include "xt_idxlist_internal.h"
#include "xt/xt_idxempty.h"
#include "xt/xt_idxvec.h"
#include "xt_idxvec_internal.h"
#include "xt/xt_idxstripes.h"
#include "xt_idxstripes_internal.h"
#include "xt_stripe_util.h"
#include "xt/xt_mpi.h"
#include "xt_idxlist_unpack.h"
#include "xt_cover.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "ensure_array_size.h"
#include "instr.h"
#include "xt_config_internal.h"

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

static void
idxstripes_delete(Xt_idxlist data);

static size_t
idxstripes_get_pack_size(Xt_idxlist data, MPI_Comm comm);

static void
idxstripes_pack(Xt_idxlist data, void *buffer, int buffer_size,
                int *position, MPI_Comm comm);

static Xt_idxlist
idxstripes_copy(Xt_idxlist idxlist);

static Xt_idxlist
idxstripes_sorted_copy(Xt_idxlist idxlist, Xt_config config);

static void
idxstripes_get_indices(Xt_idxlist idxlist, Xt_int *indices);

static const Xt_int *
idxstripes_get_indices_const(Xt_idxlist idxlist);

static int
idxstripes_get_num_index_stripes(Xt_idxlist idxlist);

static void
idxstripes_get_index_stripes(Xt_idxlist idxlist,
                             struct Xt_stripe *restrict stripes,
                             size_t num_stripes_alloc);

static int
idxstripes_get_index_at_position(Xt_idxlist idxlist, int position,
                                 Xt_int * index);

static int
idxstripes_get_indices_at_positions(Xt_idxlist idxlist, const int *positions,
                                    int num, Xt_int *index,
                                    Xt_int undef_idx);
static int
idxstripes_get_pos_exts_of_index_stripes(
  Xt_idxlist idxlist,
  int num_stripes, const struct Xt_stripe stripes[num_stripes],
  int *num_ext, struct Xt_pos_ext **pos_ext,
  int single_match_only,
  Xt_config config);

static int
idxstripes_get_position_of_index(Xt_idxlist idxlist, Xt_int index,
                                 int * position);

static int
idxstripes_get_position_of_index_off(Xt_idxlist idxlist, Xt_int index,
                                     int * position, int offset);

static Xt_int
idxstripes_get_min_index(Xt_idxlist idxlist);

static Xt_int
idxstripes_get_max_index(Xt_idxlist idxlist);

static int
idxstripes_get_sorting(Xt_idxlist idxlist);

static const struct xt_idxlist_vtable idxstripes_vtable = {
  .delete                       = idxstripes_delete,
  .get_pack_size                = idxstripes_get_pack_size,
  .pack                         = idxstripes_pack,
  .copy                         = idxstripes_copy,
  .sorted_copy			= idxstripes_sorted_copy,
  .get_indices                  = idxstripes_get_indices,
  .get_indices_const            = idxstripes_get_indices_const,
  .get_num_index_stripes	= idxstripes_get_num_index_stripes,
  .get_index_stripes            = idxstripes_get_index_stripes,
  .get_index_at_position        = idxstripes_get_index_at_position,
  .get_indices_at_positions     = idxstripes_get_indices_at_positions,
  .get_position_of_index        = idxstripes_get_position_of_index,
  .get_positions_of_indices     = NULL,
  .get_pos_exts_of_index_stripes = idxstripes_get_pos_exts_of_index_stripes,
  .get_position_of_index_off    = idxstripes_get_position_of_index_off,
  .get_positions_of_indices_off = NULL,
  .get_min_index                = idxstripes_get_min_index,
  .get_max_index                = idxstripes_get_max_index,
  .get_sorting			= idxstripes_get_sorting,
  .get_bounding_box             = NULL,
  .idxlist_pack_code            = STRIPES,
};

static MPI_Datatype stripe_dt;

#if __GNUC__ >= 11
__attribute__ ((access (none, 1)))
int MPI_Get_address(XT_MPI_GET_ADDRESS_LOCATION_CONST void *location, MPI_Aint *address);
#endif

void
xt_idxstripes_initialize(void)
{
  struct Xt_stripe stripe;

  MPI_Aint base_address, nstrides_address;

  MPI_Get_address(&stripe, &base_address);
  MPI_Get_address(&stripe.nstrides, &nstrides_address);

  enum { num_stripe_dt_elems = 2 };
  static const int block_lengths[num_stripe_dt_elems] = {2,1};
  MPI_Aint displacements[num_stripe_dt_elems]
    = { 0, nstrides_address - base_address };
  static const MPI_Datatype types[num_stripe_dt_elems]
    = { Xt_int_dt, MPI_INT };
  MPI_Datatype stripe_dt_unaligned;

#ifdef XT_MPI_TYPE_CREATE_STRUCT_CONST_ARRAY_ARGS
  xt_mpi_call(MPI_Type_create_struct(num_stripe_dt_elems,
                                     block_lengths, displacements, types,
                                     &stripe_dt_unaligned), Xt_default_comm);
#else
  xt_mpi_call(MPI_Type_create_struct(num_stripe_dt_elems,
                                     (int *)block_lengths, displacements,
                                     (MPI_Datatype *)types,
                                     &stripe_dt_unaligned), Xt_default_comm);
#endif
  xt_mpi_call(MPI_Type_create_resized(stripe_dt_unaligned, 0,
                                      (MPI_Aint)sizeof(stripe),
                                      &stripe_dt), Xt_default_comm);
  xt_mpi_call(MPI_Type_free(&stripe_dt_unaligned), Xt_default_comm);
  xt_mpi_call(MPI_Type_commit(&stripe_dt), Xt_default_comm);
}

void
xt_idxstripes_finalize(void)
{
  xt_mpi_call(MPI_Type_free(&stripe_dt), Xt_default_comm);
}

typedef struct Xt_idxstripes_ *Xt_idxstripes;

struct Xt_stripes_sort {
  struct Xt_minmax range;        /* minimal and maximal position of stripe */
  int position;                  /* position of stripe this range
                                  * corresponds to (permutation
                                  * obtained from sorting) */
  int inv_position;              /* position in stripes_sort when
                                  * indexed with unsorted index */
};

/** structure of flags:
 *
 * * bits 0-1 bits encode sorting
 *
 * * bit 2 is set iff any stripes overlap each other, i.e.
 *   xt_stripes_overlap(stripes[i], stripes[j]) is true for some stripes
 *   i and j which does not necessarily imply they also share indices
 * * bit 3 is set iff at least one stripe has stride 0
 * * bit 4 is set iff any index is part of the stripes more than once
 */
enum idxstripes_flag_bits {
  /* todo: exploit newly introduced sorting information to speed up
   * index position lookups */
  stripes_sort_bit = 0,
  stripes_do_overlap_bit = 2,
  stripes_some_have_zero_stride_bit,
  stripes_intersect_bit,
};

enum idxstripes_flag_mask {
  stripes_do_overlap_mask = 1 << stripes_do_overlap_bit,
  stripes_some_have_zero_stride_mask = 1 << stripes_some_have_zero_stride_bit,
  stripes_intersect_mask = 1 << stripes_intersect_bit,
};

struct Xt_idxstripes_ {

  struct Xt_idxlist_ parent;

  struct Xt_stripe *stripes;
  Xt_int min_index, max_index;
  int num_stripes;
  unsigned flags;

  Xt_int *index_array_cache;
  struct Xt_stripes_sort stripes_sort[];
};

#define SORT_TYPE struct Xt_stripes_sort
#define SORT_TYPE_SUFFIX stripes_sort
#define XT_SORTFUNC_DECL static
#define SORT_TYPE_CMP_LT(a,b,...) ((a).range.min < (b).range.min)
#define SORT_TYPE_CMP_LE(a,b,...) ((a).range.min <= (b).range.min)
#define SORT_TYPE_CMP_EQ(a,b,...) ((a).range.min == (b).range.min)
#include "xt_mergesort_base.h"
#undef SORT_TYPE
#undef SORT_TYPE_SUFFIX
#undef SORT_TYPE_CMP_LT
#undef SORT_TYPE_CMP_LE
#undef SORT_TYPE_CMP_EQ
#undef XT_SORTFUNC_DECL

static Xt_idxlist
idxstripes_aggregate(Xt_idxstripes idxstripes,
                     const char *caller);

Xt_idxlist
xt_idxstripes_congeal(struct Xt_stripes_alloc stripes_alloc)
{
  return idxstripes_aggregate(stripes_alloc.idxstripes, __func__);
}


static Xt_idxlist
idxstripes_aggregate(Xt_idxstripes idxstripes,
                     const char *caller)
{
  const struct Xt_stripe *restrict stripes = idxstripes->stripes;
  struct Xt_stripes_sort *restrict stripes_sort = idxstripes->stripes_sort;
  int ntrans_up=0, ntrans_dn=0;
  Xt_int min, max, prev_end;
  long long num_indices;
  {
    int nstrides = stripes[0].nstrides;
    Xt_int stride = stripes[0].stride;
    prev_end = (Xt_int)(stripes[0].start
                        + stride * (nstrides-1));
    ntrans_up += ((nstrides - 1) & ~((stride > 0)-1));
    ntrans_dn += ((nstrides - 1) & ~((stride < 0)-1));
    num_indices = (long long)nstrides;
    struct Xt_minmax stripe_range = xt_stripe2minmax(stripes[0]);
    stripes_sort[0].range = stripe_range;
    stripes_sort[0].position = 0;
    min = stripe_range.min;
    max = stripe_range.max;
  }
  size_t num_stripes = (size_t)idxstripes->num_stripes;
  int sign_err = stripes[0].nstrides;
  unsigned have_zero_stride = stripes[0].stride == 0,
    stripes_intersect = have_zero_stride & (stripes[0].nstrides > 1);
  for (size_t i = 1; i < num_stripes; ++i) {
    struct Xt_minmax stripe_range = xt_stripe2minmax(stripes[i]);
    int nstrides = stripes[i].nstrides;
    Xt_int start = stripes[i].start, stride = stripes[i].stride;
    stripes_sort[i].range = stripe_range;
    stripes_sort[i].position = (int)i;
    min = MIN(stripe_range.min, min);
    max = MAX(stripe_range.max, max);
    num_indices += (long long)nstrides;
    sign_err |= nstrides;
    have_zero_stride |= stride == 0;
    stripes_intersect |= (stride == 0 && nstrides > 1);
    ntrans_up += (start > prev_end) + ((nstrides - 1) & ~((stride > 0)-1));
    ntrans_dn += (start < prev_end) + ((nstrides - 1) & ~((stride < 0)-1));
    prev_end = (Xt_int)(start + stride * (nstrides-1));
  }
  /* test sign bit */
  if (sign_err < 0) {
    static const char template[] = "ERROR: %s called with invalid stripes";
    size_t buf_size = sizeof (template) - 2 + strlen(caller);
    char *msg = xmalloc(buf_size);
    snprintf(msg, buf_size, template, caller);
    die(msg);
  }
  if (num_indices > 0) {
    assert(num_indices <= INT_MAX);
    xt_mergesort_stripes_sort(stripes_sort, num_stripes);
    stripes_sort[stripes_sort[0].position].inv_position = 0;
    unsigned stripes_do_overlap = 0;
    for (size_t i = 1; i < num_stripes; ++i) {
      stripes_do_overlap
        |= stripes_sort[i - 1].range.max >= stripes_sort[i].range.min;
      stripes_sort[stripes_sort[i].position].inv_position = (int)i;
    }
    stripes_intersect = stripes_intersect
      || (stripes_do_overlap
          && xt_stripes_detect_duplicate(
            num_stripes, stripes, (struct Xt_minmax){ .min=min, .max=max }));
    idxstripes->flags = XT_SORT_FLAGS(ntrans_up, ntrans_dn) << stripes_sort_bit
      | stripes_do_overlap << stripes_do_overlap_bit
      | have_zero_stride << stripes_some_have_zero_stride_bit
      | stripes_intersect << stripes_intersect_bit;
    idxstripes->min_index = min;
    idxstripes->max_index = max;
    idxstripes->index_array_cache = NULL;
    Xt_idxlist_init(&idxstripes->parent, &idxstripes_vtable, (int)num_indices);
    return (Xt_idxlist)idxstripes;
  } else {
    free(idxstripes);
    return xt_idxempty_new();
  }
}

static struct Xt_stripes_alloc
idxstripes_alloc(size_t num_stripes);

struct Xt_stripes_alloc
xt_idxstripes_alloc(size_t num_stripes)
{
  return idxstripes_alloc(num_stripes);
}

static struct Xt_stripes_alloc
idxstripes_alloc(size_t num_stripes)
{
  size_t header_size = ((sizeof (struct Xt_idxstripes_)
                         + (sizeof (struct Xt_stripes_sort)
                            * (size_t)num_stripes)
                         + sizeof (struct Xt_stripe) - 1)
                        / sizeof (struct Xt_stripe))
    * sizeof (struct Xt_stripe),
    body_size = sizeof (struct Xt_stripe) * (size_t)num_stripes;
  Xt_idxstripes idxstripes = xmalloc(header_size + body_size);
  struct Xt_stripe *stripes_assign
    = (struct Xt_stripe *)(void *)((unsigned char *)idxstripes + header_size);
  assert(num_stripes <= INT_MAX);
  idxstripes->num_stripes = (int)num_stripes;
  idxstripes->stripes = stripes_assign;
  return (struct Xt_stripes_alloc){ .idxstripes = idxstripes,
                                    .stripes = stripes_assign };
}


Xt_idxlist
xt_idxstripes_new(struct Xt_stripe const * stripes, int num_stripes) {
  INSTR_DEF(instr,"xt_idxstripes_new")
  INSTR_START(instr);
  // ensure that yaxt is initialized
  assert(xt_initialized());

  Xt_idxlist result;

  if (num_stripes > 0) {
    struct Xt_stripes_alloc stripes_alloc
      = idxstripes_alloc((size_t)num_stripes);
    stripes_alloc.idxstripes->num_stripes
        = (int)xt_stripes_merge_copy((size_t)num_stripes,
                                     stripes_alloc.stripes,
                                     stripes, false);
    result = idxstripes_aggregate(stripes_alloc.idxstripes, __func__);
  } else
    result = xt_idxempty_new();
  INSTR_STOP(instr);
  return result;
}

Xt_idxlist
xt_idxstripes_from_idxlist_new(Xt_idxlist idxlist_src) {
  // ensure that yaxt is initialized
  assert(xt_initialized());
  Xt_idxlist result;
  if (idxlist_src->num_indices > 0) {
    size_t num_stripes = (size_t)xt_idxlist_get_num_index_stripes(idxlist_src);
    struct Xt_stripes_alloc stripes_alloc
      = idxstripes_alloc(num_stripes);
    xt_idxlist_get_index_stripes_keep_buf(idxlist_src, stripes_alloc.stripes,
                                          num_stripes);
    result = idxstripes_aggregate(stripes_alloc.idxstripes, __func__);
  } else
    result = xt_idxempty_new();
  return result;
}



Xt_idxlist
xt_idxstripes_prealloc_new(const struct Xt_stripe *stripes, int num_stripes)
{
  Xt_idxlist result;

  if (num_stripes > 0) {
    size_t header_size = ((sizeof (struct Xt_idxstripes_)
                           + ((size_t)num_stripes
                              * sizeof (struct Xt_stripes_sort))
                           + sizeof (struct Xt_stripe) - 1)
                          / sizeof (struct Xt_stripe))
      * sizeof (struct Xt_stripe);
    Xt_idxstripes idxstripes = xmalloc(header_size);
    idxstripes->num_stripes = num_stripes;
    /* this assignment is safe because no methods modify the pointee */
    idxstripes->stripes = (struct Xt_stripe *)(intptr_t)stripes;
    result = idxstripes_aggregate(idxstripes, __func__);
  } else
    result = xt_idxempty_new();
  return result;
}


static void
idxstripes_delete(Xt_idxlist data) {

  if (data == NULL) return;

  Xt_idxstripes stripes = (Xt_idxstripes)data;

  free(stripes->index_array_cache);
  free(stripes);
}

static size_t
idxstripes_get_pack_size(Xt_idxlist data, MPI_Comm comm) {

  Xt_idxstripes stripes = (Xt_idxstripes)data;

  int size_header, size_stripes =  0;

  xt_mpi_call(MPI_Pack_size(2, MPI_INT, comm, &size_header), comm);
  if (stripes->num_stripes)
    xt_mpi_call(MPI_Pack_size(stripes->num_stripes, stripe_dt, comm,
                              &size_stripes), comm);

  return (size_t)size_header + (size_t)size_stripes;
}

static void
idxstripes_pack(Xt_idxlist data, void *buffer, int buffer_size,
                int *position, MPI_Comm comm) {
  INSTR_DEF(instr,"idxstripes_pack")
  INSTR_START(instr);

  assert(data);
  Xt_idxstripes stripes = (Xt_idxstripes)data;
  int header[2] = { STRIPES, stripes->num_stripes };

  xt_mpi_call(MPI_Pack(header, 2, MPI_INT, buffer,
                       buffer_size, position, comm), comm);
  int num_stripes = stripes->num_stripes;
  if (num_stripes)
    xt_mpi_call(MPI_Pack(stripes->stripes, num_stripes, stripe_dt,
                         buffer, buffer_size, position, comm), comm);
  INSTR_STOP(instr);
}

Xt_idxlist xt_idxstripes_unpack(void *buffer, int buffer_size, int *position,
                                MPI_Comm comm) {

  INSTR_DEF(instr,"xt_idxstripes_unpack")
  INSTR_START(instr);

  int num_stripes;
  xt_mpi_call(MPI_Unpack(buffer, buffer_size, position,
                         &num_stripes, 1, MPI_INT, comm), comm);

  Xt_idxlist result;
  assert(num_stripes);
  struct Xt_stripes_alloc stripes_alloc =
    idxstripes_alloc((size_t)num_stripes);

  xt_mpi_call(MPI_Unpack(buffer, buffer_size, position,
                         stripes_alloc.stripes,
                         num_stripes, stripe_dt, comm),comm);
  result = idxstripes_aggregate(stripes_alloc.idxstripes, __func__);
  INSTR_STOP(instr);
  return result;
}

#define XT_IDXSTRIPES_POS_EXT_MAP_COUNT
#include "xt_idxstripes_pos_ext_map.h"
#undef XT_IDXSTRIPES_POS_EXT_MAP_COUNT
#include "xt_idxstripes_pos_ext_map.h"

static Xt_idxlist
compute_intersection_fallback(Xt_idxstripes idxstripes_src,
                              Xt_idxstripes idxstripes_dst,
                              Xt_config config)
{
  INSTR_DEF(instr,"compute_intersection_fallback")
  INSTR_START(instr);
  Xt_idxlist intersection;
  int num_ext;
  struct Xt_pos_ext *pos_exts;

  int num_unmatched =
    idxstripes_get_pos_exts_of_index_stripes(
      &idxstripes_src->parent,
      idxstripes_dst->num_stripes,
      idxstripes_dst->stripes,
      &num_ext, &pos_exts, false,
      config);
  (void)num_unmatched;
  if (num_ext > 0) {
    /* 1. count index stripes correponding to pos_exts */
    size_t num_result_stripes
      = get_mapped_stripes_count((size_t)num_ext, pos_exts, idxstripes_src);
    /* 2. allocate result idxstripes */
    struct Xt_stripes_alloc intersection_alloc
      = idxstripes_alloc(num_result_stripes);
    /* 3. compute result idxstripes */
    map_ext2stripes((size_t)num_ext, pos_exts, idxstripes_src,
                    intersection_alloc.stripes);
    /* 4. congeal */
    Xt_idxlist intersection_unsorted
      = idxstripes_aggregate(intersection_alloc.idxstripes, __func__);
    if (XT_CONFIG_GET_FORCE_NOSORT(config)
        || idxstripes_get_sorting(intersection_unsorted) == 1) {
      intersection = intersection_unsorted;
    } else {
      intersection = idxstripes_sorted_copy(intersection_unsorted, config);
      xt_idxlist_delete(intersection_unsorted);
    }
    free(pos_exts);
  } else {
    intersection = xt_idxempty_new();
  }
  INSTR_STOP(instr);
  return intersection;
}


struct extended_gcd {
  Xt_int gcd, u, v;
};

/* computes egcd of two positive integers a and b such that
   egcd.gcd == egcd.u * a + egcd.v * b */
static inline struct extended_gcd extended_gcd(Xt_int a, Xt_int b) {
  Xt_int t = 1, u = 1, v = 0, s = 0;
  while (b>0)
  {
    Xt_int q = (Xt_int)(a / b);
    Xt_int prev_a = a;
    a = b;
    b = (Xt_int)(prev_a - q * b);
    Xt_int prev_u = u;
    u = s;
    s = (Xt_int)(prev_u - q * s);
    Xt_int prev_v = v;
    v = t;
    t = (Xt_int)(prev_v - q * t);
  }
  return (struct extended_gcd){ .gcd = a, .u = u, .v = v };
}

/* This implementation uses the method outlined in
 * James M. Stichnoth,
 * Efficient Compilation of Array Statements for Private Memory Multicomputers
 * February, 1993  CMU-CS-93-109
 */
static struct Xt_stripe
get_stripe_intersection(struct Xt_stripe stripe_a,
                        struct Xt_stripe stripe_b) {

  INSTR_DEF(instr,"get_stripe_intersection")
  INSTR_START(instr);

  struct Xt_bounded_stripe {
    Xt_int min, max, stride, representative;
  };

  struct Xt_bounded_stripe bsa, bsb, bsi;

  Xt_int stride_zero_mask_a = (stripe_a.stride != 0) - 1,
    stride_zero_mask_b = (stripe_b.stride != 0) - 1;
  {
    Xt_int mask = Xt_isign_mask(stripe_a.stride);
    bsa.min = (Xt_int)(stripe_a.start
                       + (mask & (stripe_a.stride * (stripe_a.nstrides - 1))));
    bsa.max = (Xt_int)(stripe_a.start
                       + (~mask & (stripe_a.stride * (stripe_a.nstrides - 1))));
  }
  bsa.representative = stripe_a.start;
  {
    Xt_int mask = Xt_isign_mask(stripe_b.stride);
    bsb.min = (Xt_int)(stripe_b.start
                       + (mask & (stripe_b.stride * (stripe_b.nstrides - 1))));
    bsb.max = (Xt_int)(stripe_b.start
                       + (~mask & (stripe_b.stride * (stripe_b.nstrides - 1))));
  }
  bsb.representative = stripe_b.start;

  bsa.stride = (Xt_int)((stripe_a.stride & ~stride_zero_mask_a)
                            | (stride_zero_mask_a & 1));
  bsb.stride = (Xt_int)((stripe_b.stride & ~stride_zero_mask_b)
                        | (stride_zero_mask_b & 1));

  /* FIXME: adjust the one further away from 0 */
  /* adjust second representative to minimize difference to first representative */
  Xt_int abs_stride_b = XT_INT_ABS(bsb.stride);
#ifdef XT_LONG
  Xt_long start_diff = (Xt_long)stripe_a.start - stripe_b.start;
  bsb.representative
    = (Xt_int)(bsb.representative
               + (start_diff/abs_stride_b
                  + (start_diff%abs_stride_b > abs_stride_b/2))
               * abs_stride_b);
  start_diff = (Xt_long)bsb.representative - bsa.representative;
#else
  Xt_long start_diff = xiisub(stripe_a.start, stripe_b.start);
  Xt_idiv start_diff_abs_stride_b_div
    = xlidivu(xlabs(start_diff), (Xt_uint)abs_stride_b);
  bsb.representative
    = (Xt_int)(bsb.representative
               + (start_diff_abs_stride_b_div.quot
                  + start_diff_abs_stride_b_div.rem > abs_stride_b/2)
               * abs_stride_b);
  start_diff = xiisub(bsb.representative, bsa.representative);
#endif
  Xt_int abs_stride_a = XT_INT_ABS(bsa.stride);
  struct extended_gcd eg = extended_gcd(abs_stride_a, abs_stride_b);
  bsi.min = MAX(bsa.min, bsb.min);
  bsi.max = MIN(bsa.max, bsb.max);
  /* we need to temporarily represent the product
   * (bsb.representative - bsa.representative) * eg.u * abs_stride_a
   * so we figure out first, if that product would overflow Xt_int
   */
  int stride_a_bits = xt_int_bits - xinlz((Xt_uint)abs_stride_a),
    stride_bits = xt_int_bits - xinlz((Xt_uint)abs_stride_b) + stride_a_bits,
    product_bits;
  if (eg.u && bsb.representative != bsa.representative)
    product_bits =
      xt_int_bits - xinlz((Xt_uint)XT_INT_ABS(eg.u))
#ifdef XT_LONG
      + xt_int_bits - xinlz((Xt_uint)xlabs(start_diff))
#else
      + xt_int_bits - xinlz(xlabs(start_diff).lo)
#endif
      + stride_a_bits;
  else
    product_bits = 0;
  int strides_mask, nstrides;
  struct Xt_stripe intersection;
  if (product_bits < xt_int_bits && stride_bits < xt_int_bits) {
    Xt_int temp_stride = (Xt_int)((abs_stride_a * bsb.stride)/eg.gcd);
    bsi.stride = (Xt_int)temp_stride;
    Xt_int min_rep;
    {
      Xt_int some_rep
        = (Xt_int)(bsa.representative
                   + ((bsb.representative - bsa.representative) * eg.u
                      * abs_stride_a / eg.gcd));
      /* compute minimal bsi representative >= bsi.min */
      Xt_int abs_bsi_stride = XT_INT_ABS(temp_stride),
        r_diff = (Xt_int)(bsi.min - some_rep),
        steps = (Xt_int)(r_diff / abs_bsi_stride);
      steps = (Xt_int)(steps + (steps * abs_bsi_stride < r_diff));
      min_rep = (Xt_int)(some_rep + steps * abs_bsi_stride);
      bsi.representative = (Xt_int)min_rep;
    }
    nstrides = (int)((bsi.max - min_rep)/temp_stride + llsign(temp_stride));
    int even_divide
      = ((((bsb.representative - bsa.representative) % eg.gcd) == 0)
         & (bsi.stride == temp_stride || abs(nstrides) == 1));
    /* requires two's complement integers */
    strides_mask = ~(((even_divide) & (bsi.min <= bsi.max)
                      & (min_rep <= bsi.max) & (min_rep >= bsi.min)) - 1);
    Xt_int max_rep
      = (Xt_int)(min_rep + (nstrides - llsign(temp_stride)) * bsi.stride);
    intersection.start = (Xt_int)((bsa.stride >= 0) ? min_rep : max_rep);
  } else {
#ifdef XT_LONG
#define ABS(x) (((x) < 0) ? -(x) : (x))
    Xt_long temp_stride = xiimul(abs_stride_a, bsb.stride)/eg.gcd;
    bsi.stride = (Xt_int)temp_stride;
    Xt_long min_rep;
    /* did the computation of bsi.stride remain in Xt_int range? */
    bool bsi_stride_in_range = xl_is_in_xt_int_range(temp_stride);
    /* did the computation of min_rep remain in Xt_int range? */
    bool min_rep_in_range;
    {
      Xt_long some_rep;
      if (bsb.representative != bsa.representative) {
        /* computation might generate huge intermediary values, needing 3
         * times as many bits as Xt_int has, so we use long multiplication */
        Xt_long temp_1 = start_diff * eg.u;
        Xt_tword temp_2 = xlimulu(temp_1, (Xt_uint)abs_stride_a);
        /* some_rep will be representable as Xt_long
         * if dividing by eg.gcd results in a value in the range of Xt_long...*/
        min_rep_in_range
          /* ... that means either all bits in temp_2.hi are
           * identical to the high bit of temp_2.midlo ... */
          = (temp_2.hi
             == (Xt_uint)(Xt_isign_mask((Xt_int)(temp_2.midlo >> xt_int_bits))))
          /* or the aggregate value of temp_2 divided by
           * XT_LONG_MAX+1 is less than eg.gcd, instead of actually
           * dividing, we can simply shift the uppermost bit of
           * temp_2.midlo into temp_2.hi left-shifted by one */
          | (XT_INT_ABS((Xt_int)(temp_2.hi << 1) | ((Xt_long)temp_2.midlo < 0)) < eg.gcd);
        /* compute any representative */
        some_rep
          = (Xt_long)bsa.representative + ((Xt_long)temp_2.midlo / eg.gcd);
      } else {
        min_rep_in_range = true;
        some_rep = bsa.representative;
      }
      /* compute minimal bsi representative >= bsi.min */
      Xt_long abs_bsi_stride = ABS(temp_stride),
        r_diff = bsi.min - some_rep,
        steps = r_diff / abs_bsi_stride;
      bool steps_in_range = xl_is_in_xt_int_range(steps);
      steps = steps + (steps * abs_bsi_stride < r_diff);
      min_rep = some_rep;
      if (steps_in_range)
        min_rep += steps * abs_bsi_stride;
      min_rep_in_range &= xl_is_in_xt_int_range(min_rep);
      bsi.representative = (Xt_int)min_rep;
    }
    int min_rep_mask = ~((int)min_rep_in_range - 1);
    nstrides = (int)((((bsi.max - min_rep)/temp_stride)+xlsign(temp_stride))
                         & min_rep_mask)
      | (~min_rep_mask & (bsb.representative == bsa.representative));
    int even_divide
      = ((start_diff%eg.gcd == 0) & (bsi_stride_in_range || abs(nstrides) == 1));
    /* requires two's complement integers */
    strides_mask = ~(((even_divide) & (bsi.min <= bsi.max)
                      & (min_rep <= bsi.max) & (min_rep >= bsi.min)) - 1);
    Xt_int max_rep
      = (Xt_int)(min_rep + (nstrides - xlsign(temp_stride)) * bsi.stride);
    intersection.start = (Xt_int)((bsa.stride >= 0) ? min_rep : max_rep);
#else
    Xt_long temp_stride = xlldiv(xiimul(abs_stride_a, bsb.stride),
                                 xi2l(eg.gcd)).quot;
    bsi.stride = (Xt_int)temp_stride.lo;
    /* minimal value that is reachable by both stripe starts plus some
     * multiple of the respective stride */
    Xt_long min_rep;
    /* did the computation of bsi.stride remain in Xt_int range? */
    bool bsi_stride_in_range = xl_is_in_xt_int_range(temp_stride);
    /* did the computation of min_rep remain in Xt_int range? */
    bool min_rep_in_range;
    {
      Xt_long some_rep;
      if (bsb.representative != bsa.representative) {
        /* computation might generate huge intermediary values */
        Xt_long temp_1
          = xiimul(bsb.representative - bsa.representative, eg.u);
        Xt_tword temp_2 = xlimul(temp_1, abs_stride_a);
        min_rep_in_range
          = (temp_2.hi
             == (Xt_uint)(Xt_isign_mask((Xt_int)(temp_2.mid))))
          | (XT_INT_ABS((Xt_int)(temp_2.hi << 1) | ((Xt_int)temp_2.mid < 0)) < eg.gcd);
        Xt_long t2 = (Xt_long){.hi = temp_2.mid, .lo = temp_2.lo };
        Xt_ldiv temp_3;
        if (eg.gcd == 1) {
          temp_3.quot = t2;
          temp_3.rem = (Xt_long){.hi = 0, .lo = 0 };
        } else if (xlicmp_lt(t2, eg.gcd)) {
          temp_3.quot = (Xt_long){ 0, 0 };
          temp_3.rem = t2;
        } else
          temp_3 = xlldiv((Xt_long){.hi = temp_2.mid, .lo = temp_2.lo },
                          xi2l(eg.gcd));
        some_rep = xliadd(temp_3.quot, bsa.representative);
      } else {
        min_rep_in_range = true;
        some_rep = xi2l(bsa.representative);
      }
      /* compute minimal bsi representative >= bsi.min */
      Xt_long abs_bsi_stride = xlabs(temp_stride),
        r_diff = xilsub(bsi.min, some_rep);
      Xt_ldiv steps = xlldiv(r_diff, abs_bsi_stride);
      bool steps_in_range = xl_is_in_xt_int_range(steps.quot);
      steps.quot = xlinc(steps.quot, (bool)((Xt_int)steps.rem.hi >= 0
                                            && steps.rem.lo != 0));
      if (steps_in_range)
        min_rep = xlladd(some_rep,
                         xiimul((Xt_int)steps.quot.lo,
                                (Xt_int)abs_bsi_stride.lo));
      else
        min_rep = some_rep;
      min_rep_in_range &= xl_is_in_xt_int_range(min_rep);
      bsi.representative = (Xt_int)min_rep.lo;
    }
    if (min_rep_in_range)
      nstrides = (int)xlldiv(xilsub(bsi.max, min_rep), temp_stride).quot.lo + xlsign(temp_stride);
    else
      nstrides = bsb.representative == bsa.representative;
    int even_divide
      = ((((bsb.representative - bsa.representative) % eg.gcd) == 0)
         & (bsi_stride_in_range || abs(nstrides) == 1));
    /* requires two's complement integers */
    strides_mask = ~(((even_divide) & (bsi.min <= bsi.max)
                      & xlicmp_le(min_rep, bsi.max) & xlicmp_ge(min_rep, bsi.min)) - 1);
    Xt_int max_rep
      = (Xt_int)(min_rep.lo
                 + (Xt_uint)((nstrides - xlsign(temp_stride)) * bsi.stride));
    intersection.start = ((bsa.stride >= 0) ? (Xt_int)min_rep.lo : max_rep);
#endif
  }
  intersection.stride
    = (Xt_int)(Xt_isign(bsa.stride) * XT_INT_ABS(bsi.stride));
  intersection.nstrides
    = (abs(nstrides) & strides_mask
       & ~((int)stride_zero_mask_a & (int)stride_zero_mask_b))
    | (stripe_b.nstrides & (int)stride_zero_mask_a & (int)stride_zero_mask_b);
  INSTR_STOP(instr);
  return intersection;
}

// this routine only works for idxstripes where !stripes_do_overlap
static Xt_idxlist
idxstripes_compute_intersection(Xt_idxstripes idxstripes_src,
                                Xt_idxstripes idxstripes_dst) {

  INSTR_DEF(instr,"idxstripes_compute_intersection")
  INSTR_START(instr);

  struct Xt_stripe *restrict inter_stripes = NULL;
  size_t num_inter_stripes = 0;
  size_t inter_stripes_array_size = 0;

  const struct Xt_stripes_sort *restrict src_stripes_sort
    = idxstripes_src->stripes_sort,
    *restrict dst_stripes_sort = idxstripes_dst->stripes_sort;
  const struct Xt_stripe *restrict stripes_src = idxstripes_src->stripes,
    *restrict stripes_dst = idxstripes_dst->stripes;

  size_t i_src = 0, i_dst = 0;
  size_t num_stripes_src = (size_t)idxstripes_src->num_stripes,
    num_stripes_dst = (size_t)idxstripes_dst->num_stripes;

  while ((i_src < num_stripes_src) &
         (i_dst < num_stripes_dst)) {

    while (i_src < num_stripes_src &&
           src_stripes_sort[i_src].range.max
           < dst_stripes_sort[i_dst].range.min) ++i_src;

    if ( i_src >= num_stripes_src ) break;

    while (i_dst < num_stripes_dst &&
           dst_stripes_sort[i_dst].range.max
           < src_stripes_sort[i_src].range.min) ++i_dst;

    if ( i_dst >= num_stripes_dst ) break;

    if ((src_stripes_sort[i_src].range.min
         <= dst_stripes_sort[i_dst].range.max)
        & (src_stripes_sort[i_src].range.max
           >= dst_stripes_sort[i_dst].range.min)) {
      ENSURE_ARRAY_SIZE(inter_stripes, inter_stripes_array_size,
                        num_inter_stripes+1);

      struct Xt_stripe intersection_stripe;
      inter_stripes[num_inter_stripes] = intersection_stripe =
        get_stripe_intersection(stripes_src[src_stripes_sort[i_src].position],
                                stripes_dst[dst_stripes_sort[i_dst].position]);
      num_inter_stripes += intersection_stripe.nstrides > 0;
    }

    if (dst_stripes_sort[i_dst].range.max
        < src_stripes_sort[i_src].range.max)
      i_dst++;
    else
      i_src++;
  }

  if (num_inter_stripes) {
    /* invert negative and merge consecutive stripes */
    struct Xt_stripe prev_stripe = inter_stripes[0];
    if (prev_stripe.stride < 0) {
      prev_stripe.start
        = (Xt_int)(prev_stripe.start
                   + (prev_stripe.stride * (Xt_int)(prev_stripe.nstrides - 1)));
      prev_stripe.stride = (Xt_int)-prev_stripe.stride;
    }
    inter_stripes[0] = prev_stripe;
    size_t j = 0;
    for (size_t i = 1; i < num_inter_stripes; ++i) {
      struct Xt_stripe stripe = inter_stripes[i];
      if (stripe.stride < 0) {
        stripe.start = (Xt_int)(stripe.start
                                + stripe.stride * (Xt_int)(stripe.nstrides - 1));
        stripe.stride = (Xt_int)-stripe.stride;
      }
      if ((stripe.stride == prev_stripe.stride)
          & (stripe.start
             == (prev_stripe.start
                 + prev_stripe.stride * (Xt_int)prev_stripe.nstrides)))
      {
        prev_stripe.nstrides += stripe.nstrides;
        inter_stripes[j].nstrides = prev_stripe.nstrides;
      }
      else
      {
        inter_stripes[++j] = stripe;
        prev_stripe = stripe;
      }
    }
    num_inter_stripes = j + 1;
  }

  Xt_idxlist inter = xt_idxstripes_new(inter_stripes, (int)num_inter_stripes);

  free(inter_stripes);

  INSTR_STOP(instr);
  return inter;
}

Xt_idxlist
xt_idxstripes_get_intersection(Xt_idxlist idxlist_src, Xt_idxlist idxlist_dst,
                               Xt_config config)
{
  // both lists are index stripes:
  Xt_idxstripes idxstripes_src = (Xt_idxstripes)idxlist_src,
    idxstripes_dst = (Xt_idxstripes)idxlist_dst;

  if ((idxstripes_src->flags | idxstripes_dst->flags)
      & stripes_do_overlap_mask) {
    return compute_intersection_fallback(idxstripes_src, idxstripes_dst,
                                         config);
  } else
    return idxstripes_compute_intersection(idxstripes_src, idxstripes_dst);
}

static inline bool
stripe_contains_index(struct Xt_stripe stripe, Xt_int idx)
{
  Xt_int start = stripe.start, stride = stripe.stride;
  if (idx >= start && stride > 0) {
    return idx < start + stride * stripe.nstrides
      && (idx - start) % stride == 0;
  } else if (idx <= start && stride < 0) {
    return idx > start + stride * stripe.nstrides
      && (start - idx) % stride == 0;
  } else {
    return idx == start;
  }
}

PPM_DSO_INTERNAL Xt_idxlist
xt_idxstripes_get_idxvec_intersection(Xt_idxlist idxlist_src,
                                      Xt_idxlist idxlist_dst,
                                      Xt_config XT_UNUSED(config))
{
  Xt_idxstripes idxstripes_src = (Xt_idxstripes)idxlist_src;
  const struct Xt_stripe *stripes_src = idxstripes_src->stripes;
  assert(idxlist_dst->vtable->idxlist_pack_code == VECTOR
         && idxlist_src->vtable == &idxstripes_vtable);
  const Xt_int *indices = xt_idxlist_get_indices_const(idxlist_dst);
  size_t ndst = (size_t)idxlist_dst->num_indices;
  size_t nsrc_stripes = (size_t)idxstripes_src->num_stripes,
    num_found_indices = 0;
  enum {
    size_t_bits = sizeof (size_t) * CHAR_BIT,
  };
  size_t *found
    = xcalloc(((ndst + size_t_bits - 1)/ size_t_bits),
              sizeof (size_t));
  for (size_t i = 0; i < ndst; ++i) {
    Xt_int dst_idx = indices[i];
    for (size_t j = 0; j < nsrc_stripes; ++j)
      if (stripe_contains_index(stripes_src[j], dst_idx)) {
        num_found_indices++;
        found[i/size_t_bits] |= (size_t)1 << (i%size_t_bits);
        break;
      }
  }
  Xt_idxlist intersection;
  if (num_found_indices) {
    struct Xt_vec_alloc vec_alloc = xt_idxvec_alloc((int)num_found_indices);
    Xt_int *restrict found_indices = vec_alloc.vector;
    for (size_t i = 0, inter_ofs = 0; i < ndst; ++i) {
      size_t fofs = i / size_t_bits;
      if (found[fofs]&((size_t)1 << (i%size_t_bits))) {
        found_indices[inter_ofs] = indices[i];
        ++inter_ofs;
      }
    }
    intersection = xt_idxvec_congeal(vec_alloc);
  } else
    intersection = xt_idxempty_new();
  free(found);
  return intersection;
}

static Xt_idxlist
idxstripes_copy(Xt_idxlist idxlist) {

  Xt_idxstripes stripes = (Xt_idxstripes)idxlist;

  return xt_idxstripes_new(stripes->stripes, stripes->num_stripes);
}

static inline Xt_int
stripe_min(struct Xt_stripe stripe, int pfx_remove)
{
  Xt_int min = stripe.start;
  Xt_int stride = stripe.stride;
  int nstrides = stripe.nstrides;
  Xt_int adj = (Xt_int)(stride >= 0 ? pfx_remove : nstrides - 1 - pfx_remove);
  return (Xt_int)(min + stride * adj);
}

static int
stripe_min_cmp_lt(const struct Xt_stripe a, const struct Xt_stripe b,
                  const int pfx_remove_a, const int pfx_remove_b)
{
  /* emptied stripes compare minimal to move them to end of array */
  if (a.nstrides - pfx_remove_a == 0)
    return true;
  if (b.nstrides - pfx_remove_b == 0)
    return false;
  Xt_int amin = stripe_min(a, pfx_remove_a),
    bmin = stripe_min(b, pfx_remove_b);
  if (amin == bmin) {
    amin = XT_INT_ABS(a.stride);
    bmin  = XT_INT_ABS(b.stride);
  }
  return amin > bmin;
}

#define SORT_TYPE int
#define SORT_TYPE_SUFFIX stripe_by_min
#define SORT_TYPE_CMP_LT(u,v,i,j) \
  (stripe_min_cmp_lt(stripes[(u)], stripes[(v)], \
                     pfx_removed[(u)], pfx_removed[(v)]))
#define XT_SORT_EXTRA_ARGS_DECL , int *restrict pfx_removed, \
    const struct Xt_stripe stripes[]
#define XT_SORT_EXTRA_ARGS_PASS , pfx_removed, stripes
#define XT_SORTFUNC_DECL static
#include "xt_heapsort_base.h"


static Xt_int
retrieve_min(size_t num_src_stripes,
             const struct Xt_stripe *src_stripes,
             int *src_permutation,
             int *pfx_removed)
{
  size_t min_stripe_pos = (size_t)src_permutation[0];
  struct Xt_stripe minstripe = src_stripes[min_stripe_pos];
  Xt_int minidx;
  if (minstripe.stride >= 0) {
    minidx = (Xt_int)(minstripe.start
                      + minstripe.stride * pfx_removed[min_stripe_pos]);
  } else {                      /* minstripe.stride < 0 */
    minidx = (Xt_int)(minstripe.start
                      + (minstripe.nstrides - 1 - pfx_removed[min_stripe_pos])
                      * minstripe.stride);
  }
  ++pfx_removed[min_stripe_pos];
  xt_heapify_stripe_by_min(src_permutation, num_src_stripes,
                           0, pfx_removed, src_stripes);
  return minidx;
}

Xt_idxlist
xt_idxstripes_sort_new(size_t num_src_stripes,
                       const struct Xt_stripe src_stripes[],
                       Xt_config config)
{
  (void)config;
  int *restrict src_permutation
    = xmalloc(num_src_stripes * 2 * sizeof (*src_permutation)),
    *restrict pfx_removed = src_permutation + num_src_stripes;
  size_t indices_left = 0;
  for (size_t i = 0; i < num_src_stripes; ++i) {
    src_permutation[i] = (int)i;
    pfx_removed[i] = 0;
    indices_left += (size_t)src_stripes[i].nstrides;
  }
  /* this only needs a heapsort actually for the following selection */
  xt_heapsort_stripe_by_min(src_permutation, num_src_stripes,
                            pfx_removed, src_stripes);
  /* select stripe with minimal element */
  struct Xt_stripe *dst_stripes
    = xmalloc(sizeof (*dst_stripes) * num_src_stripes);
  size_t num_dst_stripes = 0, size_dst_stripes = num_src_stripes;
  {
    Xt_int minidx = retrieve_min(num_src_stripes, src_stripes,
                                 src_permutation, pfx_removed);
    dst_stripes[0] = (struct Xt_stripe){ .start = minidx, .stride = 1,
                                         .nstrides = 1 };
    num_dst_stripes = 1;
    --indices_left;
  }
  for (; indices_left != 0; --indices_left) {
    /* get lowest index */
    Xt_int minidx = retrieve_min(num_src_stripes, src_stripes,
                                 src_permutation, pfx_removed);
    /* insert lowest index into destination */
    size_t i = num_dst_stripes - 1;
    if (dst_stripes[i].nstrides == 1) {
      dst_stripes[i].stride = (Xt_int)(minidx - dst_stripes[i].start);
      dst_stripes[i].nstrides = 2;
    } else if (minidx
               == (Xt_int)(
                 dst_stripes[i].start
                 + dst_stripes[i].stride * dst_stripes[i].nstrides)) {
      ++dst_stripes[i].nstrides;
    } else {
      ENSURE_ARRAY_SIZE(dst_stripes, size_dst_stripes, num_dst_stripes+1);
      dst_stripes[num_dst_stripes] = (struct Xt_stripe){ .start = minidx,
                                                         .stride = 1,
                                                         .nstrides = 1 };
      ++num_dst_stripes;
    }
  }
  /* \todo: optimize to create destination in-place */
  Xt_idxlist copy = xt_idxstripes_new(dst_stripes, (int)num_dst_stripes);
  free(dst_stripes);
  free(src_permutation);
  return copy;
}

static Xt_idxlist
idxstripes_order_invert(Xt_idxlist idxlist)
{
  Xt_idxstripes src = (Xt_idxstripes)idxlist;
  size_t num_stripes = (size_t)src->num_stripes;
  struct Xt_stripes_alloc stripes_alloc
    = idxstripes_alloc(num_stripes);
  struct Xt_stripe *restrict inverted_stripes = stripes_alloc.stripes;
  const struct Xt_stripe *src_stripes = src->stripes;
  for (size_t i = 0; i < num_stripes; ++i) {
    Xt_int src_start = src_stripes[num_stripes-1-i].start,
      src_stride = src_stripes[num_stripes-1-i].stride;
    int nstrides = src_stripes[num_stripes-1-i].nstrides;
    inverted_stripes[i] = (struct Xt_stripe){
      .start = (Xt_int)(src_start + src_stride * (nstrides-1)),
      .stride = (Xt_int)-src_stride,
      .nstrides = nstrides,
    };
  }
  return idxstripes_aggregate(stripes_alloc.idxstripes, __func__);
}

static Xt_idxlist
idxstripes_sorted_copy(Xt_idxlist idxlist, Xt_config config)
{
  Xt_idxstripes src = (Xt_idxstripes)idxlist;
  int sort_flags = (int)((src->flags >> stripes_sort_bit) & sort_mask);
  sort_flags -= sort_flags < 3;
  Xt_idxlist sorted_idxlist;
  if (sort_flags == 1)
    sorted_idxlist = idxstripes_copy(idxlist);
  else if (sort_flags == -1)
    sorted_idxlist = idxstripes_order_invert(idxlist);
  else
    sorted_idxlist = xt_idxstripes_sort_new((size_t)src->num_stripes,
                                            src->stripes, config);
  return sorted_idxlist;
}

static void
idxstripes_get_indices(Xt_idxlist idxlist, Xt_int *indices) {
  INSTR_DEF(instr,"idxstripes_get_indices")
  INSTR_START(instr);

  Xt_idxstripes idxstripes = (Xt_idxstripes)idxlist;
  const struct Xt_stripe *restrict stripes = idxstripes->stripes;
  size_t num_stripes = (size_t)idxstripes->num_stripes, ofs = 0;
  for (size_t i = 0; i < num_stripes; ofs += (size_t)stripes[i].nstrides, ++i )
    for (size_t j = 0; j < (size_t)stripes[i].nstrides; ++j)
      indices[ofs+j] = (Xt_int)(stripes[i].start
                                + (Xt_int)j * stripes[i].stride);

  INSTR_STOP(instr);
}

static Xt_int const*
idxstripes_get_indices_const(Xt_idxlist idxlist) {

  Xt_idxstripes idxstripes = (Xt_idxstripes)idxlist;

  if (idxstripes->index_array_cache) return idxstripes->index_array_cache;

  int num_indices = idxlist->num_indices;
  Xt_int *tmp_index_array
    = xmalloc((size_t)num_indices * sizeof (*tmp_index_array));
  idxstripes_get_indices(idxlist, tmp_index_array);
  return idxstripes->index_array_cache = tmp_index_array;
}

static void
idxstripes_get_index_stripes(Xt_idxlist idxlist,
                             struct Xt_stripe *restrict stripes,
                             size_t num_stripes_alloc)
{
#ifdef NDEBUG
  (void)num_stripes_alloc;
#endif
  INSTR_DEF(instr,"idxstripes_get_index_stripes")
  INSTR_START(instr);
  Xt_idxstripes idxstripes = (Xt_idxstripes)idxlist;
  assert((size_t)idxstripes->num_stripes <= num_stripes_alloc);
  memcpy(stripes, idxstripes->stripes,
         (size_t)idxstripes->num_stripes * sizeof (stripes[0]));
  INSTR_STOP(instr);
}

const struct Xt_stripe *
xt_idxstripes_get_index_stripes_const(Xt_idxlist idxlist)
{
  assert(idxlist->vtable == &idxstripes_vtable);
  Xt_idxstripes idxstripes = (Xt_idxstripes)idxlist;
  return idxstripes->stripes;
}

static int
idxstripes_get_num_index_stripes(Xt_idxlist idxlist)
{
  Xt_idxstripes idxstripes = (Xt_idxstripes)idxlist;
  return idxstripes->num_stripes;
}

size_t
xt_idxstripes_get_num_index_stripes(Xt_idxlist idxlist)
{
  assert(idxlist->vtable == &idxstripes_vtable);
  Xt_idxstripes idxstripes = (Xt_idxstripes)idxlist;
  return (size_t)idxstripes->num_stripes;
}

static int
idxstripes_get_index_at_position(Xt_idxlist idxlist, int position,
                                 Xt_int * index) {

  INSTR_DEF(instr,"idxstripes_get_index_at_position")
  INSTR_START(instr);

  int retval;

  if (position >= 0 && position < idxlist->num_indices) {
    Xt_idxstripes idxstripes = (Xt_idxstripes)idxlist;
    const struct Xt_stripe *restrict stripes = idxstripes->stripes;
    size_t num_stripes = (size_t)idxstripes->num_stripes;
    size_t i = 0;
    while (i < num_stripes && position >= stripes[i].nstrides) {
      position-= stripes[i].nstrides;
      ++i;
    }
    *index = (Xt_int)(stripes[i].start + position * stripes[i].stride);
    retval = 0;
  } else {
    retval = 1;
  }

  INSTR_STOP(instr);
  return retval;
}


static int
idxstripes_get_indices_at_positions(Xt_idxlist idxlist, const int *positions,
                                    int num_pos, Xt_int *index,
                                    Xt_int undef_idx) {

  INSTR_DEF(instr,"idxstripes_get_indices_at_positions")
  INSTR_START(instr);

  Xt_idxstripes idxstripes = (Xt_idxstripes)idxlist;
  const struct Xt_stripe *stripes = idxstripes->stripes;

  int max_pos = idxlist->num_indices;
  int sub_pos = 0;
  /* initial guess which stripe to inspect for position */
  int istripe = 0;
  /* position range corresponding to stripe indexed by istripe */
  int stripe_start_pos = 0;
  int undef_count = 0;

  for (int ipos = 0; ipos < num_pos; ipos++) {

    int seek_pos = positions[ipos];

    if (seek_pos >= 0 && seek_pos < max_pos) {
      while (seek_pos < stripe_start_pos) {
        istripe--;
#ifndef NDEBUG
        if (istripe < 0)
          die("idxstripes_get_indices_at_positions: internal error:"
              " crossed 0-boundary");
#endif
        stripe_start_pos -= (int)stripes[istripe].nstrides;
      }

      while (seek_pos > stripe_start_pos + stripes[istripe].nstrides - 1) {
        stripe_start_pos += (int)stripes[istripe].nstrides;
        istripe++;
#ifndef NDEBUG
        if (istripe >= idxstripes->num_stripes)
          die("idxstripes_get_indices_at_positions: internal error:"
              " crossed upper boundary");
#endif
      }

      sub_pos = seek_pos - stripe_start_pos;
      index[ipos]
        = (Xt_int)(stripes[istripe].start + sub_pos * stripes[istripe].stride);
    } else {
      index[ipos] = undef_idx;
      undef_count++;
    }
  }

  INSTR_STOP(instr);

  return undef_count;
}

static int
idxstripes_get_position_of_index(Xt_idxlist idxlist, Xt_int index,
                                 int * position) {

  return idxstripes_get_position_of_index_off(idxlist, index, position, 0);
}

static int
idxstripes_get_position_of_index_off(Xt_idxlist idxlist, Xt_int index,
                                     int * position, int offset) {

  INSTR_DEF(instr,"idxstripes_get_position_of_index_off")
  INSTR_START(instr);

  int retval = 1;
  size_t offset_ = (offset < 0) ? (size_t)0 : (size_t)offset;

  Xt_idxstripes idxstripes = (Xt_idxstripes)idxlist;
  const struct Xt_stripes_sort *restrict stripes_sort
    = idxstripes->stripes_sort;
  const struct Xt_stripe *restrict stripes = idxstripes->stripes;
  if (offset_ < (size_t)idxlist->num_indices) {
    size_t num_stripes = (size_t)idxstripes->num_stripes, i;
    size_t position_offset = 0;

    for (i = 0;
         position_offset + (size_t)stripes[i].nstrides <= offset_;
         ++i)
      position_offset += (size_t)stripes[i].nstrides;

    for (; i < num_stripes;
         position_offset += (size_t)stripes[i].nstrides, ++i) {

      size_t inv_pos = (size_t)stripes_sort[i].inv_position;
      struct Xt_minmax range = stripes_sort[inv_pos].range;

      if (index > range.max || index < range.min)
        continue;

      int rel_start;
      Xt_int stride = stripes[i].stride;
      if (index != stripes[i].start) {
        Xt_int index_offset = (Xt_int)(index - stripes[i].start);
        if (index_offset % stride)
          continue;
        rel_start = (int)(index_offset/stripes[i].stride);
      } else if (stride != 0)
        rel_start = 0;
      else
        rel_start = (int)(MAX(offset_, position_offset) - position_offset);
      *position = (int)((size_t)rel_start + position_offset);
      retval = 0;
      goto fun_exit;
    }
  }

  *position = -1;

 fun_exit: ;
  INSTR_STOP(instr);
  return retval;
}

static inline void
append_ext(struct Xt_pos_ext pos_ext, struct Xt_pos_ext_vec *restrict result)
{
  size_t num_pos_exts_ = result->num_pos_ext,
    size_pos_exts_ = result->size_pos_ext;
  struct Xt_pos_ext *restrict pos_exts_ = result->pos_ext;
  if (!xt_can_merge_pos_ext(pos_exts_[num_pos_exts_ - 1], pos_ext))
  {
    if (num_pos_exts_ + 1 == size_pos_exts_)
    {
      size_pos_exts_ += 16;
      /* offsetting by 1 necessary to keep the terminator in place */
      result->pos_ext = pos_exts_ = (struct Xt_pos_ext *)
        xrealloc(pos_exts_ - 1, (size_pos_exts_ + 1) * sizeof (*pos_exts_)) + 1;
      result->size_pos_ext = size_pos_exts_;
    }
    pos_exts_[num_pos_exts_] = pos_ext;
    result->num_pos_ext = num_pos_exts_ + 1;
  }
  else {
    /* merge new ext with previous */
    pos_exts_[num_pos_exts_ - 1].size
      =  isign(pos_ext.start - pos_exts_[num_pos_exts_ - 1].start)
      * (abs(pos_exts_[num_pos_exts_ - 1].size) + abs(pos_ext.size));
  }
}

struct Xt_stripes_lookup {
  size_t num_stripes;
  const struct Xt_stripe *stripes;
  const struct Xt_stripes_sort *stripes_sort;
  int *stripes_nstrides_psum;
  size_t *stripe_candidate_bitmask;
  bool stripes_do_overlap, stripes_do_intersect;
};

static inline void
create_stripes_lookup(struct Xt_stripes_lookup *restrict db,
                      Xt_idxstripes idxstripes) {
  const struct Xt_stripe *restrict stripes;
  db->stripes = stripes = idxstripes->stripes;
  size_t num_db_stripes = (size_t)idxstripes->num_stripes;
  bool stripes_do_overlap = idxstripes->flags & stripes_do_overlap_mask,
    stripes_do_intersect = idxstripes->flags & stripes_intersect_mask;
  enum {
    size_t_bits = sizeof (size_t) * CHAR_BIT,
  };
  size_t overlap_handling = stripes_do_overlap
    ? ((num_db_stripes + size_t_bits - 1)/size_t_bits+1) * sizeof (size_t) : 0,
    align_adjust = !stripes_do_overlap ? (num_db_stripes + 1) * sizeof (int):
    ((num_db_stripes + 1)*sizeof (int)+sizeof (size_t)-1)
    / sizeof (size_t) * sizeof (size_t);
  /* using num_db_stripes + 1 ensures re-aligning is always possible */
  int *restrict db_stripes_nstrides_psum
    = xmalloc(align_adjust + overlap_handling);
  if (stripes_do_overlap)
    db->stripe_candidate_bitmask
      = (void *)((unsigned char *)db_stripes_nstrides_psum+align_adjust);
  int accum = db_stripes_nstrides_psum[0] = 0;
  for (size_t j = 0; j < num_db_stripes; ++j) {
    accum += stripes[j].nstrides;
    db_stripes_nstrides_psum[j + 1] = accum;
  }
  db->stripes_sort = idxstripes->stripes_sort;
  db->num_stripes = num_db_stripes;
  db->stripes_nstrides_psum = db_stripes_nstrides_psum;
  db->stripes_do_overlap = stripes_do_overlap;
  db->stripes_do_intersect = stripes_do_intersect;
}

static inline void
destroy_stripes_lookup(struct Xt_stripes_lookup *restrict db) {
  free(db->stripes_nstrides_psum);
}

struct int_vec
{
  size_t size, num;
  union {
    int *vec;
    int e1;
  } p;
};

static inline size_t
bsearch_stripes_sort(size_t n,
                     const struct Xt_stripes_sort a[n],
                     Xt_int min_key)
{
  size_t left = 0, right = n - 1; /* avoid overflow in `(left + right)/2' */
  if ((a && n > 0)) ; else return n; /* invalid input or empty array */
  while (left < right)
  {
    /* invariant: a[left].range.min <= min_key <= a[right].range.min
     * or not in a */
    /*NOTE: *intentionally* truncate for odd sum */
    size_t m = (left + right + 1) / 2;
    if (a[m].range.min > min_key)
      right = m - 1;   /* a[left].range.min <= min_key < a[m].range.min
                        * or min_key not in a */
    else
      left = m;/* a[m].range.min <= min_key <= a[right].range.min
                * or min_key not in a */
  }
  /* assert(left == right) */
  return a[right].range.min <= min_key ? right : n;
}


static void
find_candidates(struct Xt_stripe query,
                const struct Xt_stripes_lookup *restrict db,
                struct int_vec *candidates,
                bool single_match_only,
                struct Xt_pos_ext_vec *restrict cover,
                Xt_config config)
{
  (void)cover;
  struct Xt_minmax query_minmax = xt_stripe2minmax(query);
  size_t num_db_stripes = db->num_stripes;
  const struct Xt_stripes_sort *restrict db_stripes_sort = db->stripes_sort;
  size_t start_pos;
  /* find trivial match first if sensible */
  if (db->stripes_do_overlap) {
    if ((start_pos
         = bsearch_stripes_sort(num_db_stripes, db_stripes_sort,
                                query_minmax.min))
        != num_db_stripes) {
      int stripes_start_pos = db_stripes_sort[start_pos].position;
      if (xt_stripes_eq(db->stripes[stripes_start_pos],
                         query)
          && (!single_match_only
              || (!db->stripes_do_intersect
                  && xt_cover_search(
                    cover,
                    (struct Xt_pos_range){
                      .start = db->stripes_nstrides_psum[stripes_start_pos],
                      .end = db->stripes_nstrides_psum[stripes_start_pos]
                      + query.nstrides -1 }, 0) == SIZE_MAX
                ))) {
        candidates->num = 1;
        if (candidates->size == 0)
          candidates->p.e1 = stripes_start_pos;
        else
          candidates->p.vec[0] = stripes_start_pos;
        return;
      }
    }
  }
  start_pos
    = bsearch_stripes_sort(num_db_stripes, db_stripes_sort, query_minmax.max);
  if (start_pos != num_db_stripes) {
    assert(db_stripes_sort[start_pos].range.min <= query_minmax.max);
    size_t end_pos = start_pos;
    while (end_pos < num_db_stripes
           && db_stripes_sort[end_pos].range.min <= query_minmax.max)
      ++end_pos;
    /* find all overlaps (which is more complicated if
     * non-overlapping isn't guaranteed) */
    size_t num_candidates;
    if (!db->stripes_do_overlap)
    {
      while (start_pos > 0
             && (db_stripes_sort[start_pos - 1].range.max >= query_minmax.min))
        --start_pos;
      num_candidates = end_pos - start_pos;
      if (candidates->size == 0 && num_candidates == 1) {
        candidates->p.e1 = db_stripes_sort[start_pos].position;
        candidates->num = 1;
        return;
      }
      if (candidates->size < num_candidates) {
        int *prevp = (candidates->size == 0) ? NULL : candidates->p.vec;
        candidates->p.vec = xrealloc(prevp,
                                     num_candidates
                                     * sizeof (candidates->p.vec[0]));
        candidates->size = num_candidates;
      }
      candidates->num = num_candidates;
      int *restrict candidates_vec = candidates->p.vec;
      for (size_t i = 0; i < num_candidates; ++i)
        candidates_vec[i] = db_stripes_sort[start_pos + i].position;
    }
    else
    {
      num_candidates = 0;
      size_t min_candidate = start_pos;
      const struct Xt_stripe *stripes = db->stripes;
      Xt_int query_start = query.start, query_stride = XT_INT_ABS(query.stride);
#ifdef XT_USE_FAST_DIVISIBLE_TEST
      struct xt_fast_div_coeff stride_div_test_coeff;
      stride_div_test_coeff = get_fast_div_coeff((Xt_uint)query_stride);
#endif
      size_t *bm_store = db->stripe_candidate_bitmask;
      size_t i = end_pos - 1;
      enum {
        size_t_bits = sizeof (size_t) * CHAR_BIT,
      };
      size_t pred_accum = 0;
      if (i >= size_t_bits) {
        if ((i+1) % size_t_bits) {
          for (; (i+1) % size_t_bits; --i)
          {
            size_t predicate
              = db_stripes_sort[i].range.max >= query_minmax.min;
            size_t pos = (size_t)db_stripes_sort[i].position;
            predicate
              &= (XT_INT_ABS(stripes[pos].stride) != query_stride
                  || is_divisible(query_stride, stride_div_test_coeff,
                                  (Xt_int)(query_start - stripes[pos].start)));
            num_candidates += predicate;
            size_t predicate_mask = predicate - 1;
            min_candidate = (min_candidate & predicate_mask)
              | (i & ~predicate_mask);
            pred_accum = (pred_accum << 1) | predicate;
          }
          bm_store[(i+1)/size_t_bits] = pred_accum;
          pred_accum = 0;
        }
        assert(i % size_t_bits == size_t_bits - 1);
        for (; i >= size_t_bits;) {
          for (size_t j = 0; j < size_t_bits; ++j, --i) {
            size_t predicate
              = db_stripes_sort[i].range.max >= query_minmax.min;
            size_t pos = (size_t)db_stripes_sort[i].position;
            predicate
              &= (XT_INT_ABS(stripes[pos].stride) != query_stride
                  || is_divisible(query_stride, stride_div_test_coeff,
                                  (Xt_int)(query_start - stripes[pos].start)));
            num_candidates += predicate;
            size_t predicate_mask = predicate - 1;
            min_candidate = (min_candidate & predicate_mask)
              | (i & ~predicate_mask);
            pred_accum = (pred_accum << 1) | predicate;
          }
          bm_store[(i+1)/size_t_bits] = pred_accum;
          pred_accum = 0;
        }
      }
      for (; i != SIZE_MAX; --i)
      {
        size_t predicate
          = db_stripes_sort[i].range.max >= query_minmax.min;
        size_t pos = (size_t)db_stripes_sort[i].position;
        predicate
          &= (XT_INT_ABS(stripes[pos].stride) != query_stride
              || is_divisible(query_stride, stride_div_test_coeff,
                              (Xt_int)(query_start - stripes[pos].start)));
        num_candidates += predicate;
        size_t predicate_mask = predicate - 1;
        min_candidate = (min_candidate & predicate_mask)
          | (i & ~predicate_mask);
        pred_accum = (pred_accum << 1) | predicate;
      }
      bm_store[0] = pred_accum;
      if (candidates->size == 0 && num_candidates == 1) {
        candidates->p.e1 = db_stripes_sort[min_candidate].position;
        candidates->num = num_candidates;
        return;
      }
      if (candidates->size < num_candidates+1) {
        int *prevp = (candidates->size == 0) ? NULL : candidates->p.vec;
        candidates->p.vec = xrealloc(prevp,
                                     (num_candidates + 1)
                                     * sizeof (candidates->p.vec[0]));
        candidates->size = num_candidates + 1;
      }
      candidates->num = num_candidates;
      int *restrict candidates_vec = candidates->p.vec;
      size_t j = 0,
        next_bm_mult = (min_candidate + size_t_bits - 1)/size_t_bits*size_t_bits,
        predicates = (bm_store[min_candidate/size_t_bits]
                      >> (min_candidate % size_t_bits));
      if (num_candidates == 1) {
        end_pos = min_candidate+1;
        assert(predicates&1);
      }
      for (i = min_candidate; i < end_pos; next_bm_mult+=size_t_bits) {
        size_t ulim = MIN(end_pos, next_bm_mult);
        for (; i < ulim; ++i) {
          size_t predicate = predicates&1;
          predicates >>= 1;
          if (predicate)
            candidates_vec[j] = db_stripes_sort[i].position;
          j += predicate;
        }
        predicates = bm_store[i/size_t_bits];
      }
      assert(j == num_candidates);
    }
    config->sort_funcs->sort_int(candidates->p.vec, num_candidates);
  } else
    candidates->num = 0;
}

static struct Xt_idxstripes_ *
expand_zero_stripes(size_t num_stripes,
                    const struct Xt_stripe *restrict stripes)
{
  size_t expansion = 0;
  for (size_t i = 0; i < num_stripes; ++i) {
    expansion += stripes[i].stride == 0 ? (size_t)(stripes[i].nstrides - 1) : 0;
  }
  struct Xt_stripe *restrict expanded_stripes
    = xmalloc((num_stripes + expansion) * sizeof (expanded_stripes[0]));
  size_t j = 0;
  for (size_t i = 0; i < num_stripes; ++i) {
    struct Xt_stripe stripe = stripes[i];
    if (stripe.stride == 0) {
      for (size_t k = 0; k < (size_t)stripe.nstrides; ++k)
        expanded_stripes[j + k] = (struct Xt_stripe){ .start = stripe.start,
                                                      .stride = 1,
                                                      .nstrides = 1 };
      j += (size_t)stripe.nstrides;
    } else {
      expanded_stripes[j] = stripe;
      ++j;
    }
  }
  return (struct Xt_idxstripes_ *)
    xt_idxstripes_prealloc_new(expanded_stripes,
                               (int)(num_stripes + expansion));
}


static size_t
idxstripes_get_pos_exts_of_index_stripe(
  struct Xt_stripe query,
  const struct Xt_stripes_lookup *restrict db,
  struct Xt_pos_ext_vec *restrict result,
  struct Xt_pos_ext_vec *restrict cover,
  bool single_match_only,
  size_t num_candidates,
  int *restrict candidates);

int
idxstripes_get_pos_exts_of_index_stripes(
  Xt_idxlist idxlist,
  int num_stripes,
  const struct Xt_stripe stripes[num_stripes],
  int *num_ext,
  struct Xt_pos_ext **pos_exts,
  int single_match_only,
  Xt_config config)
{
  size_t unmatched = 0;
  struct Xt_pos_ext_vec result;
  result.num_pos_ext = 0;
  struct Xt_idxstripes_ *restrict idxstripes = (struct Xt_idxstripes_ *)idxlist;
  if (num_stripes > 0)
  {
    if (idxstripes->flags & stripes_some_have_zero_stride_mask)
      idxstripes = expand_zero_stripes((size_t)idxstripes->num_stripes,
                                       idxstripes->stripes);
    result.size_pos_ext = (size_t)MAX(MIN(idxstripes->num_stripes,
                                          num_stripes), 8);
    result.pos_ext = xmalloc((result.size_pos_ext + 1)
                             * sizeof (*result.pos_ext));
    /* put non-concatenable *terminator* at offset -1 */
    result.pos_ext[0] = (struct Xt_pos_ext){.start = INT_MIN, .size = 0 };
    ++result.pos_ext;
    struct Xt_pos_ext_vec cover;
    xt_cover_start(&cover, result.size_pos_ext);
    struct Xt_stripes_lookup stripes_db;
    create_stripes_lookup(&stripes_db, idxstripes);
    struct int_vec candidates = { .size = 0, .p.vec = NULL };
    for (size_t i = 0; i < (size_t)num_stripes; ++i) {
      struct Xt_stripe query = stripes[i];
      int j = query.nstrides;
      query.nstrides = ((query.stride != 0) | (query.nstrides == 0))
        ? query.nstrides : 1;
      query.stride = (Xt_int)((query.stride != 0 && query.nstrides != 1)
                              ? query.stride : (Xt_int)1);
      find_candidates(query, &stripes_db, &candidates,
                      single_match_only, &cover, config);
      if (candidates.num) {
        do {
          unmatched += idxstripes_get_pos_exts_of_index_stripe(
            query, &stripes_db, &result, &cover, single_match_only != 0,
            candidates.num, candidates.size != 0
            ? candidates.p.vec : &candidates.p.e1);
        } while ((stripes[i].stride == 0) & (--j > 0));
      }
    }
    if (candidates.size)
      free(candidates.p.vec);
    --(result.pos_ext);
    memmove(result.pos_ext, result.pos_ext + 1,
            sizeof (*result.pos_ext) * result.num_pos_ext);
    *pos_exts = xrealloc(result.pos_ext,
                         sizeof (*result.pos_ext) * result.num_pos_ext);
    destroy_stripes_lookup(&stripes_db);
    xt_cover_finish(&cover);
    if ((struct Xt_idxstripes_ *)idxlist != idxstripes) {
      free(idxstripes->stripes);
      xt_idxlist_delete((Xt_idxlist)idxstripes);
    }
  }
  *num_ext = (int)result.num_pos_ext;
  return (int)unmatched;
}

static size_t
conditional_pos_ext_insert(struct Xt_stripe query,
                           struct Xt_pos_ext pos_ext2add,
                           const struct Xt_stripes_lookup *restrict db,
                           struct Xt_pos_ext_vec *restrict result,
                           struct Xt_pos_ext_vec *restrict cover,
                           size_t num_candidates,
                           int *restrict candidates);

struct unmatched_tail
{
  size_t unmatched;
  struct Xt_stripe query_tail;
};

static struct unmatched_tail
idxstripes_complex_get_pos_exts_of_index_stripe(
  struct Xt_stripe query,
  const struct Xt_stripes_lookup *restrict stripes_lookup,
  struct Xt_pos_ext_vec *restrict result,
  struct Xt_pos_ext_vec *restrict cover,
  bool single_match_only,
  size_t num_candidates,
  int *restrict candidates);

static size_t
pos_ext_insert(struct Xt_stripe query,
               struct Xt_pos_ext pos_ext2add,
               const struct Xt_stripes_lookup *stripes_lookup,
               struct Xt_pos_ext_vec *restrict result,
               struct Xt_pos_ext_vec *restrict cover,
               bool single_match_only,
               size_t num_candidates,
               int *restrict candidates)
{
  size_t unmatched = 0;
  if (single_match_only)
    unmatched += conditional_pos_ext_insert(
      query, pos_ext2add, stripes_lookup, result, cover,
      num_candidates, candidates);
  else
    append_ext(pos_ext2add, result);
  return unmatched;
}



size_t
idxstripes_get_pos_exts_of_index_stripe(
  struct Xt_stripe query,
  const struct Xt_stripes_lookup *restrict db,
  struct Xt_pos_ext_vec *restrict result,
  struct Xt_pos_ext_vec *restrict cover,
  bool single_match_only,
  size_t num_candidates,
  int *restrict candidates)
{
  size_t unmatched = 0;
  struct Xt_minmax query_minmax = xt_stripe2minmax(query);
  const struct Xt_stripe *restrict db_stripes = db->stripes;
  const int *restrict db_stripes_nstrides_psum = db->stripes_nstrides_psum;
  const struct Xt_stripes_sort *restrict db_stripes_sort = db->stripes_sort;
  for (size_t j = 0; j < num_candidates; ++j) {
    size_t unsort_pos = (size_t)candidates[j];
    size_t sort_pos = (size_t)db_stripes_sort[unsort_pos].inv_position;
    if ((query_minmax.min <= db_stripes_sort[sort_pos].range.max)
        & (query_minmax.max >= db_stripes_sort[sort_pos].range. min))
      ;
    else
      continue;
    struct Xt_stripe db_stripe = db_stripes[unsort_pos];
    Xt_int stride = query.stride;
    /* determine if db_stripe and query can be easily aligned */
    if ((stride > 0)
        & (stride == db_stripe.stride)
        & ((query.start - db_stripe.start) % stride == 0)) {
      /* divide query into skipped, matching and left-over parts,
       * where skipped and left-over are non-matching */
      Xt_int overlap_start = MAX(query.start, db_stripe.start);
      /* handle skipped query part */
      int skipLen = (int)((overlap_start - query.start) / stride);
      if (skipLen)
      {
        struct Xt_stripe query_head = {
          .start = query.start,
          .stride = skipLen > 1 ? stride : (Xt_int)1,
          .nstrides = skipLen };
        unmatched
          += idxstripes_get_pos_exts_of_index_stripe(
            query_head, db, result, cover, single_match_only,
            num_candidates - j - 1, candidates + j + 1);
        query.start = (Xt_int)(query.start + stride * (Xt_int)skipLen);
        query.nstrides -= skipLen;
        query.stride = query.nstrides > 1 ? query.stride : (Xt_int)1;
      }
      int db_stripe_skip
        = (int)((overlap_start - db_stripe.start) / stride);
      int overlap_nstrides
        = imin(query.nstrides, db_stripe.nstrides - db_stripe_skip);

      unmatched +=
        pos_ext_insert(
          (struct Xt_stripe){ .start = query.start,
              .stride = overlap_nstrides > 1 ? query.stride : (Xt_int)1,
              .nstrides = overlap_nstrides},
          (struct Xt_pos_ext){ .start
              = db_stripes_nstrides_psum[unsort_pos] + db_stripe_skip,
              .size = overlap_nstrides
              }, db, result, cover, single_match_only, num_candidates - j - 1,
          candidates + j + 1);

      if (!(query.nstrides -= overlap_nstrides))
        goto search_finished;
      else {
        query.start = (Xt_int)(overlap_start + stride * overlap_nstrides);
        query.stride = query.nstrides > 1 ? query.stride : (Xt_int)1;
        query_minmax = xt_stripe2minmax(query);
        continue;
      }
    }
    else if ((stride != 0)
             & (stride == db_stripe.stride)
             & ((query.start - db_stripe.start) % stride != 0)) {
      /* skip stripe that cannnot match at any position */
    } else {
      /* Handle complex overlap */
      struct unmatched_tail search_result =
        idxstripes_complex_get_pos_exts_of_index_stripe(
          query, db, result, cover, single_match_only,
          num_candidates - j, candidates + j);
      unmatched += search_result.unmatched;
      if (!search_result.query_tail.nstrides)
        goto search_finished;
      else {
        query = search_result.query_tail;
        query_minmax = xt_stripe2minmax(query);
      }
    }
  }
  unmatched += (size_t)query.nstrides;
  /* query wasn't found, add all indices to unmatched */
search_finished:
  return unmatched;
}

/* todo: add indices into pos_exts which are sorted by end/first of ranges */
static size_t
conditional_pos_ext_insert(struct Xt_stripe query,
                           struct Xt_pos_ext pos_ext2add,
                           const struct Xt_stripes_lookup *restrict db,
                           struct Xt_pos_ext_vec *restrict result,
                           struct Xt_pos_ext_vec *restrict cover,
                           size_t num_candidates,
                           int *restrict candidates)
{
  /* single_match_only is true => never re-match positions */
  size_t unmatched = 0;
  Xt_int stride = query.stride;
  if (pos_ext2add.size == -1)
    pos_ext2add.size = 1;
  int querySizeMaskNeg = isign_mask(pos_ext2add.size);
tail_search:
  ;
  int pos_ext2add_s = pos_ext2add.start
    + (querySizeMaskNeg & (pos_ext2add.size + 1)),
    pos_ext2add_e = pos_ext2add.start
    + (~querySizeMaskNeg & (pos_ext2add.size - 1));
  struct Xt_pos_range query_range
    = { .start = pos_ext2add_s, .end = pos_ext2add_e };
  /* does overlap exist? */
  size_t overlap_pos =
    xt_cover_insert_or_overlap(cover, query_range, 0);
  if (overlap_pos == SIZE_MAX) {
    /* remaining extent does not overlap any existing one */
    append_ext(pos_ext2add, result);
  } else {
    struct Xt_pos_ext *restrict pos_exts_ = cover->pos_ext;
    int dbSizeMaskNeg = isign_mask(pos_exts_[overlap_pos].size),
      db_s = pos_exts_[overlap_pos].start
      + (dbSizeMaskNeg & (pos_exts_[overlap_pos].size + 1)),
      db_e = pos_exts_[overlap_pos].start
      + (~dbSizeMaskNeg & (pos_exts_[overlap_pos].size - 1));
    /* determine length of overlap parts */
    int lowQuerySkip = db_s - pos_ext2add_s;
    int lowDbSkip = -lowQuerySkip;
    lowQuerySkip = (int)((unsigned)(lowQuerySkip + abs(lowQuerySkip))/2);
    lowDbSkip = (int)((unsigned)(lowDbSkip + abs(lowDbSkip))/2);
    int overlapLen = MIN(db_e - db_s - lowDbSkip + 1,
                         abs(pos_ext2add.size) - lowQuerySkip);
    int highQuerySkip = abs(pos_ext2add.size) - lowQuerySkip - overlapLen;
    /* then adjust lengths to direction of overlap (from
     * perspective of pos_ext2add */
    int querySkipLen = (~querySizeMaskNeg & lowQuerySkip)
      | (querySizeMaskNeg & -highQuerySkip),
      queryTailLen = (querySizeMaskNeg & -lowQuerySkip)
      | (~querySizeMaskNeg & highQuerySkip);
    if (querySkipLen)
    {
      int absQuerySkipLen = abs(querySkipLen);
      struct Xt_stripe query_skip = {
        .start = query.start,
        .stride = absQuerySkipLen > 1 ? stride : (Xt_int)1,
        .nstrides = absQuerySkipLen,
      };
      struct Xt_pos_ext pos_ext2add_skip = {
        .start = pos_ext2add.start,
        .size = querySkipLen
      };
      unmatched
        += conditional_pos_ext_insert(
          query_skip, pos_ext2add_skip, db,
          result, cover, num_candidates, candidates);
      pos_exts_ = result->pos_ext;
      query.start = (Xt_int)(query.start
                             + stride * (Xt_int)absQuerySkipLen);
      query.nstrides -= absQuerySkipLen;
      query.stride = query.nstrides > 1 ? query.stride : (Xt_int)1;
      pos_ext2add.start += querySkipLen;
      pos_ext2add.size -= querySkipLen;
    }
    /* head part of (remaining) query matches part of already inserted index
     * range */
    struct Xt_stripe query_head = {
      .start = query.start,
      .stride = abs(overlapLen) > 1 ? stride : (Xt_int)1,
      .nstrides = abs(overlapLen) };
    /* restart search for overlapping part within following ranges */
    unmatched
      += idxstripes_get_pos_exts_of_index_stripe(
        query_head, db, result, cover, true,
        num_candidates, candidates);
    pos_exts_ = result->pos_ext;
    if (queryTailLen) {
      /* shorten query accordingly */
      query.nstrides -= abs(overlapLen);
      query.start = (Xt_int)(query.start + stride * (Xt_int)abs(overlapLen));
      query.stride = query.nstrides > 1 ? query.stride : (Xt_int)1;
      int directedOverlapLen = (~querySizeMaskNeg & overlapLen)
        | (querySizeMaskNeg & -overlapLen);
      pos_ext2add.start += directedOverlapLen;
      pos_ext2add.size -= directedOverlapLen;
      goto tail_search;
    }
    /* whole range handled, return */
  }
  return unmatched;
}

static Xt_int
idxstripes_get_min_index(Xt_idxlist idxlist) {
  Xt_idxstripes idxstripes = (Xt_idxstripes)idxlist;
  return idxstripes->min_index;
}

static Xt_int
idxstripes_get_max_index(Xt_idxlist idxlist) {
  Xt_idxstripes idxstripes = (Xt_idxstripes)idxlist;
  return idxstripes->max_index;
}

static int
idxstripes_get_sorting(Xt_idxlist idxlist)
{
  unsigned sort_flags = (((Xt_idxstripes)idxlist)->flags >> stripes_sort_bit) & sort_mask;
  return (int)sort_flags-(sort_flags < 3);
}


/* unfortunately, Cray cc 8.[56].x miscompile the following function,
 * hence optimizations must be disabled. Feel free to dig deeper into
 * the problem, but it's fixed in 8.7. */
#if defined _CRAYC && _RELEASE_MAJOR == 8 && _RELEASE_MINOR >= 5 && _RELEASE_MINOR < 7
#pragma _CRI noopt
#endif
static struct unmatched_tail
idxstripes_complex_get_pos_exts_of_index_stripe(
  struct Xt_stripe query,
  const struct Xt_stripes_lookup *restrict db,
  struct Xt_pos_ext_vec *restrict result,
  struct Xt_pos_ext_vec *restrict cover,
  bool single_match_only,
  size_t num_candidates,
  int *restrict candidates)
{
  size_t unmatched = 0;
  const struct Xt_stripe *restrict db_stripes = db->stripes;
  size_t db_stripe_pos = (size_t)candidates[0];
  struct Xt_stripe overlap = get_stripe_intersection(query,
                                                     db_stripes[db_stripe_pos]);
  if (overlap.nstrides == 0)
    return (struct unmatched_tail){ .unmatched = 0, .query_tail = query};

  int skipped = (int)((overlap.start - query.start)/query.stride);
  if (skipped)
  {
    unmatched += idxstripes_get_pos_exts_of_index_stripe(
      (struct Xt_stripe){ .start = query.start,
          .stride = skipped > 1 ? query.stride : (Xt_int)1,
          .nstrides = skipped}, db, result, cover,
      single_match_only, num_candidates - 1, candidates + 1);
    query.start = (Xt_int)(query.start + skipped * query.stride);
    query.nstrides -= skipped;
    query.stride = query.nstrides > 1 ? query.stride : 1;
  }
  /* Since overlap.nstrides > 0, overlap.start always matches, but
   * depending on the stride the remaining overlapping indices might or might
   * not be consecutive in the query, in the latter case intervening
   * parts need to be searched for too.
   * Stripes of length 1 are naturally consecutive.
   */
  int db_stripe_skip
    = (int)((overlap.start - db_stripes[db_stripe_pos].start)
            / db_stripes[db_stripe_pos].stride);
  int db_pos = db->stripes_nstrides_psum[db_stripe_pos] + db_stripe_skip;
  if (((overlap.stride == query.stride)
       & (overlap.stride == db_stripes[db_stripe_pos].stride))
      | (overlap.nstrides == 1))
  {
    unmatched += pos_ext_insert(overlap, (struct Xt_pos_ext){
        .start = db_pos,
          .size = overlap.nstrides },
      db, result, cover, single_match_only, num_candidates - 1, candidates + 1);
    query.nstrides -= overlap.nstrides;
    query.start = (Xt_int)(query.start + overlap.nstrides * query.stride);
    query.stride = query.nstrides > 1 ? query.stride : (Xt_int)1;
  }
  else if ((overlap.stride == query.stride)
           & (overlap.stride == -db_stripes[db_stripe_pos].stride))
  {
    /* all parts of the overlap can be used directly,
       but are inversely sorted in db_stripe */
    unmatched += pos_ext_insert(overlap, (struct Xt_pos_ext){
        .start = db_pos, .size = -overlap.nstrides },
      db, result, cover, single_match_only,
      num_candidates - 1, candidates + 1);
    query.nstrides -= overlap.nstrides;
    query.start = (Xt_int)(query.start + overlap.nstrides * query.stride);
    query.stride = query.nstrides > 1 ? query.stride : (Xt_int)1;
  }
  else if (overlap.stride == query.stride)
  {
    /* all parts of the overlap can be used but are non-consecutive
     * in db_stripe */
    int db_step = (int)(overlap.stride/db_stripes[db_stripe_pos].stride);
    /* todo: try to keep (prefix of) stripe together if the
     * corresponding positions cannot be inserted anyway */
    for (int i = 0; i < overlap.nstrides; ++i, db_pos += db_step)
      unmatched += pos_ext_insert(
        (struct Xt_stripe){ (Xt_int)(overlap.start + i*overlap.stride), 1, 1 },
        (struct Xt_pos_ext){ .start = db_pos, .size = 1 },
        db, result, cover, single_match_only,
        num_candidates - 1, candidates + 1);
    query.nstrides -= overlap.nstrides;
    query.start = (Xt_int)(query.start + overlap.nstrides * query.stride);
    query.stride = query.nstrides > 1 ? query.stride : (Xt_int)1;
  }
  else /* overlap.stride != query.stride => overlap.stride > query.stride */
  {
    /*
     * query.start = (Xt_int)(query.start
     *                        + overlap.stride * overlap.nstrides);
     */
    int stride_step = (int)(overlap.stride / query.stride);
    int db_step = (int)(overlap.stride/db_stripes[db_stripe_pos].stride);
    do {
      struct Xt_stripe consecutive_overlap = { .start = query.start,
                                               .stride = 1,
                                               .nstrides = 1 },
        intervening = { .start = (Xt_int)(query.start + query.stride),
                        .stride = query.stride,
                        .nstrides = MIN(query.nstrides - 1, stride_step - 1) };
      /* split off start index, then handle intervening */
      unmatched +=
        pos_ext_insert(consecutive_overlap, (struct Xt_pos_ext){
            .start = db_pos, .size = 1 },
          db, result, cover, single_match_only,
          num_candidates - 1, candidates + 1);
      db_pos += db_step;
      if (--query.nstrides > 0) {
        unmatched +=
          idxstripes_get_pos_exts_of_index_stripe(
            intervening, db, result, cover, single_match_only,
            num_candidates - 1, candidates + 1);
        query.nstrides -= intervening.nstrides;
      }
      query.start = (Xt_int)(query.start + query.stride * stride_step);
      query.stride = query.nstrides > 1 ? query.stride : (Xt_int)1;
      overlap.start = (Xt_int)(overlap.start + overlap.stride);
    } while (--overlap.nstrides);
  }
  return (struct unmatched_tail){ .unmatched = unmatched, .query_tail = query};
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
