/**
 * @file xt_idxsection_stripe_intersection.h
 * Contains code parts independent of whether section
 * is destination or source side of intersection computation.
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
#define XT_TOKEN_PASTE2_(a,b) a##b
#ifdef NUM_DIMENSIONS
#  ifdef XT_IDXSECTION_STRIPES_ISECT_SINGLE_MATCH_ONLY
#    define xt_idxsection_get_idxstripes_intersection__(ndims)       \
  XT_TOKEN_PASTE2_(xt_idxsection_get_idxstripes_intersection_sm_, ndims)
#  else
#    define xt_idxsection_get_idxstripes_intersection__(ndims)       \
  XT_TOKEN_PASTE2_(xt_idxsection_get_idxstripes_intersection_, ndims)
#  endif
#  define NUM_DIMENSIONS_ NUM_DIMENSIONS
#else
#  ifdef XT_IDXSECTION_STRIPES_ISECT_SINGLE_MATCH_ONLY
#    define xt_idxsection_get_idxstripes_intersection__(ndims)  \
  xt_idxsection_get_idxstripes_intersection_sm_
#  else
#    define xt_idxsection_get_idxstripes_intersection__(ndims)    \
  xt_idxsection_get_idxstripes_intersection_
#  endif
#  define NUM_DIMENSIONS_ num_dimensions
#endif
static Xt_idxlist
xt_idxsection_get_idxstripes_intersection__(NUM_DIMENSIONS)
  (Xt_idxsection idxsection, Xt_idxlist idxstripes_list, Xt_config config)
{
  // intersection between an idxsection and a set of stripes
  size_t num_stripes
    = xt_idxstripes_get_num_index_stripes(idxstripes_list);

  const struct Xt_stripe *stripes
    = xt_idxstripes_get_index_stripes_const(idxstripes_list);
#ifdef NUM_DIMENSIONS
#define Xt_idxsection__(ndims) XT_TOKEN_PASTE2_(Xt_idxsection_,ndims)
  struct Xt_idxsection__(NUM_DIMENSIONS) {

    struct Xt_idxlist_ parent;

    Xt_int *index_array_cache;

    Xt_int global_start_index;
    Xt_int local_start_index;
    Xt_int min_index_cache;
    Xt_int max_index_cache;
    int ndim;
    unsigned flags;
    struct dim_desc dims[NUM_DIMENSIONS];
  };
  struct Xt_idxsection__(NUM_DIMENSIONS) h;
  idxsection_init_sorted_copy(idxsection, (Xt_idxsection)&h);
  struct dim_desc *restrict dims = h.dims;
#else
  size_t num_dimensions = (size_t)idxsection->ndim;
  Xt_idxsection sssp
    = xmalloc(sizeof (struct Xt_idxsection_)
              + sizeof (struct dim_desc) * num_dimensions);
  idxsection_init_sorted_copy(idxsection, sssp);
  struct dim_desc *dims = sssp->dims;
#endif
enum {
  size_t_bits = sizeof (size_t) * CHAR_BIT,
};
#ifndef XT_IDXSECTION_STRIPES_ISECT_SINGLE_MATCH_ONLY
size_t *restrict isect_stripes_bm
  = xcalloc((num_stripes+size_t_bits-1)/size_t_bits, sizeof (size_t));
#else
size_t num_indices_in_section
  = (size_t)xt_idxlist_get_num_indices(&idxsection->parent),
  section_bm_size = (num_indices_in_section+size_t_bits-1)/size_t_bits,
  *restrict pos_used_bm = xcalloc(section_bm_size
                         + (num_stripes+size_t_bits-1)/size_t_bits,
                         sizeof (size_t)),
  *restrict isect_stripes_bm = pos_used_bm + section_bm_size;
#endif
size_t num_result_stripes_total = 0;
{
#define XT_IDXSECTION_STRIPES_ISECT_CREATE_STRIPE_MASK
  /* start with non-appendable terminator */
  Xt_int accum_start = XT_INT_MIN;
  int buf_nstrides = 0;
  Xt_int buf_stride = -1;
  Xt_int index_continuation = XT_INT_MIN-1;
  size_t isect_mask = 0;
#include "xt_idxsection_stripe_iterate.h"
  if (num_stripes % size_t_bits)
    isect_stripes_bm[num_stripes/size_t_bits] = isect_mask;
#undef XT_IDXSECTION_STRIPES_ISECT_CREATE_STRIPE_MASK
}
Xt_idxlist result;
if (num_result_stripes_total) {
#ifdef XT_IDXSECTION_STRIPES_ISECT_SINGLE_MATCH_ONLY
  memset(pos_used_bm, 0, (num_indices_in_section+size_t_bits-1)/size_t_bits *
         sizeof (size_t));
#endif
  struct Xt_stripes_alloc result_alloc
    = xt_idxstripes_alloc(num_result_stripes_total);
  struct Xt_stripe *restrict result_stripes = result_alloc.stripes;
  size_t num_result_stripes = (size_t)-1;
  Xt_int index_continuation = XT_INT_MIN-1, accum_start = XT_INT_MIN;
  int buf_nstrides = 0;
  Xt_int buf_stride = -1;
#include "xt_idxsection_stripe_iterate.h"
  result_stripes[num_result_stripes].stride = buf_stride;
  result_stripes[num_result_stripes].nstrides = buf_nstrides;
  assert(num_result_stripes_total == num_result_stripes+1);
  Xt_idxlist initial_result = xt_idxstripes_congeal(result_alloc);
  if (XT_CONFIG_GET_FORCE_NOSORT(config)
      || xt_idxlist_get_sorting(initial_result) == 1)
    result = initial_result;
  else {
    result = xt_idxlist_sorted_copy(initial_result);
    xt_idxlist_delete(initial_result);
  }
} else
  result = xt_idxempty_new();
#undef NUM_DIMENSIONS_
#ifndef NUM_DIMENSIONS
free(sssp);
#endif
#ifdef XT_IDXSECTION_STRIPES_ISECT_SINGLE_MATCH_ONLY
free(pos_used_bm);
#else
free(isect_stripes_bm);
#endif
  return result;
}
#undef xt_idxsection_get_idxstripes_intersection__
#ifdef NUM_DIMENSIONS
#undef Xt_idxsection__
#endif
#undef XT_TOKEN_PASTE2_
/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
