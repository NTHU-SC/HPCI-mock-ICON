/**
 * @file mergesort.c
 *
 * @copyright Copyright  (C)  2012 Jörg Behrens <behrens@dkrz.de>
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

#include <stdlib.h>
#include <stdio.h>

#include "xt/mergesort.h"
#include "xt/xt_sort.h"
#include "core/ppm_xfuncs.h"

#define SORT_TYPE idxpos_type
#define SORT_TYPE_SUFFIX idxpos
#define SORT_TYPE_CMP_LT(a,b,...) ((a).idx < (b).idx || ((a).idx == (b).idx && (a).pos < (b).pos))
#define SORT_TYPE_CMP_LE(a,b,...) ((a).idx < (b).idx || ((a).idx == (b).idx && (a).pos <= (b).pos))
#define SORT_TYPE_CMP_EQ(a,b,...) ((a).idx == (b).idx && (a).pos == (b).pos)
#include "xt_mergesort_base.h"


void
xt_mergesort_xt_int_permutation(Xt_int *a, size_t n, int *restrict permutation);


void xt_mergesort_index(Xt_int *restrict v_idx, int n,
                        int *restrict v_pos, int reset_pos) {

  /*  Initialization */
  int *v_pos_orig = v_pos;

  if (!v_pos) v_pos = xmalloc((size_t)n  * sizeof(v_pos[0]));

  if (v_pos != v_pos_orig || reset_pos)
    xt_assign_id_map_int((size_t)n, v_pos, 0);

  xt_mergesort_xt_int_permutation(v_idx, (size_t)n, v_pos);
  if (v_pos != v_pos_orig) free(v_pos);
}

#undef SORT_TYPE
#undef SORT_TYPE_SUFFIX
#undef SORT_TYPE_CMP_LT
#undef SORT_TYPE_CMP_LE
#undef SORT_TYPE_CMP_EQ
#define SORT_TYPE int
#define SORT_TYPE_SUFFIX int
#define SORT_TYPE_CMP_LT(a,b,...) ((a) < (b))
#define SORT_TYPE_CMP_LE(a,b,...) ((a) <= (b))
#define SORT_TYPE_CMP_EQ(a,b,...) ((a) == (b))
#include "xt_mergesort_base.h"

#undef SORT_TYPE
#undef SORT_TYPE_SUFFIX
#undef SORT_TYPE_CMP_LT
#undef SORT_TYPE_CMP_LE
#undef SORT_TYPE_CMP_EQ
#define SORT_TYPE Xt_int
#define SORT_TYPE_SUFFIX xt_int
#define SORT_TYPE_CMP_LT(a,b,...) ((a) < (b))
#define SORT_TYPE_CMP_LE(a,b,...) ((a) <= (b))
#define SORT_TYPE_CMP_EQ(a,b,...) ((a) == (b))
#include "xt_mergesort_base.h"

#undef SORT_TYPE
#undef SORT_TYPE_SUFFIX
#undef SORT_TYPE_CMP_LT
#undef SORT_TYPE_CMP_LE
#undef SORT_TYPE_CMP_EQ
#define SORT_TYPE Xt_int
#define SORT_TYPE_SUFFIX xt_int_permutation
#define SORT_TYPE_CMP_LT(u,v,i,j) \
  ((u) < (v) || ((u) == (v) && permutation_v[(i)] < permutation_v[(j)]))
#define SORT_TYPE_CMP_LE(u,v,i,j) \
  ((u) < (v) || ((u) == (v) && permutation_v[(i)] <= permutation_v[(j)]))
#define SORT_TYPE_CMP_EQ(u,v,i,j) \
  ((u) == (v) && permutation_v[(i)] == permutation_v[(j)])
#define XT_SORT_EXTRA_ARGS_DECL , int *restrict permutation
#define XT_SORT_EXTRA_ARGS_PASS , permutation
#define XT_SORT_EXTRA_ARGS_INNER_DECL , int *restrict permutation_v, int *restrict permutation_w
#define XT_SORT_EXTRA_ARGS_INNER_PASS(a,b) , TOKEN_PASTE(permutation,a), \
    TOKEN_PASTE(permutation,b)
#define XT_SORT_EXTRA_ALLOC_SIZE (sizeof(*permutation) * (n+1))
static inline void *
roundPtr(void *p, size_t mult)
{
  uintptr_t ip = (uintptr_t)p;
  return (void *)(((ip + mult - 1) / mult) * mult);
}

#define XT_SORT_EXTRA_ALLOC_DECL int *permutation_v = permutation, \
    *permutation_w = roundPtr(XT_SORT_EXTRA_ALLOC, sizeof (int))
#define XT_SORT_EXTRA_ARGS_SWAP(i,j) do {                               \
    size_t i_ = (size_t)(i), j_ = (size_t)(j);                          \
    int tp = permutation_v[i_];                                         \
    permutation_v[i_] = permutation_v[j_];                              \
    permutation_v[j_] = tp;                                             \
    (void)permutation_w;                                                \
  } while (0)
#define XT_SORT_ASSIGN(a, p, b, q) do {         \
    (a)[(p)] = (b)[(q)];                        \
    TOKEN_PASTE(permutation,a)[(p)] = TOKEN_PASTE(permutation,b)[(q)]; \
  } while (0)
#include "xt_mergesort_base.h"

#undef SORT_TYPE
#undef SORT_TYPE_SUFFIX
#undef SORT_TYPE_CMP_LT
#undef SORT_TYPE_CMP_LE
#undef SORT_TYPE_CMP_EQ
#undef XT_SORT_EXTRA_ARGS_DECL
#undef XT_SORT_EXTRA_ARGS_PASS
#undef XT_SORT_EXTRA_ARGS_INNER_DECL
#undef XT_SORT_EXTRA_ARGS_INNER_PASS
#undef XT_SORT_EXTRA_ARGS_SWAP
#undef XT_SORT_EXTRA_ALLOC_SIZE
#undef XT_SORT_EXTRA_ALLOC_DECL
#undef XT_SORT_ASSIGN
#define SORT_TYPE int
#define SORT_TYPE_SUFFIX int_permutation
#define SORT_TYPE_CMP_LT(u,v,i,j) ((u) < (v) || (u == v && permutation_v[(i)] < permutation_v[(j)]))
#define SORT_TYPE_CMP_LE(u,v,i,j) ((u) < (v) || (u == v && permutation_v[(i)] <= permutation_v[(j)]))
#define SORT_TYPE_CMP_EQ(u,v,i,j) ((u) == (v) && permutation_v[(i)] == permutation_v[(j)])
#define XT_SORT_EXTRA_ARGS_DECL , int *restrict permutation
#define XT_SORT_EXTRA_ARGS_PASS , permutation
#define XT_SORT_EXTRA_ARGS_INNER_DECL , int *restrict permutation_v, int *restrict permutation_w
#define XT_SORT_EXTRA_ARGS_INNER_PASS(a,b) , TOKEN_PASTE(permutation,a), \
    TOKEN_PASTE(permutation,b)
#define XT_SORT_EXTRA_ALLOC_SIZE (sizeof(*permutation) * (n+1))
#define XT_SORT_EXTRA_ALLOC_DECL int *permutation_v = permutation, \
    *permutation_w = XT_SORT_EXTRA_ALLOC
#define XT_SORT_EXTRA_ARGS_SWAP(i,j) do {                               \
    size_t i_ = (size_t)(i), j_ = (size_t)(j);                          \
    int tp = permutation_v[i_];                                         \
    permutation_v[i_] = permutation_v[j_];                              \
    permutation_v[j_] = tp;                                             \
    (void)permutation_w;                                                \
  } while (0)
#define XT_SORT_ASSIGN(a, p, b, q) do {         \
    (a)[(p)] = (b)[(q)];                        \
    TOKEN_PASTE(permutation,a)[(p)] = TOKEN_PASTE(permutation,b)[(q)]; \
  } while (0)
#include "xt_mergesort_base.h"

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
