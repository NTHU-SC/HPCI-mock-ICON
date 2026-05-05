/**
 * @file xt_mergesort_base.h
 * @brief macros to create mergesort implementations, 4 way top-down method
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
#ifndef XT_MERGESORT_BASE_H
#define XT_MERGESORT_BASE_H
#define TOKEN_PASTE(a,b) a##_##b
#define NAME_COMPOSE(a,b) TOKEN_PASTE(a,b)
#endif

#ifndef SORT_TYPE
#error "must define type to sort on"
#endif

#ifndef SORT_TYPE_SUFFIX
#error "must define suffix for type to name functions"
#endif

#ifndef SORT_TYPE_CMP_LT
#error "must define macro to compare SORT_TYPE for less than relation"
#endif

#ifndef SORT_TYPE_CMP_LE
#error "must define macro to compare SORT_TYPE for less than or equal relation"
#endif

#ifndef SORT_TYPE_CMP_EQ
#error "must define macro to compare SORT_TYPE for equality"
#endif

#ifndef XT_SORTFUNC_DECL
#define XT_SORTFUNC_DECL
#define XT_SORTFUNC_DECL_UNDEF
#endif

#ifndef XT_SORT_EXTRA_ALLOC_SIZE
#define XT_SORT_EXTRA_ALLOC_SIZE 0
#define XT_SORT_EXTRA_ALLOC_SIZE_UNDEF
#endif

#ifndef XT_SORT_EXTRA_ALLOC_DECL
#define XT_SORT_EXTRA_ALLOC_DECL
#define XT_SORT_EXTRA_ALLOC_DECL_UNDEF
#endif

#ifndef XT_SORT_EXTRA_ARGS_DECL
/* these declarations are appended to parameters, defaults to nothing */
#define XT_SORT_EXTRA_ARGS_DECL
#define XT_SORT_EXTRA_ARGS_DECL_UNDEF
#endif

#ifndef XT_SORT_EXTRA_ARGS_PASS
/* determines what is passed to the parameters declared in
 * XT_SORT_EXTRA_ARGS_DECL, defaults to nothing */
#define XT_SORT_EXTRA_ARGS_PASS
#define XT_SORT_EXTRA_ARGS_PASS_UNDEF
#endif

#ifndef XT_SORT_EXTRA_ARGS_INNER_DECL
/* these declarations are appended to inner implementation parameters,
   defaults to XT_SORT_EXTRA_ARGS_DECL */
#define XT_SORT_EXTRA_ARGS_INNER_DECL XT_SORT_EXTRA_ARGS_DECL
#define XT_SORT_EXTRA_ARGS_INNER_DECL_UNDEF
#endif

#ifndef XT_SORT_EXTRA_ARGS_INNER_PASS
/* determines what is passed to the parameters declared in
 * XT_SORT_EXTRA_ARGS_INNER_DECL, the macro takes the two arrays as
 * parameters, definition defaults to XT_SORT_EXTRA_ARGS_PASS */
#define XT_SORT_EXTRA_ARGS_INNER_PASS(a,b) XT_SORT_EXTRA_ARGS_PASS
#define XT_SORT_EXTRA_ARGS_INNER_PASS_UNDEF
#endif

#ifndef XT_SORT_ASSIGN
#define XT_SORT_ASSIGN(a, i, b, j) (a)[(i)] = (b)[(j)]
#define XT_SORT_ASSIGN_UNDEF
#endif

#ifndef XT_SORT_EXTRA_ARGS_SWAP
#define XT_SORT_EXTRA_ARGS_SWAP(i,j)
#define XT_SORT_EXTRA_ARGS_SWAP_UNDEF
#endif

#define XT_MERGESORT NAME_COMPOSE(xt_mergesort,SORT_TYPE_SUFFIX)
#define XT_MERGESORT_INNER NAME_COMPOSE(xt_mergesort_i,SORT_TYPE_SUFFIX)
#define XT_INSERTIONSORT NAME_COMPOSE(xt_insertionsort, SORT_TYPE_SUFFIX)
#define XT_MERGE NAME_COMPOSE(xt_merge, SORT_TYPE_SUFFIX)

#ifndef SWAP
#define SWAP(i,j) do {                                \
    SORT_TYPE t = v[i]; v[i] = v[j]; v[j] = t;        \
    XT_SORT_EXTRA_ARGS_SWAP(i, j);                    \
  } while (0)
#else
#define XT_SORT_SWAP_DEF
#endif

static inline void
XT_INSERTIONSORT(SORT_TYPE *restrict v, size_t start, size_t end
                 XT_SORT_EXTRA_ARGS_INNER_DECL)
{
  for (size_t m = start+1; m < end; ++m)
    for (size_t l = m; l > start && SORT_TYPE_CMP_LT(v[l], v[l-1],l,l-1); --l)
      SWAP(l, l-1);
}

static void
XT_MERGE(SORT_TYPE *restrict v, SORT_TYPE *restrict w,
         size_t start, size_t mid, size_t end
         XT_SORT_EXTRA_ARGS_INNER_DECL)
{
  size_t p = start, q = mid;
  while(p < mid && q < end) {
    if (SORT_TYPE_CMP_LE(v[p], v[q], p, q)) {
      XT_SORT_ASSIGN(w, start, v, p);
      ++p;
    } else {
      XT_SORT_ASSIGN(w, start, v, q);
      ++q;
    }
    start++;
  }
  for (; p < mid; ++p, ++start) {
    XT_SORT_ASSIGN(w, start, v, p);
  }
  for (; q < end; ++q, ++start) {
    XT_SORT_ASSIGN(w, start, v, q);
  }
}

static void
XT_MERGESORT_INNER(SORT_TYPE *restrict v, SORT_TYPE *restrict w,
                   size_t start, size_t end
                   XT_SORT_EXTRA_ARGS_INNER_DECL)
{
  size_t n = end - start;
  if (n<9) {
    XT_INSERTIONSORT(v, start, end XT_SORT_EXTRA_ARGS_INNER_PASS(v,w));
  } else {
    // compute 4 ranges for sub-sort
    size_t ub1 = start + n/4;
    size_t ub2 = start + 2*n/4;
    size_t ub3 = start + 3*n/4;

    XT_MERGESORT_INNER(v, w, start, ub1 XT_SORT_EXTRA_ARGS_INNER_PASS(v,w));
    XT_MERGESORT_INNER(v, w, ub1, ub2 XT_SORT_EXTRA_ARGS_INNER_PASS(v,w));
    XT_MERGESORT_INNER(v, w, ub2, ub3 XT_SORT_EXTRA_ARGS_INNER_PASS(v,w));
    XT_MERGESORT_INNER(v, w, ub3, end XT_SORT_EXTRA_ARGS_INNER_PASS(v,w));

    // 2 x 2-way merge v -> w
    XT_MERGE(v, w, start, ub1, ub2 XT_SORT_EXTRA_ARGS_INNER_PASS(v,w));
    XT_MERGE(v, w, ub2, ub3, end XT_SORT_EXTRA_ARGS_INNER_PASS(v,w));
    // final merge: w -> v
    XT_MERGE(w, v, start, ub2, end XT_SORT_EXTRA_ARGS_INNER_PASS(w,v));
  }
}

XT_SORTFUNC_DECL
void XT_MERGESORT(SORT_TYPE *restrict a, size_t n XT_SORT_EXTRA_ARGS_DECL)
{
  SORT_TYPE *v = a,
    *w = xmalloc(n * sizeof (*a) + XT_SORT_EXTRA_ALLOC_SIZE);
#define XT_SORT_EXTRA_ALLOC ((void *)(w+n))
  XT_SORT_EXTRA_ALLOC_DECL;
  XT_MERGESORT_INNER(v, w, 0, n XT_SORT_EXTRA_ARGS_INNER_PASS(v,w));
#undef XT_SORT_EXTRA_ALLOC
  free(w);
}


#ifndef XT_SORT_SWAP_DEF
#undef SWAP
#else
#undef XT_SORT_SWAP_DEF
#endif

#undef XT_MERGE
#undef XT_INSERTIONSORT
#undef XT_MERGESORT_INNER
#undef XT_MERGESORT

#ifdef XT_SORT_EXTRA_ARGS_SWAP_UNDEF
#undef XT_SORT_EXTRA_ARGS_SWAP
#undef XT_SORT_EXTRA_ARGS_SWAP_UNDEF
#endif

#ifdef XT_SORT_ASSIGN_UNDEF
#undef XT_SORT_ASSIGN
#undef XT_SORT_ASSIGN_UNDEF
#endif

#ifdef XT_SORT_EXTRA_ARGS_INNER_DECL_UNDEF
#undef XT_SORT_EXTRA_ARGS_INNER_DECL
#undef XT_SORT_EXTRA_ARGS_INNER_DECL_UNDEF
#endif

#ifdef XT_SORT_EXTRA_ARGS_INNER_PASS_UNDEF
#undef XT_SORT_EXTRA_ARGS_INNER_PASS
#undef XT_SORT_EXTRA_ARGS_INNER_PASS_UNDEF
#endif

#ifdef XT_SORT_EXTRA_ARGS_DECL_UNDEF
#undef XT_SORT_EXTRA_ARGS_DECL
#undef XT_SORT_EXTRA_ARGS_DECL_UNDEF
#endif

#ifdef XT_SORT_EXTRA_ARGS_PASS_UNDEF
#undef XT_SORT_EXTRA_ARGS_PASS
#undef XT_SORT_EXTRA_ARGS_PASS_UNDEF
#endif

#ifdef XT_SORT_EXTRA_ALLOC_DECL_UNDEF
#undef XT_SORT_EXTRA_ALLOC_DECL
#undef XT_SORT_EXTRA_ALLOC_DECL_UNDEF
#endif

#ifdef XT_SORT_EXTRA_ALLOC_SIZE_UNDEF
#undef XT_SORT_EXTRA_ALLOC_SIZE
#undef XT_SORT_EXTRA_ALLOC_SIZE_UNDEF
#endif

#ifdef XT_SORTFUNC_DECL_UNDEF
#undef XT_SORTFUNC_DECL
#undef XT_SORTFUNC_DECL_UNDEF
#endif

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
