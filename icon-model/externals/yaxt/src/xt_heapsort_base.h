/**
 * @file xt_heapsort_base.h
 * @brief macros to create heapsort implementations
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

#ifndef XT_SORTFUNC_DECL
#define XT_SORTFUNC_DECL
#define XT_SORTFUNC_DECL_UNDEF
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

#ifndef XT_SORT_ASSIGN
#define XT_SORT_ASSIGN(a, i, b, j) (a)[(i)] = (b)[(j)]
#define XT_SORT_ASSIGN_UNDEF
#endif

#ifndef XT_SORT_EXTRA_ARGS_SWAP
#define XT_SORT_EXTRA_ARGS_SWAP(i,j)
#define XT_SORT_EXTRA_ARGS_SWAP_UNDEF
#endif

#define XT_HEAPSORT NAME_COMPOSE(xt_heapsort, SORT_TYPE_SUFFIX)
#define XT_HEAPIFY NAME_COMPOSE(xt_heapify, SORT_TYPE_SUFFIX)

#ifndef SWAP
#define SWAP(i,j) do {                                \
    SORT_TYPE t = v[i]; v[i] = v[j]; v[j] = t;        \
    XT_SORT_EXTRA_ARGS_SWAP(i, j);                    \
  } while (0)
#else
#define XT_SORT_SWAP_DEF
#endif

static inline size_t
left(size_t i)
{
  return 2*i + 1;
}

static inline size_t
right(size_t i)
{
  return 2*i + 2;
}

static inline size_t
parent(size_t i)
{
  return (i - 1) / 2;
}

XT_SORTFUNC_DECL void
XT_HEAPIFY(SORT_TYPE *restrict v, size_t n, size_t i
           XT_SORT_EXTRA_ARGS_DECL)
{
  assert(i < n);
  do {
    size_t l = left(i),
      r = right(i), largest;
    if (l < n && SORT_TYPE_CMP_LT(v[i], v[l], i, l)) {
      largest = l;
    } else
      largest = i;
    if (r < n && SORT_TYPE_CMP_LT(v[largest], v[r], largest, r))
      largest = r;
    if (largest == i) break;
    SWAP(i, largest);
    i = largest;
  } while(1);
}



XT_SORTFUNC_DECL void
XT_HEAPSORT(SORT_TYPE heap[], size_t n
            XT_SORT_EXTRA_ARGS_DECL)
{
  for (size_t i = n/2; i--;)
    XT_HEAPIFY(heap, n, i XT_SORT_EXTRA_ARGS_PASS);
}


#ifndef XT_SORT_SWAP_DEF
#undef SWAP
#else
#undef XT_SORT_SWAP_DEF
#endif

#undef XT_HEAPIFY
#undef XT_HEAPSORT

#ifdef XT_SORT_EXTRA_ARGS_SWAP_UNDEF
#undef XT_SORT_EXTRA_ARGS_SWAP
#undef XT_SORT_EXTRA_ARGS_SWAP_UNDEF
#endif

#ifdef XT_SORT_ASSIGN_UNDEF
#undef XT_SORT_ASSIGN
#undef XT_SORT_ASSIGN_UNDEF
#endif

#ifdef XT_SORT_EXTRA_ARGS_DECL_UNDEF
#undef XT_SORT_EXTRA_ARGS_DECL
#undef XT_SORT_EXTRA_ARGS_DECL_UNDEF
#endif

#ifdef XT_SORT_EXTRA_ARGS_PASS_UNDEF
#undef XT_SORT_EXTRA_ARGS_PASS
#undef XT_SORT_EXTRA_ARGS_PASS_UNDEF
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
