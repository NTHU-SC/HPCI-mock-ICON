/**
 * @file xt_mpi.c
 *
 * @copyright Copyright  (C)  2013 Jörg Behrens <behrens@dkrz.de>
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
#include <inttypes.h>
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_core.h"
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"

#if ! (HAVE_DECL___BUILTIN_CTZL || HAVE_DECL___BUILTIN_CLZL)       \
  && (HAVE_DECL___LZCNT && SIZEOF_LONG == SIZEOF_INT               \
      || HAVE_DECL___LZCNT64 && SIZEOF_LONG == 8 && CHAR_BIT == 8)
#include <intrin.h>
#endif

//taken from http://beige.ucs.indiana.edu/I590/node85.html
void xt_mpi_error(int error_code, MPI_Comm comm) {
  int rank;
  MPI_Comm_rank(comm, &rank);

  char error_string[MPI_MAX_ERROR_STRING];
  int length_of_error_string, error_class;

  MPI_Error_class(error_code, &error_class);
  MPI_Error_string(error_class, error_string, &length_of_error_string);
  fprintf(stderr, "%3d: %s\n", rank, error_string);
  MPI_Error_string(error_code, error_string, &length_of_error_string);
  fprintf(stderr, "%3d: %s\n", rank, error_string);
  MPI_Abort(comm, error_code);
}


size_t
xt_disp2ext_count(size_t disp_len, const int *disp)
{
  if (!disp_len) return 0;
  size_t i = 0;
  int cur_stride = 1, cur_size = 1;
  int last_disp = disp[0];
  for (size_t p = 1; p < disp_len; ++p) {
    int new_disp = disp[p];
    int new_stride = new_disp - last_disp;
    if (cur_size == 1) {
      cur_stride = new_stride;
      cur_size = 2;
    } else if (new_stride == cur_stride) {
      // cur_size >= 2:
      cur_size++;
    } else if (cur_size > 2 || (cur_size == 2 && cur_stride == 1) ) {
      // we accept small contiguous vectors (nstrides==2, stride==1)
      i++;
      cur_stride = 1;
      cur_size = 1;
    } else { // cur_size == 2, next offset doesn't match current stride
      // break up trivial vec:
      i++;
      cur_size = 2;
      cur_stride = new_stride;
    }
    last_disp = new_disp;
  }
  // tail cases:
  if (cur_size > 2 || (cur_size == 2 && cur_stride == 1)) {
    i++;
  } else if (cur_size == 2) {
    i+=2;
  } else { // cur_size == 1
    i++;
  }

  return i;
}

size_t
xt_disp2ext(size_t disp_len, const int *disp,
            struct Xt_offset_ext *restrict v)
{
  if (disp_len<1) return 0;

  int cur_start = disp[0], cur_stride = 1, cur_size = 1;
  int last_disp = cur_start;
  size_t i = 0;
  for (size_t p = 1; p < disp_len; ++p) {
    int new_disp = disp[p];
    int new_stride = new_disp - last_disp;
    if (cur_size == 1) {
      cur_stride = new_stride;
      cur_size = 2;
    } else if (new_stride == cur_stride) {
      // cur_size >= 2:
      cur_size++;
    } else if (cur_size > 2 || (cur_size == 2 && cur_stride == 1) ) {
      // we accept small contiguous vectors (nstrides==2, stride==1)
      v[i] = (struct Xt_offset_ext){ .start = cur_start, .stride = cur_stride,
                                     .size = cur_size };
      i++;
      cur_start = new_disp;
      cur_stride = 1;
      cur_size = 1;
    } else { // cur_size == 2, next offset doesn't match current stride
      // break up trivial vec:
      v[i].start = cur_start;
      v[i].size = 1;
      v[i].stride = 1;
      i++;
      cur_start += cur_stride;
      cur_size = 2;
      cur_stride = new_stride;
    }
    last_disp = new_disp;
  }
  // tail cases:
  if (cur_size > 2 || (cur_size == 2 && cur_stride == 1)) {
    v[i] = (struct Xt_offset_ext){ .start = cur_start, .stride = cur_stride,
                                   .size = cur_size };
    i++;
  } else if (cur_size == 2) {
    v[i].start = cur_start;
    v[i].size = 1;
    v[i].stride = 1;
    i++;
    v[i].start = cur_start + cur_stride;
    v[i].size = 1;
    v[i].stride = 1;
    i++;
  } else { // cur_size == 1
    v[i].start = cur_start;
    v[i].size = 1;
    v[i].stride = 1;
    i++;
  }

  return i;
}

/* functions to handle optimizations on communicators */
static int xt_mpi_comm_internal_keyval = MPI_KEYVAL_INVALID;

typedef unsigned long used_map_elem;

enum {
  used_map_elem_bits = sizeof (used_map_elem) * CHAR_BIT,
};

struct xt_mpi_comm_internal_attr {
  int refcount;
  unsigned used_map_size;
  used_map_elem used_map[];
};

static int
xt_mpi_comm_internal_keyval_copy(
  MPI_Comm XT_UNUSED(oldcomm), int XT_UNUSED(keyval),
  void *XT_UNUSED(extra_state), void *XT_UNUSED(attribute_val_in),
  void *attribute_val_out, int *flag)
{
  struct xt_mpi_comm_internal_attr *new_comm_attr
    = malloc(sizeof (struct xt_mpi_comm_internal_attr)
             + sizeof (used_map_elem));
  int retval;
  if (new_comm_attr)
  {
    new_comm_attr->refcount = 1;
    new_comm_attr->used_map_size = 1;
    new_comm_attr->used_map[0] = 1U;
    *(void **)attribute_val_out = new_comm_attr;
    *flag = 1;
    retval = MPI_SUCCESS;
  } else {
    *flag = 0;
    retval = MPI_ERR_NO_MEM;
  }
  return retval;
}

static int
xt_mpi_comm_internal_keyval_delete(
  MPI_Comm XT_UNUSED(comm), int XT_UNUSED(comm_keyval),
  void *attribute_val, void *XT_UNUSED(extra_state))
{
  free(attribute_val);
  return MPI_SUCCESS;
}

static int xt_mpi_tag_ub_val;

void
xt_mpi_init(void) {
  assert(xt_mpi_comm_internal_keyval == MPI_KEYVAL_INVALID);
  xt_mpi_call(MPI_Comm_create_keyval(xt_mpi_comm_internal_keyval_copy,
                                     xt_mpi_comm_internal_keyval_delete,
                                     &xt_mpi_comm_internal_keyval, NULL),
              Xt_default_comm);
  void *attr;
  int flag;
  xt_mpi_call(MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &attr, &flag),
              MPI_COMM_WORLD);
  assert(flag);
  xt_mpi_tag_ub_val = *(int *)attr;
}

void
xt_mpi_finalize(void) {
  assert(xt_mpi_comm_internal_keyval != MPI_KEYVAL_INVALID);
  xt_mpi_call(MPI_Comm_free_keyval(&xt_mpi_comm_internal_keyval),
              Xt_default_comm);
}

static struct xt_mpi_comm_internal_attr *
xt_mpi_comm_get_internal_attr(MPI_Comm comm)
{
  int attr_found;
  void *attr_val;
  assert(xt_mpi_comm_internal_keyval != MPI_KEYVAL_INVALID);
  xt_mpi_call(MPI_Comm_get_attr(comm, xt_mpi_comm_internal_keyval,
                                &attr_val, &attr_found),
              comm);
  return attr_found ? attr_val : NULL;
}

#if HAVE_DECL___BUILTIN_CTZL
#define ctzl(v) (__builtin_ctzl(v))
#elif HAVE_DECL___BUILTIN_CLZL                                     \
  || HAVE_DECL___LZCNT && SIZEOF_LONG == SIZEOF_INT                \
  || HAVE_DECL___LZCNT64 && SIZEOF_LONG == 8 && CHAR_BIT == 8
static inline int
ctzl(unsigned long v) {
  enum {
    ulong_bits = sizeof (unsigned long) * CHAR_BIT,
  };
  /* clear all but lowest 1 bit */
  v = v & ~(v - 1);
  int c = ulong_bits - 1 - (int)
#if HAVE_DECL___BUILTIN_CTZL
    __builtin_clzl(v)
#elif HAVE_DECL___LZCNT && SIZEOF_LONG == SIZEOF_INT
    __lzcnt(v)
#else
    __lzcnt64(v)
#endif
    ;
  return c;
}
#else
static inline int
ctzl(unsigned long v) {
  enum {
    ulong_bits = sizeof (unsigned long) * CHAR_BIT,
  };
  // c will be the number of zero bits on the right
  unsigned int c = ulong_bits;
  v &= (unsigned long)-(long)v;
  if (v) c--;
#if SIZEOF_UNSIGNED_LONG * CHAR_BIT == 64
  if (v & UINT64_C(0x00000000ffffffff)) c -= 32;
  if (v & UINT64_C(0x0000ffff0000ffff)) c -= 16;
  if (v & UINT64_C(0x00ff00ff00ff00ff)) c -= 8;
  if (v & UINT64_C(0x0f0f0f0f0f0f0f0f)) c -= 4;
  if (v & UINT64_C(0x3333333333333333)) c -= 2;
  if (v & UINT64_C(0x5555555555555555)) c -= 1;
#elif SIZEOF_UNSIGNED_LONG * CHAR_BIT == 32
  if (v & 0x0000FFFFUL) c -= 16;
  if (v & 0x00FF00FFUL) c -= 8;
  if (v & 0x0F0F0F0FUL) c -= 4;
  if (v & 0x33333333UL) c -= 2;
  if (v & 0x55555555UL) c -= 1;
#else
  error "Unexpected size of long.\n"
#endif
  return (int)c;
}
#endif

MPI_Comm
xt_mpi_comm_smart_dup(MPI_Comm comm, int *tag_offset)
{
  MPI_Comm comm_dest;
  struct xt_mpi_comm_internal_attr *comm_xt_attr_val
    = xt_mpi_comm_get_internal_attr(comm);
  size_t position = 0;
  int refcount = comm_xt_attr_val ? comm_xt_attr_val->refcount : 0;
  if (comm_xt_attr_val
      && (refcount + 1) < xt_mpi_tag_ub_val / xt_mpi_num_tags) {
    comm_dest = comm;
    comm_xt_attr_val->refcount = ++refcount;
    size_t used_map_size = comm_xt_attr_val->used_map_size;
    while (position < used_map_size
           && comm_xt_attr_val->used_map[position] == ~(used_map_elem)0)
      ++position;
    if (position >= used_map_size) {
      /* sadly, we need to recreate the value to enlarge it */
      struct xt_mpi_comm_internal_attr *new_comm_xt_attr_val
        = xmalloc(sizeof (*new_comm_xt_attr_val)
                  + (used_map_size + 1) * sizeof (used_map_elem));
      new_comm_xt_attr_val->refcount = refcount;
      new_comm_xt_attr_val->used_map_size = (unsigned)(used_map_size + 1);
      for (size_t i = 0; i < used_map_size; ++i)
        new_comm_xt_attr_val->used_map[i] = comm_xt_attr_val->used_map[i];
      new_comm_xt_attr_val->used_map[used_map_size] = 1U;
      position *= used_map_elem_bits;
      assert(xt_mpi_comm_internal_keyval != MPI_KEYVAL_INVALID);
      xt_mpi_call(MPI_Comm_set_attr(comm_dest, xt_mpi_comm_internal_keyval,
                                    new_comm_xt_attr_val), comm_dest);
    } else {
      /* not all bits are set, find first unset position and insert */
      used_map_elem used_map_entry = comm_xt_attr_val->used_map[position],
        unset_lsb = ~used_map_entry & (used_map_entry + 1),
        bit_pos = (used_map_elem)ctzl(unset_lsb);
      comm_xt_attr_val->used_map[position] = used_map_entry | unset_lsb;
      position = position * used_map_elem_bits + (size_t)bit_pos;
    }
  } else {
    struct xt_mpi_comm_internal_attr *comm_attr
      = xmalloc(sizeof (*comm_attr) + sizeof (used_map_elem));
    comm_attr->refcount = 1;
    comm_attr->used_map_size = 1;
    comm_attr->used_map[0] = 1U;
    xt_mpi_call(MPI_Comm_dup(comm, &comm_dest), comm);
    assert(xt_mpi_comm_internal_keyval != MPI_KEYVAL_INVALID);
    xt_mpi_call(MPI_Comm_set_attr(comm_dest, xt_mpi_comm_internal_keyval,
                                  comm_attr), comm_dest);
  }
  *tag_offset = (int)(position * xt_mpi_num_tags);
  return comm_dest;
}

void
xt_mpi_comm_smart_dedup(MPI_Comm *comm, int tag_offset)
{
  struct xt_mpi_comm_internal_attr *comm_xt_attr_val
    = xt_mpi_comm_get_internal_attr(*comm);
  int refcount = comm_xt_attr_val ? --(comm_xt_attr_val->refcount) : 0;
  if (refcount < 1) {
    xt_mpi_call(MPI_Comm_free(comm), MPI_COMM_WORLD);
    *comm = MPI_COMM_NULL;
  } else {
    size_t position = (size_t)tag_offset / xt_mpi_num_tags,
      map_elem = position / used_map_elem_bits,
      in_elem_bit = position % used_map_elem_bits;
    comm_xt_attr_val->used_map[map_elem] &= ~((used_map_elem)1 << in_elem_bit);
  }
}

void
xt_mpi_comm_mark_exclusive(MPI_Comm comm) {
  struct xt_mpi_comm_internal_attr *comm_attr
    = xmalloc(sizeof (*comm_attr) + sizeof (used_map_elem));
  comm_attr->refcount = 1;
  comm_attr->used_map_size = 1;
  comm_attr->used_map[0] = 1U;
  assert(xt_mpi_comm_internal_keyval != MPI_KEYVAL_INVALID);
  xt_mpi_call(MPI_Comm_set_attr(comm, xt_mpi_comm_internal_keyval,
                                comm_attr), comm);
}

bool
xt_mpi_test_some(int *restrict num_req,
                 MPI_Request *restrict req,
                 int *restrict ops_completed, MPI_Comm comm)
{
  int done_count;
  size_t num_req_ = (size_t)*num_req;

#if __GNUC__ >= 11 && __GNUC__ <= 13
  /* GCC 11 has no means to specify that the special value pointer
   * MPI_STATUSES_IGNORE does not need to point to something of size > 0 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif
  xt_mpi_call(MPI_Testsome(*num_req, req, &done_count, ops_completed,
                           MPI_STATUSES_IGNORE), comm);
#if __GNUC__ >= 11 && __GNUC__ <= 13
#pragma GCC diagnostic pop
#endif

  if (done_count != MPI_UNDEFINED) {
    if (num_req_ > (size_t)done_count) {
      for (size_t i = 0, j = num_req_;
           i < (size_t)done_count && j >= num_req_ - (size_t)done_count;
           ++i)
        if (ops_completed[i] < (int)num_req_ - done_count) {
          while (req[--j] == MPI_REQUEST_NULL);
          req[ops_completed[i]] = req[j];
        }
      num_req_ -= (size_t)done_count;
    }
    else
      num_req_ = 0;
  }
  *num_req = (int)num_req_;
  return num_req_ == 0;
}

#ifdef _OPENMP
bool
xt_mpi_test_some_mt(int *restrict num_req,
                    MPI_Request *restrict req,
                    int *restrict ops_completed, MPI_Comm comm)
{
  int done_count;
  size_t num_req_ = (size_t)*num_req;

  size_t num_threads = (size_t)omp_get_num_threads(),
    tid = (size_t)omp_get_thread_num();
  size_t start_req = (num_req_ * tid) / num_threads,
    nreq_ = (num_req_ * (tid+1)) / num_threads - start_req;

  for (size_t i = start_req; i < start_req + nreq_; ++i)
    ops_completed[i] = -1;
#if __GNUC__ >= 11 && __GNUC__ <= 13
  /* GCC 11 has no means to specify that the special value pointer
   * MPI_STATUSES_IGNORE does not need to point to something of size > 0 */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#pragma GCC diagnostic ignored "-Wstringop-overread"
#endif
  xt_mpi_call(MPI_Testsome((int)nreq_, req+start_req, &done_count,
                           ops_completed+start_req, MPI_STATUSES_IGNORE), comm);
#if __GNUC__ >= 11 && __GNUC__ <= 13
#pragma GCC diagnostic pop
#endif
  if (done_count == MPI_UNDEFINED)
    done_count = 0;
#pragma omp barrier
#pragma omp atomic
  *num_req -= done_count;
#pragma omp barrier
  done_count = (int)num_req_ - *num_req;
#pragma omp single
  {
    if (num_req_ > (size_t)done_count) {
      for (size_t i = 0, j = 0; i < num_req_; ++i)
        if (req[i] != MPI_REQUEST_NULL)
          req[j++] = req[i];
    }
    *num_req = (int)num_req_ - done_count;
  }
  num_req_ -= (size_t)done_count;
  return num_req_ == 0;
}
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
