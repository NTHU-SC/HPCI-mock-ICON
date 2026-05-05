/**
 * @file xt_stripe.c
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
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "xt/xt_core.h"
#include "xt/xt_stripe.h"
#include "xt_stripe_util.h"

#include "xt_stripe_util.h"
#include "core/ppm_xfuncs.h"
#include "instr.h"
#include "ensure_array_size.h"

void xt_convert_indices_to_stripes(const Xt_int *restrict indices,
                                   int num_indices,
                                   struct Xt_stripe **stripes,
                                   int * num_stripes)
{
  *stripes = NULL;
  *num_stripes = 0;
  xt_convert_indices_to_stripes_keep_buf(indices, num_indices,
                                         stripes, num_stripes);
}

size_t
xt_indices_count_stripes(size_t num_indices,
                         const Xt_int indices[num_indices])
{
  size_t num_temp_stripes = 0;
  if (num_indices > 0) {
    size_t i = 0;

    while(i < (size_t)num_indices) {
      ++num_temp_stripes;

      size_t j = 1;
      Xt_int stride = 1, stripe_base_index = indices[i];
      if (i + j < (size_t)num_indices) {
        stride = (Xt_int)(indices[i + 1] - stripe_base_index);
        do {
          ++j;
        } while ((i + j) < (size_t)num_indices
                 && stripe_base_index + (Xt_int)j * stride == indices[i + j]);
      }
      j-= ((i + j + 1 < (size_t)num_indices)
           && ((indices[i + j - 1] == indices[i + j] - 1)
               & (indices[i + j] == indices[i + j + 1] - 1)));
      i += j;
    }
  }
  return num_temp_stripes;
}


void xt_convert_indices_to_stripes_keep_buf(const Xt_int *restrict indices,
                                            int num_indices,
                                            struct Xt_stripe **stripes,
                                            int * num_stripes) {

  INSTR_DEF(instr,"xt_idxstripes_convert_to_stripes")
  INSTR_START(instr);

  size_t num_indices_ = (size_t)(num_indices > 0 ? num_indices : 0);
  size_t num_stripes_alloc = xt_indices_count_stripes(num_indices_, indices);
  struct Xt_stripe *restrict temp_stripes
    = *stripes = xrealloc(*stripes, num_stripes_alloc * sizeof(**stripes));
  size_t num_stripes_
    = xt_convert_indices_to_stripes_buf(num_indices_, indices,
                                        num_stripes_alloc, temp_stripes);
  *num_stripes = (int)num_stripes_;
  INSTR_STOP(instr);
}

size_t
xt_convert_indices_to_stripes_buf(size_t num_indices,
                                  const Xt_int *restrict indices,
                                  size_t num_stripes_alloc,
                                  struct Xt_stripe *stripes)
{
  (void)num_stripes_alloc;
  size_t num_stripes = 0;
  if (num_indices > 0) {
    size_t i = 0;

    while(i < num_indices) {
      size_t j = 1;

      Xt_int stride = 1;
      if (i + j < num_indices) {
        stride = (Xt_int)(indices[i + 1] - indices[i]);
        do {
          ++j;
        } while ((i + j) < num_indices
                 && indices[i] + (Xt_int)j * stride == indices[i + j]);
      }
      j-= ((i + j + 1 < num_indices)
           && ((indices[i + j - 1] == indices[i + j] - 1)
               & (indices[i + j] == indices[i + j + 1] - 1)));

      stripes[num_stripes++]
        = (struct Xt_stripe){ .start = indices[i], .stride = stride,
                              .nstrides = (int)j };
      i += j;
    }
  }

  assert(num_stripes <= num_stripes_alloc);
  return num_stripes;
}

size_t
xt_stripes_merge_copy(size_t num_stripes,
                      struct Xt_stripe *stripes_dst,
                      const struct Xt_stripe *stripes_src,
                      bool lookback)
{
  size_t skip = 1;
  if (num_stripes) {
    if (lookback) {
      Xt_int stride = stripes_src[0].stride,
        prev_stride = stripes_dst[-1].stride,
        start = stripes_src[0].start,
        prev_start = stripes_dst[-1].start;
      if (stride == prev_stride
          && start == prev_start + stride * (Xt_int)stripes_dst[-1].nstrides) {
        /* merge perfectly aligned stripes */
        stripes_dst[-1].nstrides += stripes_src[0].nstrides;
        ++skip;
        goto copy_loop;
      }
    }
    stripes_dst[0] = stripes_src[0];
  copy_loop:
    if (num_stripes > 1)
      for (size_t i = 1; i < num_stripes; ++i) {
        Xt_int stride = stripes_src[i].stride,
          prev_stride = stripes_dst[i-skip].stride,
          start = stripes_src[i].start,
          prev_start = stripes_dst[i-skip].start;
        if (stride == prev_stride
            && start == prev_start + stride * (Xt_int)stripes_dst[i-skip].nstrides) {
          /* merge perfectly aligned stripes */
          stripes_dst[i-skip].nstrides += stripes_src[i].nstrides;
          ++skip;
        } else
          stripes_dst[i-skip+1] = stripes_src[i];
      }
  }
  return num_stripes - (skip - 1);
}

struct Xt_stripe_summary
xt_summarize_stripes(size_t num_stripes, const struct Xt_stripe stripes[num_stripes])
{
  int ntrans_up=0, ntrans_dn=0;
  long long num_indices = 0;
  if (num_stripes > 0) {
    Xt_int prev_end;
    {
      int nstrides = stripes[0].nstrides;
      Xt_int stride = stripes[0].stride;
      prev_end = (Xt_int)(stripes[0].start
                          + stride * (nstrides-1));
      ntrans_up += ((nstrides - 1) & ~((stride > 0)-1));
      ntrans_dn += ((nstrides - 1) & ~((stride < 0)-1));
      num_indices += nstrides;
    }
    for (size_t i = 1; i < num_stripes; ++i) {
      int nstrides = stripes[i].nstrides;
      Xt_int start = stripes[i].start, stride = stripes[i].stride;
      num_indices += nstrides;
      ntrans_up += (start > prev_end) + ((nstrides - 1) & ~((stride > 0)-1));
      ntrans_dn += (start < prev_end) + ((nstrides - 1) & ~((stride < 0)-1));
      prev_end = (Xt_int)(start + stride * (nstrides-1));
    }
  }
  return (struct Xt_stripe_summary){
    .num_indices = num_indices,
    .flags = XT_SORT_FLAGS(ntrans_up, ntrans_dn),
  };
}

#define MIN(X, Y)  ((X) < (Y) ? (X) : (Y))

bool
xt_stripes_detect_duplicate(size_t num_stripes,
                            const struct Xt_stripe stripes[num_stripes],
                            struct Xt_minmax index_range)
{
  bool contains_duplicate = false;
  if (num_stripes) {
    Xt_int ofs = index_range.min;
    enum {
      block_max = 1 << 24,
      size_t_bits = sizeof (size_t) * CHAR_BIT,
    };
    size_t bm_size = ((size_t)MIN(block_max, index_range.max - ofs + 1)
                      + size_t_bits - 1 ) / size_t_bits;
    size_t *occupancy_bm = xcalloc(bm_size, sizeof (size_t));
    size_t dup = 0;
    do {
      Xt_int current_range_size
        = (Xt_int)MIN(block_max, index_range.max - ofs + 1);
      struct Xt_minmax blk_range = (struct Xt_minmax){
        ofs, (Xt_int)(ofs + current_range_size - 1)
      };
      for (size_t i = 0; i < num_stripes; ++i) {
        struct Xt_stripe stripe_i = stripes[i];
        if (stripe_i.stride < 0) {
          stripe_i.start = (Xt_int)(stripe_i.start
                                    + stripe_i.stride * (stripe_i.nstrides-1));
          stripe_i.stride = XT_INT_ABS(stripe_i.stride);
        } else if (stripe_i.stride == 0) {
          if (stripe_i.nstrides > 1) {
            free(occupancy_bm);
            return true;
          }
          stripe_i.stride = 1;
        }
        struct Xt_minmax stripe_range = (struct Xt_minmax){
          stripe_i.start, (Xt_int)(stripe_i.start
                                    + stripe_i.stride * (stripe_i.nstrides-1))
        };
        if (stripe_range.max >= blk_range.min
            && stripe_range.min <= blk_range.max)
        {
          Xt_int fstep
            = (Xt_int)((Xt_doz(ofs, stripe_i.start) + stripe_i.stride - 1)
                       / stripe_i.stride);
          for (Xt_int index
                 = (Xt_int)(stripe_i.start + fstep * stripe_i.stride);
               index <= blk_range.max && fstep < stripe_i.nstrides;
               ++fstep, index = (Xt_int)(index + stripe_i.stride))
          {
            size_t mask = (size_t)1 << ((size_t)(index - ofs)%size_t_bits),
              mask_ofs = (size_t)(index - ofs)/size_t_bits;
            dup |= mask & occupancy_bm[mask_ofs];
            occupancy_bm[mask_ofs] |= mask;
          }
        }
      }
      Xt_uint remaining = (Xt_uint)((Xt_uint)index_range.max - (Xt_uint)ofs);
      if (dup || remaining >= (Xt_uint)block_max)
        break;
      ofs = (Xt_int)(ofs + block_max);
      memset(occupancy_bm, 0, bm_size * sizeof (size_t));
    } while (true);
    free(occupancy_bm);
    contains_duplicate = dup;
  }
  return contains_duplicate;
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
