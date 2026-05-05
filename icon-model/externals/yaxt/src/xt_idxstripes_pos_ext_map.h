/**
 * @file xt_idxstripes_pos_ext_map.h
 *
 * @copyright Copyright  (C)  2024 Jörg Behrens <behrens@dkrz.de>
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

/* this code maps and optionally appends to an array of stripes */
static
#ifdef XT_IDXSTRIPES_POS_EXT_MAP_COUNT
#define get_mapped_stripes get_mapped_stripes_count
size_t
#else
#define get_mapped_stripes map_ext2stripes
void
#endif
get_mapped_stripes(size_t num_ext,
                   struct Xt_pos_ext *restrict pos_exts,
                   Xt_idxstripes idxstripes_src
#ifndef XT_IDXSTRIPES_POS_EXT_MAP_COUNT
                   , struct Xt_stripe *restrict result_stripes
#endif
)
{
  size_t num_result_stripes = 0;
  size_t stripes_src_ofs = 0, pos = 0;
  struct Xt_stripe result_stripe = {
    .start = XT_INT_MIN, .stride = -1, .nstrides = 0 };
  Xt_int result_stripe_next = XT_INT_MIN - 1;
  const struct Xt_stripe *src_stripes = idxstripes_src->stripes;
  for (size_t i = 0; i < num_ext; ++i) {
    size_t start_pos = (size_t)pos_exts[i].start;
    int size = pos_exts[i].size;
    if (start_pos < pos) {
      do {
        --stripes_src_ofs;
        pos -= (size_t)src_stripes[stripes_src_ofs].nstrides;
      } while (start_pos < pos);
    } else {
      while (start_pos >= pos + (size_t)src_stripes[stripes_src_ofs].nstrides) {
        pos += (size_t)src_stripes[stripes_src_ofs].nstrides;
        ++stripes_src_ofs;
      }
    }
    assert(stripes_src_ofs < (size_t)idxstripes_src->num_stripes);
    /* add overlap with current stripe to result set */
    if (size > 0) {
    build_incr_overlap:;
      struct Xt_stripe overlapping_stripe = src_stripes[stripes_src_ofs];
      size_t remaining_nstrides
        = pos + (size_t)overlapping_stripe.nstrides - start_pos;
      int intersection_size = MIN((int)remaining_nstrides, size);
      Xt_int intersection_stride = overlapping_stripe.stride,
        intersection_stripe_start = (Xt_int)(
          overlapping_stripe.start +
          (Xt_int)((size_t)overlapping_stripe.nstrides - remaining_nstrides)
          * intersection_stride);
      if ((intersection_stride == result_stripe.stride
           || intersection_size == 1)
          && intersection_stripe_start == result_stripe_next) {
        result_stripe.nstrides += intersection_size;
        intersection_stride = result_stripe.stride;
      } else if (result_stripe.nstrides == 1
                 && result_stripe.start
                 == intersection_stripe_start - intersection_stride)
      {
        result_stripe.nstrides += intersection_size;
        result_stripe.stride = intersection_stride;
      } else if (result_stripe.nstrides == 1
                 && intersection_size == 1) {
        result_stripe.nstrides += intersection_size;
        result_stripe.stride
          = intersection_stride
          = (Xt_int)(intersection_stripe_start - result_stripe.start);
      } else {
#ifndef XT_IDXSTRIPES_POS_EXT_MAP_COUNT
        if (num_result_stripes)
          result_stripes[num_result_stripes-1] = result_stripe;
#endif
        ++num_result_stripes;
        result_stripe = (struct Xt_stripe){
          .start = intersection_stripe_start,
          .stride = intersection_stride,
          .nstrides = intersection_size };
      }
      result_stripe_next = (Xt_int)(intersection_stripe_start
                             + intersection_stride * intersection_size);
      if (size -= intersection_size) {
        start_pos += (size_t)intersection_size;
        pos += (size_t)overlapping_stripe.nstrides;
        ++stripes_src_ofs;
        goto build_incr_overlap;
      }
    } else {
    build_decr_overlap:;
      struct Xt_stripe overlapping_stripe = src_stripes[stripes_src_ofs];
      size_t remaining_nstrides
        = start_pos - pos;
      int intersection_size = MIN((int)remaining_nstrides, abs(size));
      Xt_int intersection_stride = (Xt_int)-overlapping_stripe.stride,
        intersection_stripe_start = (Xt_int)(
          overlapping_stripe.start +
          (Xt_int)remaining_nstrides * (-intersection_stride));
      if ((intersection_stride == result_stripe.stride
           || intersection_size == 1)
          && intersection_stripe_start == result_stripe_next) {
        result_stripe.nstrides += intersection_size;
        intersection_stride = result_stripe.stride;
      } else if (result_stripe.nstrides == 1
                 && result_stripe.start
                 == intersection_stripe_start - intersection_stride)
      {
        result_stripe.nstrides += intersection_size;
        result_stripe.stride = intersection_stride;
      } else if (result_stripe.nstrides == 1
                 && intersection_size == 1) {
        result_stripe.nstrides += intersection_size;
        result_stripe.stride
          = intersection_stride
          = (Xt_int)(intersection_stripe_start - result_stripe.start);
      } else {
#ifndef XT_IDXSTRIPES_POS_EXT_MAP_COUNT
        if (num_result_stripes)
          result_stripes[num_result_stripes-1] = result_stripe;
#endif
        ++num_result_stripes;
        result_stripe = (struct Xt_stripe){
          .start = intersection_stripe_start,
          .stride = intersection_stride,
          .nstrides = intersection_size };
      }
      result_stripe_next = (Xt_int)(intersection_stripe_start
                             + intersection_stride * intersection_size);
      if (size += intersection_size) {
        start_pos -= (size_t)intersection_size;
        --stripes_src_ofs;
        pos -= (size_t)src_stripes[stripes_src_ofs].nstrides;
        goto build_decr_overlap;
      }
    }
  }
#ifdef XT_IDXSTRIPES_POS_EXT_MAP_COUNT
  return num_result_stripes;
#else
  result_stripes[num_result_stripes-1] = result_stripe;
#endif
}

#undef get_mapped_stripes

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
