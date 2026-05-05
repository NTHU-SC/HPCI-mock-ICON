/**
 * @file xt_xmap_dist_dir_bucket_gen_internal.h
 *
 * @brief Default bucket generator for creation of distributed directories.
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
#ifndef XT_XMAP_DIST_DIR_BUCKET_GEN_INTERNAL_H
#define XT_XMAP_DIST_DIR_BUCKET_GEN_INTERNAL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdbool.h>

#include "core/ppm_visibility.h"
#include "xt/xt_idxlist.h"
#include "xt/xt_xmap_dist_dir_bucket_gen.h"
#include "xt/xt_xmap_dist_dir_bucket_gen2.h"

struct Xt_xmdd_bucket_gen;

struct Xt_xmdd_bucket_gen_comms_f {
  MPI_Fint intra_comm, inter_comm,
    tag_offset_intra, tag_offset_inter;
};

typedef int (*Xt_xmdd_bucket_gen_init_state_internal)(
  void *gen_state,
  Xt_idxlist src_idxlist,
  Xt_idxlist dst_idxlist,
  Xt_config config,
  const struct Xt_xmdd_bucket_gen_comms *comms,
  void *init_params,
  const struct Xt_xmdd_bucket_gen_ *gen);

typedef int (*Xt_xmdd_bucket_gen_init_state_f)(
  void *gen_state,
  Xt_idxlist src_idxlist,
  Xt_idxlist dst_idxlist,
  Xt_config config,
  struct Xt_xmdd_bucket_gen_comms_f *comms,
  void *init_params);

struct Xt_xmdd_bucket_gen_ {
  /** The init function sets up the generator state. */
  Xt_xmdd_bucket_gen_init_state_internal init;
  /** The destroy function clean up the generator state. Can be zero
   * if no cleaning is needed. */
  Xt_xmdd_bucket_gen_destroy_state destroy;
  Xt_xmdd_bucket_gen_get_intersect_max_num get_intersect_max_num;
  /** The next function returns the next bucket and corresponding rank
   * (ranks can be skipped when the intersection will be empty
   * anyway).
   * Any previously returned buckets become invalid. */
  Xt_xmdd_bucket_gen_next next;
  /** gen_state_size is the size of the generator state
   *
   * The distributed directory will provide memory aligned to
   * pointer variables of this size. */
  size_t gen_state_size;
  Xt_xmdd_bucket_gen_init_state_f init_f;
  void *init_params;
};

/**
 * This is the default implementation
 */
extern const struct Xt_xmdd_bucket_gen_ Xt_xmdd_cycl_stripe_bucket_gen_desc;

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
