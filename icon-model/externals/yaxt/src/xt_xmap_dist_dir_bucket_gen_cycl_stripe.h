/**
 * @file xt_xmap_dist_dir_bucket_gen_cycl_stripe.h
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
#ifndef XT_XMAP_DIST_DIR_BUCKET_GEN_CYCL_STRIPE_H
#define XT_XMAP_DIST_DIR_BUCKET_GEN_CYCL_STRIPE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdbool.h>

#include "core/ppm_visibility.h"
#include "xt/xt_idxlist.h"
#include "xt/xt_xmap_dist_dir_bucket_gen.h"

struct bucket_params {
  Xt_int global_interval, local_interval,
    local_index_range_lbound, local_index_range_ubound;
  int most_recent_rank_generated, max_rank_generated;
};

struct Xt_xmdd_bucket_gen_cycl_stripe_state {
  /* pointer to array of stripes holding the data for last bucket
   * generated for dist_dir_rank.  Must point to a pointer allocated
   * by malloc or related and will be resized if necessary. */
  struct Xt_stripe *stripes;
  /* number of stripes that can be stored in \a stripes */
  size_t stripes_array_size;
  Xt_idxlist last_list_generated;
  struct bucket_params dst, src;
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
