/**
 * @file xt_xmap_dist_dir_bucket_gen.c
 *
 * @brief Implements class hiding different bucket generators
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>

#include "xt/xt_xmap_dist_dir_bucket_gen.h"
#include "xt_xmap_dist_dir_bucket_gen_internal.h"
#include "core/ppm_xfuncs.h"

Xt_xmdd_bucket_gen
xt_xmdd_bucket_gen_new(void)
{
  return xcalloc(sizeof (struct Xt_xmdd_bucket_gen_), 1);
}

void
xt_xmdd_bucket_gen_delete(Xt_xmdd_bucket_gen gen)
{
  free(gen);
}

void
xt_xmdd_bucket_gen_define_interface(
  Xt_xmdd_bucket_gen gen,
  Xt_xmdd_bucket_gen_init_state init,
  Xt_xmdd_bucket_gen_destroy_state destroy,
  Xt_xmdd_bucket_gen_get_intersect_max_num get_intersect_max_num,
  Xt_xmdd_bucket_gen_next next,
  size_t gen_state_size,
  void *init_params)
{
  gen->init_f = (Xt_xmdd_bucket_gen_init_state_f)0;
  gen->init = (Xt_xmdd_bucket_gen_init_state_internal)(void (*)(void))init;
  gen->destroy = destroy;
  gen->get_intersect_max_num = get_intersect_max_num;
  gen->next = next;
  gen->gen_state_size = gen_state_size;
  gen->init_params = init_params;
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
