/**
 * @file xt_config_internal.h
 * @brief implementation of configuration object
 *
 * @copyright Copyright  (C)  2020 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 *
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
#ifndef XT_CONFIG_INTERNAL_H
#define XT_CONFIG_INTERNAL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdint.h>

#include <mpi.h>

#include "core/ppm_visibility.h"
#include "xt/xt_redist.h"
#include "xt/xt_config.h"
#include "xt/sort_common.h"
#include "xt/xt_xmap_dist_dir_bucket_gen.h"
#include "xt_exchanger.h"

enum xt_config_flags {
  exch_no_dt_dup = 1 << 0,
  xt_mthread_mode_bit_ofs = 1,
  xt_mthread_mode_num_bits = 1,
  xt_mthread_mode_mask = 1 << xt_mthread_mode_bit_ofs,
  xt_force_nosort_bit_ofs = xt_mthread_mode_bit_ofs + xt_mthread_mode_num_bits,
  /* if set, rather return unsorted data than incur costly data duplication */
  xt_force_nosort = 1 << xt_force_nosort_bit_ofs,
  xt_force_xmap_striping_bit_ofs = xt_force_nosort_bit_ofs + 1,
  xt_force_xmap_striping_num_bits = 2,
  xt_force_xmap_striping_mask = 3 << xt_force_xmap_striping_bit_ofs,
};

struct Xt_sort_algo_funcptr
{
  void (*sort_int)(int *a, size_t n);
  void (*sort_xt_int)(Xt_int *a, size_t n);
  void (*sort_index)(Xt_int *restrict a, int n, int *restrict idx,
                     int reset_index);
  void (*sort_idxpos)(idxpos_type *v, size_t n);
  void (*sort_xt_int_permutation)(Xt_int a[], size_t n, int permutation[]);
  void (*sort_int_permutation)(int a[], size_t n, int permutation[]);
};


struct Xt_config_ {
  /**
   * constructor to use when creating the exchanger of a redist
   */
  Xt_exchanger_new exchanger_new;
  /**
   * function pointers to implement sort algorithms */
  const struct Xt_sort_algo_funcptr *sort_funcs;
  /**
   * description of bucket generator */
  const struct Xt_xmdd_bucket_gen_ *xmdd_bucket_gen;
  /**
   * pointer to exchanger team share data
   */
  void *exchanger_team_share;
  /**
   * automatically compress index lists of vector type at this size
   * into another representation to save on computation/memory overall
   */
  int idxv_cnv_size;
  /**
   * binary combination of xt_config_flags */
  uint32_t flags;
};

extern struct Xt_config_ xt_default_config;

PPM_DSO_INTERNAL void
xt_config_defaults_init(void);

int
xt_sort_algo_id_by_name(const char *name);

/**
 * Get appropriate exchanger constructor.
 *
 * @param config configuration object
 * @param comm communicator to use the constructor with
 * @returns configured exchanger constructor, or a fallback if the
 * configured constructor does not apply to \a comm.
 */
PPM_DSO_INTERNAL Xt_exchanger_new
xt_config_get_exchange_new_by_comm(Xt_config config, MPI_Comm comm);

#define XT_CONFIG_GET_FORCE_NOSORT(config) \
  (((config)->flags & xt_force_nosort) != UINT32_C(0))
#define XT_CONFIG_SET_FORCE_NOSORT(config) \
  do { (config)->flags |= (uint32_t)xt_force_nosort; } while (0)
#define XT_CONFIG_UNSET_FORCE_NOSORT(config) \
  do { (config)->flags &= ~(uint32_t)xt_force_nosort; } while (0)
#define XT_CONFIG_SET_FORCE_NOSORT_BIT(config, val)                     \
  do { (config)->flags = ((config)->flags &                             \
                          ~(uint32_t)xt_force_nosort)                   \
      | ((uint32_t)(val != 0) << xt_force_nosort_bit_ofs); } while (0)

#define XT_CONFIG_GET_XMAP_STRIPING(config) \
  (((config)->flags >> xt_force_xmap_striping_bit_ofs) & 3U)
#define XT_CONFIG_SET_XMAP_STRIPING(config, v)                      \
  do { (config)->flags                                                  \
      = ((config)->flags & ~(uint32_t)xt_force_xmap_striping_mask)  \
      | (uint32_t)((v&3) << xt_force_xmap_striping_bit_ofs); }      \
  while (0)

#define XT_CONFIG_BUCKET_DESTROY(config, bucket_gen_state)         \
  do { if ((config)->xmdd_bucket_gen->destroy)                  \
      (config)->xmdd_bucket_gen->destroy((bucket_gen_state));   \
  } while (0)

#endif  /* XT_CONFIG_INTERNAL_H */
/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
