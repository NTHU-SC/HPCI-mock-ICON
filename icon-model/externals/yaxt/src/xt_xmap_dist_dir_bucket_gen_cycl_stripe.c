/**
 * @file xt_xmap_dist_dir_bucket_gen_cycl_stripe.c
 *
 * @brief Implementation of default bucket generator for the creation of
 * distributed directories.
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

#include <string.h>
#include <assert.h>

#include <mpi.h>

#include "xt/xt_xmap_dist_dir_bucket_gen.h"
#include "xt_xmap_dist_dir_bucket_gen_cycl_stripe.h"
#include "xt_xmap_dist_dir_common.h"
#include "xt_idxlist_internal.h"
#include "xt_idxstripes_internal.h"
#include "xt_idxlist_unpack.h"
#include "xt_config_internal.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_mpi.h"
#include "xt_mpi_internal.h"
#include "xt_arithmetic_util.h"
#include "ensure_array_size.h"
#include "xt_xmap_dist_dir_bucket_gen_internal.h"

static int
xt_xmdd_bucket_gen_cycl_stripe_init(
  void *gen_state_,
  Xt_idxlist src_idxlist,
  Xt_idxlist dst_idxlist,
  Xt_config config,
  struct Xt_xmdd_bucket_gen_comms *comms,
  void *init_params);

static struct Xt_com_list
xt_xmdd_cycl_stripe_get_next_bucket(void *gen_state_, int type);

static void
xt_xmdd_bucket_gen_cycl_stripe_destroy(void *gen_state);

static int
xt_xmdd_bucket_gen_cycl_stripe_get_intersect_max_num(void *gen_state, int type);

const struct Xt_xmdd_bucket_gen_ Xt_xmdd_cycl_stripe_bucket_gen_desc
= { .init
    = (Xt_xmdd_bucket_gen_init_state_internal)
    (void (*)(void))xt_xmdd_bucket_gen_cycl_stripe_init,
    .destroy = xt_xmdd_bucket_gen_cycl_stripe_destroy,
    .get_intersect_max_num = xt_xmdd_bucket_gen_cycl_stripe_get_intersect_max_num,
    .next = xt_xmdd_cycl_stripe_get_next_bucket,
    .gen_state_size = sizeof (struct Xt_xmdd_bucket_gen_cycl_stripe_state),
};


static inline Xt_int
get_intracomm_dist_dir_global_interval_size(
  Xt_idxlist src, Xt_idxlist dst,
  bool *stripify, MPI_Comm comm, int comm_size,
  Xt_config config)
{
  unsigned long long local_vals[2], global_sums[2];

  unsigned num_indices_src = (unsigned)xt_idxlist_get_num_indices(src);
  local_vals[0] = num_indices_src;
  local_vals[1] = xt_idxlist_is_stripe_conversion_profitable_(src, config)
    || xt_idxlist_is_stripe_conversion_profitable_(dst, config);

  xt_mpi_call(MPI_Allreduce(local_vals, global_sums, 2,
                            MPI_UNSIGNED_LONG_LONG, MPI_SUM, comm), comm);

  *stripify = global_sums[1] > 0;
  return (Xt_int)(MAX(((global_sums[0] + (unsigned)comm_size - 1)
                      / (unsigned)comm_size), 1) * (unsigned)comm_size);
}

static inline Xt_int get_min_idxlist_index(Xt_idxlist a, Xt_idxlist b) {

  int num_a = xt_idxlist_get_num_indices(a),
    num_b = xt_idxlist_get_num_indices(b);
  Xt_int min_index_a = num_a ? xt_idxlist_get_min_index(a) : XT_INT_MAX,
    min_index_b = num_b ? xt_idxlist_get_min_index(b) : XT_INT_MAX,
    min_index = (Xt_int)MIN(min_index_a, min_index_b);
  return min_index;
}

static inline Xt_int get_max_idxlist_index(Xt_idxlist a, Xt_idxlist b) {

  int num_a = xt_idxlist_get_num_indices(a),
    num_b = xt_idxlist_get_num_indices(b);
  Xt_int max_index_a = num_a ? xt_idxlist_get_max_index(a) : XT_INT_MIN,
    max_index_b = num_b ? xt_idxlist_get_max_index(b) : XT_INT_MIN,
    max_index = MAX(max_index_a, max_index_b);
  return max_index;
}

static struct bucket_params
get_intracomm_bucket_params(
  Xt_int global_interval,
  Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
  int comm_size)
{
  Xt_int local_interval = (Xt_int)(global_interval / comm_size);
  Xt_int local_index_range_lbound
    = get_min_idxlist_index(src_idxlist, dst_idxlist);
  Xt_int local_index_range_ubound
    = get_max_idxlist_index(src_idxlist, dst_idxlist);
  size_t first_overlapping_bucket = 0;
  /* is it impossible for early buckets to overlap our lists? */
  if (local_index_range_lbound >= 0
      && (local_index_range_ubound < global_interval)) {
    first_overlapping_bucket
      = (size_t)(local_index_range_lbound / local_interval);
  }
  /* is it impossible for later ranks to overlap our lists? */
  size_t start_of_non_overlapping_bucket_suffix
    = (size_t)(((long long)local_index_range_ubound + local_interval - 1)
               / local_interval) + 1;
  if (local_index_range_lbound < 0
      || start_of_non_overlapping_bucket_suffix > (size_t)comm_size)
    start_of_non_overlapping_bucket_suffix = (size_t)comm_size;
  /*
   * size_t max_num_intersect
   *     = start_of_non_overlapping_bucket_suffix - first_overlapping_bucket;
   */
  return (struct bucket_params){
    .global_interval = global_interval,
    .local_interval = local_interval,
    .local_index_range_lbound = local_index_range_lbound,
    .local_index_range_ubound = local_index_range_ubound,
    .most_recent_rank_generated = (int)first_overlapping_bucket - 1,
    .max_rank_generated = (int)(start_of_non_overlapping_bucket_suffix-1),
  };
}

/* unfortunately GCC 11 cannot handle the literal constants used for
 * MPI_STATUSES_IGNORE by MPICH */
#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif

/* interval_size[0] and interval_size[1] are the global interval
 * size for the local and remote group */

static inline void
get_intercomm_dist_dir_global_interval_size(
  Xt_idxlist src, Xt_idxlist dst,
  bool *stripify, Xt_int interval_size[2],
  struct Xt_xmdd_bucket_gen_comms *comms,
  int comm_size, int remote_size,
  Xt_config config)
{
  /* global_sums[0] and [1] refer to the local and remote group of
   * intercommunicator inter_comm */
  unsigned long long local_vals[2], global_sums[2][2];

  unsigned num_indices_src = (unsigned)xt_idxlist_get_num_indices(src);
  local_vals[0] = num_indices_src;
  local_vals[1] = (num_indices_src > (unsigned)config->idxv_cnv_size)
    || (xt_idxlist_get_num_indices(dst) > config->idxv_cnv_size);

  xt_mpi_call(MPI_Allreduce(local_vals, global_sums[0], 2,
                            MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                            comms->intra_comm),
              comms->intra_comm);
  /* instead of sendrecv one might use hand-programmed multi-casts
   * sending to each rank in a range from the remote group and
   * receiving from the first rank in that group,
   * the better choice probably depends on the asymmetry of the group
   * sizes, i.e. use bcast from a very small to a very large group
   * and few sends from a large to a small group */
  int comm_rank;
  xt_mpi_call(MPI_Comm_rank(comms->inter_comm, &comm_rank),
              comms->inter_comm);
  if (comm_rank == 0) {
    int tag = comms->tag_offset_inter + xt_mpi_tag_xmap_dist_dir_src_send;
    xt_mpi_call(MPI_Sendrecv(global_sums[0], 2, MPI_UNSIGNED_LONG_LONG, 0, tag,
                             global_sums[1], 2, MPI_UNSIGNED_LONG_LONG, 0, tag,
                             comms->inter_comm, MPI_STATUS_IGNORE),
                comms->inter_comm);
  }
  xt_mpi_call(MPI_Bcast(global_sums[1], 2, MPI_UNSIGNED_LONG_LONG,
                        0, comms->intra_comm), comms->intra_comm);
  *stripify = (global_sums[0][1] > 0 || global_sums[1][1] > 0);
  interval_size[0]
    = (Xt_int)(((global_sums[0][0] + (unsigned)comm_size - 1)
                / (unsigned)comm_size) * (unsigned)comm_size);
  interval_size[1]
    = (Xt_int)(((global_sums[1][0] + (unsigned)remote_size - 1)
                / (unsigned)remote_size) * (unsigned)remote_size);
}

static inline Xt_int
get_intercomm_min_idxlist_index(Xt_idxlist l)
{
  int num_idx = xt_idxlist_get_num_indices(l);
  Xt_int min_index = num_idx ? xt_idxlist_get_min_index(l) : XT_INT_MAX;
  return min_index;
}

static inline Xt_int
get_intercomm_max_idxlist_index(Xt_idxlist l)
{
  int num_idx = xt_idxlist_get_num_indices(l);
  Xt_int max_index = num_idx ? xt_idxlist_get_max_index(l) : XT_INT_MIN;
  return max_index;
}


static struct bucket_params
get_intercomm_bucket_params(Xt_idxlist idxlist,
                            Xt_int global_interval, int comm_size)
{
  /* guard vs. comm_size being larger than number of indices */
  Xt_int local_interval = MAX((Xt_int)1, (Xt_int)(global_interval / comm_size));
  Xt_int local_index_range_lbound = get_intercomm_min_idxlist_index(idxlist);
  Xt_int local_index_range_ubound = get_intercomm_max_idxlist_index(idxlist);
  size_t first_overlapping_bucket = 0;
  /* is it impossible for early buckets to overlap our lists? */
  if (local_index_range_lbound >= 0
      && (local_index_range_ubound < global_interval)) {
    first_overlapping_bucket
      = (size_t)(local_index_range_lbound / local_interval);
  }
  /* is it impossible for later ranks to overlap our lists? */
  size_t start_of_non_overlapping_bucket_suffix
    = (size_t)(((long long)local_index_range_ubound + local_interval - 1)
               / local_interval) + 1;
  if (local_index_range_lbound < 0
      || start_of_non_overlapping_bucket_suffix > (size_t)comm_size)
    start_of_non_overlapping_bucket_suffix = (size_t)comm_size;
  return (struct bucket_params){
    .global_interval = global_interval,
    .local_interval = local_interval,
    .local_index_range_lbound = local_index_range_lbound,
    .local_index_range_ubound = local_index_range_ubound,
    .most_recent_rank_generated = (int)first_overlapping_bucket - 1,
    .max_rank_generated = (int)(start_of_non_overlapping_bucket_suffix-1),
  };
}

static int
xt_xmdd_bucket_gen_cycl_stripe_init(
  void *gen_state_,
  Xt_idxlist src_idxlist,
  Xt_idxlist dst_idxlist,
  Xt_config config,
  struct Xt_xmdd_bucket_gen_comms *comms,
  void *init_params)
{
  (void)init_params;
  struct Xt_xmdd_bucket_gen_cycl_stripe_state *gen_state = gen_state_;
  bool stripify = false;
  int comm_size;
  xt_mpi_call(MPI_Comm_size(comms->intra_comm, &comm_size), comms->intra_comm);
  if (comms->inter_comm == MPI_COMM_NULL) {
    /* global_interval is a multiple of comm_size and has a size of at least
       comm_size */
    Xt_int global_interval_size
      = get_intracomm_dist_dir_global_interval_size(
        src_idxlist, dst_idxlist,
        &stripify, comms->intra_comm, comm_size, config);
    gen_state->dst
      = get_intracomm_bucket_params(global_interval_size,
                                    src_idxlist, dst_idxlist, comm_size);
  } else {
    Xt_int global_interval_size[2];
    int remote_size;
    xt_mpi_call(MPI_Comm_remote_size(comms->inter_comm, &remote_size),
                comms->inter_comm);
    get_intercomm_dist_dir_global_interval_size(
      src_idxlist, dst_idxlist,
      &stripify, global_interval_size,
      comms, comm_size, remote_size, config);
    gen_state->src
      = get_intercomm_bucket_params(src_idxlist, global_interval_size[0],
                                    comm_size);
    gen_state->dst
      = get_intercomm_bucket_params(dst_idxlist, global_interval_size[1],
                                    remote_size);
  }
  gen_state->last_list_generated = NULL;
  gen_state->stripes = NULL;
  gen_state->stripes_array_size = 0;
  return stripify;
}

static void
xt_xmdd_bucket_gen_cycl_stripe_destroy(void *gen_state_)
{
  struct Xt_xmdd_bucket_gen_cycl_stripe_state *gen_state = gen_state_;
  if (gen_state->last_list_generated)
    xt_idxlist_delete(gen_state->last_list_generated);
  free(gen_state->stripes);
}

static int
xt_xmdd_bucket_gen_cycl_stripe_get_intersect_max_num(void *gen_state_, int type)
{
  int max_num_intersect = 0;
  struct Xt_xmdd_bucket_gen_cycl_stripe_state *gen_state = gen_state_;
  struct bucket_params *params = type == Xt_dist_dir_bucket_gen_type_send
    ? &gen_state->src : &gen_state->dst;
  int max_rank_generated = params->max_rank_generated,
    most_recent_rank_generated = params->most_recent_rank_generated;
  max_num_intersect =
    params->local_index_range_ubound
    >= params->local_index_range_lbound
    && max_rank_generated >= most_recent_rank_generated
    ? max_rank_generated - most_recent_rank_generated : 0;
  return max_num_intersect;
}


/**
 * \brief generates the buckets of the distributed directory
 *
 * The buckets of the distributed directory are computed as follows:
 * - compute sum of the sizes of all src_idxlist's and round to next
 *   multiple of comm_size -> global_interval
 * - local_interval_size = global_interval / comm_size
 * - local_interval_start = rank * local_interval_size
 * - local_interval_end = (rank + 1) * local_interval_size
 * - local_interval = $[local_interval_start,local_interval_end)$
 * - bucket of each rank is the set of stripes that is defined as:
 *      stripe[i].start = start + rank*local_interval_start + i*global_interval
 *      stripe[i].stride = 1;
 *      stripe[i].nstrides = (int)local_interval;
 *      with i in [0,num_stripes)
 *      bucket[rank] = xt_idxstripes_new(stripes, num_stripes);
 * - num_stripes and start are choosen such that:
 *      start + local_interval
 *         > get_min_idxlist_index(src_idxlist, dst_idxlist)
 *      (which implies for rank 0
 *       stripe[0].start <= get_min_idxlist_index(src_idxlist, dst_idxlist))
 *   where
 *       start = global_interval * k,
 *       k is minimal
 *   and
 *      start + num_stripes*global_interval
 *         > get_max_idxlist_index(src_idxlist, dst_idxlist);
 *   i.e. local_interval is replicated modulo global_interval for the
 *   whole range of indices.
 *
 * \
 * \param[in] bucket_params the parameters
 * \param[in,out] gen_state bucket generator state
 * \param[in] dist_dir_rank rank for which to compute bucket
 * \return newly created bucket and corresponding rank. The bucket must not be used after
 * either the next call to xt_xmap_dist_dir_get_bucket or
 * xt_xmdd_bucket_gen_cycl_stripe_destroy for the same gen_state.
 */
static Xt_idxlist
xt_xmap_dist_dir_get_bucket(struct Xt_xmdd_bucket_gen_cycl_stripe_state *gen_state,
                            struct bucket_params *bucket_params,
                            int dist_dir_rank)
{
  Xt_int global_interval = bucket_params->global_interval;
  Xt_int local_interval = bucket_params->local_interval;
  Xt_int local_index_range_lbound = bucket_params->local_index_range_lbound;
  Xt_int local_index_range_ubound = bucket_params->local_index_range_ubound;
  int num_stripes = 0;

  /* find first index in bucket of dist_dir_rank
     <= local_index_range_lbound */
  Xt_int start = (Xt_int)(0 + dist_dir_rank * local_interval);
  {
    long long start_correction
      = (long long)local_index_range_lbound - (long long)start;
    Xt_int corr_steps
      = (Xt_int)((start_correction
                  - (llsign_mask(start_correction)
                     & (long long)(global_interval - 1)))
                 / (long long)global_interval);
    start = (Xt_int)(start + corr_steps * global_interval);
  }
  /* next find last stripe in bucket of dist_dir_rank
   * <= local_index_range_ubound */
  Xt_int end
    = (Xt_int)(start
               + (((long long)local_index_range_ubound - (long long)start)
                  / global_interval) * global_interval);
  Xt_int use_start_stripe
    = (Xt_int)(start + local_interval > local_index_range_lbound);
  num_stripes = (int)(((long long)end - (long long)start)/global_interval)
    + (int)use_start_stripe;
  start = (Xt_int)(start
                   + ((Xt_int)((Xt_uint)use_start_stripe - (Xt_uint)1)
                      & global_interval));
  if (!num_stripes)
    return NULL;

  struct Xt_stripe *restrict stripes = gen_state->stripes;
  ENSURE_ARRAY_SIZE(stripes, gen_state->stripes_array_size, (size_t)num_stripes);
  for (int j = 0; j < num_stripes; ++j) {
    stripes[j].start = (Xt_int)(start + j * global_interval);
    stripes[j].stride = 1;
    stripes[j].nstrides = (int)local_interval;
  }
  gen_state->stripes = stripes;
  return xt_idxstripes_prealloc_new(stripes, num_stripes);
}

static struct Xt_com_list
xt_xmdd_cycl_stripe_get_next_bucket(void *gen_state_, int type)
{
  struct Xt_xmdd_bucket_gen_cycl_stripe_state *gen_state = gen_state_;
  struct Xt_com_list bucket;
  if (gen_state->last_list_generated)
    xt_idxlist_delete(gen_state->last_list_generated);
  struct bucket_params *params =
    type == Xt_dist_dir_bucket_gen_type_send ? &gen_state->src
    : &gen_state->dst;
  int next_rank = params->most_recent_rank_generated,
    max_rank = params->max_rank_generated;
  bucket.list = NULL;
  while (next_rank < max_rank
         && !(bucket.list = xt_xmap_dist_dir_get_bucket(
                gen_state, params, ++next_rank)))
    ;
  params->most_recent_rank_generated = bucket.rank = next_rank;
  gen_state->last_list_generated = bucket.list;
  return bucket;
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
