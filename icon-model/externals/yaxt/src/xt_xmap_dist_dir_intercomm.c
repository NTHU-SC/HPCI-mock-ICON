/**
 * @file xt_xmap_dist_dir_intercomm.c
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

#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>

#include <mpi.h>

#include "xt/xt_idxlist.h"
#include "xt/xt_idxlist_collection.h"
#include "xt/xt_idxvec.h"
#include "xt/xt_idxstripes.h"
#include "xt/xt_idxempty.h"
#include "xt/xt_xmap.h"
#include "xt/xt_xmap_dist_dir.h"
#include "xt/xt_xmap_dist_dir_intercomm.h"
#include "xt/xt_mpi.h"
#include "xt_arithmetic_util.h"
#include "xt_idxstripes_internal.h"
#include "xt_mpi_internal.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt/xt_xmap_intersection.h"
#include "xt_idxlist_internal.h"
#include "xt_xmap_dist_dir_common.h"
#include "xt_config_internal.h"
#include "instr.h"
#include "xt/xt_sort.h"
#include "xt_xmap_dist_dir_bucket_gen_internal.h"

enum {
  SEND_SIZE = 0,
  SEND_NUM = 1,
  SEND_SIZE_ASIZE,
};

static inline void
rank_no_send(size_t rank, int (*restrict send_size)[SEND_SIZE_ASIZE])
{
  send_size[rank][SEND_SIZE] = 0;
  send_size[rank][SEND_NUM] = 0;
}

struct mmsg_buf
{
  size_t num_msg;
  void *buffer;
};


static struct mmsg_buf
compute_and_pack_bucket_intersections(void *bucket_gen_state,
                                      int bucket_type,
                                      Xt_idxlist idxlist,
                                      int (*send_size)[SEND_SIZE_ASIZE],
                                      MPI_Comm comm, int comm_size,
                                      Xt_config config)
{
  int nosort_forced = XT_CONFIG_GET_FORCE_NOSORT(config);
  size_t send_size_filled = (size_t)-1;
  size_t max_num_intersect
    = (size_t)config->xmdd_bucket_gen->get_intersect_max_num(
      bucket_gen_state, bucket_type);
  struct mmsg_buf result = { 0, 0 };
  if (max_num_intersect) {
    Xt_idxlist idxlist_sorted
      = nosort_forced || xt_idxlist_get_sorting(idxlist) == 1
      ? idxlist
      : xt_idxlist_sorted_copy_custom(idxlist, config);

    size_t num_msg = 0;
    size_t send_buffer_size = 0;
    struct Xt_com_list *restrict sends
      = xmalloc(max_num_intersect * sizeof(*sends));
    struct Xt_com_list bucket;
    while ((bucket = config->xmdd_bucket_gen->next(
              bucket_gen_state, bucket_type)).list) {
      size_t rank;
      for (rank = send_size_filled + 1; rank < (size_t)bucket.rank; ++rank)
        rank_no_send(rank, send_size);

      Xt_idxlist isect2send
        = xt_idxlist_get_intersection(idxlist_sorted, bucket.list);
      if (xt_idxlist_get_num_indices(isect2send) > 0) {
        sends[num_msg].list = isect2send;
        sends[num_msg].rank = (int)rank;
        send_buffer_size += xt_idxlist_get_pack_size(isect2send, comm);
        /* send_size[rank][SEND_SIZE] is set below after the actual
         * pack, because MPI_Pack_size only gives an upper bound,
         * not the actually needed size */
        send_size[rank][SEND_NUM] = 1;
        ++num_msg;
      } else {
        rank_no_send(rank, send_size);
        xt_idxlist_delete(isect2send);
      }
    }
    if (idxlist_sorted != idxlist)
      xt_idxlist_delete(idxlist_sorted);
    for (size_t rank = send_size_filled+1; rank < (size_t)comm_size; ++rank)
      rank_no_send(rank, send_size);

    unsigned char *send_buffer = xmalloc(send_buffer_size);
    size_t ofs = 0;
    for (size_t i = 0; i < num_msg; ++i) {
      int position = 0;
      xt_idxlist_pack(sends[i].list, send_buffer + ofs,
                      (int)(send_buffer_size-ofs), &position, comm);
      send_size[sends[i].rank][SEND_SIZE] = position;
      ofs += (size_t)position;
      xt_idxlist_delete(sends[i].list);
    }

    free(sends);
    result.num_msg = num_msg;
    result.buffer = send_buffer;
  } else {
    memset(send_size, 0, (size_t)comm_size * sizeof (*send_size));
    result.num_msg = 0;
    result.buffer = NULL;
  }

  return result;
}


static void
compress_sizes(int (*restrict sizes)[SEND_SIZE_ASIZE], int comm_size,
               struct Xt_xmdd_txstat *tx_stat, int *counts)
{
  size_t tx_num = 0, size_sum = 0;
  for (size_t i = 0; i < (size_t)comm_size; ++i)
    if (sizes[i][SEND_SIZE]) {
      int tx_size = sizes[i][SEND_SIZE];
      size_sum += (size_t)tx_size;
      sizes[tx_num][SEND_SIZE] = tx_size;
      if (counts) counts[tx_num] = sizes[i][SEND_NUM];
      sizes[tx_num][SEND_NUM] = (int)i;
      ++tx_num;
    }
  *tx_stat = (struct Xt_xmdd_txstat){ .bytes = size_sum, .num_msg = tx_num };
}

static void *
create_intersections(void *bucket_gen_state,
                     int bucket_type,
                     struct Xt_xmdd_txstat tx_stat[2],
                     int recv_size[][SEND_SIZE_ASIZE],
                     int send_size[][SEND_SIZE_ASIZE],
                     Xt_idxlist idxlist,
                     MPI_Comm comm, int comm_size, Xt_config config)
{
  struct mmsg_buf ddr
    = compute_and_pack_bucket_intersections(
      bucket_gen_state, bucket_type, idxlist,
      send_size, comm, comm_size, config);
  xt_mpi_call(MPI_Alltoall((int *)send_size, SEND_SIZE_ASIZE, MPI_INT,
                           (int *)recv_size, SEND_SIZE_ASIZE, MPI_INT, comm),
              comm);
  compress_sizes(recv_size, comm_size, tx_stat + 0, NULL);
  compress_sizes(send_size, comm_size, tx_stat + 1, NULL);
  assert(ddr.num_msg == tx_stat[1].num_msg);
  return ddr.buffer;
}

typedef int (*tx_fp)(void *, int, MPI_Datatype, int, int,
                     MPI_Comm, MPI_Request *);
static void
tx_intersections(size_t num_msg,
                 const int (*sizes)[SEND_SIZE_ASIZE],
                 unsigned char *buffer, MPI_Request *requests,
                 int tag, MPI_Comm comm, tx_fp tx_op)
{
  size_t ofs = 0;
  for (size_t i = 0; i < num_msg; ++i)
  {
    int rank = sizes[i][SEND_NUM], count = sizes[i][SEND_SIZE];
    xt_mpi_call(tx_op(buffer + ofs,
                      count, MPI_PACKED, rank, tag, comm, requests + i), comm);
    ofs += (size_t)count;
  }
}

static void
irecv_intersections(size_t num_msg,
                    const int (*recv_size)[SEND_SIZE_ASIZE],
                    void *recv_buffer, MPI_Request *requests,
                    int tag, MPI_Comm comm)
{
  tx_intersections(num_msg, recv_size, recv_buffer, requests, tag, comm,
                   (tx_fp)MPI_Irecv);
}

static void
isend_intersections(size_t num_msg,
                    const int (*send_size)[SEND_SIZE_ASIZE],
                    void *send_buffer, MPI_Request *requests,
                    int tag, MPI_Comm comm)
{
  tx_intersections(num_msg, send_size, send_buffer, requests, tag, comm,
                   (tx_fp)MPI_Isend);
}


static struct dist_dir *
unpack_dist_dir(struct Xt_xmdd_txstat tx_stat,
                const int (*sizes)[SEND_SIZE_ASIZE],
                void *buffer,
                MPI_Comm comm)
{
  size_t num_msg = tx_stat.num_msg, buf_size = tx_stat.bytes;
  struct dist_dir *restrict dist_dir
    = xmalloc(sizeof (*dist_dir) + sizeof (*dist_dir->entries) * num_msg);
  dist_dir->num_entries = (int)num_msg;
  int position = 0;
  for (size_t i = 0; i < num_msg; ++i)
  {
    int rank = sizes[i][SEND_NUM];
    dist_dir->entries[i].rank = rank;
    dist_dir->entries[i].list
      = xt_idxlist_unpack(buffer, (int)buf_size, &position, comm);
  }
  return dist_dir;
}

struct dd_result {
  struct dist_dir_pair dist_dirs;
  int stripify;
};

/* unfortunately GCC 11 cannot handle the literal constants used for
 * MPI_STATUSES_IGNORE by MPICH */
#if __GNUC__ >= 11 && __GNUC__ <= 13 && defined MPICH
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wstringop-overread"
#pragma GCC diagnostic ignored "-Wstringop-overflow"
#endif

static struct dd_result
generate_distributed_directories(
  Xt_idxlist src_idxlist,
  Xt_idxlist dst_idxlist,
  const struct Xt_xmdd_bucket_gen_comms *comms,
  int remote_size, int comm_size,
  Xt_config config) {

  size_t bgd_size = config->xmdd_bucket_gen->gen_state_size;
  bgd_size = (bgd_size + sizeof (void *) - 1)/sizeof (void *) * sizeof (void *);
  void *bgd[bgd_size];
  struct dd_result results;
  results.stripify
    = config->xmdd_bucket_gen->init(
      &bgd, src_idxlist, dst_idxlist, config, comms, NULL,
      config->xmdd_bucket_gen);
  int (*send_size_local)[SEND_SIZE_ASIZE]
    = xmalloc(((size_t)comm_size + (size_t)remote_size)
              * 2 * sizeof(*send_size_local)),
    (*send_size_remote)[SEND_SIZE_ASIZE] = send_size_local + comm_size,
    (*recv_size_local)[SEND_SIZE_ASIZE] = send_size_remote + remote_size,
    (*recv_size_remote)[SEND_SIZE_ASIZE] = recv_size_local + comm_size;
  struct Xt_xmdd_txstat tx_stat_local[2], tx_stat_remote[2];
  void *send_buffer_local
    = create_intersections(bgd, Xt_dist_dir_bucket_gen_type_send,
                           tx_stat_local, recv_size_local,
                           send_size_local, src_idxlist,
                           comms->intra_comm, comm_size, config);
  void *send_buffer_remote
    = create_intersections(bgd, Xt_dist_dir_bucket_gen_type_recv,
                           tx_stat_remote, recv_size_remote,
                           send_size_remote, dst_idxlist,
                           comms->inter_comm, remote_size, config);
  XT_CONFIG_BUCKET_DESTROY(config, &bgd);
  size_t num_req = tx_stat_local[0].num_msg + tx_stat_remote[0].num_msg
    + tx_stat_local[1].num_msg + tx_stat_remote[1].num_msg;
  MPI_Request *dir_init_requests
    = xmalloc(num_req * sizeof(*dir_init_requests)
              + tx_stat_local[0].bytes + tx_stat_remote[0].bytes);
  void *recv_buffer_local = dir_init_requests + num_req,
    *recv_buffer_remote = ((unsigned char *)recv_buffer_local
                           + tx_stat_local[0].bytes);
  int tag_intra = comms->tag_offset_intra
    + xt_mpi_tag_xmap_dist_dir_src_send;
  size_t req_ofs = tx_stat_local[0].num_msg;
  irecv_intersections(tx_stat_local[0].num_msg,
                      (const int (*)[SEND_SIZE_ASIZE])recv_size_local,
                      recv_buffer_local, dir_init_requests,
                      tag_intra, comms->intra_comm);
  int tag_inter = comms->tag_offset_inter
    + xt_mpi_tag_xmap_dist_dir_src_send;
  irecv_intersections(tx_stat_remote[0].num_msg,
                      (const int (*)[SEND_SIZE_ASIZE])recv_size_remote,
                      recv_buffer_remote, dir_init_requests + req_ofs,
                      tag_inter, comms->inter_comm);
  req_ofs += tx_stat_remote[0].num_msg;
  isend_intersections(tx_stat_local[1].num_msg,
                      (const int (*)[SEND_SIZE_ASIZE])send_size_local,
                      send_buffer_local, dir_init_requests + req_ofs,
                      tag_intra, comms->intra_comm);
  req_ofs += tx_stat_local[1].num_msg;
  isend_intersections(tx_stat_remote[1].num_msg,
                      (const int (*)[SEND_SIZE_ASIZE])send_size_remote,
                      send_buffer_remote, dir_init_requests + req_ofs,
                      tag_inter, comms->inter_comm);
  // wait for data transfers to complete
  xt_mpi_call(MPI_Waitall((int)num_req, dir_init_requests,
                          MPI_STATUSES_IGNORE), comms->inter_comm);
  free(send_buffer_local);
  free(send_buffer_remote);
  results.dist_dirs.src
    = unpack_dist_dir(tx_stat_local[0],
                    (const int (*)[SEND_SIZE_ASIZE])recv_size_local,
                    recv_buffer_local, comms->intra_comm);
  results.dist_dirs.dst
    = unpack_dist_dir(tx_stat_remote[0],
                      (const int (*)[SEND_SIZE_ASIZE])recv_size_remote,
                      recv_buffer_remote, comms->inter_comm);
  free(send_size_local);
  free(dir_init_requests);
  return results;
}


static size_t
send_size_from_intersections(size_t num_intersections,
                             const struct isect *restrict src_dst_intersections,
                             enum xt_xmdd_direction target,
                             MPI_Comm comm, int comm_size,
                             int (*restrict send_size_target)[SEND_SIZE_ASIZE])
{
  size_t total_send_size = 0;
  for (int i = 0; i < comm_size; ++i)
    (void)(send_size_target[i][SEND_SIZE] = 0),
      (void)(send_size_target[i][SEND_NUM] = 0);

  int rank_pack_size;
  xt_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &rank_pack_size), comm);

  for (size_t i = 0; i < num_intersections; ++i)
  {
    size_t msg_size = (size_t)rank_pack_size
      + xt_idxlist_get_pack_size(src_dst_intersections[i].idxlist, comm);
    size_t target_rank = (size_t)src_dst_intersections[i].rank[target];
    /* send_size_target[target_rank][SEND_SIZE] += msg_size; */
    ++(send_size_target[target_rank][SEND_NUM]);
    total_send_size += msg_size;
  }
  assert(total_send_size <= INT_MAX);
  return total_send_size;
}

static struct mmsg_buf
pack_dist_dirs(size_t num_intersections,
               struct isect *restrict src_dst_intersections,
               int (*send_size)[SEND_SIZE_ASIZE],
               enum xt_xmdd_direction target,
               bool isect_idxlist_delete, MPI_Comm comm, int comm_size) {

  size_t total_send_size
    = send_size_from_intersections(num_intersections,
                                   src_dst_intersections,
                                   target,
                                   comm, comm_size, send_size);

  unsigned char *send_buffer = xmalloc(total_send_size);
  qsort(src_dst_intersections, num_intersections,
        sizeof (src_dst_intersections[0]),
        target == xt_xmdd_direction_src
        ? xt_xmdd_cmp_isect_src_rank : xt_xmdd_cmp_isect_dst_rank);
  size_t ofs = 0;
  size_t num_requests
    = xt_xmap_dist_dir_pack_intersections(
      target, num_intersections, src_dst_intersections,
      isect_idxlist_delete,
      SEND_SIZE_ASIZE, SEND_SIZE, send_size,
      send_buffer, total_send_size, &ofs, comm);
  return (struct mmsg_buf){ .num_msg = num_requests,
                              .buffer = send_buffer };
}

static struct dist_dir *
unpack_dist_dir_results(struct Xt_xmdd_txstat tx_stat,
                        void *restrict recv_buffer,
                        int *restrict entry_counts,
                        MPI_Comm comm)
{
  size_t num_msg = tx_stat.num_msg;
  int buf_size = (int)tx_stat.bytes;
  int position = 0;
  size_t num_entries_sent = 0;
  for (size_t i = 0; i < num_msg; ++i)
    num_entries_sent += (size_t)entry_counts[i];
  struct dist_dir *dist_dir
    = xmalloc(sizeof (struct dist_dir)
              + (sizeof (struct Xt_com_list) * num_entries_sent));
  dist_dir->num_entries = (int)num_entries_sent;
  struct Xt_com_list *restrict entries = dist_dir->entries;
  size_t num_entries = 0;
  for (size_t i = 0; i < num_msg; ++i) {
    size_t num_entries_from_rank = (size_t)entry_counts[i];
    for (size_t j = 0; j < num_entries_from_rank; ++j) {
      xt_mpi_call(MPI_Unpack(recv_buffer, buf_size, &position,
                             &entries[num_entries].rank,
                             1, MPI_INT, comm), comm);
      entries[num_entries].list =
        xt_idxlist_unpack(recv_buffer, buf_size, &position, comm);
      ++num_entries;
    }
  }
  assert(num_entries == num_entries_sent);
  qsort(entries, num_entries_sent, sizeof(*entries), xt_com_list_rank_cmp);
  xt_xmap_dist_dir_same_rank_merge(&dist_dir);
  return dist_dir;
}


static struct dd_result
exchange_idxlists(Xt_idxlist src_idxlist,
                  Xt_idxlist dst_idxlist,
                  const struct Xt_xmdd_bucket_gen_comms *comms,
                  Xt_config config) {

  int comm_size, remote_size;
  xt_mpi_call(MPI_Comm_size(comms->inter_comm, &comm_size),
              comms->inter_comm);
  xt_mpi_call(MPI_Comm_remote_size(comms->inter_comm, &remote_size),
              comms->inter_comm);

  struct dd_result bucket_isects
    = generate_distributed_directories(src_idxlist, dst_idxlist, comms,
                                       remote_size, comm_size,
                                       config);


  int (*send_size_local)[SEND_SIZE_ASIZE]
    = xmalloc(((size_t)comm_size + (size_t)remote_size)
              * 2U * sizeof(*send_size_local)),
    (*recv_size_local)[SEND_SIZE_ASIZE] = send_size_local + comm_size,
    (*send_size_remote)[SEND_SIZE_ASIZE] = recv_size_local + comm_size,
    (*recv_size_remote)[SEND_SIZE_ASIZE] = send_size_remote + remote_size;

  /* match the source and destination entries in the local distributed
   * directories... */
  struct isect *src_dst_intersections;
  size_t num_intersections
    = xt_xmap_dist_dir_match_src_dst(bucket_isects.dist_dirs.src,
                                     bucket_isects.dist_dirs.dst,
                                     &src_dst_intersections, config);
  xt_xmdd_free_dist_dirs(bucket_isects.dist_dirs);
  /* ... and pack the results into a sendable format */
  struct mmsg_buf dd_local, dd_remote;
  dd_local
    = pack_dist_dirs(num_intersections, src_dst_intersections, send_size_local,
                     xt_xmdd_direction_src, false, comms->intra_comm,
                     comm_size);
  dd_remote
    = pack_dist_dirs(num_intersections, src_dst_intersections, send_size_remote,
                     xt_xmdd_direction_dst, true, comms->inter_comm,
                     remote_size);
  free(src_dst_intersections);

  // get the data size the local process will receive from other processes
  xt_mpi_call(MPI_Alltoall((int *)send_size_local, SEND_SIZE_ASIZE, MPI_INT,
                           (int *)recv_size_local, SEND_SIZE_ASIZE, MPI_INT,
                           comms->intra_comm), comms->intra_comm);
  xt_mpi_call(MPI_Alltoall((int *)send_size_remote, SEND_SIZE_ASIZE, MPI_INT,
                           (int *)recv_size_remote, SEND_SIZE_ASIZE, MPI_INT,
                           comms->inter_comm), comms->inter_comm);

  struct Xt_xmdd_txstat tx_stat_local[2], tx_stat_remote[2];
  int *isect_counts_recv_local
    = xmalloc(((size_t)comm_size + (size_t)remote_size) * sizeof (int)),
    *isect_counts_recv_remote = isect_counts_recv_local + comm_size;
  compress_sizes(send_size_local, comm_size, tx_stat_local+1, NULL);
  compress_sizes(recv_size_local, comm_size, tx_stat_local+0,
                 isect_counts_recv_local);
  compress_sizes(send_size_remote, remote_size, tx_stat_remote+1, NULL);
  compress_sizes(recv_size_remote, remote_size, tx_stat_remote+0,
                 isect_counts_recv_remote);
  assert(tx_stat_local[1].num_msg == dd_local.num_msg
         && tx_stat_remote[1].num_msg == dd_remote.num_msg);
  size_t num_requests
    = dd_local.num_msg + dd_remote.num_msg
    + tx_stat_local[0].num_msg + tx_stat_remote[0].num_msg;
  assert(num_requests <= INT_MAX);
  MPI_Request *requests
    = xmalloc(num_requests * sizeof(*requests)
              + tx_stat_local[0].bytes + tx_stat_remote[0].bytes);
  void *recv_buf_local = requests + num_requests,
    *recv_buf_remote = (unsigned char *)recv_buf_local + tx_stat_local[0].bytes;
  size_t req_ofs = tx_stat_local[0].num_msg;
  int tag_intra = comms->tag_offset_intra + xt_mpi_tag_xmap_dist_dir_src_send;
  irecv_intersections(tx_stat_local[0].num_msg,
                      (const int (*)[SEND_SIZE_ASIZE])recv_size_local,
                      recv_buf_local, requests, tag_intra, comms->intra_comm);
  int tag_inter = comms->tag_offset_inter + xt_mpi_tag_xmap_dist_dir_src_send;
  irecv_intersections(tx_stat_remote[0].num_msg,
                      (const int (*)[SEND_SIZE_ASIZE])recv_size_remote,
                      recv_buf_remote, requests+req_ofs, tag_inter,
                      comms->inter_comm);
  req_ofs += tx_stat_remote[0].num_msg;
  isend_intersections(tx_stat_local[1].num_msg,
                      (const int (*)[SEND_SIZE_ASIZE])send_size_local,
                      dd_local.buffer, requests+req_ofs, tag_intra,
                      comms->intra_comm);
  req_ofs += tx_stat_local[1].num_msg;
  isend_intersections(tx_stat_remote[1].num_msg,
                      (const int (*)[SEND_SIZE_ASIZE])send_size_remote,
                      dd_remote.buffer, requests+req_ofs, tag_inter,
                      comms->inter_comm);
  xt_mpi_call(MPI_Waitall((int)num_requests, requests, MPI_STATUSES_IGNORE),
              comms->inter_comm);
  free(dd_local.buffer);
  free(dd_remote.buffer);
  free(send_size_local);

  struct dd_result results;
  results.stripify = bucket_isects.stripify;
  results.dist_dirs.src
    = unpack_dist_dir_results(tx_stat_local[0], recv_buf_local,
                              isect_counts_recv_local, comms->intra_comm);
  results.dist_dirs.dst
    = unpack_dist_dir_results(tx_stat_remote[0], recv_buf_remote,
                              isect_counts_recv_remote, comms->inter_comm);
  free(requests);
  free(isect_counts_recv_local);
  return results;
}



Xt_xmap
xt_xmap_dist_dir_intercomm_custom_new(Xt_idxlist src_idxlist,
                                      Xt_idxlist dst_idxlist,
                                      MPI_Comm inter_comm,
                                      MPI_Comm intra_comm,
                                      Xt_config config)
{
  INSTR_DEF(this_instr,"xt_xmap_dist_dir_intercomm_new")
  INSTR_START(this_instr);

  // ensure that yaxt is initialized
  assert(xt_initialized());

  struct Xt_xmdd_bucket_gen_comms comms;
  comms.intra_comm = xt_mpi_comm_smart_dup(intra_comm,
                                           &comms.tag_offset_intra);
  comms.inter_comm = xt_mpi_comm_smart_dup(inter_comm,
                                           &comms.tag_offset_inter);
  struct dd_result results
    = exchange_idxlists(src_idxlist, dst_idxlist,
                        &comms, config);

  int stripify = XT_CONFIG_GET_XMAP_STRIPING(config);
  if (stripify == 2)
    stripify = results.stripify;
  Xt_xmap (*xmap_new)(int num_src_intersections,
                      const struct Xt_com_list *src_com,
                      int num_dst_intersections,
                      const struct Xt_com_list *dst_com,
                      Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
                      MPI_Comm comm)
    = stripify ? xt_xmap_intersection_ext_new : xt_xmap_intersection_new;

  Xt_xmap xmap
    = xmap_new(results.dist_dirs.src->num_entries,
               results.dist_dirs.src->entries,
               results.dist_dirs.dst->num_entries,
               results.dist_dirs.dst->entries,
               src_idxlist, dst_idxlist, comms.inter_comm);

  xt_mpi_comm_smart_dedup(&comms.inter_comm, comms.tag_offset_inter);
  xt_mpi_comm_smart_dedup(&comms.intra_comm, comms.tag_offset_intra);

  xt_xmdd_free_dist_dirs(results.dist_dirs);
  INSTR_STOP(this_instr);
  return xmap;
}

Xt_xmap
xt_xmap_dist_dir_intercomm_new(Xt_idxlist src_idxlist, Xt_idxlist dst_idxlist,
                               MPI_Comm inter_comm, MPI_Comm intra_comm)
{
  return xt_xmap_dist_dir_intercomm_custom_new(
    src_idxlist, dst_idxlist, inter_comm, intra_comm, &xt_default_config);
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
