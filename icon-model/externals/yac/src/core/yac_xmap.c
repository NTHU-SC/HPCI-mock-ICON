// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <yaxt.h>

#include "yac_xmap.h"
#include "ppm/ppm_xfuncs.h"
#include "yac_assert.h"
#include "yac_mpi_internal.h"

struct yac_xmap_msg {
  int count;
  const int * array_of_blocklengths;
  const int * array_of_displacements;
  int rank;
};

struct yac_xmap_ {
  MPI_Comm comm;
  int num_send_msg;
  int num_recv_msg;
  struct yac_xmap_msg msgs[];
};

static struct yac_xmap_msg parse_transfer_pos(
  int const * transfer_pos, size_t num_transfer_pos) {

  int count = 0;
  int * blocklengths = xmalloc(num_transfer_pos * sizeof(*blocklengths));
  int * displacements = xmalloc(num_transfer_pos * sizeof(*displacements));

  int prev_transfer_pos = transfer_pos[0];
  for (size_t i = 0; i < num_transfer_pos; ++i) {
    int curr_transfer_pos = transfer_pos[i];
    if (curr_transfer_pos != prev_transfer_pos + 1) {
      blocklengths[count] = 1;
      displacements[count] = curr_transfer_pos;
      YAC_ASSERT(
        count != INT_MAX,
        "ERROR(parse_transfer_pos): number of elements exceeds INT_MAX");
      ++count;
    } else {
      blocklengths[count-1]++;
    }
    prev_transfer_pos = curr_transfer_pos;
  }

  return
    (struct yac_xmap_msg)
      {.count = count,
       .array_of_blocklengths =
          xrealloc(blocklengths, (size_t)count * sizeof(*blocklengths)),
       .array_of_displacements =
          xrealloc(displacements, (size_t)count * sizeof(*displacements)),
       .rank = -1};
}

yac_xmap yac_xmap_from_point_infos(
  struct remote_point_infos * point_infos, size_t count, MPI_Comm comm) {

  int comm_size;
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  size_t * sendcounts, * recvcounts, * sdispls, * rdispls;
  yac_get_comm_buffers(
    1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);

  for (size_t i = 0; i < count; ++i) {
    struct remote_point_info * curr_info =
      (point_infos[i].count > 1)?
        (point_infos[i].data.multi):
        (&(point_infos[i].data.single));
    sendcounts[curr_info->rank]++;
  }

  yac_generate_alltoallv_args(
    1, sendcounts, recvcounts, sdispls, rdispls, comm);
  size_t num_src_msg = 0, num_dst_msg = 0;
  for (int i = 0; i < comm_size; ++i) {
    num_src_msg += (recvcounts[i] > 0);
    num_dst_msg += (sendcounts[i] > 0);
  }

  size_t recv_count =
    rdispls[comm_size-1] + recvcounts[comm_size-1];

  int * pos_buffer =
    xmalloc((recv_count + 2 * count) * sizeof(*pos_buffer));
  int * src_pos_buffer = pos_buffer;
  int * dst_pos_buffer = pos_buffer + recv_count;
  int * send_pos_buffer = pos_buffer + recv_count + count;

  // pack send buffer
  for (size_t i = 0; i < count; ++i) {
    struct remote_point_info * curr_info =
      (point_infos[i].count > 1)?
        (point_infos[i].data.multi):
        (&(point_infos[i].data.single));
    size_t pos = sdispls[curr_info->rank+1]++;
    dst_pos_buffer[pos] = i;
    send_pos_buffer[pos] = (int)(curr_info->orig_pos);
  }

  // redistribute positions of requested data
  yac_alltoallv_int_p2p(
    send_pos_buffer, sendcounts, sdispls,
    src_pos_buffer, recvcounts, rdispls, comm,
    "yac_xmap_from_point_infos", __LINE__);

  yac_xmap xmap =
    xmalloc(
      sizeof(*xmap) +
      ((size_t)(num_src_msg + num_dst_msg)) * sizeof(xmap->msgs[0]));

  xmap->comm = comm;
  xmap->num_send_msg = num_src_msg;
  xmap->num_recv_msg = num_dst_msg;
  struct yac_xmap_msg * send_msgs = xmap->msgs;
  struct yac_xmap_msg * recv_msgs = &(xmap->msgs[num_src_msg]);

  // set transfer_pos pointers and transfer_pos counts in com_pos's
  num_src_msg = 0;
  num_dst_msg = 0;
  for (int i = 0; i < comm_size; ++i) {
    if (recvcounts[i] > 0) {
      send_msgs[num_src_msg] =
        parse_transfer_pos(src_pos_buffer, recvcounts[i]);
      send_msgs[num_src_msg].rank = i;
      src_pos_buffer += recvcounts[i];
      ++num_src_msg;
    }
    if (sendcounts[i] > 0) {
      recv_msgs[num_dst_msg] =
        parse_transfer_pos(dst_pos_buffer, sendcounts[i]);
      recv_msgs[num_dst_msg].rank = i;
      dst_pos_buffer += sendcounts[i];
      ++num_dst_msg;
    }
  }
  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);

  free(pos_buffer);

  return xmap;
}

static void xt_redist_msg_free(
  struct Xt_redist_msg * msgs, size_t count, MPI_Comm comm) {
  for (size_t i = 0; i < count; ++i) {
    MPI_Datatype * dt = &(msgs[i].datatype);
    if (*dt != MPI_DATATYPE_NULL) yac_mpi_call(MPI_Type_free(dt), comm);
  }
  free(msgs);
}

Xt_redist yac_xmap_generate_redist(yac_xmap xmap, MPI_Datatype base_type) {

  size_t total_num_msg = (size_t)(xmap->num_send_msg + xmap->num_recv_msg);

  struct Xt_redist_msg * msgs_buffer =
    xmalloc(total_num_msg * sizeof(*msgs_buffer));

  for (size_t i = 0; i < total_num_msg; ++i) {

    MPI_Datatype msg_datatype;
    yac_mpi_call(
      MPI_Type_indexed(
        xmap->msgs[i].count,
        xmap->msgs[i].array_of_blocklengths,
        xmap->msgs[i].array_of_displacements, base_type, &msg_datatype),
      xmap->comm);
    xt_mpi_call(MPI_Type_commit(&msg_datatype), xmap->comm);

    msgs_buffer[i].rank = xmap->msgs[i].rank;
    msgs_buffer[i].datatype = msg_datatype;
  }

  Xt_config redist_config = xt_config_new();

  Xt_redist redist =
    xt_redist_single_array_base_custom_new(
      xmap->num_send_msg, xmap->num_recv_msg,
      msgs_buffer, msgs_buffer + xmap->num_send_msg, xmap->comm, redist_config);

  xt_config_delete(redist_config);
  xt_redist_msg_free(msgs_buffer, total_num_msg, xmap->comm);

  return redist;
}

void yac_xmap_delete(yac_xmap xmap) {

  for (int i = 0; i < (xmap->num_send_msg + xmap->num_recv_msg); ++i) {
    free((void*)(xmap->msgs[i].array_of_blocklengths));
    free((void*)(xmap->msgs[i].array_of_displacements));
  }
  free(xmap);
}
