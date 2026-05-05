// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include "string.h"
#include "proc_sphere_part.h"
#include "remote_point.h"
#include "ensure_array_size.h"
#include "yac_mpi_internal.h"
#ifdef SCOREP_USER_ENABLE
#include "scorep/SCOREP_User.h"
#endif // SCOREP_USER_ENABLE

enum splicomm_tags {
  DATA_SIZE_TAG,
  DATA_TAG,
};

// WARNING: before changing this datatype ensure that the MPI datatype created
// for this still matches its data layout
struct dist_vertex {
  double coord[3];
  int grid_idx;
  yac_int global_id;
  int num_owners;
  size_t owner_offset;
};

struct comm_buffers {
  int * sendcounts;
  int * recvcounts;
  int * sdispls;
  int * rdispls;
};

struct proc_sphere_part_node;

struct proc_sphere_part_node_data {
  union {
    struct proc_sphere_part_node * node;
    int rank;
  } data;
  int is_leaf;
};
struct proc_sphere_part_node {
  struct proc_sphere_part_node_data U, T;
  double gc_norm_vector[3];
};

struct neigh_search_data {
  int * ranks;
  int num_ranks;
  struct proc_sphere_part_node * node;
};

static int get_proc_sphere_part_node_data_pack_size(
  struct proc_sphere_part_node_data node_data, MPI_Comm comm);

static int get_proc_sphere_part_node_pack_size(
  struct proc_sphere_part_node node, MPI_Comm comm) {

  int vec_pack_size;
  yac_mpi_call(
    MPI_Pack_size(3, MPI_DOUBLE, comm, &vec_pack_size), comm);

  return vec_pack_size +
         get_proc_sphere_part_node_data_pack_size(node.U, comm) +
         get_proc_sphere_part_node_data_pack_size(node.T, comm);
}

static int get_proc_sphere_part_node_data_pack_size(
  struct proc_sphere_part_node_data node_data, MPI_Comm comm) {

  int int_pack_size;
  yac_mpi_call(MPI_Pack_size(1, MPI_INT, comm, &int_pack_size), comm);

  int data_size = int_pack_size;
  if (node_data.is_leaf)
    data_size += int_pack_size;
  else
    data_size +=
      get_proc_sphere_part_node_pack_size((*node_data.data.node), comm);

  return data_size;
}

static void pack_proc_sphere_part_node_data(
  struct proc_sphere_part_node_data node_data, void * pack_buffer,
  int pack_buffer_size, int * position, MPI_Comm comm);

static void pack_proc_sphere_part_node(
  struct proc_sphere_part_node * node, void * pack_buffer,
  int pack_buffer_size, int * position, MPI_Comm comm) {

  yac_mpi_call(MPI_Pack(&(node->gc_norm_vector[0]), 3, MPI_DOUBLE, pack_buffer,
                        pack_buffer_size, position, comm), comm);
  pack_proc_sphere_part_node_data(
    node->U, pack_buffer, pack_buffer_size, position, comm);
  pack_proc_sphere_part_node_data(
    node->T, pack_buffer, pack_buffer_size, position, comm);
}

static void pack_proc_sphere_part_node_data(
  struct proc_sphere_part_node_data node_data, void * pack_buffer,
  int pack_buffer_size, int * position, MPI_Comm comm) {

  yac_mpi_call(MPI_Pack(&(node_data.is_leaf), 1, MPI_INT, pack_buffer,
                        pack_buffer_size, position, comm), comm);

  if (node_data.is_leaf)
    yac_mpi_call(MPI_Pack(&(node_data.data.rank), 1, MPI_INT, pack_buffer,
                          pack_buffer_size, position, comm), comm);
  else
    pack_proc_sphere_part_node(node_data.data.node, pack_buffer,
                               pack_buffer_size, position, comm);
}

static struct proc_sphere_part_node_data unpack_proc_sphere_part_node_data(
  void * pack_buffer, int pack_buffer_size, int * position,
  MPI_Comm comm);

static struct proc_sphere_part_node * unpack_proc_sphere_part_node(
  void * pack_buffer, int pack_buffer_size, int * position,
  MPI_Comm comm) {

  struct proc_sphere_part_node * node = xmalloc(1 * sizeof(*node));

  yac_mpi_call(
    MPI_Unpack(pack_buffer, pack_buffer_size, position,
               &(node->gc_norm_vector[0]), 3, MPI_DOUBLE, comm), comm);

  node->U = unpack_proc_sphere_part_node_data(pack_buffer, pack_buffer_size,
                                              position, comm);
  node->T = unpack_proc_sphere_part_node_data(pack_buffer, pack_buffer_size,
                                              position, comm);

  return node;
}

static struct proc_sphere_part_node_data unpack_proc_sphere_part_node_data(
  void * pack_buffer, int pack_buffer_size, int * position,
  MPI_Comm comm) {

  struct proc_sphere_part_node_data node_data;

  yac_mpi_call(
    MPI_Unpack(pack_buffer, pack_buffer_size, position,
               &(node_data.is_leaf), 1, MPI_INT, comm), comm);

  if (node_data.is_leaf)
    yac_mpi_call(
      MPI_Unpack(pack_buffer, pack_buffer_size, position,
                 &(node_data.data.rank), 1, MPI_INT, comm), comm);
  else
    node_data.data.node =
      unpack_proc_sphere_part_node(
        pack_buffer, pack_buffer_size, position, comm);

  return node_data;
}

static struct proc_sphere_part_node_data get_remote_data(
  struct proc_sphere_part_node_data local_data,
  struct yac_group_comm local_group_comm,
  struct yac_group_comm remote_group_comm) {

  int comm_rank;
  MPI_Comm comm = local_group_comm.comm;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);

  int data_size;
  void * recv_buffer = NULL;

  int order = local_group_comm.start < remote_group_comm.start;

  for (int i = 0; i < 2; ++i) {

    if (order == i) {

      yac_bcast_group(
        &data_size, 1, MPI_INT, remote_group_comm.start, local_group_comm);
      recv_buffer = xmalloc((size_t)data_size);
      yac_bcast_group(recv_buffer, data_size, MPI_PACKED,
                      remote_group_comm.start, local_group_comm);

    } else {
      if (comm_rank == local_group_comm.start) {

        // pack local_data
        int pack_buffer_size =
          get_proc_sphere_part_node_data_pack_size(local_data, comm);
        void * pack_buffer = xmalloc((size_t)pack_buffer_size);
        int position = 0;
        pack_proc_sphere_part_node_data(
          local_data, pack_buffer, pack_buffer_size, &position, comm);

        // broadcast data size to other group
        yac_bcast_group(
          &position, 1, MPI_INT, local_group_comm.start, remote_group_comm);

        // broadcast remote_data to other root
        yac_bcast_group(pack_buffer, position, MPI_PACKED,
                        local_group_comm.start, remote_group_comm);
        free(pack_buffer);
      }
    }
  }

  // unpack node data
  int position = 0;
  struct proc_sphere_part_node_data other_data =
    unpack_proc_sphere_part_node_data(
      recv_buffer, data_size, &position, comm);
  free(recv_buffer);

  return other_data;
}

static void compute_redist_recvcounts_rdispls(
  int comm_rank, int split_rank, int comm_size,
  size_t global_sizes[2], size_t (*all_bucket_sizes)[2],
  int * counts, int * displs, size_t * recv_count) {

  int color = comm_rank >= split_rank;

  int split_comm_size = split_rank;
  int split_comm_rank = comm_rank;
  if (color) {
    split_comm_rank -= split_comm_size;
    split_comm_size = comm_size - split_comm_size;
  }

  size_t global_size = global_sizes[color];
  size_t local_interval_start =
    (size_t)
      (((long long)global_size * (long long)split_comm_rank +
        (long long)(split_comm_size - 1)) /
       (long long)split_comm_size);
  size_t local_interval_end =
    (size_t)
      (((long long)global_size * (long long)(split_comm_rank+1) +
        (long long)(split_comm_size - 1)) /
      (long long)split_comm_size);

  *recv_count = (size_t)(local_interval_end - local_interval_start);

  size_t start_idx = 0;
  for (int i = 0; i < comm_size; ++i) {

    size_t next_start_idx = start_idx + all_bucket_sizes[i][color];
    size_t interval_start = MAX(start_idx, local_interval_start);
    size_t interval_end = MIN(next_start_idx, local_interval_end);

    if (interval_start < interval_end) {

      size_t count = interval_end - interval_start;
      size_t disp = interval_start - local_interval_start;

      YAC_ASSERT(
        (count <= INT_MAX) && (disp <= INT_MAX),
        "ERROR(compute_redist_recvcounts_rdispls): invalid interval")

      counts[i] = (int)count;
      displs[i] = (int)disp;
    } else {
      counts[i] = 0;
      displs[i] = 0;
    }

    start_idx = next_start_idx;
  }
}

static void compute_redist_sendcounts_sdispls(
  int comm_rank, int split_rank, int comm_size,
  size_t global_sizes[2], size_t (*all_bucket_sizes)[2],
  int * counts, int * displs) {

  size_t U_size = all_bucket_sizes[comm_rank][0];

  size_t local_interval_start[2] = {0, 0};
  size_t local_interval_end[2];
  for (int i = 0; i < comm_rank; ++i)
    for (int j = 0; j < 2; ++j)
      local_interval_start[j] += all_bucket_sizes[i][j];
  for (int j = 0; j < 2; ++j)
    local_interval_end[j] =
      local_interval_start[j] + all_bucket_sizes[comm_rank][j];

  int comm_sizes[2] = {split_rank, comm_size - split_rank};

  size_t start_idx[2] = {0,0};
  for (int i = 0; i < comm_size; ++i) {

    int color = i >= split_rank;
    size_t global_size = global_sizes[color];
    int split_comm_rank = i - (color?(split_rank):0);
    int split_comm_size = comm_sizes[color];
    size_t next_start_idx =
      (size_t)(
        ((long long)global_size * (long long)(split_comm_rank + 1) +
         (long long)(split_comm_size - 1)) /
        (long long)split_comm_size);
    size_t interval_start = MAX(start_idx[color], local_interval_start[color]);
    size_t interval_end = MIN(next_start_idx, local_interval_end[color]);

    if (interval_start < interval_end) {

      size_t count = interval_end - interval_start;
      size_t disp = interval_start - local_interval_start[color] +
                      ((color)?(U_size):(0));

      YAC_ASSERT(
        (count <= INT_MAX) && (disp <= INT_MAX),
        "ERROR(compute_redist_sendcounts_sdispls): invalid interval")

      counts[i] = (int)count;
      displs[i] = (int)disp;
    } else {
      counts[i] = 0;
      displs[i] = 0;
    }

    start_idx[color] = next_start_idx;
  }
}

static void reorder_data_remote_point_info(
  size_t * reorder_idx, size_t count, struct remote_point_info * data) {

  // this routine assumes that all entries in reorder_idx are unique and in the
  // range [0;count[

  // for all elements
  for (size_t i = 0; i < count; ++i) {

    if (reorder_idx[i] != i) {

      size_t j = i;
      struct remote_point_info temp = data[i];

      while(1) {
        { // swap(j, reorder_idx[j]}
          size_t swap = reorder_idx[j];
          reorder_idx[j] = j;
          j = swap;
        }
        if (j == i) break;
        { // swap(temp, data[j]
          struct remote_point_info swap = data[j];
          data[j] = temp;
          temp = swap;
        }
      };

      data[i] = temp;
    }
  }
}

static int compare_dist_vertices(const void * a, const void * b) {

  int ret;

  if ((ret = ((const struct dist_vertex *)a)->grid_idx -
             ((const struct dist_vertex *)b)->grid_idx)) return ret;
  if (((const struct dist_vertex *)a)->global_id != XT_INT_MAX)
    return
      (((const struct dist_vertex *)a)->global_id >
       ((const struct dist_vertex *)b)->global_id) -
      (((const struct dist_vertex *)a)->global_id <
       ((const struct dist_vertex *)b)->global_id);
  else
    return
      compare_coords(
        ((const struct dist_vertex *)a)->coord,
        ((const struct dist_vertex *)b)->coord);
}

static void remove_duplicated_vertices(
  struct dist_vertex * vertices, size_t * num_vertices,
  struct remote_point_info * owners, size_t num_owners,
  size_t ** owners_reorder_idx, size_t * owners_reorder_idx_array_size) {

  if (*num_vertices == 0) return;

  ENSURE_ARRAY_SIZE(
    *owners_reorder_idx, *owners_reorder_idx_array_size, num_owners);

  // initialise owner offset
  {
    size_t owner_offset = 0;
    for (size_t i = 0; i < *num_vertices; ++i) {
      vertices[i].owner_offset = owner_offset;
      owner_offset += (size_t)(vertices[i].num_owners);
    }
    YAC_ASSERT(
      owner_offset == num_owners,
      "ERROR(remove_duplicated_vertices): internal error "
      "(owner_offset != num_owners)");
  }

  // sort vertices by grid id and global id, if available otherwise compare
  // coordinates
  qsort(vertices, *num_vertices, sizeof(*vertices), compare_dist_vertices);
  // on GPU/NEC VE do seperate sorts for each criteria using radix sort

  size_t old_num_vertices = *num_vertices;
  size_t new_num_vertices = 0;
  struct dist_vertex dummy_vertex = {.grid_idx = INT_MAX};
  struct dist_vertex * prev_vertex = &dummy_vertex;
  struct dist_vertex * curr_vertex = vertices;
  size_t * owners_reorder_idx_ = *owners_reorder_idx;
  for (size_t i = 0, reorder_idx = 0; i < old_num_vertices;
       ++i, ++curr_vertex) {
    size_t owner_offset = curr_vertex->owner_offset;
    for (int j = 0; j < curr_vertex->num_owners;
         ++j, ++owner_offset, ++reorder_idx)
      owners_reorder_idx_[owner_offset] = reorder_idx;
    if (compare_dist_vertices(prev_vertex, curr_vertex)) {
      prev_vertex = vertices + new_num_vertices;
      ++new_num_vertices;
      if (prev_vertex != curr_vertex) *prev_vertex = *curr_vertex;
    } else {
      prev_vertex->num_owners += curr_vertex->num_owners;
    }
  }
  *num_vertices = new_num_vertices;

  // reorder owners according to new vertex order
  reorder_data_remote_point_info(owners_reorder_idx_, num_owners, owners);
}

static void redistribute_dist_vertices(
  struct dist_vertex ** vertices, size_t * num_vertices,
  struct remote_point_info ** owners, size_t * num_owners,
  size_t global_bucket_sizes[2], size_t (*all_bucket_sizes)[2],
  int split_rank, struct comm_buffers comm_buffers,
  size_t ** owners_reorder_idx, size_t * owners_reorder_idx_array_size,
  MPI_Datatype dist_vertex_dt, MPI_Datatype remote_point_info_dt,
  struct yac_group_comm group_comm) {

#ifdef SCOREP_USER_ENABLE
SCOREP_USER_REGION_DEFINE( redist_data_region )
SCOREP_USER_REGION_BEGIN(
  redist_data_region, "data redistribution", SCOREP_USER_REGION_TYPE_COMMON )
#endif

  int group_rank = yac_group_comm_get_rank(group_comm);
  int group_size = yac_group_comm_get_size(group_comm);

  // compute send and receive counts and respective ranks for data
  // redistribution
  compute_redist_sendcounts_sdispls(
    group_rank, split_rank, group_size, global_bucket_sizes, all_bucket_sizes,
    comm_buffers.sendcounts, comm_buffers.sdispls);
  size_t new_num_vertices;
  compute_redist_recvcounts_rdispls(
    group_rank, split_rank, group_size, global_bucket_sizes, all_bucket_sizes,
    comm_buffers.recvcounts, comm_buffers.rdispls, &new_num_vertices);

  // redistribute vertices
  struct dist_vertex * new_vertices =
    xmalloc(new_num_vertices * sizeof(*new_vertices));
  yac_alltoallv_p2p_group(
    *vertices, comm_buffers.sendcounts, comm_buffers.sdispls,
    new_vertices, comm_buffers.recvcounts, comm_buffers.rdispls,
    sizeof(**vertices), dist_vertex_dt, group_comm);

  // adjust comm_buffers for redistribution of owner data
  size_t new_num_owners;
  {
    size_t saccu = 0, raccu = 0, send_vertex_idx = 0, recv_vertex_idx = 0;
    for (int i = 0; i < group_size; ++i) {
      size_t send_count = 0, recv_count = 0;
      for (int j = 0; j < comm_buffers.sendcounts[i]; ++j, ++send_vertex_idx)
        send_count += (size_t)((*vertices)[send_vertex_idx].num_owners);
      for (int j = 0; j < comm_buffers.recvcounts[i]; ++j, ++recv_vertex_idx)
        recv_count += (size_t)(new_vertices[recv_vertex_idx].num_owners);
      YAC_ASSERT(
        (saccu <= INT_MAX) && (raccu <= INT_MAX),
        "ERROR(redistribute_dist_vertices): displacement exceeds INT_MAX");
      YAC_ASSERT(
        (send_count <= INT_MAX) && (recv_count <= INT_MAX),
        "ERROR(redistribute_dist_vertices): counts exceeds INT_MAX");
      comm_buffers.sendcounts[i] = (int)send_count;
      comm_buffers.recvcounts[i] = (int)recv_count;
      comm_buffers.sdispls[i] = (int)saccu;
      comm_buffers.rdispls[i] = (int)raccu;
      saccu += send_count;
      raccu += recv_count;
    }
    YAC_ASSERT(
      saccu == *num_owners,
      "ERROR(redistribute_dist_vertices): inconsistent owner count");
    new_num_owners = raccu;
  }

  // redistribute owner data
  struct remote_point_info * new_owners =
    xmalloc(new_num_owners * sizeof(*new_owners));
  yac_alltoallv_p2p_group(
    *owners, comm_buffers.sendcounts, comm_buffers.sdispls,
    new_owners, comm_buffers.recvcounts, comm_buffers.rdispls,
    sizeof(**owners), remote_point_info_dt, group_comm);

  free(*vertices);
  free(*owners);
  *vertices = new_vertices;
  *num_vertices = new_num_vertices;
  *owners = new_owners;
  *num_owners = new_num_owners;

  // check for duplicated verties
  // (identified either by global id (if available) or coordinates)
  remove_duplicated_vertices(
    *vertices, num_vertices, *owners, *num_owners,
    owners_reorder_idx, owners_reorder_idx_array_size);

#ifdef SCOREP_USER_ENABLE
SCOREP_USER_REGION_END( redist_data_region )
#endif
}

static struct proc_sphere_part_node * generate_proc_sphere_part_node_recursive(
  struct dist_vertex ** vertices, size_t * num_vertices,
  struct remote_point_info ** owners, size_t * num_owners,
  size_t (*all_bucket_sizes)[2], struct comm_buffers comm_buffers,
  size_t ** reorder_idx, size_t * reorder_idx_array_size,
  int ** list_flag_, size_t * list_flag_array_size,
  MPI_Datatype dist_vertex_dt, MPI_Datatype remote_point_info_dt,
  struct yac_group_comm group_comm, double prev_gc_norm_vector[3]) {

#ifdef SCOREP_USER_ENABLE
SCOREP_USER_REGION_DEFINE( local_balance_point_region )
SCOREP_USER_REGION_DEFINE( global_balance_point_region )
SCOREP_USER_REGION_DEFINE( splitting_region )
SCOREP_USER_REGION_DEFINE( comm_split_region )
#endif

  int group_rank = yac_group_comm_get_rank(group_comm);
  int group_size = yac_group_comm_get_size(group_comm);

  //--------------------
  // compute split plane
  //--------------------

#ifdef SCOREP_USER_ENABLE
SCOREP_USER_REGION_BEGIN(
  local_balance_point_region, "local balance point",
  SCOREP_USER_REGION_TYPE_COMMON )
#endif
  // compute local balance point
  double balance_point[3] = {0.0, 0.0, 0.0};
  for (size_t i = 0; i < *num_vertices; ++i) {
    double * vertex_coord = (*vertices)[i].coord;
    for (int j = 0; j < 3; ++j) balance_point[j] += vertex_coord[j];
  }
#ifdef SCOREP_USER_ENABLE
SCOREP_USER_REGION_END( local_balance_point_region )
SCOREP_USER_REGION_BEGIN(
  global_balance_point_region, "global balance point",
  SCOREP_USER_REGION_TYPE_COMMON )
#endif

  // compute global balance point (make sure that the allreduce operation
  // generates bit-identical results on all processes)
  yac_allreduce_sum_dble(&(balance_point[0]), 3, group_comm);

  // check whether the computed balance_point is unambiguous, otherwise use
  // a point which is perpendicularly to the previous split plane
  if ((fabs(balance_point[0]) > 1e-9) ||
      (fabs(balance_point[1]) > 1e-9) ||
      (fabs(balance_point[2]) > 1e-9)) {
     normalise_vector(balance_point);
  } else {
     balance_point[0] = prev_gc_norm_vector[2];
     balance_point[1] = prev_gc_norm_vector[0];
     balance_point[2] = prev_gc_norm_vector[1];
  }

  // compute norm vector of new split plane
  double gc_norm_vector[3];
  crossproduct_kahan(balance_point, prev_gc_norm_vector, gc_norm_vector);

  // check whether the computed norm vector is unambiguous, otherwise use
  // a norm vector which is perpendicularly to the previous split plane
  if ((fabs(gc_norm_vector[0]) > 1e-9) ||
      (fabs(gc_norm_vector[1]) > 1e-9) ||
      (fabs(gc_norm_vector[2]) > 1e-9)) {
     normalise_vector(gc_norm_vector);
  } else {
     gc_norm_vector[0] = prev_gc_norm_vector[2];
     gc_norm_vector[1] = prev_gc_norm_vector[0];
     gc_norm_vector[2] = prev_gc_norm_vector[1];
  }

#ifdef SCOREP_USER_ENABLE
SCOREP_USER_REGION_END( global_balance_point_region )
#endif

  //-----------------
  // split local data
  //-----------------

#ifdef SCOREP_USER_ENABLE
SCOREP_USER_REGION_BEGIN( splitting_region, "splitting data", SCOREP_USER_REGION_TYPE_COMMON )
#endif

  ENSURE_ARRAY_SIZE(*list_flag_, *list_flag_array_size, *num_vertices);
  int * list_flag = *list_flag_;

  // angle between a vertex coord and the great circle plane:
  // acos(dot(gc_norm_vector, vertex_coord)) = angle(gc_norm_vector, vertex_coord)
  // acos(dot(gc_norm_vector, vertex_coord)) - PI/2 = angle(gc_plane, vertex_coord)
  // dot <= 0.0    -> U list
  // dot >  0.0    -> T list

  // compute for each vertex the list it belongs to
  struct dist_vertex * vertices_ = *vertices;
  size_t num_vertices_ = *num_vertices;
  YAC_OMP_PARALLEL
  {
    YAC_OMP_FOR
    for (size_t i = 0; i < num_vertices_; ++i) {
      double * curr_coordinates_xyz = vertices_[i].coord;
      double dot = curr_coordinates_xyz[0] * gc_norm_vector[0] +
                   curr_coordinates_xyz[1] * gc_norm_vector[1] +
                   curr_coordinates_xyz[2] * gc_norm_vector[2];

      // if (angle >= M_PI_2)
      list_flag[i] = dot <= 0.0;
    }
  }
  size_t U_size = 0, T_size = 0;
  for (size_t i = 0; i < num_vertices_; ++i) {
    if (list_flag[i]) ++U_size;
    else ++T_size;
  }

  // initialise owner offset
  for (size_t i = 0, owner_offset = 0; i < num_vertices_; ++i) {
    vertices_[i].owner_offset = owner_offset;
    owner_offset += (size_t)(vertices_[i].num_owners);
  }

  // The number of T-vertices among the first U-size number of vertices in the
  // array is equal to the number of U-vertices in the remaining array. These
  // have to be exchanged in order to get a sorted vertex array

  // search for all T-vertices in the U-part of the array and swap them with a
  // U-vertex in the T-part
  for (size_t i = 0, j = U_size; i < U_size; ++i) {
    // if the current vertex belongs to the T-list
    if (!list_flag[i]) {
      // search for a matching U-vertex
      for (;!list_flag[j];++j);
      struct dist_vertex temp_vertex = vertices_[i];
      vertices_[i] = vertices_[j];
      vertices_[j] = temp_vertex;
      ++j;
    }
  }

  // reorder owners accoring to new vertex order
  ENSURE_ARRAY_SIZE(*reorder_idx, *reorder_idx_array_size, *num_owners);
  size_t * reorder_idx_ = *reorder_idx;
  for (size_t i = 0, k = 0; i < num_vertices_; ++i) {
    size_t offset = vertices_[i].owner_offset;
    int num_owners = vertices_[i].num_owners;
    for (int j = 0; j < num_owners; ++j, ++k, ++offset) {
      reorder_idx_[offset] = k;
    }
  }
  reorder_data_remote_point_info(reorder_idx_, *num_owners, *owners);

#ifdef SCOREP_USER_ENABLE
SCOREP_USER_REGION_END( splitting_region )
#endif

  size_t bucket_sizes[2] = {U_size, T_size};

  // exchange local U/T sizes between all processes
  yac_allgather_size_t(
    &(bucket_sizes[0]), &(all_bucket_sizes[0][0]), 2, group_comm);

  // determine global U/T sizes
  size_t global_bucket_sizes[2] = {0, 0};
  for (int i = 0; i < group_size; ++i)
    for (int j = 0; j < 2; ++j)
      global_bucket_sizes[j] += all_bucket_sizes[i][j];
  size_t global_num_vertices =
    global_bucket_sizes[0] + global_bucket_sizes[1];

  //----------------------
  // split into two groups
  //----------------------
#ifdef SCOREP_USER_ENABLE
SCOREP_USER_REGION_BEGIN(
  comm_split_region, "creating splitcomm", SCOREP_USER_REGION_TYPE_COMMON )
#endif
  // determine processor groups
  int split_rank =
    MIN(
      (int)MAX(
        ((global_bucket_sizes[0] * (size_t)group_size +
          global_num_vertices/2) / global_num_vertices), 1),
      group_size - 1);

  // generate processor groups
  struct yac_group_comm local_group_comm, remote_group_comm;
  yac_group_comm_split(
    group_comm, split_rank, &local_group_comm, &remote_group_comm);

#ifdef SCOREP_USER_ENABLE
SCOREP_USER_REGION_END( comm_split_region )
#endif

  //------------------
  // redistribute data
  //------------------

  redistribute_dist_vertices(
    vertices, num_vertices, owners, num_owners,
    global_bucket_sizes, all_bucket_sizes, split_rank, comm_buffers,
    reorder_idx, reorder_idx_array_size, dist_vertex_dt, remote_point_info_dt,
    group_comm);

  //----------
  // recursion
  //----------

  // generate proc_sphere_part node for remaining data
  struct proc_sphere_part_node_data local_data;

  if (yac_group_comm_get_size(local_group_comm) > 1) {
    local_data.data.node =
      generate_proc_sphere_part_node_recursive(
        vertices, num_vertices, owners, num_owners,
        all_bucket_sizes, comm_buffers, reorder_idx, reorder_idx_array_size,
        list_flag_, list_flag_array_size, dist_vertex_dt, remote_point_info_dt,
        local_group_comm, gc_norm_vector);
    local_data.is_leaf = 0;
  } else {

    local_data.data.rank = yac_group_comm_get_global_rank(group_comm);
    local_data.is_leaf = 1;
  }

  // get proc_sphere_part_node_data from remote group
  struct proc_sphere_part_node_data remote_data =
    get_remote_data(local_data, local_group_comm, remote_group_comm);

  // generate node
  struct proc_sphere_part_node * node = xmalloc(1 * sizeof(*node));
  if (group_rank < split_rank) {
    node->U = local_data;
    node->T = remote_data;
  } else {
    node->U = remote_data;
    node->T = local_data;
  }
  memcpy(node->gc_norm_vector, gc_norm_vector, sizeof(gc_norm_vector));

  return node;
}

static MPI_Datatype yac_get_dist_vertex_mpi_datatype(MPI_Comm comm) {

  struct dist_vertex dummy;
  MPI_Datatype dist_vertex_dt;
  int array_of_blocklengths[] = {3, 1, 1, 1, 1};
  const MPI_Aint array_of_displacements[] =
    {(MPI_Aint)(intptr_t)(const void *)&(dummy.coord[0]) -
     (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.grid_idx) -
     (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.global_id) -
     (MPI_Aint)(intptr_t)(const void *)&dummy,
     (MPI_Aint)(intptr_t)(const void *)&(dummy.num_owners) -
     (MPI_Aint)(intptr_t)(const void *)&dummy};
  const MPI_Datatype array_of_types[] =
    {MPI_DOUBLE, MPI_INT, yac_int_dt, MPI_INT};
  yac_mpi_call(
    MPI_Type_create_struct(4, array_of_blocklengths, array_of_displacements,
                           array_of_types, &dist_vertex_dt), comm);
  return yac_create_resized(dist_vertex_dt, sizeof(dummy), comm);
}

static void inform_dist_vertex_owners(
  struct dist_vertex * dist_vertices, size_t num_dist_vertices,
  struct remote_point_info * dist_owners,
  size_t ** reorder_idx, size_t * reorder_idx_array_size,
  yac_int **global_vertex_ids[2], int * global_ids_missing,
  int **vertex_ranks[2], size_t * num_vertices, MPI_Comm comm) {

  // Here we assume that the vertices are sorted first by their grid index and
  // second by global id/coordinate (depending on availablity).
  // Additionally, the vertices should contain no duplications.

  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  // count the number of vertices per grid
  size_t num_vertices_per_grid[2] = {0, 0};
  size_t num_owners_per_grid[2] = {0, 0};
  for (size_t i = 0; i < num_dist_vertices; ++i) {
    num_vertices_per_grid[dist_vertices[i].grid_idx]++;
    num_owners_per_grid[dist_vertices[i].grid_idx] +=
      (size_t)(dist_vertices[i].num_owners);
  }

  size_t * sendcounts, * recvcounts, * sdispls, * rdispls;
  yac_get_comm_buffers(
    1, &sendcounts, &recvcounts, &sdispls, &rdispls, comm);

  ENSURE_ARRAY_SIZE(
    *reorder_idx, *reorder_idx_array_size,
    MAX((num_owners_per_grid[0] + num_vertices[0]),
        (num_owners_per_grid[1] + num_vertices[1])));
  yac_int * global_vertex_ids_buffer =
    xmalloc(
      MAX(
        (global_ids_missing[0]?(num_owners_per_grid[0] + num_vertices[0]):0),
        (global_ids_missing[1]?(num_owners_per_grid[1] + num_vertices[1]):0)) *
      sizeof(*global_vertex_ids_buffer));

  // for both grids
  for (int grid_idx = 0; grid_idx < 2; ++grid_idx) {

    struct dist_vertex * dist_grid_vertices =
      dist_vertices + ((grid_idx == 0)?0:num_vertices_per_grid[0]);
    struct remote_point_info * dist_grid_owners =
      dist_owners + ((grid_idx == 0)?0:num_owners_per_grid[0]);

    // if the user did not provide global vertex ids for the current grid
    yac_int id_offset = 0;
    if (global_ids_missing[grid_idx]) {

      YAC_ASSERT(
        num_vertices_per_grid[grid_idx] <= (size_t)XT_INT_MAX,
        "ERROR(inform_dist_vertex_owners): global_id out of bounds");

      // determine exclusive scan of sum of numbers of unique
      // coordinates on all ranks
      yac_int yac_int_num_vertices = (yac_int)num_vertices_per_grid[grid_idx];
      yac_mpi_call(MPI_Exscan(&yac_int_num_vertices, &id_offset, 1, yac_int_dt,
                              MPI_SUM, comm), comm);
      if (comm_rank == 0) id_offset = 0;

      YAC_ASSERT(
        ((size_t)id_offset + num_vertices_per_grid[grid_idx]) <=
        (size_t)XT_INT_MAX,
        "ERROR(inform_dist_vertex_owners): global_id out of bounds")
    }

    memset(sendcounts, 0, (size_t)(comm_size + 1) * sizeof(*sendcounts));

    // determine send counts for vertex information
    for (size_t i = 0; i < num_owners_per_grid[grid_idx]; ++i)
      sendcounts[dist_grid_owners[i].rank]++;
    yac_generate_alltoallv_args(
      1, sendcounts, recvcounts, sdispls, rdispls, comm);

    size_t * send_reorder_idx = *reorder_idx;
    size_t * recv_reorder_idx = *reorder_idx + num_owners_per_grid[grid_idx];
    yac_int * send_global_vertex_ids = global_vertex_ids_buffer;
    yac_int * recv_global_vertex_ids = global_vertex_ids_buffer +
                                       num_owners_per_grid[grid_idx];

    for (size_t i = 0, k = 0; i < num_vertices_per_grid[grid_idx]; ++i) {

      int num_vertex_owners = dist_grid_vertices[i].num_owners;

      for (int j = 0; j < num_vertex_owners; ++j, ++k) {

        size_t pos = sdispls[dist_grid_owners[k].rank + 1]++;
        send_reorder_idx[pos] = dist_grid_owners[k].orig_pos;
        if (global_ids_missing[grid_idx])
          send_global_vertex_ids[pos] = id_offset + (yac_int)i;
      }
    }

    // exchange reorder idx
    yac_alltoallv_p2p(
      send_reorder_idx, sendcounts, sdispls,
      recv_reorder_idx, recvcounts, rdispls,
      sizeof(*send_reorder_idx), YAC_MPI_SIZE_T, comm,
      "inform_dist_vertex_owners", __LINE__);

    // generate vertex ranks
    {
      int * curr_vertex_ranks =
        xmalloc(num_vertices[grid_idx] * sizeof(*curr_vertex_ranks));
      size_t j = 0;
      for (int rank = 0; rank < comm_size; ++rank)
        for (size_t i = 0; i < recvcounts[rank]; ++i, ++j)
          curr_vertex_ranks[recv_reorder_idx[j]] = rank;
      *(vertex_ranks[grid_idx]) = curr_vertex_ranks;
    }

    // exchange and set global ids (if not provided by the user)
    if (global_ids_missing[grid_idx]) {

      yac_alltoallv_p2p(
        send_global_vertex_ids, sendcounts, sdispls,
        recv_global_vertex_ids, recvcounts, rdispls,
        sizeof(*send_global_vertex_ids), yac_int_dt, comm,
        "inform_dist_vertex_owners", __LINE__);

      yac_int * curr_global_vertex_ids =
        xmalloc(num_vertices[grid_idx] * sizeof(*curr_global_vertex_ids));

      for (size_t i = 0; i < num_vertices[grid_idx]; ++i)
        curr_global_vertex_ids[recv_reorder_idx[i]] = recv_global_vertex_ids[i];

      *(global_vertex_ids[grid_idx]) = curr_global_vertex_ids;
    }
  }

  free(global_vertex_ids_buffer);

  yac_free_comm_buffers(sendcounts, recvcounts, sdispls, rdispls);
}

void yac_proc_sphere_part_new(
  yac_coordinate_pointer vertex_coordinates[2], size_t * num_vertices,
  struct proc_sphere_part_node ** proc_sphere_part,
  yac_int **global_vertex_ids_[2], int **vertex_ranks[2], MPI_Comm comm) {

  // check whether there are user-provided global vertex ids
  yac_int *global_vertex_ids[2] =
    {*(global_vertex_ids_[0]), *(global_vertex_ids_[1])};
  int global_ids_missing[2] =
    {(num_vertices[0] > 0) && (global_vertex_ids[0] == NULL),
     (num_vertices[1] > 0) && (global_vertex_ids[1] == NULL)};
  yac_mpi_call(
    MPI_Allreduce(MPI_IN_PLACE, global_ids_missing, 2, MPI_INT, MPI_MAX, comm),
    comm);

  size_t total_num_vertices = num_vertices[0] + num_vertices[1];

  int comm_rank, comm_size;
  yac_mpi_call(MPI_Comm_rank(comm, &comm_rank), comm);
  yac_mpi_call(MPI_Comm_size(comm, &comm_size), comm);

  double base_gc_norm_vector[3] = {0.0,0.0,1.0};

  int vertices_available = total_num_vertices > 0;
  yac_mpi_call(
    MPI_Allreduce(
      MPI_IN_PLACE, &vertices_available, 1, MPI_INT, MPI_MAX, comm), comm);

  if ((comm_size > 1) && vertices_available) {

    // generate basic owner information for all vertices
    struct remote_point_info * dist_owners =
      xmalloc(total_num_vertices * sizeof(*dist_owners));
    for (size_t i = 0; i < num_vertices[0]; ++i)
      dist_owners[i] =
        (struct remote_point_info){.rank = comm_rank, .orig_pos = i};
    for (size_t i = 0, j = num_vertices[0]; i < num_vertices[1]; ++i, ++j)
      dist_owners[j] =
        (struct remote_point_info){.rank = comm_rank, .orig_pos = i};

    // generate vertex information for all vertices
    struct dist_vertex * dist_vertices =
      xmalloc(total_num_vertices * sizeof(*dist_vertices));
    {
      size_t k = 0;
      // for both grids
      for (int i = 0; i < 2; ++i) {
        // for all vertices of the grid
        for (size_t j = 0; j < num_vertices[i]; ++j, ++k) {
          memcpy(
            dist_vertices[k].coord, vertex_coordinates[i][j],
            3 * sizeof(vertex_coordinates[i][j][0]));
          dist_vertices[k].grid_idx = i;
          dist_vertices[k].global_id =
            (global_ids_missing[i])?XT_INT_MAX:global_vertex_ids[i][j];
          dist_vertices[k].num_owners = 1;
        }
      }
    }

    // set up buffers for generation of proc_sphere_part
    size_t num_dist_vertices = total_num_vertices;
    size_t num_dist_owners = total_num_vertices;
    size_t (*all_bucket_sizes)[2] =
      xmalloc((size_t)comm_size * sizeof(*(all_bucket_sizes)));
    struct comm_buffers comm_buffers;
    comm_buffers.sendcounts =
      xmalloc(4 * (size_t)comm_size * sizeof(*(comm_buffers.sendcounts)));
    comm_buffers.recvcounts = comm_buffers.sendcounts + 1 * comm_size;
    comm_buffers.sdispls =    comm_buffers.sendcounts + 2 * comm_size;
    comm_buffers.rdispls =    comm_buffers.sendcounts + 3 * comm_size;
    size_t * reorder_idx = NULL, reorder_idx_array_size = 0;
    int * list_flag = NULL;
    size_t list_flag_array_size = 0;
    MPI_Datatype dist_vertex_dt = yac_get_dist_vertex_mpi_datatype(comm);
    MPI_Datatype remote_point_info_dt =
      yac_get_remote_point_info_mpi_datatype(comm);
    struct yac_group_comm group_comm = yac_group_comm_new(comm);

    // initial redistribute of all vertices
    // (in case one of the two grids has significantly more vertices per process
    //  than the other, this improves the load balance in the initial step)
    yac_allgather_size_t(
      (size_t[2]){total_num_vertices,0}, &(all_bucket_sizes[0][0]), 2,
      group_comm);
    size_t global_num_vertices = 0;
    for (int i = 0; i < comm_size; ++i)
      global_num_vertices += all_bucket_sizes[i][0];
    redistribute_dist_vertices(
      &dist_vertices, &num_dist_vertices, &dist_owners, &num_dist_owners,
      (size_t[2]){global_num_vertices, 0}, all_bucket_sizes, comm_size,
      comm_buffers, &reorder_idx, &reorder_idx_array_size,
      dist_vertex_dt, remote_point_info_dt, group_comm);

    // generate proc_sphere_part
    *proc_sphere_part =
      generate_proc_sphere_part_node_recursive(
        &dist_vertices, &num_dist_vertices, &dist_owners, &num_dist_owners,
        all_bucket_sizes, comm_buffers, &reorder_idx, &reorder_idx_array_size,
        &list_flag, &list_flag_array_size, dist_vertex_dt, remote_point_info_dt,
        group_comm, base_gc_norm_vector);

    // cleanup
    yac_group_comm_delete(group_comm);
    yac_mpi_call(MPI_Type_free(&remote_point_info_dt), comm);
    yac_mpi_call(MPI_Type_free(&dist_vertex_dt), comm);
    free(list_flag);
    free(comm_buffers.sendcounts);
    free(all_bucket_sizes);

    // return information about distributed vertices to original owners
    inform_dist_vertex_owners(
      dist_vertices, num_dist_vertices, dist_owners,
      &reorder_idx, &reorder_idx_array_size, global_vertex_ids_,
      global_ids_missing, vertex_ranks, num_vertices, comm);

    // cleanup
    free(reorder_idx);
    free(dist_owners);
    free(dist_vertices);

  } else {
    *proc_sphere_part = xmalloc(1 * sizeof(**proc_sphere_part));
    (*proc_sphere_part)->U.data.rank = 0;
    (*proc_sphere_part)->U.is_leaf = 1;
    (*proc_sphere_part)->T.data.rank = 0;
    (*proc_sphere_part)->T.is_leaf = 1;
    (*proc_sphere_part)->gc_norm_vector[0] = base_gc_norm_vector[0];
    (*proc_sphere_part)->gc_norm_vector[1] = base_gc_norm_vector[1];
    (*proc_sphere_part)->gc_norm_vector[2] = base_gc_norm_vector[2];

    // generate global ids and vertex ranks
    for (int grid_idx = 0; grid_idx < 2; ++grid_idx) {
      *(vertex_ranks[grid_idx]) =
        xmalloc(num_vertices[grid_idx] * sizeof(**(vertex_ranks[grid_idx])));
      for (size_t i = 0; i < num_vertices[grid_idx]; ++i)
        (*(vertex_ranks[grid_idx]))[i] = 0;
      if (global_ids_missing[grid_idx]) {
        *(global_vertex_ids_[grid_idx]) =
          xmalloc(
            num_vertices[grid_idx] * sizeof(**(global_vertex_ids_[grid_idx])));
        for (size_t i = 0; i < num_vertices[grid_idx]; ++i)
          (*(global_vertex_ids_[grid_idx]))[i] = (yac_int)i;
      }
    }
  }
}

static int is_serial_node(struct proc_sphere_part_node * node) {
  return (node->U.is_leaf) && (node->T.is_leaf) &&
         (node->U.data.rank == 0) && (node->T.data.rank == 0);
}

// #define YAC_NEC_EXPERIMENTAL
#ifdef YAC_NEC_EXPERIMENTAL
// the following code may be better for a vector machine
static void yac_proc_sphere_part_do_point_search_recursive(
  struct proc_sphere_part_node * node, yac_coordinate_pointer search_coords,
  size_t * search_idx, size_t * temp_search_idx, int * flag,
  size_t count, int * ranks) {

  for (size_t i = 0; i < count; ++i) {

    // compute cos angle between current search point and norm vector
    double dot = search_coords[search_idx[i]][0] * node->gc_norm_vector[0] +
                 search_coords[search_idx[i]][1] * node->gc_norm_vector[1] +
                 search_coords[search_idx[i]][2] * node->gc_norm_vector[2];

    // if (angle >= M_PI_2)
    flag[i] = (dot <= 0.0);
  }

  size_t u_size = 0, t_size = 0;
  for (size_t i = 0; i < count; ++i) {
    if (flag[i]) search_idx[u_size++] = search_idx[i];
    else         temp_search_idx[t_size++] = search_idx[i];
  }

  if (node->U.is_leaf) {
    size_t rank = node->U.data.rank;;
    for (size_t i = 0; i < u_size; ++i) ranks[search_idx[i]] = rank;
  } else {
    yac_proc_sphere_part_do_point_search_recursive(
      node->U.data.node, search_coords, search_idx, temp_search_idx + t_size,
      flag, u_size, ranks);
  }

  if (node->T.is_leaf) {
    size_t rank = node->T.data.rank;
    for (size_t i = 0; i < t_size; ++i) ranks[temp_search_idx[i]] = rank;
  } else {
    yac_proc_sphere_part_do_point_search_recursive(
      node->T.data.node, search_coords, temp_search_idx, search_idx + u_size,
      flag, t_size, ranks);
  }
}
#endif // YAC_NEC_EXPERIMENTAL

void yac_proc_sphere_part_do_point_search(
  struct proc_sphere_part_node * node, yac_coordinate_pointer search_coords,
  size_t count, int * ranks) {

  if (is_serial_node(node)) {
    for (size_t i = 0; i < count; ++i) ranks[i] = 0;
    return;
  }

#ifdef YAC_NEC_EXPERIMENTAL

  size_t * search_idx = xmalloc(2 * count * sizeof(*search_idx));
  for (size_t i = 0; i < count; ++i) search_idx[i] = i;
  int * flag = xmalloc(count * sizeof(*flag));
  yac_proc_sphere_part_do_point_search_recursive(
    node, search_coords, search_idx, search_idx + count, flag, count, ranks);
  free(flag);
  free(search_idx);

#else

  YAC_OMP_PARALLEL
  {
    YAC_OMP_FOR
    for (size_t i = 0; i < count; ++i) {

      double * curr_coord = search_coords[i];

      struct proc_sphere_part_node * curr_node = node;

      while (1) {

        // compute cos angle between current search point and norm vector
        double dot = curr_coord[0] * curr_node->gc_norm_vector[0] +
                     curr_coord[1] * curr_node->gc_norm_vector[1] +
                     curr_coord[2] * curr_node->gc_norm_vector[2];

        // if (angle >= M_PI_2)
        if (dot <= 0.0) {
          if (curr_node->U.is_leaf) {
            ranks[i] = curr_node->U.data.rank;
            break;
          } else {
            curr_node = curr_node->U.data.node;
            continue;
          }
        // if (angle < M_PI_2)
        } else {
          if (curr_node->T.is_leaf) {
            ranks[i] = curr_node->T.data.rank;
            break;
          } else {
            curr_node = curr_node->T.data.node;
            continue;
          }
        }
      }
    }
  }
#endif // YAC_NEC_EXPERIMENTAL
}

static void bnd_circle_search(
  struct proc_sphere_part_node * node, struct bounding_circle bnd_circle,
  int * ranks, int * rank_count) {

  double dot = bnd_circle.base_vector[0] * node->gc_norm_vector[0] +
               bnd_circle.base_vector[1] * node->gc_norm_vector[1] +
               bnd_circle.base_vector[2] * node->gc_norm_vector[2];

  // angle < M_PI_2 + bnd_circle.inc_angle
  if (dot > - bnd_circle.inc_angle.sin) {

    if (node->T.is_leaf) {

      ranks[*rank_count] = node->T.data.rank;
      ++*rank_count;

    } else {
      bnd_circle_search(node->T.data.node, bnd_circle, ranks, rank_count);
    }
  }

  // angle > M_PI_2 - bnd_circle.inc_angle
  if (dot < bnd_circle.inc_angle.sin) {

    if (node->U.is_leaf) {

      ranks[*rank_count] = node->U.data.rank;
      ++*rank_count;

    } else {
      bnd_circle_search(node->U.data.node, bnd_circle, ranks, rank_count);
    }
  }
}

static void bnd_circle_search_big_angle(
  struct proc_sphere_part_node * node, struct bounding_circle bnd_circle,
  int * ranks, int * rank_count) {

  if (node->T.is_leaf) {

    ranks[*rank_count] = node->T.data.rank;
    ++*rank_count;

  } else {
    bnd_circle_search_big_angle(
      node->T.data.node, bnd_circle, ranks, rank_count);
  }

  if (node->U.is_leaf) {

    ranks[*rank_count] = node->U.data.rank;
    ++*rank_count;

  } else {
    bnd_circle_search_big_angle(
      node->U.data.node, bnd_circle, ranks, rank_count);
  }
}

void yac_proc_sphere_part_do_bnd_circle_search(
  struct proc_sphere_part_node * node, struct bounding_circle bnd_circle,
  int * ranks, int * rank_count) {

  YAC_ASSERT(
    compare_angles(bnd_circle.inc_angle, SIN_COS_M_PI) == -1,
    "ERROR(yac_proc_sphere_part_do_bnd_circle_search): angle is >= PI")

  // special case in which the proc_sphere_part_node only contains a single
  // rank
  if (is_serial_node(node)) {
    ranks[0] = 0;
    *rank_count = 1;
  } else if (bnd_circle.inc_angle.cos <= 0.0) {
    *rank_count = 0;
    bnd_circle_search_big_angle(node, bnd_circle, ranks, rank_count);
  } else {
    *rank_count = 0;
    bnd_circle_search(node, bnd_circle, ranks, rank_count);
  }
}

static int get_leaf_ranks(struct proc_sphere_part_node * node, int * ranks) {

  int curr_size;

  if (node->U.is_leaf) {
    *ranks = node->U.data.rank;
    curr_size = 1;
  } else {
    curr_size = get_leaf_ranks(node->U.data.node, ranks);
  }
  if (node->T.is_leaf) {
    ranks[curr_size++] = node->T.data.rank;
  } else {
    curr_size += get_leaf_ranks(node->T.data.node, ranks + curr_size);
  }

  return curr_size;
}

static void get_neigh_ranks(
  struct proc_sphere_part_node * node, uint64_t * leaf_sizes, uint64_t min_size,
  uint64_t ** inner_node_sizes, int * send_flags, int * recv_flags,
  int comm_rank, struct neigh_search_data * last_valid_node) {

  int curr_node_is_valid = (**inner_node_sizes >= min_size);

  struct neigh_search_data * curr_valid_node;
  struct neigh_search_data temp_valid_node;

  int * ranks = last_valid_node->ranks;

  if (curr_node_is_valid) {
    temp_valid_node.ranks = ranks;
    temp_valid_node.num_ranks = 0;
    temp_valid_node.node = node;
    curr_valid_node = &temp_valid_node;
    last_valid_node->num_ranks = 0;
  } else {
    curr_valid_node = last_valid_node;
  }

  for (int j = 0; j < 2; ++j) {

    struct proc_sphere_part_node_data node_data = (j == 0)?(node->U):(node->T);

    if (node_data.is_leaf) {

      int rank = node_data.data.rank;

      if (leaf_sizes[rank] < min_size) {

        if (curr_valid_node->num_ranks == 0)
          curr_valid_node->num_ranks =
            get_leaf_ranks(curr_valid_node->node, ranks);

        // if the current leaf is the local process
        if (rank == comm_rank)
          for (int i = 0; i < curr_valid_node->num_ranks; ++i)
            recv_flags[ranks[i]] = 1;

        // if the process of the current leaf required data from the
        // local process
        for (int i = 0; i < curr_valid_node->num_ranks; ++i) {
          if (ranks[i] == comm_rank) {
            send_flags[rank] = 1;
            break;
          }
        }
      }

    } else {
      ++*inner_node_sizes;
      get_neigh_ranks(
        node_data.data.node, leaf_sizes, min_size, inner_node_sizes,
        send_flags, recv_flags, comm_rank, curr_valid_node);
    }
  }
}

static uint64_t determine_node_sizes(
  struct proc_sphere_part_node * node, uint64_t * leaf_sizes,
  uint64_t ** inner_node_sizes) {

  uint64_t * curr_inner_node_size = *inner_node_sizes;
  uint64_t node_size;
  if (node->U.is_leaf) {
    node_size = leaf_sizes[node->U.data.rank];
  } else {
    ++*inner_node_sizes;
    node_size =
      determine_node_sizes(node->U.data.node, leaf_sizes, inner_node_sizes);
  }
  if (node->T.is_leaf) {
    node_size += leaf_sizes[node->T.data.rank];
  } else {
    ++*inner_node_sizes;
    node_size +=
      determine_node_sizes(node->T.data.node, leaf_sizes, inner_node_sizes);
  }
  return (*curr_inner_node_size = node_size);
}

void yac_proc_sphere_part_get_neigh_ranks(
  struct proc_sphere_part_node * node, uint64_t * leaf_sizes,
  uint64_t min_size, int * send_flags, int * recv_flags,
  int comm_rank, int comm_size) {

  uint64_t * inner_node_sizes =
    xcalloc((size_t)comm_size, sizeof(*inner_node_sizes));

  uint64_t * temp_inner_node_sizes = inner_node_sizes;
  determine_node_sizes(node, leaf_sizes, &temp_inner_node_sizes);

  YAC_ASSERT(
    *inner_node_sizes >= min_size,
    "ERROR(yac_proc_sphere_part_get_neigh_ranks): sum of global leaf sizes "
    "is < min_size")

  struct neigh_search_data search_data = {
    .ranks = xmalloc((size_t)comm_size * sizeof(int)),
    .num_ranks = 0,
    .node = node
  };

  temp_inner_node_sizes = inner_node_sizes;
  get_neigh_ranks(node, leaf_sizes, min_size, &temp_inner_node_sizes,
                  send_flags, recv_flags, comm_rank, &search_data);

  send_flags[comm_rank] = 0;
  recv_flags[comm_rank] = 0;

  free(search_data.ranks);
  free(inner_node_sizes);
}

void yac_proc_sphere_part_node_delete(struct proc_sphere_part_node * node) {

  if (!(node->U.is_leaf)) yac_proc_sphere_part_node_delete(node->U.data.node);
  if (!(node->T.is_leaf)) yac_proc_sphere_part_node_delete(node->T.data.node);
  free(node);
}
