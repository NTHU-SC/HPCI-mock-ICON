// Copyright (c) 2025 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>

#include "read_exodus_grid.h"
#include "utils_common.h"
#include "io_utils.h"
#include "geometry.h"
#include "read_grid.h"
#include "yac_mpi_internal.h"

#ifdef YAC_NETCDF_ENABLED
#include <netcdf.h>

static size_t * generate_offsets(size_t N, int * counts) {

  size_t * offsets = xmalloc(N * sizeof(*offsets));
  for (size_t i = 0, accu = 0; i < N; ++i) {
    offsets[i] = accu;
    accu += (size_t)(counts[i]);
  }
  return offsets;
}

// taken from scales-ppm library
// https://www.dkrz.de/redmine/projects/scales-ppm
static inline int
partition_idx_from_element_idx(size_t element_idx, size_t num_elements,
                               int num_partitions) {

  return (int)((((unsigned long)element_idx) * ((unsigned long)num_partitions) +
                (unsigned long)num_partitions - 1) /
               ((unsigned long)num_elements));
}

void yac_read_exodus_grid_information_parallel(
  const char * filename, MPI_Comm comm, yac_coordinate_pointer * node_coords,
  yac_int ** elem_ids, yac_int ** node_ids,
  size_t * num_elem, size_t * num_nodes,
  int ** num_nodes_per_elem, int ** num_elem_per_node,
  size_t ** elem_to_node, size_t ** node_to_elem) {

  int comm_rank, comm_size;

  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &comm_size);

  int local_is_io, * io_ranks, num_io_ranks;
  yac_get_io_ranks(comm, &local_is_io, &io_ranks, &num_io_ranks);

  size_t num_global_nodes, num_global_elem, num_nod_per_elem;

  size_t read_local_start_elem = 0;
  size_t read_num_local_elems = 0;
  size_t read_local_start_node = 0;
  size_t read_num_local_nodes = 0;

  yac_coordinate_pointer read_node_coords = NULL;
  int * read_dist_elem_to_node = NULL;

  if (local_is_io) {

    unsigned long io_proc_idx = ULONG_MAX;
    for (int i = 0; (i < num_io_ranks) && (io_proc_idx == ULONG_MAX); ++i)
      if (io_ranks[i] == comm_rank)
        io_proc_idx = (unsigned long)i;

    // open file
    int ncid;
    yac_nc_open(filename, NC_NOWRITE, &ncid);

    // get number of cells and vertices
    int dim_id;
    size_t num_el_blk, num_el_in_blk1;
    yac_nc_inq_dimid(ncid, "num_nodes", &dim_id);
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dim_id, &num_global_nodes));
    yac_nc_inq_dimid(ncid, "num_elem", &dim_id);
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dim_id, &num_global_elem));
    yac_nc_inq_dimid(ncid, "num_el_blk", &dim_id);
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dim_id, &num_el_blk));
    yac_nc_inq_dimid(ncid, "num_el_in_blk1", &dim_id);
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dim_id, &num_el_in_blk1));
    yac_nc_inq_dimid(ncid, "num_nod_per_el1", &dim_id);
    YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dim_id, &num_nod_per_elem));

    YAC_ASSERT(
      num_el_blk == 1,
      "ERROR(yac_read_exodus_grid_information_parallel): "
      "reader currently only supports a single block of elements");
    YAC_ASSERT(
      num_global_elem == num_el_in_blk1,
      "ERROR(yac_read_exodus_grid_information_parallel): "
      "total number of elements and number of elements in first block "
      "of elements do not match");

    // determine local range for element and node data
    read_local_start_elem =
      ((unsigned long)num_global_elem * io_proc_idx) / (unsigned long)num_io_ranks;
    read_num_local_elems =
      ((unsigned long)num_global_elem * (io_proc_idx+1)) / (unsigned long)num_io_ranks -
      (unsigned long)read_local_start_elem;
    read_local_start_node =
      ((unsigned long)num_global_nodes * io_proc_idx) / (unsigned long)num_io_ranks;
    read_num_local_nodes =
      ((unsigned long)num_global_nodes * (io_proc_idx+1)) / (unsigned long)num_io_ranks -
      (unsigned long)read_local_start_node;

    // read basic grid data (each process its individual part)
    double * read_node_coord_x =
      xmalloc(read_num_local_nodes * sizeof(*read_node_coord_x));
    double * read_node_coord_y =
      xmalloc(read_num_local_nodes * sizeof(*read_node_coord_y));
    double * read_node_coord_z =
      xmalloc(read_num_local_nodes * sizeof(*read_node_coord_z));;
    int varid;
    yac_nc_inq_varid(ncid, "coord", &varid);
    YAC_HANDLE_ERROR(
      nc_get_vara_double(
        ncid, varid, (size_t[]){0, read_local_start_node},
        (size_t[]){1,read_num_local_nodes}, read_node_coord_x));
    YAC_HANDLE_ERROR(
      nc_get_vara_double(
        ncid, varid, (size_t[]){1, read_local_start_node},
        (size_t[]){1,read_num_local_nodes}, read_node_coord_y));
    YAC_HANDLE_ERROR(
      nc_get_vara_double(
        ncid, varid, (size_t[]){2, read_local_start_node},
        (size_t[]){1,read_num_local_nodes}, read_node_coord_z));
    read_node_coords =
      xmalloc(read_num_local_nodes * sizeof(*read_node_coords));
    for (size_t i = 0; i < read_num_local_nodes; ++i) {
      read_node_coords[i][0] = read_node_coord_x[i];
      read_node_coords[i][1] = read_node_coord_y[i];
      read_node_coords[i][2] = read_node_coord_z[i];
    }
    free(read_node_coord_z);
    free(read_node_coord_y);
    free(read_node_coord_x);

    read_dist_elem_to_node =
      xmalloc(
        read_num_local_elems * num_nod_per_elem *
        sizeof(*read_dist_elem_to_node));
    yac_nc_inq_varid(ncid, "connect1", &varid);
    YAC_HANDLE_ERROR(
      nc_get_vara_int(
        ncid, varid, (size_t[]){read_local_start_elem, 0},
        (size_t[]){read_num_local_elems, num_nod_per_elem},
        read_dist_elem_to_node));
    for (size_t i = 0; i < read_num_local_elems * num_nod_per_elem; ++i)
      read_dist_elem_to_node[i]--;

    YAC_HANDLE_ERROR(nc_close(ncid));

  } else {
    read_node_coords = xmalloc(1 * sizeof(*read_node_coords));
    read_dist_elem_to_node = xmalloc(1 * sizeof(*read_dist_elem_to_node));
  }

  free(io_ranks);

  {
    size_t tmp;
    if (comm_rank == 0) tmp = num_global_nodes;
    MPI_Bcast(&tmp, 1, YAC_MPI_SIZE_T, 0, comm);
    num_global_nodes = tmp;
    if (comm_rank == 0) tmp = num_global_elem;
    MPI_Bcast(&tmp, 1, YAC_MPI_SIZE_T, 0, comm);
    num_global_elem = tmp;
    if (comm_rank == 0) tmp = num_nod_per_elem;
    MPI_Bcast(&tmp, 1, YAC_MPI_SIZE_T, 0, comm);
    num_nod_per_elem = tmp;
  }

  // determine local range for element and node data
  size_t local_start_elem =
    ((unsigned long)num_global_elem * (unsigned long)comm_rank) /
    (unsigned long)comm_size;
  size_t num_local_elems =
    ((unsigned long)num_global_elem * ((unsigned long)comm_rank+1)) /
    (unsigned long)comm_size - (unsigned long)local_start_elem;
  size_t local_start_node =
    ((unsigned long)num_global_nodes * (unsigned long)comm_rank) /
    (unsigned long)comm_size;
  size_t num_local_nodes =
    ((unsigned long)num_global_nodes * ((unsigned long)comm_rank+1)) /
    (unsigned long)comm_size - (unsigned long)local_start_node;

  // redistribute basic element data (from io decomposition)
  int * dist_elem_to_node =
    xmalloc(num_local_elems * num_nod_per_elem * sizeof(*dist_elem_to_node));
  {
    int * send_count = xcalloc(comm_size, sizeof(*send_count));
    int * recv_count = xcalloc(comm_size, sizeof(*recv_count));

    for (size_t i = 0; i < read_num_local_elems; ++i)
      send_count[
        partition_idx_from_element_idx(
          read_local_start_elem + i, num_global_elem, comm_size)] +=
            num_nod_per_elem;

    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

    int * send_displ = xmalloc(comm_size * sizeof(*send_displ));
    int * recv_displ = xmalloc(comm_size * sizeof(*recv_displ));
    int send_accum = 0, recv_accum = 0;
    for (int i = 0; i < comm_size; ++i) {
      send_displ[i] = send_accum;
      recv_displ[i] = recv_accum;
      send_accum += send_count[i];
      recv_accum += recv_count[i];
    }

    MPI_Alltoallv(read_dist_elem_to_node, send_count, send_displ, MPI_INT,
                  dist_elem_to_node, recv_count, recv_displ, MPI_INT, comm);

    free(recv_displ);
    free(send_displ);
    free(recv_count);
    free(send_count);
    free(read_dist_elem_to_node);
  }

  // redistribute basic node data (from io decomposition)
  yac_coordinate_pointer dist_node_coords =
    xmalloc(num_local_nodes * sizeof(*dist_node_coords));
  {
    int * send_count = xcalloc(comm_size, sizeof(*send_count));
    int * recv_count = xcalloc(comm_size, sizeof(*recv_count));

    for (size_t i = 0; i < read_num_local_nodes; ++i)
      send_count[
        partition_idx_from_element_idx(
          read_local_start_node + i, num_global_nodes, comm_size)] += 3;

    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

    int * send_displ = xmalloc(comm_size * sizeof(*send_displ));
    int * recv_displ = xmalloc(comm_size * sizeof(*recv_displ));
    int send_accum = 0, recv_accum = 0;
    for (int i = 0; i < comm_size; ++i) {
      send_displ[i] = send_accum;
      recv_displ[i] = recv_accum;
      send_accum += send_count[i];
      recv_accum += recv_count[i];
    }

    MPI_Alltoallv(read_node_coords, send_count, send_displ, MPI_DOUBLE,
                  dist_node_coords, recv_count, recv_displ, MPI_DOUBLE, comm);

    free(recv_displ);
    free(send_displ);
    free(recv_count);
    free(send_count);
    free(read_node_coords);
  }

  // determine required nodes for core elements
  // in additional compute elem_to_node, node_to_elem, and num_elem_per_node
  size_t num_core_nodes;
  {
    size_t N = num_local_elems * num_nod_per_elem;
    *node_ids = xmalloc(N * sizeof(**node_ids));
    *num_elem_per_node = xmalloc(N * sizeof(**num_elem_per_node));
    *elem_to_node = xmalloc(N * sizeof(**elem_to_node));
    *node_to_elem = xmalloc(N * sizeof(**node_to_elem));
    for (size_t i = 0; i < N; ++i)
      (*node_ids)[i] = (yac_int)dist_elem_to_node[i];
    size_t * permutation = *node_to_elem;
    for (size_t i = 0; i < N; ++i) permutation[i] = i;
    yac_quicksort_index_yac_int_size_t(*node_ids, N, permutation);
    // remove duplicated core nodes and count number of elements per node
    yac_int prev_node_id = XT_INT_MAX;
    num_core_nodes = 0;
    for (size_t i = 0; i < N; ++i) {
      yac_int curr_node_id = (*node_ids)[i];
      if (prev_node_id == curr_node_id) {
        (*num_elem_per_node)[num_core_nodes-1]++;
      } else {
        (*num_elem_per_node)[num_core_nodes] = 1;
        (*node_ids)[num_core_nodes] = (prev_node_id = curr_node_id);
        ++num_core_nodes;
      }
      (*elem_to_node)[permutation[i]] = num_core_nodes-1;
      permutation[i] /= num_nod_per_elem;
    }
    *node_ids =
      xrealloc(*node_ids, num_core_nodes * sizeof(**node_ids));
    *num_elem_per_node =
      xrealloc(*num_elem_per_node,
               num_core_nodes * sizeof(**num_elem_per_node));
    free(dist_elem_to_node);
  }

  // get node coordinate data
  {
    *node_coords = xmalloc(num_core_nodes * sizeof(**node_coords));
    int * send_count = xcalloc(comm_size, sizeof(*send_count));
    int * recv_count = xcalloc(comm_size, sizeof(*recv_count));

    for (size_t i = 0; i < num_core_nodes; ++i)
      send_count[
        partition_idx_from_element_idx(
          (*node_ids)[i], num_global_nodes, comm_size)]++;

    MPI_Alltoall(send_count, 1, MPI_INT, recv_count, 1, MPI_INT, comm);

    int * send_displ = xmalloc(comm_size * sizeof(*send_displ));
    int * recv_displ = xmalloc(comm_size * sizeof(*recv_displ));
    int send_accum = 0, recv_accum = 0;
    for (int i = 0; i < comm_size; ++i) {
      send_displ[i] = send_accum;
      recv_displ[i] = recv_accum;
      send_accum += send_count[i];
      recv_accum += recv_count[i];
    }

    int num_all_local_nodes_remote = 0;
    for (int i = 0; i < comm_size; ++i)
      num_all_local_nodes_remote += recv_count[i];

    yac_int * remote_node_buffer =
      xmalloc(num_all_local_nodes_remote * sizeof(*remote_node_buffer));

    MPI_Alltoallv(
      *node_ids, send_count, send_displ, yac_int_dt,
      remote_node_buffer, recv_count, recv_displ, yac_int_dt, comm);

    yac_coordinate_pointer send_node_coords =
      xmalloc(num_all_local_nodes_remote * sizeof(*send_node_coords));

    for (int i = 0, l = 0; i < comm_size; ++i) {
      for (int j = 0; j < recv_count[i]; ++j, ++l) {
        size_t idx = (size_t)(remote_node_buffer[l]) - local_start_node;
        send_node_coords[l][0] = dist_node_coords[idx][0];
        send_node_coords[l][1] = dist_node_coords[idx][1];
        send_node_coords[l][2] = dist_node_coords[idx][2];
      }
      send_count[i] *= 3;
      recv_count[i] *= 3;
      send_displ[i] *= 3;
      recv_displ[i] *= 3;
    }

    free(remote_node_buffer);
    free(dist_node_coords);

    MPI_Alltoallv(send_node_coords, recv_count, recv_displ, MPI_DOUBLE,
                  *node_coords, send_count, send_displ, MPI_DOUBLE, comm);

    free(send_node_coords);
    free(recv_displ);
    free(send_displ);
    free(recv_count);
    free(send_count);
  }

  // generate elem ids for local partition
  *elem_ids = xmalloc(num_local_elems * sizeof(**elem_ids));
  for (size_t i = 0; i < num_local_elems; ++i)
    (*elem_ids)[i] = (yac_int)(local_start_elem + i);

  // generate num_nodes_per_elem
  *num_nodes_per_elem =
    xmalloc(num_local_elems * sizeof(**num_nodes_per_elem));
  for (size_t i = 0; i < num_local_elems; ++i)
    (*num_nodes_per_elem)[i] = num_nod_per_elem;

  *num_elem = num_local_elems;
  *num_nodes = num_core_nodes;
}

struct temp_edge {
  size_t node[2];
  size_t elem_to_edge_idx;
};

static int compare_temp_edges(void const * a, void const * b) {

  struct temp_edge const * edge_a = (struct temp_edge const *)a;
  struct temp_edge const * edge_b = (struct temp_edge const *)b;

  if (edge_a->node[0] != edge_b->node[0])
    return (edge_a->node[0] > edge_b->node[0])?1:-1;
  return (edge_a->node[1] > edge_b->node[1]) -
         (edge_a->node[1] < edge_b->node[1]);
}

static int check_pole(double * vertex) {

  return fabs(1.0 - fabs(vertex[2])) < 1e-8;
}

static int check_lon_edge(double * vertex_a_, double * vertex_b_) {

  double vertex_a[3] = {vertex_a_[0], vertex_a_[1], 0.0};
  double vertex_b[3] = {vertex_b_[0], vertex_b_[1], 0.0};

  normalise_vector(vertex_a);
  normalise_vector(vertex_b);

  return get_vector_angle(vertex_a, vertex_b) < 1e-6;
}

static int check_lat_edge(double * vertex_a, double * vertex_b) {

  return fabs(acos(vertex_a[2]) - acos(vertex_b[2])) < 1e-6;
}

struct yac_basic_grid_data yac_read_exodus_basic_grid_data_parallel(
  const char * filename, int use_ll_edges, MPI_Comm comm) {

  yac_coordinate_pointer node_coords;
  yac_int * elem_ids;
  yac_int * node_ids;
  size_t num_elem;
  size_t num_nodes;
  int * num_nodes_per_elem;
  int * num_elem_per_node;
  size_t * elem_to_node;
  size_t * node_to_elem;

  yac_read_exodus_grid_information_parallel(
    filename, comm, &node_coords, &elem_ids, &node_ids, &num_elem, &num_nodes,
    &num_nodes_per_elem, &num_elem_per_node, &elem_to_node, &node_to_elem);

  size_t num_edges;
  size_t * elem_to_edge;
  yac_size_t_2_pointer edge_to_node;
  enum yac_edge_type * edge_type;

  { // compute edge data

    // compute the maximum number of edge
    size_t max_num_edges = 0;
    for (size_t i = 0; i < num_elem; ++i)
      max_num_edges += num_nodes_per_elem[i];

    // generate temporary array containing edge information
    struct temp_edge * temp_edges =
      xmalloc(max_num_edges * sizeof(*temp_edges));
    for (size_t i = 0, offset = 0, k = 0; i < num_elem; ++i) {
      size_t * curr_elem_to_node = elem_to_node + offset;
      size_t curr_num_edges = num_nodes_per_elem[i];
      offset += curr_num_edges;
      for (size_t j = 0; j < curr_num_edges; ++j, ++k) {
        int order =
          curr_elem_to_node[j] > curr_elem_to_node[(j+1)%curr_num_edges];
        temp_edges[k].node[order] = curr_elem_to_node[j];
        temp_edges[k].node[order^1] = curr_elem_to_node[(j+1)%curr_num_edges];
        temp_edges[k].elem_to_edge_idx = k;
      }
    }
    qsort(temp_edges, max_num_edges,
          sizeof(*temp_edges), compare_temp_edges);

    // generate elem_to_edge and edge_to_node; count total number of edges
    elem_to_edge = xmalloc(max_num_edges * sizeof(*elem_to_edge));
    num_edges = 0;
    edge_to_node = (yac_size_t_2_pointer)temp_edges;
    for (size_t i = 0, prev_indices[2] = {SIZE_MAX, SIZE_MAX};
        i < max_num_edges; ++i) {

      size_t curr_elem_to_edge_idx = temp_edges[i].elem_to_edge_idx;
      if ((prev_indices[0] != temp_edges[i].node[0]) ||
          (prev_indices[1] != temp_edges[i].node[1])) {

        prev_indices[0] = temp_edges[i].node[0];
        prev_indices[1] = temp_edges[i].node[1];
        edge_to_node[num_edges][0] = prev_indices[0];
        edge_to_node[num_edges][1] = prev_indices[1];
        ++num_edges;
      }

      elem_to_edge[curr_elem_to_edge_idx] = num_edges - 1;
    }
    edge_to_node =
      xrealloc(edge_to_node, num_edges * sizeof(*edge_to_node));

    edge_type = xmalloc(num_edges * sizeof(*edge_type));
    if (use_ll_edges) {
      for (size_t i = 0; i < num_edges; ++i) {
        double * edge_vertex_a = node_coords[edge_to_node[i][0]];
        double * edge_vertex_b = node_coords[edge_to_node[i][1]];
        int vertex_a_is_pole = check_pole(edge_vertex_a);
        int vertex_b_is_pole = check_pole(edge_vertex_b);
        int is_lon_edge =
          (vertex_a_is_pole ^ vertex_b_is_pole) ||
          (!vertex_a_is_pole && !vertex_b_is_pole &&
           check_lon_edge(edge_vertex_a, edge_vertex_b));
        int is_lat_edge =
          (vertex_a_is_pole && vertex_b_is_pole) ||
          check_lat_edge(edge_vertex_a, edge_vertex_b);
        YAC_ASSERT_F(
          is_lon_edge || is_lat_edge,
          "ERROR(yac_read_exodus_basic_grid_data_parallel): "
          "\"use_ll_edges == true\" but edge is neither lon nor lat "
          "((%lf,%lf,%lf),(%lf,%lf,%lf))",
          edge_vertex_a[0], edge_vertex_a[1], edge_vertex_a[2],
          edge_vertex_b[0], edge_vertex_b[1], edge_vertex_b[2]);
        edge_type[i] = (is_lon_edge)?YAC_LON_CIRCLE_EDGE:YAC_LAT_CIRCLE_EDGE;
      }
    } else {
      for (size_t i = 0; i < num_edges; ++i)
        edge_type[i] = YAC_GREAT_CIRCLE_EDGE;
    }
  }

  struct yac_basic_grid_data grid_data;
  grid_data.vertex_coordinates      = node_coords;
  grid_data.cell_ids                = elem_ids;
  grid_data.vertex_ids              = node_ids;
  grid_data.edge_ids                = NULL;
  grid_data.num_cells               = num_elem;
  grid_data.num_vertices            = num_nodes;
  grid_data.num_edges               = num_edges;
  grid_data.core_cell_mask          = NULL;
  grid_data.core_vertex_mask        = NULL;
  grid_data.core_edge_mask          = NULL;
  grid_data.num_vertices_per_cell   = num_nodes_per_elem;
  grid_data.num_cells_per_vertex    = num_elem_per_node;
  grid_data.cell_to_vertex          = elem_to_node;
  grid_data.cell_to_vertex_offsets  = generate_offsets(num_elem, num_nodes_per_elem);
  grid_data.cell_to_edge            = elem_to_edge;
  grid_data.cell_to_edge_offsets    = grid_data.cell_to_vertex_offsets;
  grid_data.vertex_to_cell          = node_to_elem;
  grid_data.vertex_to_cell_offsets  = generate_offsets(num_nodes, num_elem_per_node);
  grid_data.edge_to_vertex          = edge_to_node;
  grid_data.edge_type               = edge_type;
  grid_data.num_total_cells         = num_elem;
  grid_data.num_total_vertices      = num_nodes;
  grid_data.num_total_edges         = num_edges;

  return grid_data;
}

struct yac_basic_grid * yac_read_exodus_basic_grid_parallel(
  char const * filename, char const * gridname, int use_ll_edges,
  MPI_Comm comm) {

  return
    yac_basic_grid_new(
      gridname,
      yac_read_exodus_basic_grid_data_parallel(filename, use_ll_edges, comm));
}

#else

void yac_read_exodus_grid_information_parallel(
  const char * filename, MPI_Comm comm, yac_coordinate_pointer * node_coords,
  yac_int ** elem_ids, yac_int ** node_ids,
  size_t * num_elem, size_t * num_nodes,
  int ** num_nodes_per_elem, int ** num_elem_per_node,
  size_t ** elem_to_node, size_t ** node_to_elem) {

   UNUSED(filename);
   UNUSED(comm);
   UNUSED(node_coords);
   UNUSED(elem_ids);
   UNUSED(node_ids);
   UNUSED(num_elem);
   UNUSED(num_nodes);
   UNUSED(num_nodes_per_elem);
   UNUSED(num_elem_per_node);
   UNUSED(elem_to_node);
   UNUSED(node_to_elem);
   die(
     "ERROR(yac_read_exodus_grid_information_parallel): "
     "YAC is built without the NetCDF support");
}

struct yac_basic_grid_data yac_read_exodus_basic_grid_data_parallel(
  const char * filename, int use_ll_edges, MPI_Comm comm) {

   UNUSED(filename);
   UNUSED(use_ll_edges);
   UNUSED(comm);
   die(
     "ERROR(yac_read_exodus_basic_grid_data_parallel): "
     "YAC is built without the NetCDF support");

   return
    yac_generate_basic_grid_data_reg_2d(
      (size_t[]){0,0}, (int[]){0,0}, NULL, NULL);
}

struct yac_basic_grid * yac_read_exodus_basic_grid_parallel(
  char const * filename, char const * gridname, int use_ll_edges,
  MPI_Comm comm) {

   UNUSED(filename);
   UNUSED(gridname);
   UNUSED(use_ll_edges);
   UNUSED(comm);
   die(
     "ERROR(yac_read_exodus_basic_grid_parallel): "
     "YAC is built without the NetCDF support");

   return NULL;
}

#endif // YAC_NETCDF_ENABLED
