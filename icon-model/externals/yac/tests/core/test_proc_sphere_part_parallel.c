// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <mpi.h>

#include "tests.h"
#include "test_common.h"
#include "geometry.h"
#include "proc_sphere_part.h"
#include "yac_mpi.h"

/** \file test_proc_sphere_part_parallel.c
 *  \test
* This contains a test of the proc_sphere_part grid search algorithm.
*/

int main(void) {

  MPI_Init(NULL, NULL);

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  YAC_ASSERT(comm_size == 5, "ERROR wrong number of processes (has to be 5)")

  {
    // one process has no data
    double x_vertices[18] = {0,20,40,60,80,
                             100,120,140,160,180,
                             200,220,240,260,280,
                             300,320,340};
    double y_vertices[9] = {-80,-60,-40,-20,0,20,40,60,80};
    size_t local_start[5][2] = {{0,0},{7,0},{0,0},{0,3},{7,3}};
    size_t local_count[5][2] = {{10,6},{10,6},{0,0},{10,6},{10,6}};
    size_t num_vertices_a =
      local_count[comm_rank][0] * local_count[comm_rank][1];
    yac_coordinate_pointer vertices_a =
      xmalloc(num_vertices_a * sizeof(*vertices_a));

    for (size_t i = 0, k = 0; i < local_count[comm_rank][1]; ++i)
      for (size_t j = 0; j < local_count[comm_rank][0]; ++j, ++k)
        LLtoXYZ_deg(
          x_vertices[local_start[comm_rank][0]+j],
          y_vertices[local_start[comm_rank][1]+i],
          &(vertices_a[k][0]));

    yac_coordinate_pointer vertices[2] = {vertices_a, NULL};
    size_t num_vertices[2] = {num_vertices_a, 0};
    yac_int * global_ids_a = NULL, * global_ids_b = NULL;
    yac_int **global_vertex_ids[2] = {&global_ids_a, &global_ids_b};
    int * vertex_ranks_a, * vertex_ranks_b;
    int **vertex_ranks[2] = {&vertex_ranks_a, &vertex_ranks_b};
    struct proc_sphere_part_node * proc_sphere_part;
    yac_proc_sphere_part_new(
      vertices, num_vertices, &proc_sphere_part, global_vertex_ids,
      vertex_ranks, MPI_COMM_WORLD);

    free(vertex_ranks_a);
    free(vertex_ranks_b);
    free(global_ids_a);
    free(global_ids_b);
    free(vertices_a);

    yac_proc_sphere_part_node_delete(proc_sphere_part);
  }

  {
    double x_vertices[] = {  0, 10, 20, 30, 40, 50, 60, 70, 80, 90,
                           100,110,120,130,140,150,160,170,180,190,
                           200,210,220,230,240,250,260,270,280,290,
                           300,310,320,330,340,350};
    double y_vertices[] =
      {-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90};
    size_t local_start[5][2] = {{0,0},{0,5},{12,5},{24,5},{0,14}};
    size_t local_count[5][2] = {{36,5},{12,9},{12,9},{12,9},{36,5}};
    size_t num_vertices_a =
      local_count[comm_rank][0] * local_count[comm_rank][1];
    yac_coordinate_pointer vertices_a =
      xmalloc(num_vertices_a * sizeof(*vertices_a));

    for (size_t i = 0, k = 0; i < local_count[comm_rank][1]; ++i)
      for (size_t j = 0; j < local_count[comm_rank][0]; ++j, ++k)
        LLtoXYZ_deg(
          x_vertices[local_start[comm_rank][0]+j],
          y_vertices[local_start[comm_rank][1]+i],
          vertices_a[k]);

    yac_coordinate_pointer vertices[2] = {vertices_a, NULL};
    size_t num_vertices[2] = {num_vertices_a, 0};
    yac_int * global_ids_a = NULL, * global_ids_b = NULL;
    yac_int **global_vertex_ids[2] = {&global_ids_a, &global_ids_b};
    int * vertex_ranks_a, * vertex_ranks_b;
    int **vertex_ranks[2] = {&vertex_ranks_a, &vertex_ranks_b};
    struct proc_sphere_part_node * proc_sphere_part;
    yac_proc_sphere_part_new(
      vertices, num_vertices, &proc_sphere_part, global_vertex_ids,
      vertex_ranks, MPI_COMM_WORLD);

    free(vertex_ranks_a);
    free(vertex_ranks_b);
    free(global_ids_a);
    free(global_ids_b);
    free(vertices_a);

    double x_search = 45.0, y_search = 80.0;
    struct bounding_circle bnd_circle;
    LLtoXYZ_deg(x_search, y_search, bnd_circle.base_vector);
    bnd_circle.inc_angle.sin = sin(3.13);
    bnd_circle.inc_angle.cos = cos(3.13);

    int search_ranks[5], rank_count;
    yac_proc_sphere_part_do_bnd_circle_search(
      proc_sphere_part, bnd_circle, search_ranks, &rank_count);

    if (rank_count != 5)
      PUT_ERR("error in yac_proc_sphere_part_do_bnd_circle_search");

    yac_proc_sphere_part_node_delete(proc_sphere_part);
  }

  {
    double coords[5][2][3] =
      {{{1,0,0},{-1,0,0}},
       {{0,1,0}},
       {{0,-1,0}},
       {{0,0,1}},
       {{0,0,-1}}};
    size_t num_vertices_a[5] = {2,1,1,1,1};
    yac_coordinate_pointer vertices_a =
      xmalloc(num_vertices_a[comm_rank] * sizeof(*vertices_a));

    for (size_t i = 0; i < num_vertices_a[comm_rank]; ++i)
      for (int j = 0; j < 3; ++j)
        vertices_a[i][j] = coords[comm_rank][i][j];

    yac_coordinate_pointer vertices[2] = {vertices_a, NULL};
    size_t num_vertices[2] = {num_vertices_a[comm_rank], 0};
    yac_int * global_ids_a = NULL, * global_ids_b = NULL;
    yac_int **global_vertex_ids[2] = {&global_ids_a, &global_ids_b};
    int * vertex_ranks_a = NULL, * vertex_ranks_b = NULL;
    int **vertex_ranks[2] = {&vertex_ranks_a, &vertex_ranks_b};
    struct proc_sphere_part_node * proc_sphere_part;
    yac_proc_sphere_part_new(
      vertices, num_vertices, &proc_sphere_part, global_vertex_ids,
      vertex_ranks, MPI_COMM_WORLD);

    free(vertex_ranks_a);
    free(vertex_ranks_b);
    free(global_ids_a);
    free(global_ids_b);
    free(vertices_a);

    yac_proc_sphere_part_node_delete(proc_sphere_part);
  }

  {
    yac_coordinate_pointer vertices[2] = {NULL, NULL};
    size_t num_vertices[2] = {0, 0};
    yac_int * global_ids_a = NULL, * global_ids_b = NULL;
    yac_int **global_vertex_ids[2] = {&global_ids_a, &global_ids_b};
    int * vertex_ranks_a = NULL, * vertex_ranks_b = NULL;
    int **vertex_ranks[2] = {&vertex_ranks_a, &vertex_ranks_b};
    struct proc_sphere_part_node * proc_sphere_part;
    yac_proc_sphere_part_new(
      vertices, num_vertices, &proc_sphere_part, global_vertex_ids,
      vertex_ranks, MPI_COMM_WORLD);

    free(vertex_ranks_a);
    free(vertex_ranks_b);
    free(global_ids_a);
    free(global_ids_b);

    yac_proc_sphere_part_node_delete(proc_sphere_part);
  }

  {
    double x_vertices[] = {  0, 30, 60, 90, 120, 150, 180};
    double y_vertices[] = {-90,-60,-30,  0,  30,  60,  90};

    size_t local_start[2][5][2] = {{{0,0},{0,3},{3,3},{0,0},{0,0}},
                                   {{0,0},{0,0},{0,0},{3,0},{0,3}}};
    size_t local_count[2][5][2] = {{{7,4},{4,4},{4,4},{0,0},{0,0}},
                                   {{0,0},{0,0},{4,4},{4,4},{7,4}}};
    size_t num_vertices_[2];

    yac_coordinate_pointer vertices_[2];
    yac_int *global_vertex_ids_in[2];

    for (int grid_idx = 0; grid_idx < 2; ++grid_idx) {

      num_vertices_[grid_idx] =
        local_count[grid_idx][comm_rank][0] *
        local_count[grid_idx][comm_rank][1];

      vertices_[grid_idx] =
        xmalloc(num_vertices_[grid_idx] * sizeof(*(vertices_[0])));
      global_vertex_ids_in[grid_idx] =
        xmalloc(
          num_vertices_[grid_idx] * sizeof(*(global_vertex_ids_in[0])));

      for (size_t i = 0, k = 0; i < local_count[grid_idx][comm_rank][1]; ++i) {
        for (size_t j = 0; j < local_count[grid_idx][comm_rank][0]; ++j, ++k) {
          LLtoXYZ_deg(
            x_vertices[local_start[grid_idx][comm_rank][0]+j],
            y_vertices[local_start[grid_idx][comm_rank][1]+i],
            vertices_[grid_idx][k]);
          global_vertex_ids_in[grid_idx][k] =
            (local_start[grid_idx][comm_rank][1]+i) * 7 +
             local_start[grid_idx][comm_rank][0]+j;
        }
      }
    }

    int with_grid[2], with_global_ids[2];

    for (with_grid[0] = 0; with_grid[0] < 2; ++with_grid[0]) {

      for (with_global_ids[0] = 0; with_global_ids[0] < 2;
           with_global_ids[0]++) {

        for (with_grid[1] = 0; with_grid[1] < 2; ++with_grid[1]) {

          for (with_global_ids[1] = 0; with_global_ids[1] < 2;
               with_global_ids[1]++) {

            int * vertex_ranks_out[2] = {NULL, NULL};
            yac_int * global_vertex_ids_out[2] = {NULL,NULL};

            int **vertex_ranks[2];
            yac_int **global_vertex_ids[2];
            yac_coordinate_pointer vertices[2];
            size_t num_vertices[2];

            for (int grid_idx = 0; grid_idx < 2; ++grid_idx){
              vertices[grid_idx] =
                (with_grid[grid_idx] && (num_vertices_[grid_idx] > 0))?
                  vertices_[grid_idx]:NULL;
              num_vertices[grid_idx] =
                (with_grid[grid_idx]&&(num_vertices_[grid_idx] > 0))?
                  num_vertices_[grid_idx]:0;
              global_vertex_ids[grid_idx] =
                (with_global_ids[grid_idx])?
                  &(global_vertex_ids_in[grid_idx]):
                  &(global_vertex_ids_out[grid_idx]);
              vertex_ranks[grid_idx] = &(vertex_ranks_out[grid_idx]);
            }

            struct proc_sphere_part_node * proc_sphere_part;
            yac_proc_sphere_part_new(
              vertices, num_vertices, &proc_sphere_part, global_vertex_ids,
              vertex_ranks, MPI_COMM_WORLD);

            // check results
            for (int grid_idx = 0; grid_idx < 2; ++grid_idx) {

              int * search_ranks =
                xmalloc(num_vertices[grid_idx] * sizeof(*search_ranks));
              yac_proc_sphere_part_do_point_search(
                proc_sphere_part, vertices[grid_idx], num_vertices[grid_idx],
                search_ranks);
              for (size_t i = 0; i < num_vertices[grid_idx]; ++i)
                if (search_ranks[i] != vertex_ranks_out[grid_idx][i])
                  PUT_ERR("error in yac_proc_sphere_part_do_point_search");
              free(search_ranks);

              if (with_grid[grid_idx] && !with_global_ids[grid_idx]) {

                int ref_num_global_ids = 37;
                int min_global_id = INT_MAX, max_global_id = INT_MIN;
                int global_id_count[ref_num_global_ids];

                for (int i = 0; i < ref_num_global_ids; ++i)
                  global_id_count[i] = 0;

                for (size_t i = 0; i < num_vertices_[grid_idx]; ++i) {
                  if ((int)(global_vertex_ids_out[grid_idx][i]) < min_global_id)
                    min_global_id = (int)(global_vertex_ids_out[grid_idx][i]);
                  if ((int)(global_vertex_ids_out[grid_idx][i]) > max_global_id)
                    max_global_id = (int)(global_vertex_ids_out[grid_idx][i]);
                  global_id_count[global_vertex_ids_out[grid_idx][i]]++;
                }

                MPI_Allreduce(
                  MPI_IN_PLACE, &min_global_id, 1, MPI_INT, MPI_MIN,
                  MPI_COMM_WORLD);
                MPI_Allreduce(
                  MPI_IN_PLACE, &max_global_id, 1, MPI_INT, MPI_MAX,
                  MPI_COMM_WORLD);
                MPI_Allreduce(
                  MPI_IN_PLACE, global_id_count, ref_num_global_ids,
                  MPI_INT, MPI_SUM, MPI_COMM_WORLD);

                if (min_global_id != 0)
                  PUT_ERR("error in yac_proc_sphere_part_new (min_global_id)");
                if (max_global_id != ref_num_global_ids-1)
                  PUT_ERR("error in yac_proc_sphere_part_new (min_global_id)");

                int ref_num_owners_per_vertex[2][5][7*4] =
                  {{{7,7,7,7,7,7,7,
                     1,1,1,1,1,1,1,
                     1,1,1,1,1,1,1,
                     2,2,2,3,2,2,2},
                    {2,2,2,3, 1,1,1,2, 1,1,1,2, 8,8,8,8},
                    {3,2,2,2, 2,1,1,1, 2,1,1,1, 8,8,8,8},
                    {-1},{-1}},
                   {{-1},{-1},
                    {8,8,8,8, 1,1,1,2, 1,1,1,2, 2,2,2,3},
                    {8,8,8,8, 2,1,1,1, 2,1,1,1, 3,2,2,2},
                    {2,2,2,3,2,2,2,
                    1,1,1,1,1,1,1,
                    1,1,1,1,1,1,1,
                    7,7,7,7,7,7,7}}};
                for (size_t i = 0; i < num_vertices_[grid_idx]; ++i)
                  if (ref_num_owners_per_vertex[grid_idx][comm_rank][i] !=
                      global_id_count[global_vertex_ids_out[grid_idx][i]])
                    PUT_ERR(
                      "error in yac_proc_sphere_part_new (global_id count)");
              }
            }

            for (int grid_idx = 0; grid_idx < 2; ++grid_idx) {
              free(vertex_ranks_out[grid_idx]);
              if (!with_global_ids[grid_idx])
                free(global_vertex_ids_out[grid_idx]);
            }

            yac_proc_sphere_part_node_delete(proc_sphere_part);

          } // with_global_ids_b
        } // with_grid_b
      } // with_global_ids_a
    } // with_grid_a

    for (int grid_idx = 0; grid_idx < 2; ++grid_idx) {
      free(vertices_[grid_idx]);
      free(global_vertex_ids_in[grid_idx]);
    }
  }

  MPI_Finalize();

  return TEST_EXIT_CODE;
}
