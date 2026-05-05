// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <mpi.h>
#include <yaxt.h>

#include "tests.h"
#include "yac_xmap.h"

/** \file test_yac_xmap.c
 *  \test
 * This example tests the yac_xmap functionality.
 */

int main (void) {

  MPI_Init(NULL, NULL);

  xt_initialize(MPI_COMM_WORLD);

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Barrier(MPI_COMM_WORLD);

  enum {NUM_PROCS = 3};

  if (comm_size != 3) {
    PUT_ERR("ERROR: wrong number of processes");
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  { // test xmap without any exchanges
    struct remote_point_infos * point_infos[NUM_PROCS] = {NULL, NULL, NULL};
    size_t count[NUM_PROCS] = {0, 0, 0};

    yac_xmap xmap =
      yac_xmap_from_point_infos(
        point_infos[comm_rank], count[comm_rank], MPI_COMM_WORLD);

    Xt_redist redist = yac_xmap_generate_redist(xmap, MPI_INT);

    int src_data[1] = {0}, tgt_data[1] = {-1};
    xt_redist_s_exchange1(redist, src_data, tgt_data);

    int ref_tgt_data[NUM_PROCS][1] = {{-1},{-1},{-1}};

    for (size_t i = 0; i < count[comm_rank]; ++i)
      if (ref_tgt_data[comm_rank][i] != tgt_data[i])
        PUT_ERR("error in empty exchange");

    xt_redist_delete(redist);
    yac_xmap_delete(xmap);
  }

  { // test xmap with bcast operation
    struct remote_point_infos point_infos[NUM_PROCS][1] =
      {{(struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 0, .orig_pos = 0}}},
       {(struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 0, .orig_pos = 0}}},
       {(struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 0, .orig_pos = 0}}}};
    size_t count[NUM_PROCS] = {1, 1, 1};

    yac_xmap xmap =
      yac_xmap_from_point_infos(
        point_infos[comm_rank], count[comm_rank], MPI_COMM_WORLD);

    Xt_redist redist = yac_xmap_generate_redist(xmap, MPI_INT);

    int src_data[1] = {0}, tgt_data[1] = {-1};
    xt_redist_s_exchange1(redist, src_data, tgt_data);

    int ref_tgt_data[NUM_PROCS][1] = {{0},{0},{0}};

    for (size_t i = 0; i < count[comm_rank]; ++i)
      if (ref_tgt_data[comm_rank][i] != tgt_data[i])
        PUT_ERR("error in bcast exchange");

    xt_redist_delete(redist);
    yac_xmap_delete(xmap);
  }

  { // test xmap with gather operation
    struct remote_point_infos point_infos[NUM_PROCS][3] =
      {{(struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 0, .orig_pos = 0}},
        (struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 1, .orig_pos = 0}},
        (struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 2, .orig_pos = 0}}},
       {(struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = -1, .orig_pos = 999}}},
        {(struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = -1, .orig_pos = 999}}}};
    size_t count[NUM_PROCS] = {3, 0, 0};

    yac_xmap xmap =
      yac_xmap_from_point_infos(
        point_infos[comm_rank], count[comm_rank], MPI_COMM_WORLD);

    Xt_redist redist = yac_xmap_generate_redist(xmap, MPI_INT);

    int src_data[1] = {comm_rank}, tgt_data[3] = {-1, -1, -1};
    xt_redist_s_exchange1(redist, src_data, tgt_data);

    int ref_tgt_data[NUM_PROCS][3] = {{0, 1, 2},{-1},{-1}};

    for (size_t i = 0; i < count[comm_rank]; ++i)
      if (ref_tgt_data[comm_rank][i] != tgt_data[i])
        PUT_ERR("error in gather exchange");

    xt_redist_delete(redist);
    yac_xmap_delete(xmap);
  }

  { // test xmap with allgather operation
    struct remote_point_infos point_infos[NUM_PROCS][3] =
      {{(struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 0, .orig_pos = 0}},
        (struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 1, .orig_pos = 0}},
        (struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 2, .orig_pos = 0}}},
       {(struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 0, .orig_pos = 0}},
        (struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 1, .orig_pos = 0}},
        (struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 2, .orig_pos = 0}}},
        {(struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 0, .orig_pos = 0}},
        (struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 1, .orig_pos = 0}},
        (struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 2, .orig_pos = 0}}}};
    size_t count[NUM_PROCS] = {3, 3, 3};

    yac_xmap xmap =
      yac_xmap_from_point_infos(
        point_infos[comm_rank], count[comm_rank], MPI_COMM_WORLD);

    Xt_redist redist = yac_xmap_generate_redist(xmap, MPI_INT);

    int src_data[1] = {comm_rank}, tgt_data[3] = {-1, -1, -1};
    xt_redist_s_exchange1(redist, src_data, tgt_data);

    int ref_tgt_data[NUM_PROCS][3] = {{0, 1, 2},{0, 1, 2},{0, 1, 2}};

    for (size_t i = 0; i < count[comm_rank]; ++i)
      if (ref_tgt_data[comm_rank][i] != tgt_data[i])
        PUT_ERR("error in allgather exchange");

    xt_redist_delete(redist);
    yac_xmap_delete(xmap);
  }

  { // test xmap with all2all operation
    struct remote_point_infos point_infos[NUM_PROCS][3] =
      {{(struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 0, .orig_pos = 0}},
        (struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 1, .orig_pos = 0}},
        (struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 2, .orig_pos = 0}}},
       {(struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 0, .orig_pos = 1}},
        (struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 1, .orig_pos = 1}},
        (struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 2, .orig_pos = 1}}},
        {(struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 0, .orig_pos = 2}},
        (struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 1, .orig_pos = 2}},
        (struct remote_point_infos){
           .count = 1,
             .data.single =
               (struct remote_point_info){
                 .rank = 2, .orig_pos = 2}}}};
    size_t count[NUM_PROCS] = {3, 3, 3};

    yac_xmap xmap =
      yac_xmap_from_point_infos(
        point_infos[comm_rank], count[comm_rank], MPI_COMM_WORLD);

    Xt_redist redist = yac_xmap_generate_redist(xmap, MPI_INT);

    int src_data[3] = {comm_rank * NUM_PROCS + 0,
                       comm_rank * NUM_PROCS + 1,
                       comm_rank * NUM_PROCS + 2}, tgt_data[3] = {-1, -1, -1};
    xt_redist_s_exchange1(redist, src_data, tgt_data);

    int ref_tgt_data[NUM_PROCS][3] = {{0, 3, 6},{1, 4, 7},{2, 5, 8}};

    for (size_t i = 0; i < count[comm_rank]; ++i)
      if (ref_tgt_data[comm_rank][i] != tgt_data[i])
        PUT_ERR("error in all2all exchange");

    xt_redist_delete(redist);
    yac_xmap_delete(xmap);
  }

  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}
