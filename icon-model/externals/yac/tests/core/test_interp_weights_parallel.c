// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <yaxt.h>
#include <netcdf.h>
#include <string.h>

#include "tests.h"
#include "test_common.h"
#include "geometry.h"
#include "read_icon_grid.h"
#include "interp_method_internal.h"
#include "interp_method_file.h"
#include "interp_weights_internal.h"
#include "yac_mpi.h"
#include "dist_grid_utils.h"
#include "io_utils.h"
#include "interpolation_exchange.h"

/** \file test_interp_weights_parallel.c
 *  \test
 * This contains some examples on how to use interp_weights.
 */

enum weight_type {
  NONE = 0,
  DIRECT = 1,
  SUM = 2,
  WSUM = 3,
  FIXED = 4,
};

#define FALLBACK_VALUE (-1.0)
#define FIXED_VALUE (999.0)
#define TOL (1e-6)

static void utest_get_basic_grid_data(
  char * filename, size_t * num_cells, size_t * num_vertices,
  size_t * num_edges);
static int utest_check_fixed_results(
  struct yac_basic_grid_data * grid, int is_src,
  struct yac_interp_weights * weights,
  enum yac_interp_weights_reorder_type reorder_type,
  double * field_data);
static int utest_check_direct_results(
  struct yac_basic_grid_data * grid, int is_src,
  struct yac_interp_weights * weights,
  enum yac_interp_weights_reorder_type reorder_type, size_t collection_size,
  double *** src_data, double *** src_frac_masks, double ** tgt_data);
static int utest_check_results_ref(
  int is_src, struct yac_interp_weights * weights,
  enum yac_interp_weights_reorder_type reorder_type, double * src_data,
  double * tgt_data, double * ref_tgt_data, size_t tgt_size,
  int is_active_src, int is_active_tgt);

static void utest_get_basic_weight_file_info(
  char const * weight_file_name, int * contains_links, int * contains_fixed);
static void on_existing_abort_handler(
  MPI_Comm comm, char const * msg, char const * source, int line);
static char const * weight_file_name_on_existing =
  "test_interp_weights_parallel_on_existing.nc";

static void utest_interpolation_execute_raw(
  struct yac_interpolation_exchange * interpolation_exchange,
  struct yac_interp_weights_data interp_weights_data,
  size_t collection_size, double *** src_fields, double ** tgt_field);

static void utest_interpolation_execute_frac_raw(
  struct yac_interpolation_exchange * interpolation_exchange,
  struct yac_interp_weights_data interp_weights_data,
  size_t collection_size, double *** src_fields, double *** src_frac_masks,
  double ** tgt_field);

static struct yac_basic_grid_data
  utest_generate_dummy_grid_data(
    size_t num_cells, size_t global_num_cells, size_t num_cells_offset);

static char const * grid_names[2] = {"grid_a", "grid_b"};

int main(int argc, char** argv) {

  MPI_Init(NULL, NULL);

  xt_initialize(MPI_COMM_WORLD);

  int comm_rank, comm_size;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

  if (comm_size != 5) {
    PUT_ERR("This test requires 5 processes\n");
    xt_finalize();
    MPI_Finalize();
    return TEST_EXIT_CODE;
  }

  enum yac_interp_weights_reorder_type reorder_types[] = {
    YAC_MAPPING_ON_SRC,
    YAC_MAPPING_ON_TGT
  };
  size_t reorder_type_count =
    sizeof(reorder_types) / sizeof(reorder_types[0]);

  { // test with 2 target and 3 source processes
    // * number source and target points is identical
    // * first target process owns all even target cells and the other
    //   all odd target cells
    // * each source process owns one third of the source points
    //   (they are consecutive)
    // * each process generates a part of the weights for a consecutive part
    //   of the target points
    // * checks interpolation generation for cases, where target status is
    //   deactivated on a subset of target processes
    int is_src = comm_rank < 3;

    enum {
      SRC_PROC_COUNT = 3,
      TGT_PROC_COUNT = 2,
    };
    enum {
      GLOBAL_SIZE =
        SRC_PROC_COUNT * TGT_PROC_COUNT *
        (SRC_PROC_COUNT + TGT_PROC_COUNT) * 2,
      LOCAL_SIZE =
        SRC_PROC_COUNT * TGT_PROC_COUNT * 2,
    };

    double wsum_weights[3*LOCAL_SIZE];
    double sum_weights[3*LOCAL_SIZE];
    double direct_weights[LOCAL_SIZE];
    size_t wsum_num_src_per_tgt[LOCAL_SIZE];
    size_t direct_num_src_per_tgt[LOCAL_SIZE];
    struct remote_points * default_tgts = xmalloc(1 * sizeof(*default_tgts));
    default_tgts->data = xmalloc(LOCAL_SIZE * sizeof(*default_tgts->data));
    struct remote_points none_tgts;
    struct remote_point * wsum_srcs =
      xmalloc(3 * LOCAL_SIZE * sizeof(*wsum_srcs));
    struct remote_point * direct_srcs =
      xmalloc(LOCAL_SIZE * sizeof(*direct_srcs));

    default_tgts->count = LOCAL_SIZE;
    for (size_t i = 0; i < LOCAL_SIZE; ++i) {
      for (size_t j = 0; j < 3; ++j) {
        wsum_weights[3 * i + j] = (double)(j+1);
        sum_weights[3 * i + j] = 1.0;
      }
      direct_weights[i] = 1.0;
      wsum_num_src_per_tgt[i] = 3;
      direct_num_src_per_tgt[i] = 1;

      yac_int tgt_global_id = (yac_int)(comm_rank * LOCAL_SIZE + i);
      default_tgts->data[i].global_id = tgt_global_id;
      default_tgts->data[i].data.count = 1;
      default_tgts->data[i].data.data.single.rank =
        (int)(tgt_global_id & 1) + (int)SRC_PROC_COUNT;
      default_tgts->data[i].data.data.single.orig_pos =
        (uint64_t)(tgt_global_id / 2);

      for (size_t j = 0; j < 3; ++j) {
        yac_int src_global_id =
          (tgt_global_id + (yac_int)j)%((yac_int)GLOBAL_SIZE);
        wsum_srcs[3 * i + j].global_id = src_global_id;
        wsum_srcs[3 * i + j].data.count = 1;
        wsum_srcs[3 * i + j].data.data.single.rank =
          (int)src_global_id % (int)SRC_PROC_COUNT;
        wsum_srcs[3 * i + j].data.data.single.orig_pos =
          (uint64_t)(src_global_id / SRC_PROC_COUNT);
      }

      yac_int src_global_id = tgt_global_id;
      direct_srcs[i].global_id = src_global_id;
      direct_srcs[i].data.count = 1;
      direct_srcs[i].data.data.single.rank =
        (int)src_global_id % (int)SRC_PROC_COUNT;
      direct_srcs[i].data.data.single.orig_pos =
        (uint64_t)(src_global_id / SRC_PROC_COUNT);
    }
    none_tgts.count = 0;
    none_tgts.data = 0;

    double * out_data =
      (is_src)?xmalloc((GLOBAL_SIZE / SRC_PROC_COUNT) * sizeof(*out_data)):NULL;
    if (is_src)
      for (size_t i = 0; i < GLOBAL_SIZE / SRC_PROC_COUNT; ++i)
        out_data[i] = (double)((yac_int)comm_rank + i * SRC_PROC_COUNT);
    double * in_data =
      (is_src)?NULL:xmalloc((GLOBAL_SIZE / TGT_PROC_COUNT) * sizeof(*in_data));
    double * ref_in_data =
      (is_src)?NULL:xmalloc((GLOBAL_SIZE / TGT_PROC_COUNT) * sizeof(*ref_in_data));

    struct yac_basic_grid_data grid_data =
      utest_generate_dummy_grid_data(
        GLOBAL_SIZE / ((is_src)?SRC_PROC_COUNT:TGT_PROC_COUNT), GLOBAL_SIZE,
        (size_t)comm_rank - ((is_src)?0:SRC_PROC_COUNT));

    struct yac_basic_grid * src_grid, * tgt_grid;
    if (is_src) {
      src_grid = yac_basic_grid_new(grid_names[0], grid_data);
      tgt_grid = yac_basic_grid_empty_new(grid_names[1]);
    } else {
      src_grid = yac_basic_grid_empty_new(grid_names[0]);
      tgt_grid = yac_basic_grid_new(grid_names[1], grid_data);
    }

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(src_grid, tgt_grid, MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(
        grid_pair, yac_basic_grid_get_name(src_grid),
        yac_basic_grid_get_name(tgt_grid),
        num_src_fields, src_fields, tgt_field);

    for (size_t i = 0; i < (1 << (2 * (SRC_PROC_COUNT + TGT_PROC_COUNT))); ++i) {

      enum weight_type weight_type[SRC_PROC_COUNT + TGT_PROC_COUNT];
      int max_weight_type = 0;

      for (size_t j = 0; j < SRC_PROC_COUNT + TGT_PROC_COUNT; ++j) {
        weight_type[j] = (enum weight_type)((i >> (3 * j)) & 3);
        if ((int)(weight_type[j]) > max_weight_type)
          max_weight_type = (int)(weight_type[j]);
      }

      struct remote_points * tgts;
      size_t * num_src_per_tgt;
      struct remote_point * srcs;
      double * w;

      switch (weight_type[comm_rank]) {
        case(NONE):
        default:
          tgts = &none_tgts;
          num_src_per_tgt = NULL;
          srcs = NULL;
          w = NULL;
          break;
        case(DIRECT):
          tgts = default_tgts;
          num_src_per_tgt = direct_num_src_per_tgt;
          srcs = direct_srcs;
          w = direct_weights;
          break;
        case(SUM):
          tgts = default_tgts;
          num_src_per_tgt = wsum_num_src_per_tgt;
          srcs = wsum_srcs;
          w = sum_weights;
          break;
        case(WSUM):
          tgts = default_tgts;
          num_src_per_tgt = wsum_num_src_per_tgt;
          srcs = wsum_srcs;
          w = wsum_weights;
          break;
      }

      for (int j = MAX(max_weight_type, DIRECT); j <= (int)WSUM; ++j) {

        struct yac_interp_weights * weights =
          yac_interp_weights_new(
            MPI_COMM_WORLD, YAC_LOC_CELL,
            (enum yac_location[]){YAC_LOC_CELL}, 1);

        switch (j) {
          case(NONE):
          default:
          break;
          case(DIRECT):
            yac_interp_weights_add_direct(weights, tgts, srcs);
            break;
          case(SUM):
            yac_interp_weights_add_sum(weights, tgts, num_src_per_tgt, srcs);
            break;
          case(WSUM):
            yac_interp_weights_add_wsum(
              weights, tgts, num_src_per_tgt, srcs, w);
            break;
        }

        enum {ACTIVE_TARGET_CONFIG_COUNT = 4};
        int active_target_configs[ACTIVE_TARGET_CONFIG_COUNT][TGT_PROC_COUNT] =
          {{0,0},{0,1},{1,0},{1,1}};

        for (size_t tgt_config_idx = 0;
             tgt_config_idx < ACTIVE_TARGET_CONFIG_COUNT; ++tgt_config_idx) {

          int is_active_src = 0;
          int is_active_tgt = 0;

          if (is_src) {

            is_active_src = 1;

          } else {

            int tgt_rank = (yac_int)(comm_rank - (int)SRC_PROC_COUNT);
            is_active_tgt = active_target_configs[tgt_config_idx][tgt_rank];

            for (size_t k = 0; k < GLOBAL_SIZE / TGT_PROC_COUNT; ++k) {
              yac_int tgt_global_id = (yac_int)tgt_rank + k * TGT_PROC_COUNT;
              // which rank generates the weights
              int weight_rank = (int)(tgt_global_id / (yac_int)LOCAL_SIZE);
              double ref_value;
              switch (weight_type[weight_rank]) {
                case(NONE):
                default:
                  ref_value = FALLBACK_VALUE;
                  break;
                case(DIRECT):
                  ref_value = (double)tgt_global_id;
                  break;
                case(SUM):
                  ref_value = (double)(((tgt_global_id + 0) % GLOBAL_SIZE) +
                                       ((tgt_global_id + 1) % GLOBAL_SIZE) +
                                       ((tgt_global_id + 2) % GLOBAL_SIZE));
                  break;
                case(WSUM):
                  ref_value = (double)(1 * ((tgt_global_id + 0) % GLOBAL_SIZE) +
                                       2 * ((tgt_global_id + 1) % GLOBAL_SIZE) +
                                       3 * ((tgt_global_id + 2) % GLOBAL_SIZE));
                  break;
              };
              ref_in_data[k] = is_active_tgt?ref_value:FALLBACK_VALUE;
              in_data[k] = FALLBACK_VALUE;
            }
          }

          // checking with different reorder types and collection sizes
          for (size_t reorder_type = 0; reorder_type < reorder_type_count;
              ++reorder_type)
            if (utest_check_results_ref(
                  is_src, weights, reorder_types[reorder_type],
                  out_data, in_data, ref_in_data, GLOBAL_SIZE / TGT_PROC_COUNT,
                  is_active_src, is_active_tgt))
              PUT_ERR("ERROR(yac_interp_weights_add_*): "
                      "invalid interpolation result");
        }

        MPI_Comm weights_comm =
          yac_interp_weights_get_comm(weights);

        int compare_result;
        MPI_Comm_compare(weights_comm, MPI_COMM_WORLD, &compare_result);
        if ((compare_result != MPI_IDENT) &&
            (compare_result != MPI_CONGRUENT))
          PUT_ERR("ERROR(yac_interp_weights_get_comm): wrong communicator");

        yac_interp_weights_delete(weights);
      }
    }
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(tgt_grid);
    yac_basic_grid_delete(src_grid);
    free(ref_in_data);
    free(in_data);
    free(out_data);
    free(direct_srcs);
    free(wsum_srcs);
    free(default_tgts->data);
    free(default_tgts);
  }

  { // test with 2 target and 3 source processes
    // * checks interpolation generation for cases, where source status is
    //   deactivated on a subset of source processes
    int is_src = comm_rank < 3;

    enum {
      SRC_PROC_COUNT = 3,
      TGT_PROC_COUNT = 2,
    };
    enum {
      SRC_NUM_POINTS = 3,
      TGT_NUM_POINTS = 1,
    };

    int src_ranks[SRC_PROC_COUNT] = {0, 1, 2};
    int tgt_ranks[TGT_PROC_COUNT] = {3, 4};
    double wsum_weights[SRC_NUM_POINTS];
    double sum_weights[SRC_NUM_POINTS];
    double direct_weights[TGT_NUM_POINTS];
    size_t wsum_num_src_per_tgt[TGT_NUM_POINTS];
    size_t direct_num_src_per_tgt[TGT_NUM_POINTS];
    struct remote_points * default_tgts =
      xmalloc(TGT_NUM_POINTS * sizeof(*default_tgts));
    default_tgts->data = xmalloc(TGT_NUM_POINTS * sizeof(*default_tgts->data));
    struct remote_points none_tgts;
    struct remote_point * wsum_srcs =
      xmalloc(SRC_NUM_POINTS * sizeof(*wsum_srcs));
    struct remote_point * direct_srcs =
      xmalloc(TGT_NUM_POINTS * sizeof(*direct_srcs));

    default_tgts->count = TGT_NUM_POINTS;
    for (size_t tgt_point_idx = 0; tgt_point_idx < TGT_NUM_POINTS;
         ++tgt_point_idx) {
      for (size_t src_point_idx = 0; src_point_idx < SRC_NUM_POINTS;
           ++src_point_idx) {
        wsum_weights[SRC_NUM_POINTS * tgt_point_idx + src_point_idx] =
          (double)(src_point_idx+1);
        sum_weights[SRC_NUM_POINTS * tgt_point_idx + src_point_idx] = 1.0;
      }
      direct_weights[tgt_point_idx] = 1.0;
      wsum_num_src_per_tgt[tgt_point_idx] = SRC_NUM_POINTS;
      direct_num_src_per_tgt[tgt_point_idx] = 1;

      yac_int tgt_global_id = (yac_int)tgt_point_idx;
      default_tgts->data[tgt_point_idx].global_id = tgt_global_id;
      default_tgts->data[tgt_point_idx].data.count = TGT_PROC_COUNT;
      default_tgts->data[tgt_point_idx].data.data.multi =
        xmalloc(
          TGT_PROC_COUNT *
          sizeof(*(default_tgts->data[tgt_point_idx].data.data.multi)));
      for (size_t tgt_rank_idx = 0; tgt_rank_idx < TGT_PROC_COUNT;
           ++tgt_rank_idx) {
        default_tgts->data[tgt_point_idx].data.data.multi[tgt_rank_idx].
          rank = tgt_ranks[tgt_rank_idx];
        default_tgts->data[tgt_point_idx].data.data.multi[tgt_rank_idx].
          orig_pos = tgt_point_idx;
      }

      for (size_t src_point_idx = 0; src_point_idx < SRC_NUM_POINTS;
           ++src_point_idx) {
        yac_int src_global_id = (yac_int)src_point_idx;
        wsum_srcs[SRC_NUM_POINTS * tgt_point_idx + src_point_idx].
          global_id = src_global_id;
        wsum_srcs[SRC_NUM_POINTS * tgt_point_idx + src_point_idx].
          data.count = SRC_PROC_COUNT;
        wsum_srcs[SRC_NUM_POINTS * tgt_point_idx + src_point_idx].
          data.data.multi =
            xmalloc(SRC_PROC_COUNT * sizeof(*(wsum_srcs->data.data.multi)));
        for (size_t src_rank_idx = 0; src_rank_idx < SRC_PROC_COUNT;
             ++src_rank_idx) {
          wsum_srcs[SRC_NUM_POINTS * tgt_point_idx + src_point_idx].
            data.data.multi[src_rank_idx].rank = src_ranks[src_rank_idx];
          wsum_srcs[SRC_NUM_POINTS * tgt_point_idx + src_point_idx].
            data.data.multi[src_rank_idx].orig_pos =  (uint64_t)src_point_idx;
        }
      }

      yac_int src_global_id = tgt_global_id;
      direct_srcs[tgt_point_idx].global_id = src_global_id;
      direct_srcs[tgt_point_idx].data.count = SRC_PROC_COUNT;
      direct_srcs[tgt_point_idx].data.data.multi =
        xmalloc(SRC_PROC_COUNT * sizeof(*(direct_srcs->data.data.multi)));
      for (size_t src_rank_idx = 0; src_rank_idx < SRC_PROC_COUNT;
           ++src_rank_idx) {
        direct_srcs[tgt_point_idx].data.data.multi[src_rank_idx].rank =
          src_ranks[src_rank_idx];
        direct_srcs[tgt_point_idx].data.data.multi[src_rank_idx].orig_pos = 0;
      }
    }
    none_tgts.count = 0;
    none_tgts.data = 0;

    double * out_data = (is_src)?xmalloc(SRC_NUM_POINTS * sizeof(*out_data)):NULL;
    if (is_src)
      for (size_t i = 0; i < SRC_NUM_POINTS; ++i)
        out_data[i] = (double)((yac_int)comm_rank + i * SRC_PROC_COUNT);
    double * in_data =
      (is_src)?NULL:xmalloc(TGT_NUM_POINTS * sizeof(*in_data));
    double * ref_in_data =
      (is_src)?NULL:xmalloc(TGT_NUM_POINTS * sizeof(*ref_in_data));

    struct yac_basic_grid_data grid_data =
      utest_generate_dummy_grid_data(
        is_src?SRC_NUM_POINTS:TGT_NUM_POINTS,
        is_src?SRC_NUM_POINTS:TGT_NUM_POINTS, 0);

    struct yac_basic_grid * src_grid, * tgt_grid;
    if (is_src) {
      src_grid = yac_basic_grid_new(grid_names[0], grid_data);
      tgt_grid = yac_basic_grid_empty_new(grid_names[1]);
    } else {
      src_grid = yac_basic_grid_empty_new(grid_names[0]);
      tgt_grid = yac_basic_grid_new(grid_names[1], grid_data);
    }

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(src_grid, tgt_grid, MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(
        grid_pair, yac_basic_grid_get_name(src_grid),
        yac_basic_grid_get_name(tgt_grid),
        num_src_fields, src_fields, tgt_field);

    enum weight_type weight_types[] = {NONE, FIXED, DIRECT, SUM, WSUM};
    enum {NUM_WEIGHT_TYPES = sizeof(weight_types) / sizeof(weight_types[0])};

    for (size_t weight_type_idx = 0; weight_type_idx < NUM_WEIGHT_TYPES;
         ++weight_type_idx) {

      enum weight_type weight_type = weight_types[weight_type_idx];

      struct remote_points * tgts;
      size_t * num_src_per_tgt;
      struct remote_point * srcs;
      double * w;

      struct yac_interp_weights * weights =
        yac_interp_weights_new(
          MPI_COMM_WORLD, YAC_LOC_CELL, (enum yac_location[]){YAC_LOC_CELL}, 1);

      switch (weight_type) {
        case(NONE):
        default:
          tgts = &none_tgts;
          num_src_per_tgt = NULL;
          srcs = NULL;
          w = NULL;
          break;
        case (FIXED):
          tgts = default_tgts;
          yac_interp_weights_add_fixed(weights, tgts, FIXED_VALUE);
          break;
        case(DIRECT):
          tgts = default_tgts;
          num_src_per_tgt = direct_num_src_per_tgt;
          srcs = direct_srcs;
          w = direct_weights;
          yac_interp_weights_add_direct(weights, tgts, srcs);
          break;
        case(SUM):
          tgts = default_tgts;
          num_src_per_tgt = wsum_num_src_per_tgt;
          srcs = wsum_srcs;
          w = sum_weights;
          yac_interp_weights_add_sum(weights, tgts, num_src_per_tgt, srcs);
          break;
        case(WSUM):
          tgts = default_tgts;
          num_src_per_tgt = wsum_num_src_per_tgt;
          srcs = wsum_srcs;
          w = wsum_weights;
          yac_interp_weights_add_wsum(
            weights, tgts, num_src_per_tgt, srcs, w);
          break;
      }

      for (int active_source_rank = 0; active_source_rank < SRC_PROC_COUNT;
           ++active_source_rank) {

        int is_active_src = 0;
        int is_active_tgt = 0;

        if (is_src) {

          is_active_src = active_source_rank == comm_rank;

        } else {

          is_active_tgt = 1;

          for (size_t tgt_point_idx = 0; tgt_point_idx < TGT_NUM_POINTS;
               ++tgt_point_idx) {
            switch (weight_type) {
              case(NONE):
              default:
                ref_in_data[tgt_point_idx] = FALLBACK_VALUE;
                break;
              case(FIXED):
                ref_in_data[tgt_point_idx] = FIXED_VALUE;
                break;
              case(DIRECT):
                ref_in_data[tgt_point_idx] =
                  (double)(tgt_point_idx + active_source_rank);
                break;
              case(SUM):
                ref_in_data[tgt_point_idx] =
                  (double)(
                    (tgt_point_idx + active_source_rank + 0 * SRC_NUM_POINTS) +
                    (tgt_point_idx + active_source_rank + 1 * SRC_NUM_POINTS) +
                    (tgt_point_idx + active_source_rank + 2 * SRC_NUM_POINTS));
                break;
              case(WSUM):
                ref_in_data[tgt_point_idx] =
                  (double)(
                    1 * (tgt_point_idx + active_source_rank + 0 * SRC_NUM_POINTS) +
                    2 * (tgt_point_idx + active_source_rank + 1 * SRC_NUM_POINTS) +
                    3 * (tgt_point_idx + active_source_rank + 2 * SRC_NUM_POINTS));
                break;
            };
            in_data[tgt_point_idx] = FALLBACK_VALUE;
          }
        }

        // checking with different reorder types and collection sizes
        for (size_t reorder_type = 0; reorder_type < reorder_type_count;
            ++reorder_type)
          if (utest_check_results_ref(
                is_src, weights, reorder_types[reorder_type],
                out_data, in_data, ref_in_data, TGT_NUM_POINTS,
                is_active_src, is_active_tgt))
            PUT_ERR("ERROR(yac_interp_weights_add_*): "
                    "invalid interpolation result");
      }

      MPI_Comm weights_comm =
        yac_interp_weights_get_comm(weights);

      int compare_result;
      MPI_Comm_compare(weights_comm, MPI_COMM_WORLD, &compare_result);
      if ((compare_result != MPI_IDENT) &&
          (compare_result != MPI_CONGRUENT))
        PUT_ERR("ERROR(yac_interp_weights_get_comm): wrong communicator");

      yac_interp_weights_delete(weights);
    }

    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(tgt_grid);
    yac_basic_grid_delete(src_grid);
    free(ref_in_data);
    free(in_data);
    free(out_data);
    for (size_t tgt_point_idx = 0; tgt_point_idx < TGT_NUM_POINTS;
         ++tgt_point_idx) {
      free(direct_srcs[tgt_point_idx].data.data.multi);
      free(default_tgts->data[tgt_point_idx].data.data.multi);
      for (size_t src_point_idx = 0; src_point_idx < SRC_NUM_POINTS;
            ++src_point_idx) {
        free(
          wsum_srcs[
            SRC_NUM_POINTS * tgt_point_idx + src_point_idx].data.data.multi);
      }
    }
    free(direct_srcs);
    free(wsum_srcs);
    free(default_tgts->data);
    free(default_tgts);
  }

  {
    int is_src = comm_rank < (comm_size / 2);

    MPI_Comm local_grid_comm;
    int local_grid_comm_size;
    MPI_Comm_split(MPI_COMM_WORLD, is_src, 0, &local_grid_comm);
    MPI_Comm_size(local_grid_comm, &local_grid_comm_size);

    set_even_io_rank_list(local_grid_comm);

    setenv("YAC_YAXT_EXCHANGER", "irecv_isend", 1);

    if (argc != 2) {
      PUT_ERR("ERROR: missing grid file directory");
      xt_finalize();
      MPI_Finalize();
      return TEST_EXIT_CODE;
    }

    char * filenames[2];
    char * grid_filenames[] =
      {"icon_grid_0030_R02B03_G.nc", "icon_grid_0043_R02B04_G.nc"};
    for (int i = 0; i < 2; ++i)
      filenames[i] =
        strcat(
          strcpy(
            malloc(strlen(argv[1]) + strlen(grid_filenames[i]) + 2), argv[1]),
          grid_filenames[i]);

    struct yac_basic_grid_data grid_data =
      yac_read_icon_basic_grid_data_parallel(
        filenames[is_src], local_grid_comm);
    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[is_src], grid_data),
       yac_basic_grid_empty_new(grid_names[is_src^1])};

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(
        grid_pair, grid_names[1], grid_names[0], num_src_fields,
        src_fields, tgt_field);

    struct yac_interp_weights * weights =
      yac_interp_weights_new(
        MPI_COMM_WORLD, YAC_LOC_CELL, (enum yac_location[]){YAC_LOC_CELL}, 1);

    size_t num_tgt_cells, num_tgt_vertices, num_tgt_edges;
    utest_get_basic_grid_data(
      filenames[0], &num_tgt_cells, &num_tgt_vertices, &num_tgt_edges);
    size_t num_src_cells, num_src_vertices, num_src_edges;
    utest_get_basic_grid_data(
      filenames[1], &num_src_cells, &num_src_vertices, &num_src_edges);

    {
      size_t count;
      size_t * tgt_points;
      yac_interp_grid_get_tgt_points(interp_grid, &tgt_points, &count);
      yac_int * tgt_global_ids = xmalloc(count * sizeof(*tgt_global_ids));
      size_t * size_t_buffer = xmalloc(2 * count * sizeof(*size_t_buffer));
      size_t * odd_tgt_points = size_t_buffer;
      size_t * even_tgt_points = size_t_buffer + count;
      yac_interp_grid_get_tgt_global_ids(
        interp_grid, tgt_points, count, tgt_global_ids);

      size_t odd_count = 0;
      size_t even_count = 0;
      for (size_t i = 0; i < count; ++i) {
        if (tgt_global_ids[i] & 1) odd_tgt_points[odd_count++] = tgt_points[i];
        else even_tgt_points[even_count++] = tgt_points[i];
      }

      struct remote_points odd_tgts = {
        .data =
          yac_interp_grid_get_tgt_remote_points(
            interp_grid, odd_tgt_points, odd_count),
        .count = odd_count};
      struct remote_points even_tgts = {
        .data =
          yac_interp_grid_get_tgt_remote_points(
            interp_grid, even_tgt_points, even_count),
        .count = even_count};

      yac_interp_weights_add_fixed(weights, &odd_tgts, -1.0);
      yac_interp_weights_add_fixed(weights, &even_tgts, -2.0);

      free(odd_tgts.data);
      free(even_tgts.data);
      free(size_t_buffer);
      free(tgt_global_ids);
      free(tgt_points);
    }

    double * field_data =
      xmalloc(grid_data.num_cells * sizeof(*field_data));

    for (size_t reorder_type = 0; reorder_type < reorder_type_count;
         ++reorder_type)
      if (utest_check_fixed_results(
            &grid_data, is_src, weights,
            reorder_types[reorder_type], field_data))
        PUT_ERR("ERROR(yac_interp_weights_add_fixed): "
                "invalid interpolation result");

    // check writing and reading
    {
      char const * weight_file_name =
        "test_interp_weights_parallel_fixed.nc";
      // write weights to file
      yac_interp_weights_write_to_file(
        weights, weight_file_name, grid_names[1], grid_names[0],
        num_src_cells, num_tgt_cells, YAC_WEIGHT_FILE_ERROR);

      // read weights
      struct interp_method * method_stack[2] = {
        yac_interp_method_file_new(
          weight_file_name, YAC_INTERP_FILE_ON_MISSING_FILE_DEFAULT,
          YAC_INTERP_FILE_ON_SUCCESS_DEFAULT), NULL};
      struct yac_interp_weights * weights_from_file =
        yac_interp_method_do_search(method_stack, interp_grid);
      yac_interp_method_delete(method_stack);

      if (utest_check_fixed_results(
            &grid_data, is_src, weights_from_file, YAC_MAPPING_ON_SRC, field_data))
        PUT_ERR("ERROR(yac_interp_weights_add_fixed): "
                "invalid interpolation result (read from file)");
      if (comm_rank == 0) unlink(weight_file_name);
      yac_interp_weights_delete(weights_from_file);
    }

    free(field_data);

    yac_interp_weights_delete(weights);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(grids[1]);
    yac_basic_grid_delete(grids[0]);
    free(filenames[1]);
    free(filenames[0]);
    MPI_Comm_free(&local_grid_comm);

    unsetenv("YAC_YAXT_EXCHANGER");
  }

  {
    int is_src = comm_rank < (comm_size / 2);

    MPI_Comm local_grid_comm;
    int local_grid_comm_size;
    MPI_Comm_split(MPI_COMM_WORLD, is_src, 0, &local_grid_comm);
    MPI_Comm_size(local_grid_comm, &local_grid_comm_size);

    set_even_io_rank_list(local_grid_comm);

    if (argc != 2) {
      PUT_ERR("ERROR: missing grid file directory");
      xt_finalize();
      MPI_Finalize();
      return TEST_EXIT_CODE;
    }

    char * filenames[2];
    char * grid_filenames[] =
      {"icon_grid_0030_R02B03_G.nc", "icon_grid_0043_R02B04_G.nc"};
    for (int i = 0; i < 2; ++i)
      filenames[i] =
        strcat(
          strcpy(
            malloc(strlen(argv[1]) + strlen(grid_filenames[i]) + 2), argv[1]),
          grid_filenames[i]);

    struct yac_basic_grid_data grid_data =
      yac_read_icon_basic_grid_data_parallel(
        filenames[is_src], local_grid_comm);
    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[is_src], grid_data),
       yac_basic_grid_empty_new(grid_names[is_src^1])};

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(
        grid_pair, grid_names[1], grid_names[0], num_src_fields,
        src_fields, tgt_field);

    struct yac_interp_weights * weights =
      yac_interp_weights_new(
        MPI_COMM_WORLD, YAC_LOC_CELL, (enum yac_location[]){YAC_LOC_CELL}, 1);

    size_t num_tgt_cells, num_tgt_vertices, num_tgt_edges;
    utest_get_basic_grid_data(
      filenames[0], &num_tgt_cells, &num_tgt_vertices, &num_tgt_edges);

    {
      // get all local target points
      size_t count;
      size_t * tgt_points;
      yac_interp_grid_get_tgt_points(interp_grid, &tgt_points, &count);

      // get global ids of local target points
      yac_int * tgt_global_ids = xmalloc(count * sizeof(*tgt_global_ids));
      yac_interp_grid_get_tgt_global_ids(
        interp_grid, tgt_points, count, tgt_global_ids);

      // get local ids of matching source points
      size_t * src_points = xmalloc(count * sizeof(*src_points));
      yac_interp_grid_src_global_to_local(
        interp_grid, 0, tgt_global_ids, count, src_points);
      free(tgt_global_ids);

      struct remote_points tgts = {
        .data =
          yac_interp_grid_get_tgt_remote_points(
            interp_grid, tgt_points, count),
        .count = count};

      struct remote_point * srcs =
        yac_interp_grid_get_src_remote_points(
          interp_grid, 0, src_points, count);

      yac_interp_weights_add_direct(weights, &tgts, srcs);

      free(tgts.data);
      free(tgt_points);
      free(srcs);
      free(src_points);
    }

    double *** src_data, *** src_frac_masks, ** tgt_data;

    if (is_src) {
      src_data = xmalloc(4 * sizeof(*src_data));
      for (size_t i = 0; i < 4; ++i)
        src_data[i] = xmalloc(1 * sizeof(**src_data));

      for (size_t collection_idx = 0; collection_idx < 4;
           ++collection_idx) {
        double * src_field =
          xmalloc(grid_data.num_cells * sizeof(*src_field));
        for (size_t i = 0; i < grid_data.num_cells; ++i)
          src_field[i] =
            (grid_data.core_cell_mask[i])?
              (double)(grid_data.cell_ids[i]):(-1.0);
        src_data[collection_idx][0] = src_field;
      }
      src_frac_masks = xmalloc(4 * sizeof(*src_frac_masks));
      for (size_t i = 0; i < 4; ++i)
        src_frac_masks[i] = xmalloc(1 * sizeof(**src_frac_masks));

      for (size_t collection_idx = 0; collection_idx < 4;
           ++collection_idx) {
        double * src_frac_mask =
          xmalloc(grid_data.num_cells * sizeof(*src_frac_mask));
        for (size_t i = 0; i < grid_data.num_cells; ++i)
          src_frac_mask[i] = (grid_data.cell_ids[i]&1)?1.0:0.0;
        src_frac_masks[collection_idx][0] = src_frac_mask;
      }
      tgt_data = NULL;
    } else {
      tgt_data = xmalloc(4 * sizeof(*tgt_data));
      for (size_t collection_idx = 0; collection_idx < 4;
           ++collection_idx)
        tgt_data[collection_idx] =
          xmalloc(grid_data.num_cells * sizeof(**tgt_data));
      src_data = NULL;
      src_frac_masks = NULL;
    }

    // checking with different reorder types and collection sizes
    for (size_t reorder_type = 0; reorder_type < reorder_type_count;
         ++reorder_type)
      for (size_t collection_size = 1; collection_size <= 4;
           collection_size *= 2)
          if (utest_check_direct_results(
                &grid_data, is_src, weights, reorder_types[reorder_type],
                collection_size, src_data, src_frac_masks, tgt_data))
            PUT_ERR("ERROR(yac_interp_weights_add_direct): "
                    "invalid interpolation result");

    // check writing and reading
    {
      char const * weight_file_name =
        "test_interp_weights_parallel_direct.nc";
      // write weights to file
      yac_interp_weights_write_to_file(
        weights, weight_file_name, grid_names[1], grid_names[0], 0, 0,
        YAC_WEIGHT_FILE_ERROR);

      // read weights
      struct interp_method * method_stack[2] = {
        yac_interp_method_file_new(
          weight_file_name, YAC_INTERP_FILE_ON_MISSING_FILE_DEFAULT,
          YAC_INTERP_FILE_ON_SUCCESS_DEFAULT), NULL};
      struct yac_interp_weights * weights_from_file =
        yac_interp_method_do_search(method_stack, interp_grid);
      yac_interp_method_delete(method_stack);

      // check results
      if (utest_check_direct_results(
            &grid_data, is_src, weights, YAC_MAPPING_ON_SRC, 1,
            src_data, src_frac_masks, tgt_data))
        PUT_ERR("ERROR(yac_interp_weights_add_direct): "
                "invalid interpolation result (read from file)");
      if (comm_rank == 0) unlink(weight_file_name);
      yac_interp_weights_delete(weights_from_file);
    }

    if (is_src) {
      for (size_t collection_idx = 0; collection_idx < 4; ++collection_idx) {
        free(src_data[collection_idx][0]);
        free(src_data[collection_idx]);
        free(src_frac_masks[collection_idx][0]);
        free(src_frac_masks[collection_idx]);
      }
      free(src_data);
      free(src_frac_masks);
    } else {
      for (size_t collection_idx = 0; collection_idx < 4; ++collection_idx)
        free(tgt_data[collection_idx]);
      free(tgt_data);
    }

    yac_interp_weights_delete(weights);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(grids[1]);
    yac_basic_grid_delete(grids[0]);
    free(filenames[1]);
    free(filenames[0]);
    MPI_Comm_free(&local_grid_comm);
  }

  { // check what happens when not all processes receive data
    int is_src = comm_rank < (comm_size / 2);

    MPI_Comm local_grid_comm;
    int local_grid_comm_size;
    MPI_Comm_split(MPI_COMM_WORLD, is_src, 0, &local_grid_comm);
    MPI_Comm_size(local_grid_comm, &local_grid_comm_size);

    set_even_io_rank_list(local_grid_comm);

    if (argc != 2) {
      PUT_ERR("ERROR: missing grid file directory");
      xt_finalize();
      MPI_Finalize();
      return TEST_EXIT_CODE;
    }

    char * filenames[2];
    char * grid_filenames[] =
      {"icon_grid_0030_R02B03_G.nc", "icon_grid_0043_R02B04_G.nc"};
    for (int i = 0; i < 2; ++i)
      filenames[i] =
        strcat(
          strcpy(
            malloc(strlen(argv[1]) + strlen(grid_filenames[i]) + 2), argv[1]),
          grid_filenames[i]);

    struct yac_basic_grid_data grid_data =
      yac_read_icon_basic_grid_data_parallel(
        filenames[is_src], local_grid_comm);
    struct yac_basic_grid * grids[2] =
      {yac_basic_grid_new(grid_names[is_src], grid_data),
       yac_basic_grid_empty_new(grid_names[is_src^1])};

    struct yac_dist_grid_pair * grid_pair =
      yac_dist_grid_pair_new(grids[0], grids[1], MPI_COMM_WORLD);

    struct yac_interp_field src_fields[] =
      {{.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX}};
    size_t num_src_fields = sizeof(src_fields) / sizeof(src_fields[0]);
    struct yac_interp_field tgt_field =
      {.location = YAC_LOC_CELL, .coordinates_idx = SIZE_MAX, .masks_idx = SIZE_MAX};

    struct yac_interp_grid * interp_grid =
      yac_interp_grid_new(
        grid_pair, grid_names[1], grid_names[0], num_src_fields,
        src_fields, tgt_field);

    struct yac_interp_weights * weights =
      yac_interp_weights_new(
        MPI_COMM_WORLD, YAC_LOC_CELL, (enum yac_location[]){YAC_LOC_CELL}, 1);

    size_t num_tgt_cells, num_tgt_vertices, num_tgt_edges;
    utest_get_basic_grid_data(
      filenames[0], &num_tgt_cells, &num_tgt_vertices, &num_tgt_edges);

    if (comm_rank == 0) {

      struct remote_points tgts = {
        .data =
          yac_interp_grid_get_tgt_remote_points(
            interp_grid, (size_t[]){0}, 1), .count = 1};

      struct remote_point * srcs =
        yac_interp_grid_get_src_remote_points(
          interp_grid, 0, (size_t[]){0}, 1);

      yac_interp_weights_add_direct(weights, &tgts, srcs);

      free(tgts.data);
      free(srcs);
    } else {
      struct remote_points tgts = {.data = NULL, .count = 0};

      struct remote_point * srcs =
        yac_interp_grid_get_src_remote_points(interp_grid, 0, NULL, 0);

      yac_interp_weights_add_direct(weights, &tgts, srcs);

      free(tgts.data);
      free(srcs);
    }

    for (size_t i = 1; i <= 16; i *= 2) {

      struct yac_interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, YAC_MAPPING_ON_SRC, i,
          YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0, NULL, 1, 1);

      yac_interpolation_delete(interpolation);
    }

    yac_interp_weights_delete(weights);
    yac_interp_grid_delete(interp_grid);
    yac_dist_grid_pair_delete(grid_pair);
    yac_basic_grid_delete(grids[1]);
    yac_basic_grid_delete(grids[0]);
    free(filenames[1]);
    free(filenames[0]);
    MPI_Comm_free(&local_grid_comm);
  }

  { // test yac_interp_weights_wcopy_weights

    struct yac_interp_weights * weights =
      yac_interp_weights_new(
        MPI_COMM_WORLD, YAC_LOC_CELL, (enum yac_location[]){YAC_LOC_CELL}, 1);

    // predefined stencils:
    // on rank 0:
    //   global_id=0 : fixed value 1.0
    //   global_id=1 : fixed value 2.0
    // on rank 1:
    //   global_id=2 : fixed value 3.0
    //   global_id=3 : direct from src 0
    // on rank 2:
    //   global_id=4 : direct from src 1
    //   global_id=5 : sum from src 2 and 3
    // on rank 3:
    //   global_id=6 : sum from src 1, 2, and 3
    //   global_id=7 : wsum from src 0.5*1 and 0.5*3
    // on rank 4:
    //   global_id=8 : wsum from src 0.2*3, 0.5*5, and 0.4*4
    switch (comm_rank) {
      default:
      case (0): {
        {
          struct remote_point remote_tgt_points[] =
            {{.global_id = 0,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}}};
          struct remote_points tgts = {
            .data = remote_tgt_points,
            .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
          };

          double fixed_value = 1.0;
          yac_interp_weights_add_fixed(weights, &tgts, fixed_value);
        }
        {
          struct remote_point remote_tgt_points[] =
            {{.global_id = 1,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 1}}}};
          struct remote_points tgts = {
            .data = remote_tgt_points,
            .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
          };

          double fixed_value = 2.0;
          yac_interp_weights_add_fixed(weights, &tgts, fixed_value);
        }
        break;
      }
      case (1): {
        {
          struct remote_point remote_tgt_points[] =
            {{.global_id = 2,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 2}}}};
          struct remote_points tgts = {
            .data = remote_tgt_points,
            .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
          };

          double fixed_value = 3.0;
          yac_interp_weights_add_fixed(weights, &tgts, fixed_value);
        }
        {
          struct remote_point remote_tgt_points[] =
            {{.global_id = 3,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 3}}}};
          struct remote_points tgts = {
            .data = remote_tgt_points,
            .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
          };
          struct remote_point srcs[] =
            {{.global_id = 0,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}}};
          yac_interp_weights_add_direct(weights, &tgts, srcs);
        }
        break;
      }
      case (2): {
        {
          struct remote_point remote_tgt_points[] =
            {{.global_id = 4,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 4}}}};
          struct remote_points tgts = {
            .data = remote_tgt_points,
            .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
          };
          struct remote_point srcs[] =
            {{.global_id = 1,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 1}}}};
          yac_interp_weights_add_direct(weights, &tgts, srcs);
        }
        {
          struct remote_point remote_tgt_points[] =
            {{.global_id = 5,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 5}}}};
          struct remote_points tgts = {
            .data = remote_tgt_points,
            .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
          };
          size_t num_src_per_tgt[] = {2};
          struct remote_point srcs[] =
            {{.global_id = 2,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 2}}},
             {.global_id = 3,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 3}}}};
          yac_interp_weights_add_sum(weights, &tgts, num_src_per_tgt, srcs);
        }
        break;
      }
      case (3): {
        {
          struct remote_point remote_tgt_points[] =
            {{.global_id = 6,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 6}}}};
          struct remote_points tgts = {
            .data = remote_tgt_points,
            .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
          };
          size_t num_src_per_tgt[] = {3};
          struct remote_point srcs[] =
            {{.global_id = 1,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 1}}},
             {.global_id = 2,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 2}}},
             {.global_id = 3,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 3}}}};
          yac_interp_weights_add_sum(weights, &tgts, num_src_per_tgt, srcs);
        }
        {
          struct remote_point remote_tgt_points[] =
            {{.global_id = 7,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 7}}}};
          struct remote_points tgts = {
            .data = remote_tgt_points,
            .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
          };
          size_t num_src_per_tgt[] = {2};
          struct remote_point srcs[] =
            {{.global_id = 1,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 1}}},
             {.global_id = 3,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 3}}}};
          double w[] = {0.5, 0.5};
          yac_interp_weights_add_wsum(weights, &tgts, num_src_per_tgt, srcs, w);
        }
        break;
      }
      case (4): {
        {
          struct remote_point remote_tgt_points[] =
            {{.global_id = 8,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 8}}}};
          struct remote_points tgts = {
            .data = remote_tgt_points,
            .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
          };
          size_t num_src_per_tgt[] = {3};
          struct remote_point srcs[] =
            {{.global_id = 3,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 3}}},
             {.global_id = 5,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 5}}},
             {.global_id = 4,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 4}}}};
          double w[] = {0.2, 0.5, 0.4};
          yac_interp_weights_add_wsum(weights, &tgts, num_src_per_tgt, srcs, w);
        }
        break;
      }
    };

    // copied stencils:
    // on rank 0:
    //   global_id 10: 1.0 * (fixed_value 1.0 from global_id 0)
    //   global_id 11: 4.0 * (fixed_value 2.0 from global_id 1) +
    //                 0.5 * (fixed_value 3.0 from global_id 2)
    // on rank 1:
    //   global_id 12: 0.3 * (direct from global_id 3)
    //   global_id 13: 0.5 * (direct from global_id 4) +
    //                 0.1 * (direct from global_id 3)
    // on rank 2:
    //   global_id 14: 0.2 * (sum from global_id 5)
    //   global_id 15: 0.1 * (direct from global_id 3) +
    //                 0.9 * (sum from global_id 6)
    //   global_id 16: 0.3 * (direct from global_id 4) +
    //                 0.3 * (sum from global_id 6) +
    //                 0.4 * (sum from global_id 5)
    // on rank 3:
    //   global_id 17: 0.1 * (wsum from global_id 7)
    //   global_id 18: 0.2 * (wsum from global_id 7) +
    //                 0.8 * (direct from global_id 3)
    //   global_id 19: 0.1 * (direct from global_id 4) +
    //                 0.2 * (sum from global_id 5) +
    //                 0.3 * (wsum from global_id 8)
    //   global_id 20: 1.0 * (direct from global_id 3) +
    //                 1.0 * (sum from global_id 5)
    //   global_id 21: 0.5 * (direct from global_id 3)
    //   global_id 22: 0.9 * (fixed_value 1.0 from global_id 0)
    switch (comm_rank) {
      case(0): {
        struct remote_point remote_tgt_points[] =
          {{.global_id = 10,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 10}}},
          {.global_id = 11,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 11}}}
          };
        struct remote_points tgts = {
          .data = remote_tgt_points,
          .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
        };

        size_t num_stencils_per_tgt[] = {1,2};
        size_t stencil_indices[] = {0,1,0};
        int stencil_ranks[] = {0,0,1};
        double w[] = {1.0,4.0,0.5};

        yac_interp_weights_wcopy_weights(
          weights, &tgts, num_stencils_per_tgt, stencil_indices,
          stencil_ranks, w);

        break;
      }
      case(1): {
        struct remote_point remote_tgt_points[] =
          {{.global_id = 12,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 12}}},
           {.global_id = 13,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 13}}}
          };
        struct remote_points tgts = {
          .data = remote_tgt_points,
          .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
        };

        size_t num_stencils_per_tgt[] = {1,2};
        size_t stencil_indices[] = {1,0,1};
        int stencil_ranks[] = {1,2,1};
        double w[] = {0.3,0.5,0.1};

        yac_interp_weights_wcopy_weights(
          weights, &tgts, num_stencils_per_tgt, stencil_indices,
          stencil_ranks, w);

        break;
      }
      case(2): {
        struct remote_point remote_tgt_points[] =
          {{.global_id = 14,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 14}}},
           {.global_id = 15,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 15}}},
           {.global_id = 16,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 16}}}
          };
        struct remote_points tgts = {
          .data = remote_tgt_points,
          .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
        };

        size_t num_stencils_per_tgt[] = {1,2,3};
        size_t stencil_indices[] = {1,1,0,0,0,1};
        int stencil_ranks[] = {2,1,3,2,3,2};
        double w[] = {0.2,0.1,0.9,0.3,0.3,0.4};

        yac_interp_weights_wcopy_weights(
          weights, &tgts, num_stencils_per_tgt, stencil_indices,
          stencil_ranks, w);

        break;
      }
      case(3): {
        struct remote_point remote_tgt_points[] =
          {{.global_id = 17,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 17}}},
           {.global_id = 18,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 18}}},
           {.global_id = 19,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 19}}},
           {.global_id = 20,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 20}}},
           {.global_id = 21,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 21}}},
           {.global_id = 22,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 22}}}};
        struct remote_points tgts = {
          .data = remote_tgt_points,
          .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
        };

        size_t num_stencils_per_tgt[] = {1,2,3,2,1,1};
        size_t stencil_indices[] = {1,1,1,0,1,0,1,1,1,0};
        int stencil_ranks[] = {3,3,1,2,2,4,1,2,1,0};
        double w[] = {0.1,0.2,0.8,0.1,0.2,0.3,1.0,1.0,0.5,0.9};

        yac_interp_weights_wcopy_weights(
          weights, &tgts, num_stencils_per_tgt, stencil_indices,
          stencil_ranks, w);

        break;
      }
      default: {

        yac_interp_weights_wcopy_weights(
          weights, NULL, NULL, NULL, NULL, NULL);
        break;
      }
    }

    struct yac_interpolation * interpolation =
      yac_interp_weights_get_interpolation(
        weights, YAC_MAPPING_ON_SRC, 1,
        YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0, NULL, 1, 1);

    if (comm_rank == 0) {

      double src_field_data[] = {1,2,3,4,5,6,7,8,9,10},
             tgt_field_data[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                                 -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                                 -1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
             ref_tgt_field_data[] =
              {1.0,2.0,3.0,1.0,2.0,3.0+4.0,2.0+3.0+4.0,0.5*2.0+0.5*4.0,
                 0.2*4.0+0.5*6.0+0.4*5.0,-1,
               1.0,9.5,0.3*(1.0),0.5*(2.0)+0.1*(1.0),0.2*(3.0+4.0),
                 0.1*(1.0)+0.9*(2.0+3.0+4.0),
                 0.3*(2.0)+0.3*(2.0+3.0+4.0)+0.4*(3.0+4.0),
                 0.1*(0.5*2.0 + 0.5*4.0),0.2*(0.5*2.0+0.5*4.0)+0.8*(1.0),
                 0.1*(2.0)+0.2*(3.0+4.0)+0.3*(0.2*4.0+0.5*6.0+0.4*5.0),
                 1.0*1.0+1.0*(3.0+4.0),0.5*1.0,0.9,-1,-1,-1,-1,-1,-1,-1};
      double * src_field = src_field_data,
             * tgt_field = tgt_field_data;
      double ** src_fields = &src_field;

      yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

      if (sizeof(tgt_field_data) != sizeof(ref_tgt_field_data))
        PUT_ERR("invalid reference data");

      for (size_t i = 0;
           i < sizeof(tgt_field_data)/sizeof(tgt_field_data[0]); ++i)
        if (fabs(tgt_field_data[i] - ref_tgt_field_data[i]) > TOL)
          PUT_ERR("ERROR in yac_interp_weights_wcopy_weights");

    } else {

      yac_interpolation_execute(interpolation, NULL, NULL);
    }

    yac_interpolation_delete(interpolation);

    yac_interp_weights_delete(weights);
  }

  { // test yac_interp_weights_wcopy_weights, sources are distributed among
    // processes and partially available on multiple ones

    struct yac_interp_weights * weights =
      yac_interp_weights_new(
        MPI_COMM_WORLD, YAC_LOC_CELL, (enum yac_location[]){YAC_LOC_CELL}, 1);

    // predefined stencils:
    // on rank 0:
    //   global_id=0 : sum from src 0, 1, and 2
    //   global_id=1 : sum from src 1, 3, and 4
    // on rank 1:
    //   global_id=2 : sum from src 2, 3, and 4
    //   global_id=3 : sum from src 0, 4, and 6
    // on rank 2:
    //   global_id=4 : wsum from src 0.3*0, 0.4*1, and 0.5*2
    //   global_id=5 : wsum from src 0.1*1, 0.1*3, and 0.1*4
    // on rank 3:
    //   global_id=6 : wsum from src 1.0*2, 1.0*3, and 1.0*4
    //   global_id=7 : wsum from src 1.0*0, 2.0*4, and 3.0*6
    switch (comm_rank) {
      case (0): {
        struct remote_point remote_tgt_points[] =
          {{.global_id = 0,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}},
           {.global_id = 1,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 1}}}};
        struct remote_points tgts = {
          .data = remote_tgt_points,
          .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
        };
        size_t num_src_per_tgt[] = {3,3};
        struct remote_point srcs[] =
          {{.global_id = 0,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}},
           {.global_id = 1,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 1}}},
           {.global_id = 2,
            .data =
              {.count = 2,
               .data.multi =
                 (struct remote_point_info[])
                  {{.rank = 0, .orig_pos = 2},{.rank = 1, .orig_pos = 0}}}},
           {.global_id = 1,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 1}}},
           {.global_id = 3,
            .data = {.count = 1, .data.single = {.rank = 1, .orig_pos = 1}}},
           {.global_id = 4,
            .data =
              {.count = 2,
               .data.multi =
                 (struct remote_point_info[])
                  {{.rank = 1, .orig_pos = 2},{.rank = 2, .orig_pos = 0}}}}};
        yac_interp_weights_add_sum(weights, &tgts, num_src_per_tgt, srcs);
        break;
      }
      case (1): {
        struct remote_point remote_tgt_points[] =
          {{.global_id = 2,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 2}}},
           {.global_id = 3,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 3}}}};
        struct remote_points tgts = {
          .data = remote_tgt_points,
          .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
        };
        size_t num_src_per_tgt[] = {3,3};
        struct remote_point srcs[] =
          {{.global_id = 2,
            .data =
              {.count = 2,
               .data.multi =
                 (struct remote_point_info[])
                  {{.rank = 0, .orig_pos = 2},{.rank = 1, .orig_pos = 0}}}},
           {.global_id = 3,
            .data = {.count = 1, .data.single = {.rank = 1, .orig_pos = 1}}},
           {.global_id = 4,
            .data =
              {.count = 2,
               .data.multi =
                 (struct remote_point_info[])
                  {{.rank = 1, .orig_pos = 2},{.rank = 2, .orig_pos = 0}}}},
           {.global_id = 0,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}},
           {.global_id = 4,
            .data =
              {.count = 2,
               .data.multi =
                 (struct remote_point_info[])
                  {{.rank = 1, .orig_pos = 2},{.rank = 2, .orig_pos = 0}}}},
           {.global_id = 6,
            .data = {.count = 1, .data.single = {.rank = 2, .orig_pos = 2}}}};
        yac_interp_weights_add_sum(weights, &tgts, num_src_per_tgt, srcs);
        break;
      }
      case (2): {
        struct remote_point remote_tgt_points[] =
          {{.global_id = 4,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 4}}},
           {.global_id = 5,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 5}}}};
        struct remote_points tgts = {
          .data = remote_tgt_points,
          .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
        };
        size_t num_src_per_tgt[] = {3,3};
        struct remote_point srcs[] =
          {{.global_id = 0,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}},
           {.global_id = 1,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 1}}},
           {.global_id = 2,
            .data =
              {.count = 2,
               .data.multi =
                 (struct remote_point_info[])
                  {{.rank = 0, .orig_pos = 2},{.rank = 1, .orig_pos = 0}}}},
           {.global_id = 1,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 1}}},
           {.global_id = 3,
            .data = {.count = 1, .data.single = {.rank = 1, .orig_pos = 1}}},
           {.global_id = 4,
            .data =
              {.count = 2,
               .data.multi =
                 (struct remote_point_info[])
                  {{.rank = 1, .orig_pos = 2},{.rank = 2, .orig_pos = 0}}}}};
          double w[] = {0.3,0.4,0.5, 0.1,0.1,0.1};
          yac_interp_weights_add_wsum(weights, &tgts, num_src_per_tgt, srcs, w);
        break;
      }
      case (3): {
        struct remote_point remote_tgt_points[] =
          {{.global_id = 6,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 6}}},
           {.global_id = 7,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 7}}}};
        struct remote_points tgts = {
          .data = remote_tgt_points,
          .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
        };
        size_t num_src_per_tgt[] = {3,3};
        struct remote_point srcs[] =
          {{.global_id = 2,
            .data =
              {.count = 2,
               .data.multi =
                 (struct remote_point_info[])
                  {{.rank = 0, .orig_pos = 2},{.rank = 1, .orig_pos = 0}}}},
           {.global_id = 3,
            .data = {.count = 1, .data.single = {.rank = 1, .orig_pos = 1}}},
           {.global_id = 4,
            .data =
              {.count = 2,
               .data.multi =
                 (struct remote_point_info[])
                  {{.rank = 1, .orig_pos = 2},{.rank = 2, .orig_pos = 0}}}},
           {.global_id = 0,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}},
           {.global_id = 4,
            .data =
              {.count = 2,
               .data.multi =
                 (struct remote_point_info[])
                  {{.rank = 1, .orig_pos = 2},{.rank = 2, .orig_pos = 0}}}},
           {.global_id = 6,
            .data = {.count = 1, .data.single = {.rank = 2, .orig_pos = 2}}}};
          double w[] = {1.0,1.0,1.0, 1.0,2.0,3.0};
          yac_interp_weights_add_wsum(weights, &tgts, num_src_per_tgt, srcs, w);
        break;
      }
      default: break;
    };

    // copied stencils:
    // on rank 0:
    //   global_id 9: 1.0 * (sum from global_id 0) +
    //                1.0 * (sum from global_id 2)
    //   global_id 9: 0.1 * (sum from global_id 0) +
    //                0.2 * (sum from global_id 1) +
    //                0.3 * (sum from global_id 2) +
    //                0.4 * (sum from global_id 3)
    //   global_id 10: 0.2 * (sum from global_id 2) +
    //                 0.2 * (sum from global_id 3) +
    //                 0.2 * (wsum from global_id 4) +
    //                 0.2 * (wsum from global_id 5)
    //   global_id 11: 0.25 * (wsum from global_id 5) +
    //                 0.25 * (wsum from global_id 6) +
    //                 0.25 * (wsum from global_id 7) +
    //                 0.25 * (wsum from global_id 8)
    switch (comm_rank) {
      case(0): {
        struct remote_point remote_tgt_points[] =
          {{.global_id = 8,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 8}}},
          {.global_id = 9,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 9}}},
          {.global_id = 10,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 10}}},
          {.global_id = 11,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 11}}}};
        struct remote_points tgts = {
          .data = remote_tgt_points,
          .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
        };

        size_t num_stencils_per_tgt[] = {2,4,4,4};
        size_t stencil_indices[] = {0,0, 0,1,0,1, 0,1,0,1, 0,1,0,1};
        int stencil_ranks[] = {0,1, 0,0,1,1, 1,1,2,2, 2,2,3,3};
        double w[] =
          {1.0,1.0, 0.1,0.2,0.3,0.4, 0.2,0.2,0.2,0.2, 0.25,0.25,0.25,0.25};

        yac_interp_weights_wcopy_weights(
          weights, &tgts, num_stencils_per_tgt, stencil_indices,
          stencil_ranks, w);

        break;
      }
      default: {

        yac_interp_weights_wcopy_weights(
          weights, NULL, NULL, NULL, NULL, NULL);
        break;
      }
    }
    for (size_t i = 0; i < reorder_type_count; ++i) {

      struct yac_interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_types[i], 1,
          YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0, NULL, 1, 1);

      double src_field_data[3][3] = {{1,2,3},{3,4,5},{5,6,7}};
      double * src_field[5] =
        {src_field_data[0],src_field_data[1],src_field_data[2],NULL,NULL};
      double ** src_fields = &(src_field[comm_rank]);
      double tgt_field_data[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
      double * tgt_field = tgt_field_data;

      yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

      if (comm_rank == 0) {

        double ref_tgt_field_data[] =
          {1.0+2.0+3.0,
           2.0+4.0+5.0,
           3.0+4.0+5.0,
           1.0+5.0+7.0,
           0.3*1.0+0.4*2.0+0.5*3.0,
           0.1*2.0+0.1*4.0+0.1*5.0,
           1.0*3.0+1.0*4.0+1.0*5.0,
           1.0*1.0+2.0*5.0+3.0*7.0,
           1.0*(1.0+2.0+3.0)+1.0*(3.0+4.0+5.0),
           0.1*(1.0+2.0+3.0)+0.2*(2.0+4.0+5.0)+
             0.3*(3.0+4.0+5.)+0.4*(1.0+5.0+7.0),
           0.2*(3.0+4.0+5.0)+0.2*(1.0+5.0+7.0)+
             0.2*(0.3*1.0+0.4*2.0+0.5*3.0)+0.2*(0.1*2.0+0.1*4.0+0.1*5.0),
           0.25*(0.3*1.0+0.4*2.0+0.5*3.0)+0.25*(0.1*2.0+0.1*4.0+0.1*5.0)+
             0.25*(1.0*3.0+1.0*4.0+1.0*5.0)+0.25*(1.0*1.0+2.0*5.0+3.0*7.0)};

        if (sizeof(tgt_field_data) != sizeof(ref_tgt_field_data))
          PUT_ERR("invalid reference data");

        for (size_t i = 0;
            i < sizeof(tgt_field_data)/sizeof(tgt_field_data[0]); ++i)
          if (fabs(tgt_field_data[i] - ref_tgt_field_data[i]) > TOL)
            PUT_ERR("ERROR in yac_interp_weights_wcopy_weights");
      }

      yac_interpolation_delete(interpolation);
    }

    yac_interp_weights_delete(weights);
  }

  { // test yac_interp_weights_wcopy_weights, with multiple stencils
    // containing the same sources

    struct yac_interp_weights * weights =
      yac_interp_weights_new(
        MPI_COMM_WORLD, YAC_LOC_CELL, (enum yac_location[]){YAC_LOC_CELL}, 1);

    if (comm_rank == 0) {

      // predefined stencil
      {
        struct remote_point remote_tgt_points[] =
          {{.global_id = 0,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}}};
        struct remote_points tgts = {
          .data = remote_tgt_points,
          .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
        };
        size_t num_src_per_tgt[] = {1};
        struct remote_point srcs[] =
          {{.global_id = 0,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}}};
        double w[] = {1.0};
        yac_interp_weights_add_wsum(weights, &tgts, num_src_per_tgt, srcs, w);
      }

      // copy stencil
      {
        struct remote_point remote_tgt_points[] =
          {{.global_id = 1,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 1}}}};
        struct remote_points tgts = {
          .data = remote_tgt_points,
          .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
        };

#define N (20)
        size_t num_stencils_per_tgt[] = {1 + 9 * N};
        size_t stencil_indices[1 + 9 * N];
        int stencil_ranks[1 + 9 * N];
        double w[1 + 9 * N];

        for (size_t i = 0; i < (1 + 9 * N); ++i) {
          stencil_indices[i] = 0;
          stencil_ranks[i] = 0;
        }
        for (int i = 0, k = 0; i < N; ++i)
          for (int j = 0; j < 9; ++j, ++k)
            w[k] = pow(10.0, (double)(-i));
        w[9 * N] = pow(10.0, (double)(-N+1));

        // copy stencil
        yac_interp_weights_wcopy_weights(
          weights, &tgts, num_stencils_per_tgt, stencil_indices,
          stencil_ranks, w);
      }
    } else {
      yac_interp_weights_wcopy_weights(
          weights, NULL, NULL, NULL, NULL, NULL);
    }
    for (size_t i = 0; i < reorder_type_count; ++i) {

      struct yac_interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_types[i], 1,
          YAC_FRAC_MASK_NO_VALUE, 1.0, 0.0, NULL, 1, 1);

      double src_field_data[1][1] = {{1}};
      double * src_field[5] = {src_field_data[0],NULL,NULL,NULL,NULL};
      double ** src_fields = &(src_field[comm_rank]);
      double tgt_field_data[] = {-1,-1};
      double * tgt_field = tgt_field_data;

      yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

      if (comm_rank == 0) {

        double ref_tgt_field_data[] = {1.0, 10.0};

        if (sizeof(tgt_field_data) != sizeof(ref_tgt_field_data))
          PUT_ERR("invalid reference data");

        for (size_t i = 0;
            i < sizeof(tgt_field_data)/sizeof(tgt_field_data[0]); ++i)
          if (fabs(tgt_field_data[i] - ref_tgt_field_data[i]) > 0.0)
            PUT_ERR("ERROR in yac_interp_weights_wcopy_weights");
      }

      yac_interpolation_delete(interpolation);
    }

    yac_interp_weights_delete(weights);
  }

  { // test weights using multiple source fields

    struct yac_interp_weights * weights =
      yac_interp_weights_new(
        MPI_COMM_WORLD, YAC_LOC_CELL,
        (enum yac_location[]){YAC_LOC_CELL, YAC_LOC_CORNER, YAC_LOC_EDGE}, 3);

    switch (comm_rank) {
      case (0): {
        // tgt global_id: 0
        //   src_field 0:
        //     0.3 * src global_id 0
        //  src field 1:
        //     0.3 * src global_id 0
        //  src field 2:
        //     0.4 * src global_id 0
        // tgt global_id: 2
        //   src_field 0:
        //     0.1 * src global_id 1
        //     0.1 * src global_id 2
        //  src field 1:
        //     0.2 * src global_id 1
        //  src field 2:
        //     0.3 * src global_id 3
        //     0.3 * src global_id 2
        //     0.3 * src global_id 1
        // tgt global_id: 3
        //   src_field 0:
        //     1.0 * src global_id 5
        //     1.0 * src global_id 2
        //  src field 1:
        //  src field 2:
        //     1.0 * src global_id 4
        struct remote_point remote_tgt_points[] =
          {{.global_id = 0,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}},
           {.global_id = 2,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 2}}},
           {.global_id = 3,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 3}}}
          };
        struct remote_points tgts = {
          .data = remote_tgt_points,
          .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
        };
        size_t num_src_per_field_per_tgt[3][3] = {{1,1,1},{2,1,3},{2,0,1}};
        struct remote_point srcs_per_field_data[3][5] =
          {{{.global_id = 0,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}},
            {.global_id = 1,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 1}}},
            {.global_id = 2,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 2}}},
            {.global_id = 5,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 5}}},
            {.global_id = 2,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 2}}}},
           {{.global_id = 0,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}},
            {.global_id = 1,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 1}}}},
           {{.global_id = 0,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}},
            {.global_id = 3,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 3}}},
            {.global_id = 2,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 2}}},
            {.global_id = 1,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 1}}},
            {.global_id = 4,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 4}}}}};
        struct remote_point * srcs_per_field[] =
        {srcs_per_field_data[0],
         srcs_per_field_data[1],
         srcs_per_field_data[2]};
        double w[] = {0.3,0.3,1.0, 0.1,0.1,0.2,0.3,0.3,0.3, 1.0,1.0,1.0};
        yac_interp_weights_add_wsum_mf(
          weights, &tgts, &num_src_per_field_per_tgt[0][0],
          srcs_per_field, w, 3);
        break;
      }
      case (1): {
        // tgt global_id: 4
        //   src_field 0:
        //     src global_id 0
        //  src field 1:
        //     src global_id 0
        //  src field 2:
        //     src global_id 0
        // tgt global_id: 5
        //   src_field 0:
        //     src global_id 1
        //     src global_id 2
        //  src field 1:
        //     src global_id 1
        //  src field 2:
        //     src global_id 3
        //     src global_id 2
        //     src global_id 1
        // tgt global_id: 6
        //   src_field 0:
        //     src global_id 5
        //     src global_id 2
        //  src field 1:
        //  src field 2:
        //     src global_id 4
        struct remote_point remote_tgt_points[] =
          {{.global_id = 4,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 4}}},
           {.global_id = 5,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 5}}},
           {.global_id = 6,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 6}}}
          };
        struct remote_points tgts = {
          .data = remote_tgt_points,
          .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
        };
        size_t num_src_per_field_per_tgt[3][3] = {{1,1,1},{2,1,3},{2,0,1}};
        struct remote_point srcs_per_field_data[3][5] =
          {{{.global_id = 0,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}},
            {.global_id = 1,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 1}}},
            {.global_id = 2,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 2}}},
            {.global_id = 5,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 5}}},
            {.global_id = 2,
              .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 2}}}},
           {{.global_id = 0,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}},
            {.global_id = 1,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 1}}}},
           {{.global_id = 0,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}},
            {.global_id = 3,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 3}}},
            {.global_id = 2,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 2}}},
            {.global_id = 1,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 1}}},
            {.global_id = 4,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 4}}}}};
        struct remote_point * srcs_per_field[] =
        {srcs_per_field_data[0],
         srcs_per_field_data[1],
         srcs_per_field_data[2]};
        yac_interp_weights_add_sum_mf(
          weights, &tgts, &num_src_per_field_per_tgt[0][0], srcs_per_field, 3);
        break;
      }
      case (2): {
        // tgt global_id: 7
        //   src_field 0:
        //     src global_id 0
        // tgt global_id: 8
        //   src_field 2:
        //     src global_id 3
        // tgt global_id: 9
        //   src_field 1:
        //     src global_id 2
        // tgt global_id: 10
        //   src_field 2:
        //     src global_id 2
        // tgt global_id: 11
        //   src_field 2:
        //     src global_id 3
        struct remote_point remote_tgt_points[] =
          {{.global_id = 7,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 7}}},
           {.global_id = 8,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 8}}},
           {.global_id = 9,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 9}}},
           {.global_id = 10,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 10}}},
           {.global_id = 11,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 11}}}};
        struct remote_points tgts = {
          .data = remote_tgt_points,
          .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
        };
        size_t src_field_indices[] = {0,2,1,2,2};
        struct remote_point srcs_per_field_data[3][3] =
          {{{.global_id = 0,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}}},
           {{.global_id = 2,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 2}}}},
           {{.global_id = 3,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 3}}},
            {.global_id = 2,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 2}}},
            {.global_id = 3,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 3}}}}};
        struct remote_point * srcs_per_field[] =
        {srcs_per_field_data[0],
         srcs_per_field_data[1],
         srcs_per_field_data[2]};
        yac_interp_weights_add_direct_mf(
          weights, &tgts, src_field_indices, srcs_per_field, 3);
        break;
      }
      default: break;
    };


    for (size_t i = 0; i < reorder_type_count; ++i) {

      size_t collection_size = 1;
      double frac_mask_fallback_value = YAC_FRAC_MASK_NO_VALUE;
      double scaling_factor = 1.0;
      double scaling_summand = 0.0;
      char const * yaxt_exchanger_name = NULL;
      struct yac_interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, reorder_types[i], collection_size,
          frac_mask_fallback_value, scaling_factor, scaling_summand,
          yaxt_exchanger_name, 1, 1);
      struct yac_interpolation_exchange * interpolation_exchange;
      struct yac_interp_weights_data interp_weights_data;
      yac_interp_weights_get_interpolation_raw(
        weights, collection_size, frac_mask_fallback_value,
        scaling_factor, scaling_summand, yaxt_exchanger_name,
        &interpolation_exchange, &interp_weights_data, 1, 1);

      if (comm_rank == 0) {

        double src_field_data[3][10] = {{1,2,3,4,5,6,7,8,9,10},
                                        {10,20,30,40,50,60,70,80,90,100},
                                        {100,200,300,400,500,600,700,800,900,1000}};
        double tgt_field_data[] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
                                   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1};

        double ref_tgt_field_data[] =
          {0.3*1.0+0.3*10.0+1.0*100.0,
           -1,
           0.1*2.0+0.1*3.0+0.2*20.0+0.3*400.0+0.3*300.0+0.3*200.0,
           1.0*6.0+1.0*3.0+1.0*500.0,
           1.0+10.0+100.0,
           2.0+3.0+20.0+400.0+300.0+200.0,
           6.0+3.0+500.0,
           1.0,
           400.0,
           30.0,
           300.0,
           400.0,
           -1,-1,-1,-1,-1,-1,-1,-1};
        double * src_field[3] =
          {src_field_data[0], src_field_data[1], src_field_data[2]};
        double * tgt_field = tgt_field_data;
        double ** src_fields = src_field;

        yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

        if (sizeof(tgt_field_data) != sizeof(ref_tgt_field_data))
          PUT_ERR("invalid reference data");

        for (size_t i = 0;
             i < sizeof(tgt_field_data)/sizeof(tgt_field_data[0]); ++i)
          if (fabs(tgt_field_data[i] - ref_tgt_field_data[i]) > TOL)
            PUT_ERR("ERROR in multi source field support");

        for (size_t i = 0;
             i < sizeof(tgt_field_data) / sizeof(tgt_field_data[0]); ++i)
          tgt_field_data[i] = -1;

        utest_interpolation_execute_raw(
          interpolation_exchange, interp_weights_data, collection_size,
          &src_fields, &tgt_field);

        if (sizeof(tgt_field_data) != sizeof(ref_tgt_field_data))
          PUT_ERR("invalid reference data");

        for (size_t i = 0;
             i < sizeof(tgt_field_data)/sizeof(tgt_field_data[0]); ++i)
          if (fabs(tgt_field_data[i] - ref_tgt_field_data[i]) > TOL)
            PUT_ERR("ERROR in multi source field support raw");

      } else {

        yac_interpolation_execute(interpolation, NULL, NULL);

        utest_interpolation_execute_raw(
          interpolation_exchange, interp_weights_data, collection_size,
          NULL, NULL);
      }

      yac_interp_weights_data_free(interp_weights_data);
      yac_interpolation_exchange_delete(
        interpolation_exchange, "direct_mf");
      yac_interpolation_delete(interpolation);
    }

    yac_interp_weights_delete(weights);
  }

  { // test interfaces for multiple source fields but only using
    // a single source field

    struct yac_interp_weights * weights =
      yac_interp_weights_new(
        MPI_COMM_WORLD, YAC_LOC_CELL, (enum yac_location[]){YAC_LOC_CELL}, 1);

    switch (comm_rank) {
      case (0): {
        // tgt global_id: 2
        //   src_field 0:
        //     src global_id 2
        struct remote_point remote_tgt_points[] =
          {{.global_id = 2,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 2}}}};
        struct remote_points tgts = {
          .data = remote_tgt_points,
          .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
        };
        size_t src_field_indices[] = {0,0,0};
        struct remote_point srcs_per_field_data[1][1] =
          {{{.global_id = 2,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 2}}}}};
        struct remote_point * srcs_per_field[] = {srcs_per_field_data[0]};
        yac_interp_weights_add_direct_mf(
          weights, &tgts, src_field_indices, srcs_per_field, 1);
        break;
      }
      case (1): {
        // tgt global_id: 4
        //   src_field 0:
        //     src global_id 0
        //     src global_id 2
        struct remote_point remote_tgt_points[] =
          {{.global_id = 4,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 4}}}};
        struct remote_points tgts = {
          .data = remote_tgt_points,
          .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
        };
        size_t num_src_per_field_per_tgt[1][1] = {{2}};
        struct remote_point srcs_per_field_data[1][2] =
          {{{.global_id = 0,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}},
            {.global_id = 2,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 2}}}}};
        struct remote_point * srcs_per_field[] = {srcs_per_field_data[0]};
        yac_interp_weights_add_sum_mf(
          weights, &tgts, &num_src_per_field_per_tgt[0][0], srcs_per_field, 1);
        break;
      }
      case (2): {
        // tgt global_id: 6
        //   src_field 0:
        //     0.7 * src global_id 5
        //     0.3 * src global_id 2
        struct remote_point remote_tgt_points[] =
          {{.global_id = 6,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 6}}}};
        struct remote_points tgts = {
          .data = remote_tgt_points,
          .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
        };
        size_t num_src_per_field_per_tgt[1][1] = {{2}};
        struct remote_point srcs_per_field_data[1][2] =
          {{{.global_id = 5,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 5}}},
            {.global_id = 2,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 2}}}}};
        struct remote_point * srcs_per_field[] = {srcs_per_field_data[0]};
        double w[] = {0.7,0.3};
        yac_interp_weights_add_wsum_mf(
          weights, &tgts, &num_src_per_field_per_tgt[0][0],
          srcs_per_field, w, 1);
        break;
      }
      default: break;
    };

    size_t collection_size = 1;
    double frac_mask_fallback_value = YAC_FRAC_MASK_NO_VALUE;
    double scaling_factor = 1.0;
    double scaling_summand = 0.0;
    struct yac_interpolation * interpolation =
      yac_interp_weights_get_interpolation(
        weights, YAC_MAPPING_ON_SRC, collection_size,
        frac_mask_fallback_value, scaling_factor, scaling_summand,
        NULL, 1, 1);
    struct yac_interpolation_exchange * interpolation_exchange;
    struct yac_interp_weights_data interp_weights_data;
    yac_interp_weights_get_interpolation_raw(
      weights, collection_size, frac_mask_fallback_value,
      scaling_factor, scaling_summand, NULL,
      &interpolation_exchange, &interp_weights_data, 1, 1);

    if (comm_rank == 0) {

      double src_field_data[1][10] = {{1,2,3,4,5,6,7,8,9,10}};
      double tgt_field_data[] = {-1,-1,-1,-1,-1,-1,-1,-1};

      double ref_tgt_field_data[] =
        {-1,-1,3.0,-1,1.0+3.0,-1,0.7*6.0+0.3*3.0,-1};
      double * src_field[1] = {src_field_data[0]};
      double * tgt_field = tgt_field_data;
      double ** src_fields = src_field;

      yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

      if (sizeof(tgt_field_data) != sizeof(ref_tgt_field_data))
        PUT_ERR("invalid reference data");

      for (size_t i = 0;
            i < sizeof(tgt_field_data)/sizeof(tgt_field_data[0]); ++i)
        if (fabs(tgt_field_data[i] - ref_tgt_field_data[i]) > TOL)
          PUT_ERR("ERROR in multi source field support");

      for (size_t i = 0;
            i < sizeof(tgt_field_data)/sizeof(tgt_field_data[0]); ++i)
        tgt_field_data[i] = -1;
      utest_interpolation_execute_raw(
        interpolation_exchange, interp_weights_data, collection_size,
        &src_fields, &tgt_field);

      for (size_t i = 0;
            i < sizeof(tgt_field_data)/sizeof(tgt_field_data[0]); ++i)
        if (fabs(tgt_field_data[i] - ref_tgt_field_data[i]) > TOL)
          PUT_ERR("ERROR in multi source field support");

    } else {

      yac_interpolation_execute(interpolation, NULL, NULL);

      utest_interpolation_execute_raw(
        interpolation_exchange, interp_weights_data, collection_size,
        NULL, NULL);
    }

    yac_interp_weights_data_free(interp_weights_data);
    yac_interpolation_exchange_delete(
      interpolation_exchange, "interpolation raw");
    yac_interpolation_delete(interpolation);

    yac_interp_weights_delete(weights);
  }

  { // test interfaces for multiple source fields,
    // the target points are available on all processes

    struct yac_interp_weights * weights =
      yac_interp_weights_new(
        MPI_COMM_WORLD, YAC_LOC_CELL,
        (enum yac_location[]){YAC_LOC_CELL, YAC_LOC_CORNER}, 2);

    switch (comm_rank) {
      case (0): {
        // tgt global_id: 2
        //   src_field 1:
        //     src global_id 2
        struct remote_point remote_tgt_points[] =
          {{.global_id = 2,
            .data =
              {.count = 5,
               .data.multi =
                 (struct remote_point_info[]){{.rank = 0, .orig_pos = 2},
                                              {.rank = 1, .orig_pos = 2},
                                              {.rank = 2, .orig_pos = 2},
                                              {.rank = 3, .orig_pos = 2},
                                              {.rank = 4, .orig_pos = 2}}}}};
        struct remote_points tgts = {
          .data = remote_tgt_points,
          .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
        };
        size_t src_field_indices[] = {1};
        struct remote_point srcs_per_field_data[2][1] =
          {{{.global_id = -1}},
           {{.global_id = 2,
             .data = {.count = 1, .data.single = {.rank = 1, .orig_pos = 2}}}}};
        struct remote_point * srcs_per_field[] =
          {srcs_per_field_data[0], srcs_per_field_data[1]};
        yac_interp_weights_add_direct_mf(
          weights, &tgts, src_field_indices, srcs_per_field, 2);
        break;
      }
      case (1): {
        // tgt global_id: 4
        //   src_field 0:
        //     src global_id 0
        //   src_field 1:
        //     src global_id 2
        struct remote_point remote_tgt_points[] =
          {{.global_id = 4,
            .data =
              {.count = 5,
               .data.multi =
                 (struct remote_point_info[]){{.rank = 0, .orig_pos = 4},
                                              {.rank = 1, .orig_pos = 4},
                                              {.rank = 2, .orig_pos = 4},
                                              {.rank = 3, .orig_pos = 4},
                                              {.rank = 4, .orig_pos = 4}}}}};
        struct remote_points tgts = {
          .data = remote_tgt_points,
          .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
        };
        size_t num_src_per_field_per_tgt[2][1] = {{1},{1}};
        struct remote_point srcs_per_field_data[2][1] =
          {{{.global_id = 0,
             .data = {.count = 1, .data.single = {.rank = 1, .orig_pos = 0}}}},
           {{.global_id = 2,
             .data = {.count = 1, .data.single = {.rank = 1, .orig_pos = 2}}}}};
        struct remote_point * srcs_per_field[] =
          {srcs_per_field_data[0], srcs_per_field_data[1]};
        yac_interp_weights_add_sum_mf(
          weights, &tgts, &num_src_per_field_per_tgt[0][0], srcs_per_field, 2);
        break;
      }
      case (2): {
        // tgt global_id: 6
        //   src_field 0:
        //     0.7 * src global_id 5
        //   src_field 1:
        //     0.3 * src global_id 2
        struct remote_point remote_tgt_points[] =
          {{.global_id = 6,
            .data =
              {.count = 5,
               .data.multi =
                 (struct remote_point_info[]){{.rank = 0, .orig_pos = 6},
                                              {.rank = 1, .orig_pos = 6},
                                              {.rank = 2, .orig_pos = 6},
                                              {.rank = 3, .orig_pos = 6},
                                              {.rank = 4, .orig_pos = 6}}}}};
        struct remote_points tgts = {
          .data = remote_tgt_points,
          .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
        };
        size_t num_src_per_field_per_tgt[2][1] = {{1},{1}};
        struct remote_point srcs_per_field_data[2][1] =
          {{{.global_id = 5,
             .data = {.count = 1, .data.single = {.rank = 1, .orig_pos = 5}}}},
           {{.global_id = 2,
             .data = {.count = 1, .data.single = {.rank = 1, .orig_pos = 2}}}}};
        struct remote_point * srcs_per_field[] =
          {srcs_per_field_data[0],srcs_per_field_data[1]};
        double w[] = {0.7,0.3};
        yac_interp_weights_add_wsum_mf(
          weights, &tgts, &num_src_per_field_per_tgt[0][0],
          srcs_per_field, w, 2);
        break;
      }
      default: break;
    };

    size_t collection_size = 1;
    double frac_mask_fallback_value = YAC_FRAC_MASK_NO_VALUE;
    double scaling_factor = 1.0;
    double scaling_summand = 0.0;
    struct yac_interpolation * interpolation =
      yac_interp_weights_get_interpolation(
        weights, YAC_MAPPING_ON_SRC, collection_size,
        frac_mask_fallback_value, scaling_factor, scaling_summand, NULL, 1, 1);
    struct yac_interpolation_exchange * interpolation_exchange;
    struct yac_interp_weights_data interp_weights_data;
    yac_interp_weights_get_interpolation_raw(
      weights, collection_size, frac_mask_fallback_value,
      scaling_factor, scaling_summand, NULL,
      &interpolation_exchange, &interp_weights_data, 1, 1);

    double src_field_data[2][10] =
      {{1,2,3,4,5,6,7,8,9,10},{10,20,30,40,50,60,70,80,90,100}};
    double tgt_field_data[] = {-1,-1,-1,-1,-1,-1,-1,-1};

    double ref_tgt_field_data[] =
      {-1,-1,30.0,-1,1.0+30.0,-1,0.7*6.0+0.3*30.0,-1};
    double * src_field[2] = {src_field_data[0],src_field_data[1]};
    double * tgt_field = tgt_field_data;
    double ** src_fields = src_field;

    yac_interpolation_execute(
      interpolation, (comm_rank == 1)?(&src_fields):NULL, &tgt_field);

    if (sizeof(tgt_field_data) != sizeof(ref_tgt_field_data))
      PUT_ERR("invalid reference data");

    for (size_t i = 0;
          i < sizeof(tgt_field_data)/sizeof(tgt_field_data[0]); ++i)
      if (fabs(tgt_field_data[i] - ref_tgt_field_data[i]) > TOL)
        PUT_ERR("ERROR in multi source field support");

    for (size_t i = 0;
          i < sizeof(tgt_field_data)/sizeof(tgt_field_data[0]); ++i)
      tgt_field_data[i] = -1;
    utest_interpolation_execute_raw(
      interpolation_exchange, interp_weights_data, collection_size,
      (comm_rank == 1)?(&src_fields):NULL, &tgt_field);

    for (size_t i = 0;
          i < sizeof(tgt_field_data)/sizeof(tgt_field_data[0]); ++i)
      if (fabs(tgt_field_data[i] - ref_tgt_field_data[i]) > TOL)
        PUT_ERR("ERROR in multi source field support");

    yac_interp_weights_data_free(interp_weights_data);
    yac_interpolation_exchange_delete(
      interpolation_exchange, "interpolation raw");
    yac_interpolation_delete(interpolation);

    yac_interp_weights_delete(weights);
  }

  { // test with stencil that writes to multiple locations on a single process

    struct yac_interp_weights * weights =
      yac_interp_weights_new(
        MPI_COMM_WORLD, YAC_LOC_CELL,
        (enum yac_location[]){YAC_LOC_CELL, YAC_LOC_CORNER}, 2);

    switch (comm_rank) {
      case (2): {
        // tgt global_id: 0
        //   src_field 0:
        //     0.7 * src global_id 5
        //   src_field 1:
        //     0.3 * src global_id 2
        struct remote_point remote_tgt_points[] =
          {{.global_id = 1,
            .data =
              {.count = 10,
               .data.multi =
                 (struct remote_point_info[]){{.rank = 0, .orig_pos = 1},
                                              {.rank = 2, .orig_pos = 1},
                                              {.rank = 2, .orig_pos = 2},
                                              {.rank = 1, .orig_pos = 2},
                                              {.rank = 0, .orig_pos = 5},
                                              {.rank = 1, .orig_pos = 4},
                                              {.rank = 2, .orig_pos = 0},
                                              {.rank = 3, .orig_pos = 0},
                                              {.rank = 0, .orig_pos = 3},
                                              {.rank = 4, .orig_pos = 1}}}}};
        struct remote_points tgts = {
          .data = remote_tgt_points,
          .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
        };
        size_t num_src_per_field_per_tgt[2][1] = {{1},{1}};
        struct remote_point srcs_per_field_data[2][1] =
          {{{.global_id = 5,
             .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 1}}}},
           {{.global_id = 2,
             .data = {.count = 1, .data.single = {.rank = 1, .orig_pos = 2}}}}};
        struct remote_point * srcs_per_field[] =
          {srcs_per_field_data[0],srcs_per_field_data[1]};
        double w[] = {0.7,0.3};
        yac_interp_weights_add_wsum_mf(
          weights, &tgts, &num_src_per_field_per_tgt[0][0],
          srcs_per_field, w, 2);
        break;
      }
      default: break;
    };

    for (size_t i = 0; i < reorder_type_count; ++i) {

      size_t collection_size = 1;
      double frac_mask_fallback_value = YAC_FRAC_MASK_NO_VALUE;
      double scaling_factor = 1.0;
      double scaling_summand = 0.0;
      struct yac_interpolation * interpolation =
        yac_interp_weights_get_interpolation(
          weights, YAC_MAPPING_ON_SRC, collection_size,
          frac_mask_fallback_value, scaling_factor, scaling_summand,
          NULL, 1, 1);
      struct yac_interpolation_exchange * interpolation_exchange;
      struct yac_interp_weights_data interp_weights_data;
      yac_interp_weights_get_interpolation_raw(
        weights, collection_size, frac_mask_fallback_value,
        scaling_factor, scaling_summand, NULL,
        &interpolation_exchange, &interp_weights_data, 1, 1);

      double src_field_data[2][10] =
        {{1,2,3,4,5,6,7,8,9,10},{10,20,30,40,50,60,70,80,90,100}};
      double tgt_field_data[] = {-1,-1,-1,-1,-1,-1,-1,-1};

      double ref_target_value = 0.7*2.0 + 0.3*30.0;
      double ref_tgt_field_data[5][8] =
        {{-1,ref_target_value,-1,ref_target_value,-1,ref_target_value,-1,-1},
        {-1,-1,ref_target_value,-1,ref_target_value,-1,-1,-1},
        {ref_target_value,ref_target_value,ref_target_value,-1,-1,-1,-1,-1},
        {ref_target_value,-1,-1,-1,-1,-1,-1,-1},
        {-1,ref_target_value,-1,-1,-1,-1,-1,-1}};
      double * src_field[2] = {(comm_rank==0)?(src_field_data[0]):NULL,
                              (comm_rank==1)?(src_field_data[1]):NULL};
      double * tgt_field = tgt_field_data;
      double ** src_fields = src_field;

      yac_interpolation_execute(
        interpolation, &src_fields, &tgt_field);

      if (sizeof(tgt_field_data) != sizeof(ref_tgt_field_data[comm_rank]))
        PUT_ERR("invalid reference data");

      for (size_t i = 0;
            i < sizeof(tgt_field_data)/sizeof(tgt_field_data[0]); ++i)
        if (fabs(tgt_field_data[i] - ref_tgt_field_data[comm_rank][i]) > TOL)
          PUT_ERR("ERROR in multi source field support");

      for (size_t i = 0;
            i < sizeof(tgt_field_data)/sizeof(tgt_field_data[0]); ++i)
        tgt_field_data[i] = -1;
      utest_interpolation_execute_raw(
        interpolation_exchange, interp_weights_data, collection_size,
        &src_fields, &tgt_field);

      for (size_t i = 0;
            i < sizeof(tgt_field_data)/sizeof(tgt_field_data[0]); ++i)
        if (fabs(tgt_field_data[i] - ref_tgt_field_data[comm_rank][i]) > TOL)
          PUT_ERR("ERROR in multi source field support");

      yac_interp_weights_data_free(interp_weights_data);
      yac_interpolation_exchange_delete(
        interpolation_exchange, "interpolation raw");
      yac_interpolation_delete(interpolation);
    }

    yac_interp_weights_delete(weights);
  }

  { // test writing empty stencil data to file

    struct yac_interp_weights * weights =
      yac_interp_weights_new(
        MPI_COMM_WORLD, YAC_LOC_CELL, (enum yac_location[]){YAC_LOC_CELL}, 1);

    struct remote_points tgts = {
      .data = NULL,
      .count = 0,
    };
    yac_interp_weights_add_fixed(weights, &tgts, -1.0);

    // check writing and reading
    {
      char const * weight_file_name =
        "test_interp_weights_parallel_empty.nc";
      // write weights to file
      yac_interp_weights_write_to_file(
        weights, weight_file_name, grid_names[0], grid_names[1], 0, 0,
        YAC_WEIGHT_FILE_ERROR);

      MPI_Barrier(MPI_COMM_WORLD);

      if (comm_rank == 0) unlink(weight_file_name);
    }

    yac_interp_weights_delete(weights);
  }

  { // test writing empty stencil data to file

    struct yac_interp_weights * weights =
      yac_interp_weights_new(
        MPI_COMM_WORLD, YAC_LOC_CELL, (enum yac_location[]){YAC_LOC_CELL}, 1);

    struct remote_points tgts = {
      .data = NULL,
      .count = 0,
    };
    yac_interp_weights_add_direct(weights, &tgts, NULL);

    // check writing and reading
    {
      char const * weight_file_name =
        "test_interp_weights_parallel_empty.nc";
      // write weights to file
      yac_interp_weights_write_to_file(
        weights, weight_file_name, grid_names[0], grid_names[1], 0, 0,
        YAC_WEIGHT_FILE_ERROR);

      MPI_Barrier(MPI_COMM_WORLD);

      if (comm_rank == 0) unlink(weight_file_name);
    }

    yac_interp_weights_delete(weights);
  }

  { // test writing empty stencil data to file

    struct yac_interp_weights * weights =
      yac_interp_weights_new(
        MPI_COMM_WORLD, YAC_LOC_CELL, (enum yac_location[]){YAC_LOC_CELL}, 1);

    struct remote_points tgts = {
      .data = NULL,
      .count = 0,
    };
    size_t num_src_per_tgt[] = {0};
    yac_interp_weights_add_sum(weights, &tgts, num_src_per_tgt, NULL);

    // check writing and reading
    {
      char const * weight_file_name =
        "test_interp_weights_parallel_empty.nc";
      // write weights to file
      yac_interp_weights_write_to_file(
        weights, weight_file_name, grid_names[0], grid_names[1], 0, 0,
        YAC_WEIGHT_FILE_ERROR);

      MPI_Barrier(MPI_COMM_WORLD);

      if (comm_rank == 0) unlink(weight_file_name);
    }

    yac_interp_weights_delete(weights);
  }

  { // test writing empty stencil data to file

    struct yac_interp_weights * weights =
      yac_interp_weights_new(
        MPI_COMM_WORLD, YAC_LOC_CELL, (enum yac_location[]){YAC_LOC_CELL}, 1);

    struct remote_point remote_tgt_points[] =
      {{.global_id = 0,
        .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}}};
    struct remote_points tgts = {
      .data = remote_tgt_points,
      .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
    };
    size_t num_src_per_tgt[] = {0};
    yac_interp_weights_add_sum(weights, &tgts, num_src_per_tgt, NULL);

    // check writing and reading
    {
      char const * weight_file_name =
        "test_interp_weights_parallel_empty.nc";
      // write weights to file
      yac_interp_weights_write_to_file(
        weights, weight_file_name, grid_names[0], grid_names[1], 0, 0,
        YAC_WEIGHT_FILE_ERROR);

      MPI_Barrier(MPI_COMM_WORLD);

      if (comm_rank == 0) unlink(weight_file_name);
    }

    yac_interp_weights_delete(weights);
  }

  { // test writing empty stencil data to file

    struct yac_interp_weights * weights =
      yac_interp_weights_new(
        MPI_COMM_WORLD, YAC_LOC_CELL, (enum yac_location[]){YAC_LOC_CELL}, 1);

    struct remote_point remote_tgt_points[] =
      {{.global_id = 0,
        .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}}};
    struct remote_points tgts = {
      .data = remote_tgt_points,
      .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
    };
    size_t num_src_per_tgt[] = {0};
    yac_interp_weights_add_wsum(weights, &tgts, num_src_per_tgt, NULL, NULL);

    // check writing and reading
    {
      char const * weight_file_name =
        "test_interp_weights_parallel_empty.nc";
      // write weights to file
      yac_interp_weights_write_to_file(
        weights, weight_file_name, grid_names[0], grid_names[1], 0, 0,
        YAC_WEIGHT_FILE_ERROR);

      MPI_Barrier(MPI_COMM_WORLD);

      if (comm_rank == 0) unlink(weight_file_name);
    }

    yac_interp_weights_delete(weights);
  }

  { // test writing empty stencil data to file

    struct yac_interp_weights * weights =
      yac_interp_weights_new(
        MPI_COMM_WORLD, YAC_LOC_CELL, (enum yac_location[]){YAC_LOC_CELL}, 1);

    struct remote_points tgts = {
      .data = NULL,
      .count = 0,
    };
    size_t num_src_per_tgt[] = {0};
    yac_interp_weights_add_wsum(weights, &tgts, num_src_per_tgt, NULL, NULL);

    // check writing and reading
    {
      char const * weight_file_name =
        "test_interp_weights_parallel_empty.nc";
      // write weights to file
      yac_interp_weights_write_to_file(
        weights, weight_file_name, grid_names[0], grid_names[1], 0, 0,
        YAC_WEIGHT_FILE_ERROR);

      MPI_Barrier(MPI_COMM_WORLD);

      if (comm_rank == 0) unlink(weight_file_name);
    }

    yac_interp_weights_delete(weights);
  }

  { // checking for a bug that occured then writing weights to file

    struct yac_interp_weights * weights =
      yac_interp_weights_new(
        MPI_COMM_WORLD, YAC_LOC_CELL,
        (enum yac_location[]){YAC_LOC_CELL, YAC_LOC_CORNER}, 2);

    if (comm_rank == 0) {

      {
        struct remote_point tgt_points[] =
          {{.global_id = 4,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 4}}},
           {.global_id = 5,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 5}}}};
        struct remote_points tgts = {
          .data = tgt_points,
          .count = sizeof(tgt_points) / sizeof(tgt_points[0]),
        };
        struct remote_point srcs[] =
          {{.global_id = 4,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 4}}},
           {.global_id = 5,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 5}}}};
        yac_interp_weights_add_direct(weights, &tgts, srcs);
      }

      {
        struct remote_point tgt_points[] =
          {{.global_id = 0,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}},
           {.global_id = 1,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 1}}},
           {.global_id = 2,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 2}}},
           {.global_id = 3,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 3}}},
           {.global_id = 6,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 6}}},
           {.global_id = 7,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 7}}},
           {.global_id = 8,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 8}}},
           {.global_id = 9,
            .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 9}}}};
        struct remote_points tgts = {
          .data = tgt_points,
          .count = sizeof(tgt_points) / sizeof(tgt_points[0]),
        };
        yac_interp_weights_add_fixed(weights, &tgts, -1.0);
      }
    }

    char const * weight_file_name =
      "test_interp_weights_parallel_bug_check.nc";
    // write weights to file
    yac_interp_weights_write_to_file(
      weights, weight_file_name, grid_names[0], grid_names[1], 0, 0,
      YAC_WEIGHT_FILE_ERROR);
    MPI_Barrier(MPI_COMM_WORLD);
    if (comm_rank == 0) unlink(weight_file_name);

    yac_interp_weights_delete(weights);
  }



  { // test handling of existing weight file
    // (the last test is expected to call abort, so no further tests afterwards)

    char const * src_grid_name = "src_grid";
    char const * tgt_grid_name = "tgt_grid";

    struct yac_interp_weights * weights_existing =
      yac_interp_weights_new(
        MPI_COMM_WORLD, YAC_LOC_CELL, (enum yac_location[]){YAC_LOC_CELL}, 1);
    {
      struct remote_point remote_tgt_points[] =
        {{.global_id = 0,
          .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}}};
      struct remote_points tgts = {
        .data = remote_tgt_points,
        .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
      };
      double fixed_value = 1.0;
      yac_interp_weights_add_fixed(weights_existing, &tgts, fixed_value);
    }

    struct yac_interp_weights * weights =
      yac_interp_weights_new(
        MPI_COMM_WORLD, YAC_LOC_CELL, (enum yac_location[]){YAC_LOC_CELL}, 1);
    {
      struct remote_point remote_tgt_points[] =
        {{.global_id = 0,
          .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 0}}}};
      struct remote_points tgts = {
        .data = remote_tgt_points,
        .count = sizeof(remote_tgt_points) / sizeof(remote_tgt_points[0]),
      };
      size_t num_src_per_tgt[] = {1};
      struct remote_point srcs[] =
        {{.global_id = 0,
          .data = {.count = 1, .data.single = {.rank = 0, .orig_pos = 2}}}};
      yac_interp_weights_add_sum(weights, &tgts, num_src_per_tgt, srcs);
    }

    for (int with_existing_weight_file = 0; with_existing_weight_file < 2;
         ++with_existing_weight_file) {

      // the last test should be "error on existing file"
      enum yac_weight_file_on_existing on_existing[] = {
        YAC_WEIGHT_FILE_KEEP,
        YAC_WEIGHT_FILE_OVERWRITE,
        YAC_WEIGHT_FILE_ERROR};
      enum {NUM_ON_EXISTING = sizeof(on_existing) / sizeof(on_existing[0])};

      for (size_t on_existing_idx = 0; on_existing_idx < NUM_ON_EXISTING;
           ++on_existing_idx) {

        if (with_existing_weight_file)
          yac_interp_weights_write_to_file(
            weights_existing, weight_file_name_on_existing,
            src_grid_name, tgt_grid_name, 0, 0, YAC_WEIGHT_FILE_OVERWRITE);

        enum yac_weight_file_on_existing curr_on_existing =
          on_existing[on_existing_idx];

        // in case we test for an error for an existing weight file,
        // we expect the custom abort handler to be called
        if ((curr_on_existing == YAC_WEIGHT_FILE_ERROR) &&
            with_existing_weight_file)
          yac_set_abort_handler((yac_abort_func)on_existing_abort_handler);
        else
          yac_restore_default_abort_handler();

        MPI_Barrier(MPI_COMM_WORLD);
        yac_interp_weights_write_to_file(
          weights, weight_file_name_on_existing, src_grid_name, tgt_grid_name,
          0, 0, curr_on_existing);

        int contains_links = 0;
        int contains_fixed = 0;
        if (!((curr_on_existing == YAC_WEIGHT_FILE_ERROR) &&
              with_existing_weight_file))
          utest_get_basic_weight_file_info(
            weight_file_name_on_existing, &contains_links, &contains_fixed);

        switch(curr_on_existing) {
          case(YAC_WEIGHT_FILE_KEEP):
            if (with_existing_weight_file) {
              // the orig weight file only contains a fixed value but not links
              if (!contains_fixed || contains_links)
                PUT_ERR("error in yac_interp_weights_write_to_file");
            } else {
              if (contains_fixed || !contains_links)
                PUT_ERR("error in yac_interp_weights_write_to_file");
            }
            break;
          case(YAC_WEIGHT_FILE_OVERWRITE):
            // the new weight file only contains links but no fixed value
            if (contains_fixed || !contains_links)
              PUT_ERR("error in yac_interp_weights_write_to_file");
            break;
          case(YAC_WEIGHT_FILE_ERROR):
            // if there was no existing weight file
            if (!with_existing_weight_file) {
              // the new weight file only contains links but no fixed value
              if (contains_fixed || !contains_links)
                PUT_ERR("error in yac_interp_weights_write_to_file");
            } else {

              // at least on process should have called the custom abort handler
              int abort_handler_was_called = 0;
              MPI_Allreduce(
                MPI_IN_PLACE, &abort_handler_was_called, 1,
                MPI_INT, MPI_MAX, MPI_COMM_WORLD);

              if (!abort_handler_was_called)
                PUT_ERR("error in yac_interp_weights_write_to_file");

              MPI_Barrier(MPI_COMM_WORLD);
              if (comm_rank == 0) unlink(weight_file_name_on_existing);

              xt_finalize();
              MPI_Finalize();
              exit(TEST_EXIT_CODE);
            }
            break;
          default:
            PUT_ERR("internal test error");
            break;
        } // switch(curr_on_existing)

        MPI_Barrier(MPI_COMM_WORLD);
        if (comm_rank == 0) unlink(weight_file_name_on_existing);
      } // on_existing_idx
    } // with_existing_weight_file
  }

  // The test above is expected to call the abort handler, hence the code
  // should never reach this point. New tests should be added above.

  PUT_ERR("internal test error");
  xt_finalize();
  MPI_Finalize();

  return TEST_EXIT_CODE;
}

static int utest_check_fixed_results_(
  struct yac_basic_grid_data * grid, int is_src,
  struct yac_interpolation * interpolation,
  double * field_data) {

  int err_count = 0;

  if (is_src)
    for (size_t j = 0; j < grid->num_cells; ++j)
      field_data[j] = (double)1;
  else
    for (size_t j = 0; j < grid->num_cells; ++j)
      field_data[j] = (double)-3;

  double ** src_fields, * tgt_field;
  if (is_src) {
    src_fields = &field_data;
    tgt_field = NULL;
  } else {
    src_fields = NULL;
    tgt_field = field_data;
  }

  yac_interpolation_execute(interpolation, &src_fields, &tgt_field);

  if (!is_src) {
    for (size_t j = 0; j < grid->num_cells; ++j) {
      if (grid->core_cell_mask[j]) {
        if (grid->cell_ids[j] & 1) {
          if (tgt_field[j] != -1.0) ++err_count;
        } else {
          if (tgt_field[j] != -2.0) ++err_count;
        }
      } else {
        if (tgt_field[j] != -3.0) ++err_count;
      }
    }
  }

  return err_count;
}

static double ** allocate_src_data_raw(
  size_t collection_size, struct yac_interp_weights_data interp_weights_data) {

  double ** src_data_raw;
  static double dummy;
  int with_frac_mask =
    YAC_FRAC_MASK_VALUE_IS_VALID(interp_weights_data.frac_mask_fallback_value);
  src_data_raw =
    xmalloc(
      (1 + with_frac_mask) * collection_size *
      interp_weights_data.num_src_fields * sizeof(*src_data_raw));
  {
    size_t k = 0;
    for (size_t i = 0; i < collection_size; ++i) {
      for (size_t j = 0; j < interp_weights_data.num_src_fields; ++j, ++k) {
        src_data_raw[k] =
          interp_weights_data.src_field_buffer_size[j]?
            xmalloc(
              interp_weights_data.src_field_buffer_size[j] *
              sizeof(*(src_data_raw[k]))):&dummy;
      }
    }
    if (with_frac_mask) {
      for (size_t i = 0; i < collection_size; ++i) {
        for (size_t j = 0; j < interp_weights_data.num_src_fields; ++j, ++k) {
          src_data_raw[k] =
            interp_weights_data.src_field_buffer_size[j]?
              xmalloc(
                interp_weights_data.src_field_buffer_size[j] *
                sizeof(*(src_data_raw[k]))):&dummy;
        }
      }
    }
  }
  return src_data_raw;
}

static void deallocate_src_data_raw(
  size_t collection_size, struct yac_interp_weights_data interp_weights_data,
  double ** src_data_raw) {

  int with_frac_mask =
    YAC_FRAC_MASK_VALUE_IS_VALID(interp_weights_data.frac_mask_fallback_value);
  size_t k = 0;
  for (size_t i = 0; i < collection_size; ++i)
    for (size_t j = 0; j < interp_weights_data.num_src_fields; ++j, ++k)
      if (interp_weights_data.src_field_buffer_size[j])
        free(src_data_raw[k]);
  if (with_frac_mask)
    for (size_t i = 0; i < collection_size; ++i)
      for (size_t j = 0; j < interp_weights_data.num_src_fields; ++j, ++k)
        if (interp_weights_data.src_field_buffer_size[j])
          free(src_data_raw[k]);
  free(src_data_raw);
}

#if defined __NVCOMPILER && (__NVCOMPILER_MAJOR__ <= 24)
// NVHPC has problems compiling this routine if it is inlined
#pragma noinline
#endif
static void apply_interp_weights_data(
  size_t collection_size, struct yac_interp_weights_data interp_weights_data,
  double ** src_data_raw, double ** tgt_field) {

#define FRAC_MASK_TOL (1e-12)

  int with_frac_mask =
    YAC_FRAC_MASK_VALUE_IS_VALID(interp_weights_data.frac_mask_fallback_value);

  // apply fixed values
  for (size_t i = 0; i < collection_size; ++i) {
    for (size_t j = 0, l = 0; j < interp_weights_data.num_fixed_values; ++j) {
      for (size_t k = 0; k < interp_weights_data.num_tgt_per_fixed_value[j];
          ++k, ++l) {
        tgt_field[i][interp_weights_data.tgt_idx_fixed[l]] =
          interp_weights_data.fixed_values[j];
      }
    }
  }

  // apply weights
  if (with_frac_mask) {
    for (size_t i = 0; i < collection_size; ++i) {
      size_t wgt_offset = 0;
      for (size_t j = 0; j < interp_weights_data.num_wgt_tgt; ++j) {
        double tgt_value = 0.0;
        double frac_weight_sum = 0.0;
#if defined __NVCOMPILER && (__NVCOMPILER_MAJOR__ <= 24)
// NVHPC has problems compiling the following code
#pragma novector
#endif
        for (size_t k = 0; k < interp_weights_data.num_src_per_tgt[j];
              ++k, ++wgt_offset) {
          tgt_value +=
            src_data_raw
              [i * interp_weights_data.num_src_fields +
              interp_weights_data.src_field_idx[wgt_offset]]
              [interp_weights_data.src_idx[wgt_offset]] *
            interp_weights_data.weights[wgt_offset];
          frac_weight_sum +=
            src_data_raw
              [(collection_size * interp_weights_data.num_src_fields) +
              i * interp_weights_data.num_src_fields +
              interp_weights_data.src_field_idx[wgt_offset]]
              [interp_weights_data.src_idx[wgt_offset]] *
            interp_weights_data.weights[wgt_offset];
        }
#if defined __NVCOMPILER && (__NVCOMPILER_MAJOR__ <= 24)
#pragma vector
#endif
        tgt_field[i][interp_weights_data.wgt_tgt_idx[j]] =
          (fabs(frac_weight_sum) > FRAC_MASK_TOL)?
          (interp_weights_data.scaling_factor * (tgt_value / frac_weight_sum) +
          interp_weights_data.scaling_summand):
          interp_weights_data.frac_mask_fallback_value;
      }
    }
  } else {
    for (size_t i = 0; i < collection_size; ++i) {
      size_t wgt_offset = 0;
      for (size_t j = 0; j < interp_weights_data.num_wgt_tgt; ++j) {
        double tgt_value = 0.0;
#if defined __NVCOMPILER && (__NVCOMPILER_MAJOR__ <= 24)
// NVHPC has problems compiling the following code
#pragma novector
#endif
        for (size_t k = 0; k < interp_weights_data.num_src_per_tgt[j];
              ++k, ++wgt_offset) {
          tgt_value +=
            src_data_raw
              [i * interp_weights_data.num_src_fields +
              interp_weights_data.src_field_idx[wgt_offset]]
              [interp_weights_data.src_idx[wgt_offset]] *
            interp_weights_data.weights[wgt_offset];
        }
#if defined __NVCOMPILER && (__NVCOMPILER_MAJOR__ <= 24)
#pragma vector
#endif
        tgt_field[i][interp_weights_data.wgt_tgt_idx[j]] =
          interp_weights_data.scaling_factor * tgt_value +
          interp_weights_data.scaling_summand;
      }
    }
  }
#undef FRAC_MASK_TOL
}

static void utest_interpolation_execute_frac_raw(
  struct yac_interpolation_exchange * interpolation_exchange,
  struct yac_interp_weights_data interp_weights_data,
  size_t collection_size, double *** src_fields, double *** src_frac_masks,
  double ** tgt_field) {

  int with_frac_mask =
    YAC_FRAC_MASK_VALUE_IS_VALID(interp_weights_data.frac_mask_fallback_value);

  double const * src_fields_2d[
    (1 + with_frac_mask) * collection_size *
    interp_weights_data.num_src_fields];
  double dummy;
  {
    size_t k = 0;
    for (size_t i = 0; i < collection_size; ++i)
      for (size_t j = 0; j < interp_weights_data.num_src_fields; ++j, ++k)
        src_fields_2d[k] = src_fields?src_fields[i][j]:&dummy;
    if (with_frac_mask)
      for (size_t i = 0; i < collection_size; ++i)
        for (size_t j = 0; j < interp_weights_data.num_src_fields; ++j, ++k)
          src_fields_2d[k] = src_frac_masks?src_frac_masks[i][j]:&dummy;
  }

  // allocate receive buffer
  double ** src_data_raw =
    allocate_src_data_raw(collection_size, interp_weights_data);

  // execute interpolation
  yac_interpolation_exchange_execute(
    interpolation_exchange, src_fields_2d, src_data_raw,
    "utest_interpolation_execute_raw");

  // compute target field
  apply_interp_weights_data(
    collection_size, interp_weights_data, src_data_raw, tgt_field);

  // free receive buffer
  deallocate_src_data_raw(collection_size, interp_weights_data, src_data_raw);
}

static void utest_interpolation_execute_raw(
  struct yac_interpolation_exchange * interpolation_exchange,
  struct yac_interp_weights_data interp_weights_data,
  size_t collection_size, double *** src_fields, double ** tgt_field) {

  utest_interpolation_execute_frac_raw(
    interpolation_exchange, interp_weights_data, collection_size,
    src_fields, NULL, tgt_field);
}

static void interpolation_put_frac_raw(
  struct yac_interpolation_exchange * interpolation_exchange,
  struct yac_interp_weights_data interp_weights_data,
  size_t collection_size, double *** src_fields, double *** src_frac_masks) {

  int with_frac_mask =
    YAC_FRAC_MASK_VALUE_IS_VALID(interp_weights_data.frac_mask_fallback_value);

  double const * src_fields_2d[
    (1 + with_frac_mask) * collection_size * interp_weights_data.num_src_fields];
  static double dummy;
  {
    size_t k = 0;
    for (size_t i = 0; i < collection_size; ++i)
      for (size_t j = 0; j < interp_weights_data.num_src_fields; ++j, ++k)
        src_fields_2d[k] = src_fields?src_fields[i][j]:&dummy;
    if (with_frac_mask)
      for (size_t i = 0; i < collection_size; ++i)
        for (size_t j = 0; j < interp_weights_data.num_src_fields; ++j, ++k)
          src_fields_2d[k] = src_frac_masks?src_frac_masks[i][j]:&dummy;
  }

  // send source data
  yac_interpolation_exchange_execute_put(
    interpolation_exchange, src_fields_2d, "interpolation_put_frac_raw");
}

static void interpolation_put_raw(
  struct yac_interpolation_exchange * interpolation_exchange,
  struct yac_interp_weights_data interp_weights_data,
  size_t collection_size, double *** src_fields) {

  interpolation_put_frac_raw(
    interpolation_exchange, interp_weights_data, collection_size,
    src_fields, NULL);
}

static void interpolation_get_raw(
  struct yac_interpolation_exchange * interpolation_exchange,
  struct yac_interp_weights_data interp_weights_data,
  size_t collection_size, double ** tgt_field) {

  // allocate receive buffer
  double ** src_data_raw =
    allocate_src_data_raw(collection_size, interp_weights_data);

  // receive source data
  yac_interpolation_exchange_execute_get(
    interpolation_exchange, src_data_raw, "interpolation_get_raw");

  // compute target field
  apply_interp_weights_data(
    collection_size, interp_weights_data, src_data_raw, tgt_field);

  // free receive buffer
  deallocate_src_data_raw(collection_size, interp_weights_data, src_data_raw);
}

static int utest_check_fixed_results_raw_(
  struct yac_basic_grid_data * grid, int is_src,
  struct yac_interpolation_exchange * interpolation_exchange,
  double * field_data, struct yac_interp_weights_data interp_weights_data) {

  int err_count = 0;

  if (is_src)
    for (size_t j = 0; j < grid->num_cells; ++j)
      field_data[j] = (double)1;
  else
    for (size_t j = 0; j < grid->num_cells; ++j)
      field_data[j] = (double)-3;

  double ** src_fields, * tgt_field;
  if (is_src) {
    src_fields = &field_data;
    tgt_field = NULL;
  } else {
    src_fields = NULL;
    tgt_field = field_data;
  }

  size_t collection_size = 1;
  utest_interpolation_execute_raw(
    interpolation_exchange, interp_weights_data, collection_size,
    &src_fields, &tgt_field);

  if (!is_src) {
    for (size_t j = 0; j < grid->num_cells; ++j) {
      if (grid->core_cell_mask[j]) {
        if (grid->cell_ids[j] & 1) {
          if (tgt_field[j] != -1.0) ++err_count;
        } else {
          if (tgt_field[j] != -2.0) ++err_count;
        }
      } else {
        if (tgt_field[j] != -3.0) ++err_count;
      }
    }
  }

  return err_count;
}

static int utest_check_fixed_results(
  struct yac_basic_grid_data * grid, int is_src,
  struct yac_interp_weights * weights,
  enum yac_interp_weights_reorder_type reorder_type,
  double * field_data) {

  int err_count = 0;

  size_t const collection_size = 1;
  double const frac_mask_fallback_value = YAC_FRAC_MASK_NO_VALUE;
  double const scaling_factor = 1.0;
  double const scaling_summand = 0.0;
  char const * yaxt_exchanger_name = NULL;
  {
    struct yac_interpolation * interpolation =
      yac_interp_weights_get_interpolation(
        weights, reorder_type, collection_size,
        frac_mask_fallback_value, scaling_factor, scaling_summand,
        yaxt_exchanger_name, 1, 1);
    struct yac_interpolation * interpolation_cpy =
      yac_interpolation_copy(interpolation);

    err_count +=
      utest_check_fixed_results_(grid, is_src, interpolation, field_data) +
      utest_check_fixed_results_(grid, is_src, interpolation_cpy, field_data);

    yac_interpolation_delete(interpolation_cpy);
    yac_interpolation_delete(interpolation);
  }

  {
    struct yac_interpolation_exchange * interpolation_exchange;
    struct yac_interp_weights_data interp_weights_data;
    yac_interp_weights_get_interpolation_raw(
      weights, collection_size, frac_mask_fallback_value,
      scaling_factor, scaling_summand, yaxt_exchanger_name,
      &interpolation_exchange, &interp_weights_data, 1, 1);

    err_count +=
      utest_check_fixed_results_raw_(
        grid, is_src, interpolation_exchange, field_data, interp_weights_data);

    yac_interpolation_exchange_delete(
      interpolation_exchange, "utest_check_fixed_results");
    yac_interp_weights_data_free(interp_weights_data);
  }

  return err_count;
}

static int utest_check_direct_results(
  struct yac_basic_grid_data * grid, int is_src,
  struct yac_interp_weights * weights,
  enum yac_interp_weights_reorder_type reorder_type, size_t collection_size,
  double *** src_fields, double *** src_frac_masks, double ** tgt_field) {

  int err_count = 0;

  { // test without fractional mask

    double const frac_mask_fallback_value = YAC_FRAC_MASK_NO_VALUE;
    double const scaling_factor = 1.0;
    double const scaling_summand = 0.0;
    char const * yaxt_exchanger_name = NULL;

    {
      struct yac_interpolation * interpolations[2];
      interpolations[0] =
        yac_interp_weights_get_interpolation(
          weights, reorder_type, collection_size,
          frac_mask_fallback_value, scaling_factor, scaling_summand,
          yaxt_exchanger_name, 1, 1);
      interpolations[1] = yac_interpolation_copy(interpolations[0]);

      for (int interp_idx = 0; interp_idx < 2; ++interp_idx) {

        if (!is_src)
          for (size_t collection_idx = 0; collection_idx < collection_size;
              ++collection_idx)
            for (size_t k = 0; k < grid->num_cells; ++k)
              tgt_field[collection_idx][k] = -1;

        yac_interpolation_execute(
          interpolations[interp_idx], src_fields, tgt_field);

        if (!is_src) {
          for (size_t collection_idx = 0; collection_idx < collection_size;
              ++collection_idx) {
            for (size_t j = 0; j < grid->num_cells; ++j) {
              if (grid->core_cell_mask[j]) {
                if (tgt_field[collection_idx][j] !=
                    (double)(grid->cell_ids[j])) err_count++;
              } else {
                if (tgt_field[collection_idx][j] != -1.0) err_count++;
              }
            }
          }
        }

        if (is_src) {
          yac_interpolation_execute_put(interpolations[interp_idx], src_fields);
        } else {
          for (size_t collection_idx = 0; collection_idx < collection_size;
              ++collection_idx)
            for (size_t k = 0; k < grid->num_cells; ++k)
              tgt_field[collection_idx][k] = -1;
          yac_interpolation_execute_get(interpolations[interp_idx], tgt_field);
          for (size_t collection_idx = 0; collection_idx < collection_size;
              ++collection_idx) {
            for (size_t j = 0; j < grid->num_cells; ++j) {
              if (grid->core_cell_mask[j]) {
                if (tgt_field[collection_idx][j] !=
                    (double)(grid->cell_ids[j])) err_count++;
              } else {
                if (tgt_field[collection_idx][j] != -1.0) err_count++;
              }
            }
          }
        }
        yac_interpolation_delete(interpolations[interp_idx]);
      }
    }

    { // test raw interpolation
      struct yac_interpolation_exchange * interpolation_exchange;
      struct yac_interp_weights_data interp_weights_data;
      yac_interp_weights_get_interpolation_raw(
        weights, collection_size, frac_mask_fallback_value,
        scaling_factor, scaling_summand, yaxt_exchanger_name,
        &interpolation_exchange, &interp_weights_data, 1, 1);

      { // exchange

        for (int exchange_type = 0; exchange_type < 2; ++exchange_type) {

          if (!is_src)
            for (size_t collection_idx = 0; collection_idx < collection_size;
                ++collection_idx)
              for (size_t k = 0; k < grid->num_cells; ++k)
                tgt_field[collection_idx][k] = -1;

          if (exchange_type) {
            utest_interpolation_execute_raw(
              interpolation_exchange, interp_weights_data, collection_size,
              src_fields, tgt_field);
          } else {
            if (is_src)
              interpolation_put_raw(
                interpolation_exchange, interp_weights_data, collection_size,
                src_fields);
            else
              interpolation_get_raw(
                interpolation_exchange, interp_weights_data, collection_size,
                tgt_field);
          }

          if (!is_src) {
            for (size_t collection_idx = 0; collection_idx < collection_size;
                ++collection_idx) {
              for (size_t j = 0; j < grid->num_cells; ++j) {
                if (grid->core_cell_mask[j]) {
                  if (tgt_field[collection_idx][j] !=
                      (double)(grid->cell_ids[j])) err_count++;
                } else {
                  if (tgt_field[collection_idx][j] != -1.0) err_count++;
                }
              }
            }
          }
        }
      }
      yac_interp_weights_data_free(interp_weights_data);
      yac_interpolation_exchange_delete(
        interpolation_exchange, "utest_check_direct_results");
    }
  }

  { // test with fractional mask

    double frac_mask_fallback_value = strtod("Inf", NULL);
    double scaling_factor = 1.0;
    double scaling_summand = 0.0;
    char const * yaxt_exchanger_name = NULL;

    struct yac_interpolation * interpolations[2];
    interpolations[0] =
      yac_interp_weights_get_interpolation(
        weights, reorder_type, collection_size, frac_mask_fallback_value,
        scaling_factor, scaling_summand, yaxt_exchanger_name, 1, 1);
    interpolations[1] = yac_interpolation_copy(interpolations[0]);

    for (int interp_idx = 0; interp_idx < 2; ++interp_idx) {

      if (!is_src)
        for (size_t collection_idx = 0; collection_idx < collection_size;
            ++collection_idx)
          for (size_t k = 0; k < grid->num_cells; ++k)
            tgt_field[collection_idx][k] = -1;

      yac_interpolation_execute_frac(
        interpolations[interp_idx], src_fields, src_frac_masks, tgt_field);

      if (!is_src) {
        for (size_t collection_idx = 0; collection_idx < collection_size;
            ++collection_idx) {
          for (size_t j = 0; j < grid->num_cells; ++j) {
            if (grid->core_cell_mask[j]) {
              if (grid->cell_ids[j]&1) {
                if (tgt_field[collection_idx][j] !=
                    (double)(grid->cell_ids[j])) err_count++;
              } else {
                if (!isinf(tgt_field[collection_idx][j])) err_count++;
              }
            } else {
              if (tgt_field[collection_idx][j] != -1.0) err_count++;
            }
          }
        }
      }

      if (is_src) {
        yac_interpolation_execute_put_frac(
            interpolations[interp_idx], src_fields, src_frac_masks);
      } else {
        for (size_t collection_idx = 0; collection_idx < collection_size;
            ++collection_idx)
          for (size_t k = 0; k < grid->num_cells; ++k)
            tgt_field[collection_idx][k] = -1;
        yac_interpolation_execute_get(interpolations[interp_idx], tgt_field);
        for (size_t collection_idx = 0; collection_idx < collection_size;
            ++collection_idx) {
          for (size_t j = 0; j < grid->num_cells; ++j) {
            if (grid->core_cell_mask[j]) {
              if (grid->cell_ids[j]&1) {
                if (tgt_field[collection_idx][j] !=
                    (double)(grid->cell_ids[j])) err_count++;
              } else {
                if (!isinf(tgt_field[collection_idx][j])) err_count++;
              }
            } else {
              if (tgt_field[collection_idx][j] != -1.0) err_count++;
            }
          }
        }
      }
      yac_interpolation_delete(interpolations[interp_idx]);
    }

    { // test raw interpolation
      struct yac_interpolation_exchange * interpolation_exchange;
      struct yac_interp_weights_data interp_weights_data;
      yac_interp_weights_get_interpolation_raw(
        weights, collection_size, frac_mask_fallback_value,
        scaling_factor, scaling_summand, yaxt_exchanger_name,
        &interpolation_exchange, &interp_weights_data, 1, 1);

      { // exchange

        for (int exchange_type = 0; exchange_type < 2; ++exchange_type) {

          if (!is_src)
            for (size_t collection_idx = 0; collection_idx < collection_size;
                ++collection_idx)
              for (size_t k = 0; k < grid->num_cells; ++k)
                tgt_field[collection_idx][k] = -1;

          if (exchange_type) {
            utest_interpolation_execute_frac_raw(
              interpolation_exchange, interp_weights_data, collection_size,
              src_fields, src_frac_masks, tgt_field);
          } else {
            if (is_src)
              interpolation_put_frac_raw(
                interpolation_exchange, interp_weights_data, collection_size,
                src_fields, src_frac_masks);
            else
              interpolation_get_raw(
                interpolation_exchange, interp_weights_data, collection_size,
                tgt_field);
          }

          if (!is_src) {
            for (size_t collection_idx = 0; collection_idx < collection_size;
                ++collection_idx) {
              for (size_t j = 0; j < grid->num_cells; ++j) {
                if (grid->core_cell_mask[j]) {
                  if (grid->cell_ids[j]&1) {
                    if (tgt_field[collection_idx][j] !=
                        (double)(grid->cell_ids[j])) err_count++;
                  } else {
                    if (!isinf(tgt_field[collection_idx][j])) err_count++;
                  }
                } else {
                  if (tgt_field[collection_idx][j] != -1.0) err_count++;
                }
              }
            }
          }
        }
      }
      yac_interp_weights_data_free(interp_weights_data);
      yac_interpolation_exchange_delete(
        interpolation_exchange, "utest_check_direct_results");
    }
  }

  return err_count;
}

static void utest_get_basic_grid_data(
  char * filename, size_t * num_cells, size_t * num_vertices,
  size_t * num_edges) {

  int ncid;

  yac_nc_open(filename, NC_NOWRITE, &ncid);

  int dimid;
  yac_nc_inq_dimid(ncid, "cell", &dimid);
  YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, num_cells));
  yac_nc_inq_dimid(ncid, "vertex", &dimid);
  YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, num_vertices));
  yac_nc_inq_dimid(ncid, "edge", &dimid);
  YAC_HANDLE_ERROR(nc_inq_dimlen(ncid, dimid, num_edges));

  YAC_HANDLE_ERROR(nc_close(ncid));
}

static struct yac_basic_grid_data
  utest_generate_dummy_grid_data(
    size_t local_num_cells, size_t global_num_cells, size_t num_cells_offset) {

  double coordinates_x[global_num_cells + 1];
  double coordinates_y[2] = {0.0, 1.0};
  size_t num_cells[2] = {global_num_cells, 1};
  size_t local_start[2] = {num_cells_offset, 0};
  size_t local_count[2] = {local_num_cells, 1};

  for (size_t i = 0; i <= global_num_cells; ++i) coordinates_x[i] = (double)i;

  return
    yac_generate_basic_grid_data_reg2d(
      coordinates_x, coordinates_y, num_cells, local_start, local_count, 0);
}

static int utest_check_results_ref(
  int is_src, struct yac_interp_weights * weights,
  enum yac_interp_weights_reorder_type reorder_type, double * src_data,
  double * tgt_data, double * ref_tgt_data, size_t tgt_size,
  int is_active_src, int is_active_tgt) {

  int err_count = 0;

  size_t collection_size = 1;
  double frac_mask_fallback_value = YAC_FRAC_MASK_NO_VALUE;
  double scale_factors[] = {1.0, -1.0, 0.5};
  double scale_summands[] = {0.0, -1.0, 10.0};
  char const * yaxt_exchangers[] = {NULL, "irecv_isend"};

  for (size_t scale_factor_idx = 0;
       scale_factor_idx < sizeof(scale_factors) / sizeof(scale_factors[0]);
       ++scale_factor_idx) {
    for (size_t scale_summand_idx = 0;
         scale_summand_idx < sizeof(scale_summands) / sizeof(scale_summands[0]);
         ++scale_summand_idx) {
      for (size_t yaxt_exchanger_idx = 0;
           yaxt_exchanger_idx < sizeof(yaxt_exchangers) / sizeof(yaxt_exchangers[0]);
           ++yaxt_exchanger_idx) {

        struct yac_interpolation * interpolations[2];
        interpolations[0] =
          yac_interp_weights_get_interpolation(
            weights, reorder_type, collection_size, frac_mask_fallback_value,
            scale_factors[scale_factor_idx], scale_summands[scale_summand_idx],
            yaxt_exchangers[yaxt_exchanger_idx], is_active_src, is_active_tgt);
        interpolations[1] = yac_interpolation_copy(interpolations[0]);
        struct yac_interpolation_exchange * interpolation_exchange;
        struct yac_interp_weights_data interp_weights_data;
        yac_interp_weights_get_interpolation_raw(
          weights, collection_size, frac_mask_fallback_value,
          scale_factors[scale_factor_idx], scale_summands[scale_summand_idx],
          yaxt_exchangers[yaxt_exchanger_idx],
          &interpolation_exchange, &interp_weights_data,
          is_active_src, is_active_tgt);

        for (int interp_idx = 0; interp_idx < 2; ++interp_idx) {

          double * tgt_field[1] = {tgt_data};
          double * src_field[1] = {src_data};
          double ** src_fields[1] = {src_field};

          if (!is_src)
            for (size_t k = 0; k < tgt_size; ++k) tgt_data[k] = -1;

          yac_interpolation_execute(
            interpolations[interp_idx], src_fields, tgt_field);

          if (!is_src)
            for (size_t j = 0; j < tgt_size; ++j)
              if (!(((ref_tgt_data[j] == FALLBACK_VALUE) == (tgt_data[j] == FALLBACK_VALUE)) ||
                    ((ref_tgt_data[j] == FIXED_VALUE) == (tgt_data[j] == FIXED_VALUE)) ||
                    ((ref_tgt_data[j] != FALLBACK_VALUE) && (ref_tgt_data[j] != FIXED_VALUE) &&
                     (fabs(
                       (ref_tgt_data[j] * scale_factors[scale_factor_idx] +
                        scale_summands[scale_summand_idx]) - tgt_data[j]) < TOL))))
                err_count++;

          if (is_src) {
            yac_interpolation_execute_put(
              interpolations[interp_idx], src_fields);
          } else {
            for (size_t k = 0; k < tgt_size; ++k) tgt_data[k] = -1;
            yac_interpolation_execute_get(
              interpolations[interp_idx], tgt_field);
            for (size_t j = 0; j < tgt_size; ++j)
              if (!(((ref_tgt_data[j] == FALLBACK_VALUE) == (tgt_data[j] == FALLBACK_VALUE)) ||
                    ((ref_tgt_data[j] == FIXED_VALUE) == (tgt_data[j] == FIXED_VALUE)) ||
                    ((ref_tgt_data[j] != FALLBACK_VALUE) && (ref_tgt_data[j] != FIXED_VALUE) &&
                     (fabs(
                       (ref_tgt_data[j] * scale_factors[scale_factor_idx] +
                        scale_summands[scale_summand_idx]) - tgt_data[j]) < TOL))))
                err_count++;
          }
          yac_interpolation_delete(interpolations[interp_idx]);
        }

        {
          double * tgt_field[1] = {tgt_data};
          double * src_field[1] = {src_data};
          double ** src_fields[1] = {src_field};

          if (!is_src)
            for (size_t k = 0; k < tgt_size; ++k) tgt_data[k] = -1;

          utest_interpolation_execute_raw(
            interpolation_exchange, interp_weights_data, collection_size,
            src_fields, tgt_field);

          if (!is_src)
            for (size_t j = 0; j < tgt_size; ++j)
              if (!(((ref_tgt_data[j] == FALLBACK_VALUE) == (tgt_data[j] == FALLBACK_VALUE)) ||
                    ((ref_tgt_data[j] == FIXED_VALUE) == (tgt_data[j] == FIXED_VALUE)) ||
                    ((ref_tgt_data[j] != FALLBACK_VALUE) && (ref_tgt_data[j] != FIXED_VALUE) &&
                     (fabs(
                       (ref_tgt_data[j] * scale_factors[scale_factor_idx] +
                        scale_summands[scale_summand_idx]) - tgt_data[j]) < TOL))))
                err_count++;

          yac_interp_weights_data_free(interp_weights_data);
          yac_interpolation_exchange_delete(
            interpolation_exchange, "direct_mf");
        }
      }
    }
  }

  return err_count;
}

int get_logical_attribute(int ncid, char const * attribute_name) {

  int attribute;

  size_t attlen;
  YAC_HANDLE_ERROR(nc_inq_attlen(ncid, NC_GLOBAL, attribute_name, &attlen));
  char att_text[attlen+1];
  memset(att_text, '\0', attlen+1);
  YAC_HANDLE_ERROR(nc_get_att_text(ncid, NC_GLOBAL, attribute_name, att_text));

  attribute = (strlen("TRUE") == attlen) &&
              !strncmp("TRUE", att_text, attlen);
  YAC_ASSERT_F(
    attribute ||
    ((strlen("FALSE") == attlen) &&
    !strncmp("FALSE", att_text, attlen)),
    "ERROR(get_logical_attribute): "
    " invalid global %s value %s", attribute_name, att_text);

  return attribute;
}

static void utest_get_basic_weight_file_info(
  char const * weight_file_name, int * contains_links, int * contains_fixed) {

  MPI_Barrier(MPI_COMM_WORLD);

  int ncid;

  if (!yac_file_exists(weight_file_name)) {
    PUT_ERR("error in yac_interp_weights_write_to_file");
    *contains_fixed = 0;
    *contains_links = 0;
    return;
  }

  yac_nc_open(weight_file_name, NC_NOWRITE, &ncid);
  *contains_links = get_logical_attribute(ncid, "contains_links");
  *contains_fixed = get_logical_attribute(ncid, "contains_fixed_dst");
  YAC_HANDLE_ERROR(nc_close(ncid));
}

static void on_existing_abort_handler(
  MPI_Comm comm, char const * msg, char const * source, int line) {

  UNUSED(comm);
  UNUSED(msg);
  UNUSED(source);
  UNUSED(line);

  int abort_handler_was_called = 1;
  MPI_Allreduce(
    MPI_IN_PLACE, &abort_handler_was_called, 1,
    MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  int comm_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
  MPI_Barrier(MPI_COMM_WORLD);
  if (comm_rank == 0) unlink(weight_file_name_on_existing);

  // MPI_Abort may yield non-zero error codes of mpirun hence we
  // terminate the programm gracefully
  xt_finalize();
  MPI_Finalize();
  exit(TEST_EXIT_CODE);
}
