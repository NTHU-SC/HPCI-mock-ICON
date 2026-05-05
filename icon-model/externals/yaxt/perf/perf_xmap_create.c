/**
 * @file perf_xmap_create.c
 *
 * Demonstrate performance impact of various choices in xmap creation.
 *
 * @copyright Copyright  (C)  2024 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Thomas Jahns <jahns@dkrz.de>
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

#include <assert.h>
#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <mpi.h>
#include "yaxt.h"

static void
parse_options(int argc, const char *argv[], Xt_config config);

enum {
  INDEX_GEN_COL_MAJOR,
  INDEX_GEN_ROW_MAJOR,
};

enum {
  INDEX_GEN_IDXVEC,
  INDEX_GEN_IDXSTRIPES,
  INDEX_GEN_IDXSECTION,
};

static struct {
  Xt_xmdd_bucket_gen bg;
  long num_columns;
  long nlev;
  int comm_world_size, comm_world_rank, local_size;
  int index_generation_sequence,
    index_generation_method,
    bucket_generation_sequence,
    bucket_generation_method;
  bool is_src;
} run_config = {
  .bg = NULL,
  .num_columns = 1240000, .nlev = 50,
  .index_generation_sequence = INDEX_GEN_COL_MAJOR,
  .index_generation_method = INDEX_GEN_IDXVEC,
  .bucket_generation_sequence = -1,
  .bucket_generation_method = -1,
};


static Xt_idxlist
generate_3d_list(int local_rank, int local_size, void **buf, size_t *buf_size,
                 int index_generation_method, int index_generation_sequence);

int main(int argc, const char *argv[]) {

  MPI_Init(NULL, NULL);
  xt_initialize(MPI_COMM_WORLD);
  MPI_Comm_rank(MPI_COMM_WORLD, &run_config.comm_world_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &run_config.comm_world_size);

  Xt_config conf = xt_config_new();

  parse_options(argc, argv, conf);

  void *buf = NULL;
  size_t buf_size = 0;
  Xt_idxlist idxlist_3d = generate_3d_list(
    run_config.comm_world_rank % run_config.local_size, run_config.local_size, &buf, &buf_size,
    run_config.index_generation_method, run_config.index_generation_sequence),
    idxlist_empty = xt_idxempty_new();
  free(buf);

  double tic = MPI_Wtime();

  Xt_xmap xmap =
    xt_xmap_dist_dir_custom_new(
      run_config.is_src?idxlist_3d:idxlist_empty,
      run_config.is_src?idxlist_empty:idxlist_3d,
      MPI_COMM_WORLD, conf);

  double dt = MPI_Wtime() - tic;

  Xt_xmap reo = xt_xmap_reorder_custom(xmap, XT_REORDER_RECV_UP, conf);

  xt_xmap_delete(reo);
  xt_xmap_delete(xmap);
  xt_idxlist_delete(idxlist_empty);
  xt_idxlist_delete(idxlist_3d);

  double min, max, sum;
  MPI_Reduce(&dt, &min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce(&dt, &max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce(&dt, &sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if (run_config.comm_world_rank == 0)
    fprintf(stderr, "min %e max %e avg %e\n", min, max,
            sum / (double)run_config.comm_world_size);

  xt_config_delete(conf);
  if (run_config.bg)
    xt_xmdd_bucket_gen_delete(run_config.bg);

  xt_finalize();
  MPI_Finalize();

  return EXIT_SUCCESS;
}

static inline Xt_int
start_of_rank(int rank, int num_ranks, Xt_int num_columns)
{
  return (Xt_int)(((long long)num_columns * (long long)rank) / num_ranks);
}

static inline int
rank_of_column(Xt_int column, int num_ranks, Xt_int num_columns)
{
  return (int)(((long long)column * num_ranks)/num_columns);
}

static Xt_idxlist
generate_3d_list(int local_rank, int local_size_, void **buf, size_t *buf_size,
                 int index_generation_method, int index_generation_sequence)
{
  Xt_int start_idx = start_of_rank(local_rank, local_size_,
                                   (Xt_int)run_config.num_columns);

  Xt_int next_start_idx = start_of_rank(local_rank+1, local_size_,
                                        (Xt_int)run_config.num_columns);

  Xt_int local_num_columns = (Xt_int)(next_start_idx - start_idx);

  Xt_idxlist idxlist_3d;
  if (index_generation_method == INDEX_GEN_IDXVEC) {
    Xt_int *idx3D = *buf;
    size_t needed_buf_size
      = (size_t)(local_num_columns * run_config.nlev) * sizeof(Xt_int);
    if (*buf_size < needed_buf_size) {
      *buf_size = needed_buf_size;
      *buf = idx3D = realloc(idx3D, needed_buf_size);
    }

    if (index_generation_sequence == INDEX_GEN_ROW_MAJOR) {
      for (Xt_int lev = 0, j = 0; lev < run_config.nlev; ++lev)
        for (Xt_int idx = 0; idx < local_num_columns; ++idx, ++j)
          idx3D[j] = (Xt_int)(start_idx + idx + lev * run_config.num_columns);
    }
    else if (index_generation_sequence == INDEX_GEN_COL_MAJOR) {
      for (Xt_int idx = 0, j = 0; idx < local_num_columns; ++idx)
        for (Xt_int lev = 0; lev < run_config.nlev; ++lev, ++j)
          idx3D[j] = (Xt_int)(start_idx + idx + lev * run_config.num_columns);
    }
    idxlist_3d
      = xt_idxvec_new(idx3D, (int)(run_config.nlev * local_num_columns));
  } else if (index_generation_method == INDEX_GEN_IDXSTRIPES) {
    size_t num_stripes
      = index_generation_sequence == INDEX_GEN_ROW_MAJOR
      ? (size_t)run_config.nlev
      : (size_t)local_num_columns;
    struct Xt_stripe *idx3D = *buf;
    size_t needed_buf_size = num_stripes * sizeof(*idx3D);
    if (*buf_size < needed_buf_size) {
      *buf_size = needed_buf_size;
      *buf = idx3D = realloc(idx3D, needed_buf_size);
    }

    if (index_generation_sequence == INDEX_GEN_ROW_MAJOR) {
      for (Xt_int lev = 0; lev < run_config.nlev; ++lev)
        idx3D[lev] = (struct Xt_stripe){
          .start = (Xt_int)(start_idx + lev * run_config.num_columns),
          .stride = 1,
          .nstrides = (int)local_num_columns
        };
    }
    else if (index_generation_sequence == INDEX_GEN_COL_MAJOR) {
      for (Xt_int idx = 0; idx < local_num_columns; ++idx)
        idx3D[idx] = (struct Xt_stripe){
          .start = (Xt_int)(start_idx + idx),
          .stride = (Xt_int)run_config.num_columns,
          .nstrides = (int)run_config.nlev
        };
    }
    idxlist_3d
      = xt_idxstripes_new(idx3D, (int)num_stripes);
  } else /* if (index_generation_method == INDEX_GEN_IDXSECTION) */ {
    if (index_generation_sequence == INDEX_GEN_ROW_MAJOR)
      idxlist_3d = xt_idxsection_new(
        0, 2,
        (Xt_int []){  (Xt_int)run_config.nlev, (Xt_int)run_config.num_columns },
        (int []){  (int)run_config.nlev, (int)local_num_columns },
        (Xt_int []){  (Xt_int)0, start_idx });
    else {
      abort();
      /* currently unsupported */
    }
  }
  return idxlist_3d;
}


static long
s2long(const char *sv, const char *context1, const char *context2)
{
  char *endptr;
  char msg[256+strlen(context1)+strlen(context2)];
  errno = 0;
  long v = strtol(sv, &endptr, 0);
  if ((errno == ERANGE && (v == LONG_MAX || v == LONG_MIN))
      || (errno != 0 && v == 0)) {
    sprintf(msg, "failed to parse argument value ('%s') of %s%s",
            sv, context1, context2);
    perror(msg);
    exit(1);
  } else if (endptr == sv) {
    fprintf(stderr, "malformed value ('%s') of %s%s"
            ", no digits were found\n", sv, context1, context2);
    exit(1);
  }
  /*
   * } else if (v < 1 || v > INT_MAX) {
   *   fprintf(stderr, "value of XT_CONFIG_DEFAULT_IDXVEC_AUTOCONVERT_SIZE"
   *             " environment variable (%ld) out of range [1,%d]\n",
   *             v, INT_MAX);
   *   } else
   */
  return v;
}

extern int
xt_sort_algo_id_by_name(const char *name);

static void
set_custom_bucket_gen(Xt_config config);

static void
set_list_gen_methods(int *generation_sequence,
                     int *generation_method,
                     const char *arg)
{
  size_t cur_parm_ofs = 0;
  /* split arg into comma-separated parts */
  while (arg[cur_parm_ofs]) {
    size_t cur_parm_len = strcspn(arg+cur_parm_ofs, ",");
    if (!strncmp(arg+cur_parm_ofs, "col-major", cur_parm_len))
      *generation_sequence=INDEX_GEN_COL_MAJOR;
    else if (!strncmp(arg+cur_parm_ofs, "row-major", cur_parm_len))
      *generation_sequence=INDEX_GEN_ROW_MAJOR;
    else if (!strncmp(arg+cur_parm_ofs, "idxvec", cur_parm_len))
      *generation_method=INDEX_GEN_IDXVEC;
    else if (!strncmp(arg+cur_parm_ofs, "idxstripes", cur_parm_len))
      *generation_method=INDEX_GEN_IDXSTRIPES;
    else if (!strncmp(arg+cur_parm_ofs, "idxsection", cur_parm_len))
      *generation_method=INDEX_GEN_IDXSECTION;
    else {
      fprintf(stderr, "error: unknown index generation method name: %s\n",
              arg);
      exit(1);
    }
    cur_parm_ofs += cur_parm_len
      + (arg[cur_parm_ofs+cur_parm_len] == ',');
  }
}

static void
parse_options(int argc, const char *argv[], Xt_config config)
{
  static const struct option longopts[] = {
    { .name = "index-generation", .has_arg = 1, NULL, 0 }, /* 0 */
    { .name = "num-columns", .has_arg = 1, NULL, 0 },      /* 1 */
    { .name = "nlev", .has_arg = 1, NULL, 0 },             /* 2 */
    { .name = "num-levels", .has_arg = 1, NULL, 0 },       /* 3 */
    { .name = "sort-algorithm", .has_arg = 1, NULL, 0 },   /* 4 */
    { .name = "enable-mem-saving", .has_arg = 0, NULL, 0 }, /* 5 */
    { .name = "disable-mem-saving", .has_arg = 0, NULL, 0 }, /* 6 */
    { .name = "enable-custom-bucket-generator", .has_arg = 2, NULL, 0 }, /* 7 */
    { .name = "stripe-alignment-mode", .has_arg = 1, NULL, 0 }, /* 8 */
    { .name = "disable-stripe-alignment", .has_arg = 0, NULL, 0 }, /* 9 */
    { .name = "enable-stripe-alignment", .has_arg = 0, NULL, 0 }, /* 10 */
    { 0, 0, 0, 0 },
  };
  long *opt_var[]
    = { [1] = &run_config.num_columns, [2] = &run_config.nlev, [3] = &run_config.nlev, };
  while (1) {
    int option_index = 0;
    int c = getopt_long(argc, (char * const*)(intptr_t)argv, "",
                        longopts, &option_index);
    if (c == -1)
      break;

    long align;
    switch (c) {
    case 0:
      switch (option_index) {
      case 0:
        set_list_gen_methods(&run_config.index_generation_sequence,
                             &run_config.index_generation_method,
                             optarg);
        break;
      case 1:
      case 2:
      case 3:
        *opt_var[option_index] = s2long(optarg, "command-line option --",
                                        longopts[option_index].name);
        break;
      case 4:
        {
          int algo
            = xt_sort_algo_id_by_name(optarg);
          if (algo == -1) {
            fprintf(stderr, "error: invalid sort algorithm name: %s\n",
                    optarg);
            exit(1);
          }
          xt_config_set_sort_algorithm_by_id(config, algo);
        }
        break;
      case 5:
      case 6:
        xt_config_set_mem_saving(config, option_index < 6);
        break;
      case 7:
        set_custom_bucket_gen(config);
        if (optarg)
          set_list_gen_methods(&run_config.bucket_generation_sequence,
                               &run_config.bucket_generation_method,
                               optarg);

        break;
      case 8:
        align = s2long(optarg, "command-line option --",
                       longopts[option_index].name);
        goto set_dd_stripe_align;
      case 9:
        align = 0;
        goto set_dd_stripe_align;
      case 10:
        align = 1;
      set_dd_stripe_align:
        xt_config_set_xmap_stripe_align(config, (int)align);
        break;
      }
    }
  }
  int size_hamocc;
  if (optind < argc) {
    size_hamocc = (int)s2long(argv[optind], "hamocc_size", "");
  } else {
    size_hamocc = run_config.comm_world_size/2;
  }
  run_config.is_src = run_config.comm_world_rank < (run_config.comm_world_size - size_hamocc);
  run_config.local_size = run_config.is_src
    ? run_config.comm_world_size - size_hamocc : size_hamocc;
  if (run_config.bucket_generation_sequence == -1)
    run_config.bucket_generation_sequence = run_config.index_generation_sequence;
  if (run_config.bucket_generation_method == -1)
    run_config.bucket_generation_method = run_config.index_generation_sequence;
}

struct custom_bucket_gen_state
{
  struct Xt_xmdd_bucket_gen_comms comms;
  int prev_rank_generated, last_rank_generated, intra_comm_size;
  size_t stripe_buf_size;
  void *stripe_buf;
  Xt_idxlist prev_list;
};

static int
custom_bucket_gen_init_state(
  void *gen_state,
  Xt_idxlist src_idxlist,
  Xt_idxlist dst_idxlist,
  Xt_config config,
  const struct Xt_xmdd_bucket_gen_comms *comms,
  void *init_params);

static void custom_bucket_gen_destroy_state(void *gen_state);

static int custom_bucket_gen_get_intersect_max_num(
  void *gen_state, int type);

static struct Xt_com_list custom_bucket_gen_next(
  void *gen_state, int type);

static void
set_custom_bucket_gen(Xt_config config)
{
  run_config.bg = xt_xmdd_bucket_gen_new();
  xt_xmdd_bucket_gen_define_interface(
    run_config.bg,
    custom_bucket_gen_init_state,
    custom_bucket_gen_destroy_state,
    custom_bucket_gen_get_intersect_max_num,
    custom_bucket_gen_next,
    sizeof (struct custom_bucket_gen_state),
    NULL);
  xt_config_set_xmdd_bucket_gen(config, run_config.bg);
}

static int
custom_bucket_gen_init_state(
  void *gen_state_,
  Xt_idxlist src_idxlist,
  Xt_idxlist dst_idxlist,
  Xt_config config,
  const struct Xt_xmdd_bucket_gen_comms *comms,
  void *init_params)
{
  (void)init_params;
  struct custom_bucket_gen_state *gen_state = gen_state_;
  memcpy(&gen_state->comms, comms, sizeof (gen_state->comms));
  MPI_Comm_size(comms->intra_comm, &gen_state->intra_comm_size);
  gen_state->stripe_buf_size = 0;
  gen_state->stripe_buf = NULL;
  gen_state->prev_list = NULL;
  int p = run_config.index_generation_method == INDEX_GEN_IDXSTRIPES
    || xt_idxlist_is_stripe_conversion_profitable(src_idxlist, config)
    || xt_idxlist_is_stripe_conversion_profitable(dst_idxlist, config);
  MPI_Allreduce(MPI_IN_PLACE, &p, 1, MPI_INT, MPI_LOR,
                comms->intra_comm);
  /* establish local ranges */
  int num_indices_src = xt_idxlist_get_num_indices(src_idxlist),
    num_indices_dst = xt_idxlist_get_num_indices(dst_idxlist);
  Xt_int first_column = XT_INT_MAX, last_column = XT_INT_MIN;
  if (num_indices_src) {
    Xt_int start_src = xt_idxlist_get_min_index(src_idxlist),
      end_src = (Xt_int)(xt_idxlist_get_max_index(src_idxlist)
                         - (run_config.nlev-1) * run_config.num_columns);
    if (start_src < first_column) first_column = start_src;
    if (end_src > last_column) last_column = end_src;
  }
  if (num_indices_dst) {
    Xt_int start_dst = xt_idxlist_get_min_index(dst_idxlist),
      end_dst = (Xt_int)(xt_idxlist_get_max_index(dst_idxlist)
                         - (run_config.nlev-1) * run_config.num_columns);
    if (start_dst < first_column) first_column = start_dst;
    if (end_dst > last_column) last_column = end_dst;
  }
  int start_rank = rank_of_column(first_column, gen_state->intra_comm_size,
                                  (Xt_int)run_config.num_columns),
    end_rank = rank_of_column(last_column, gen_state->intra_comm_size,
                              (Xt_int)run_config.num_columns);
  gen_state->prev_rank_generated = start_rank-1;
  gen_state->last_rank_generated = end_rank;
  return p;
}


static void custom_bucket_gen_destroy_state(void *gen_state_)
{
  struct custom_bucket_gen_state *gen_state = gen_state_;
  if (gen_state->prev_list)
    xt_idxlist_delete(gen_state->prev_list);
  free(gen_state->stripe_buf);
}


static int custom_bucket_gen_get_intersect_max_num(
  void *gen_state_, int type)
{
  struct custom_bucket_gen_state *gen_state = gen_state_;
  (void)type;
  return gen_state->last_rank_generated-gen_state->prev_rank_generated;
}


static struct Xt_com_list custom_bucket_gen_next(
  void *gen_state_, int type)
{
  struct custom_bucket_gen_state *gen_state = gen_state_;
  if (gen_state->prev_list)
    xt_idxlist_delete(gen_state->prev_list);
  int next_rank = ++gen_state->prev_rank_generated;
  Xt_idxlist bucket = NULL;
  if (next_rank <= gen_state->last_rank_generated)
    bucket = generate_3d_list(next_rank, gen_state->intra_comm_size,
                              &gen_state->stripe_buf,
                              &gen_state->stripe_buf_size,
                              run_config.bucket_generation_method,
                              run_config.bucket_generation_sequence);
  (void)type;
  gen_state->prev_list = bucket;
  return (struct Xt_com_list){ .list = bucket, .rank = next_rank };
}

/*
 * Local Variables:
 * coding: utf-8
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * license-project-url: "https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/"
 * license-default: "bsd"
 * End:
 */
