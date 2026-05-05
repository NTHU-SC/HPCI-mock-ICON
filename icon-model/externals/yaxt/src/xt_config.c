/**
 * @file xt_config.c
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
#include <errno.h>
#include <stdint.h>
#include <string.h>

#include <mpi.h>

#include <xt/xt_config.h>
#include <xt/xt_mpi.h>
#include "xt_config_internal.h"
#include "xt_exchanger_irecv_send.h"
#include "xt_exchanger_irecv_isend.h"
#include "xt_exchanger_mix_isend_irecv.h"
#include "xt_exchanger_irecv_isend_packed.h"
#include "xt_exchanger_irecv_isend_ddt_packed.h"
#include "xt_exchanger_neigh_alltoall.h"
#include "xt_idxlist_internal.h"
#include "xt/quicksort.h"
#include "xt/mergesort.h"
#include "xt/xt_xmap_dist_dir_bucket_gen.h"
#include "xt_xmap_dist_dir_bucket_gen_cycl_stripe.h"
#include "xt_xmap_dist_dir_bucket_gen_internal.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"

static const char filename[] = "xt_config.c";

Xt_config xt_config_new(void)
{
  Xt_config config = xmalloc(sizeof(*config));
  *config = xt_default_config;
  return config;
}

void xt_config_delete(Xt_config config)
{
  if (config->xmdd_bucket_gen != &Xt_xmdd_cycl_stripe_bucket_gen_desc)
    free((void *)(intptr_t)config->xmdd_bucket_gen);
  free(config);
}

static const struct {
  char name[32];
  Xt_exchanger_new f;
  int code;
} exchanger_table[] = {
  { "irecv_send",
    xt_exchanger_irecv_send_new, xt_exchanger_irecv_send },
  { "irecv_isend",
    xt_exchanger_irecv_isend_new, xt_exchanger_irecv_isend },
  { "irecv_isend_packed",
    xt_exchanger_irecv_isend_packed_new, xt_exchanger_irecv_isend_packed },
  { "irecv_isend_ddt_packed",
#ifdef XT_ENABLE_DDT_EXCHANGER
    xt_exchanger_irecv_isend_ddt_packed_new,
#else
    (Xt_exchanger_new)0,
#endif
    xt_exchanger_irecv_isend_ddt_packed },
  { "mix_irecv_isend",
    xt_exchanger_mix_isend_irecv_new, xt_exchanger_mix_isend_irecv },
  { "neigh_alltoall",
#ifdef XT_CAN_USE_MPI_NEIGHBOR_ALLTOALL
    xt_exchanger_neigh_alltoall_new,
#else
    (Xt_exchanger_new)0,
#endif
    xt_exchanger_neigh_alltoall },
};

enum {
  num_exchanger = sizeof (exchanger_table) / sizeof (exchanger_table[0]),
};

int
xt_exchanger_id_by_name(const char *name)
{
  for (size_t i = 0; i < num_exchanger; ++i)
    if (!strcmp(name, exchanger_table[i].name))
      return exchanger_table[i].code;
  return -1;
}

static inline size_t
exchanger_by_function(Xt_exchanger_new exchanger_new)
{
  for (size_t i = 0; i < num_exchanger; ++i)
    if (exchanger_table[i].f == exchanger_new)
      return i;
  return SIZE_MAX;
}


int xt_config_get_exchange_method(Xt_config config)
{
  Xt_exchanger_new exchanger_new = config->exchanger_new;
  size_t eentry = exchanger_by_function(exchanger_new);
  if (eentry != SIZE_MAX)
    return exchanger_table[eentry].code;
  static const char fmt[]
    = "error: unexpected exchanger function (%p)!";
  char buf[sizeof (fmt) + 3*sizeof(void *)];
  sprintf(buf, fmt, (void *)exchanger_new);
  Xt_abort(Xt_default_comm, buf, filename, __LINE__);
}

Xt_exchanger_new
xt_config_get_exchange_new_by_comm(Xt_config config, MPI_Comm comm)
{
  Xt_exchanger_new exchanger_new = config->exchanger_new;
#ifdef XT_CAN_USE_MPI_NEIGHBOR_ALLTOALL
  if (exchanger_new == xt_exchanger_neigh_alltoall_new) {
    int flag;
    xt_mpi_call(MPI_Comm_test_inter(comm, &flag), comm);
    if (flag)
      exchanger_new = xt_exchanger_mix_isend_irecv_new;
  }
#else
  (void)comm;
#endif
  return exchanger_new;
}

static const struct {
  char name[16];
  struct Xt_sort_algo_funcptr func;
} sort_algo_table[] = {
  { .name = "quicksort",
    .func = {
      .sort_int = xt_quicksort_int,
      .sort_xt_int = xt_quicksort_xt_int,
      .sort_index = xt_quicksort_index,
      .sort_idxpos = xt_quicksort_idxpos,
      .sort_xt_int_permutation = xt_quicksort_xt_int_permutation,
      .sort_int_permutation = xt_quicksort_int_permutation
    },
  },
  { .name = "mergesort",
    .func = {
      .sort_int = xt_mergesort_int,
      .sort_xt_int = xt_mergesort_xt_int,
      .sort_index = xt_mergesort_index,
      .sort_idxpos = xt_mergesort_idxpos,
      .sort_xt_int_permutation = xt_mergesort_xt_int_permutation,
      .sort_int_permutation = xt_mergesort_int_permutation
    },
  },
};

enum {
  num_sort_algo = sizeof (sort_algo_table) / sizeof (sort_algo_table[0]),
};

struct Xt_config_ xt_default_config = {
  .exchanger_new = xt_exchanger_mix_isend_irecv_new,
  .exchanger_team_share = NULL,
  .idxv_cnv_size = CHEAP_VECTOR_SIZE,
  .sort_funcs = &sort_algo_table[XT_QUICKSORT].func,
  .xmdd_bucket_gen = &Xt_xmdd_cycl_stripe_bucket_gen_desc,
  .flags = 2U << xt_force_xmap_striping_bit_ofs,
};


int
xt_sort_algo_id_by_name(const char *name)
{
  for (size_t i = 0; i < num_sort_algo; ++i)
    if (!strcmp(name, sort_algo_table[i].name))
      return (int)i;
  return -1;
}

static inline size_t
sort_algo_by_table(const struct Xt_sort_algo_funcptr *sort_funcs)
{
  for (size_t i = 0; i < num_sort_algo; ++i)
    if (&sort_algo_table[i].func == sort_funcs)
      return i;
  return SIZE_MAX;
}

int xt_config_get_sort_algorithm_id(Xt_config config)
{
  const struct Xt_sort_algo_funcptr *sort_funcs = config->sort_funcs;
  size_t eentry = sort_algo_by_table(sort_funcs);
  if (eentry != SIZE_MAX)
    return (int)eentry;
  static const char fmt[]
    = "error: unexpected exchanger function (%p)!";
  char buf[sizeof (fmt) + 3*sizeof(void *)];
  sprintf(buf, fmt, (const void *)sort_funcs);
  Xt_abort(Xt_default_comm, buf, filename, __LINE__);
}

void xt_config_set_sort_algorithm_by_id(Xt_config config, int algo)
{
  if (algo >= 0 && algo < num_sort_algo) {
    config->sort_funcs = &sort_algo_table[algo].func;
    return;
  }
  static const char fmt[]
    = "error: user-requested exchanger code (%d) does not exist!";
  char buf[sizeof (fmt) + 3*sizeof(int)];
  sprintf(buf, fmt, algo);
  Xt_abort(Xt_default_comm, buf, filename, __LINE__);
}

enum {
  XT_MAX_MEM_SAVING = 1,
};

void
xt_config_set_mem_saving(Xt_config config, int memconserve)
{
  XT_CONFIG_SET_FORCE_NOSORT_BIT(config, memconserve);
}

int
xt_config_get_mem_saving(Xt_config config)
{
  return XT_CONFIG_GET_FORCE_NOSORT(config);
}


Xt_xmdd_bucket_gen
xt_config_get_xmdd_bucket_gen(Xt_config config)
{
  return (Xt_xmdd_bucket_gen)(intptr_t)config->xmdd_bucket_gen;
}

void
xt_config_set_xmdd_bucket_gen(Xt_config config,
                              Xt_xmdd_bucket_gen bucket_gen_iface)
{
  if (bucket_gen_iface == &Xt_xmdd_cycl_stripe_bucket_gen_desc)
    config->xmdd_bucket_gen = bucket_gen_iface;
  else {
    Xt_xmdd_bucket_gen gen = xmalloc(sizeof (*bucket_gen_iface));
    config->xmdd_bucket_gen = gen;
    *gen = *bucket_gen_iface;
  }
}




void xt_config_set_exchange_method(Xt_config config, int method)
{
  static const char fmt[]
    = "error: user-requested exchanger code (%d) does not exist!";
  char buf[sizeof (fmt) + 3*sizeof(int)];
  const char *msg = buf;
  for (size_t i = 0; i < num_exchanger; ++i)
    if (exchanger_table[i].code == method) {
      Xt_exchanger_new exchanger_new;
      if (exchanger_table[i].f) {
        exchanger_new = exchanger_table[i].f;
      } else {
        exchanger_new = xt_default_config.exchanger_new;
        size_t default_entry = exchanger_by_function(exchanger_new);
        if (default_entry == SIZE_MAX) {
          msg = "error: invalid default exchanger constructor!";
          goto abort;
        }
        fprintf(stderr, "warning: %s exchanger unavailable, using "
                "%s instead\n",
                exchanger_table[i].name, exchanger_table[default_entry].name);
      }
      config->exchanger_new = exchanger_new;
      return;
    }
  sprintf(buf, fmt, method);
abort:
  Xt_abort(Xt_default_comm, msg, filename, __LINE__);
}

int xt_config_get_idxvec_autoconvert_size(Xt_config config)
{
  return config->idxv_cnv_size;
}

void
xt_config_set_idxvec_autoconvert_size(Xt_config config, int cnvsize)
{
  if (cnvsize > 3)
    config->idxv_cnv_size = cnvsize;
}

int
xt_config_get_redist_mthread_mode(Xt_config config)
{
  return (int)((config->flags & (uint32_t)xt_mthread_mode_mask)
               >> xt_mthread_mode_bit_ofs);
}

void
xt_config_set_redist_mthread_mode(Xt_config config, int mode)
{
  assert(mode >= XT_MT_NONE && mode <= XT_MT_OPENMP);
#ifndef _OPENMP
  if (mode == XT_MT_OPENMP)
    Xt_abort(Xt_default_comm,
             "error: automatic opening of OpenMP parallel regions requested,"
             " but OpenMP is not configured.\n", filename, __LINE__);
#else
  if (mode == XT_MT_OPENMP) {
    int thread_support_provided = MPI_THREAD_SINGLE;
    xt_mpi_call(MPI_Query_thread(&thread_support_provided), Xt_default_comm);
    if (thread_support_provided != MPI_THREAD_MULTIPLE)
      Xt_abort(Xt_default_comm,
               "error: automatic opening of OpenMP parallel regions requested,"
               "\n       but MPI is not running in thread-safe mode.\n",
               filename, __LINE__);
  }
#endif
  config->flags = (config->flags & ~(uint32_t)xt_mthread_mode_mask)
    | ((uint32_t)mode << xt_mthread_mode_bit_ofs);
}

void
xt_config_set_dist_dir_stripe_alignment(Xt_config config, int use_stripe_alignment);

void
xt_config_set_dist_dir_stripe_alignment(Xt_config config, int use_stripe_alignment)
{
  xt_config_set_xmap_stripe_align(config, use_stripe_alignment);
}

int
xt_config_get_dist_dir_stripe_alignment(Xt_config config);

int
xt_config_get_dist_dir_stripe_alignment(Xt_config config)
{
  return XT_CONFIG_GET_XMAP_STRIPING(config);
}

void
xt_config_set_xmap_stripe_align(Xt_config config, int use_stripe_align)
{
  if (use_stripe_align < 0 || use_stripe_align > 2)
    Xt_abort(Xt_default_comm,
             "error: invalid value passed to "
             "xt_config_set_xmap_stripe_align.\n", filename, __LINE__);
  XT_CONFIG_SET_XMAP_STRIPING(config, use_stripe_align);
}

int
xt_config_get_xmap_stripe_align(Xt_config config)
{
  return XT_CONFIG_GET_XMAP_STRIPING(config);
}

/**
 * \page defaults Customizable defaults
 *
 * \tableofcontents
 * YAXT relies on some defaults that usually provide an adequate
 * compromise of space and time efficiency. In cases when these fail
 * to provide, one can either override the global default with an
 * environment variable or use a custom constructor method that takes
 * an extra argument of type \ref Xt_config for exactly the place where
 * the usual default fails.
 *
 * \section default_exchanger Exchanger choice
 *
 * The internal exchanger class handles the message passing part of the
 * redist.
 *
 * The default exchanger, mix_irecv_isend, is suited for MPI
 * implementations with strong support for derived MPI datatypes and
 * efficient handling of nonblocking point-to-point communication.
 *
 * By setting the environment variable XT_CONFIG_DEFAULT_EXCHANGE_METHOD,
 * one of the other available exchangers can be selected:
 * <ul>
 *
 * <li> irecv_send uses non-blocking receives but then uses blocking
 * sends. This may be beneficial for platforms that can elide some
 * preparations for non-blocking transfers in this case but is
 * expected to have downsides as soon as actual network latencies
 * become relevant.
 *
 * <li> irecv_isend differs from the default in that it doesn't bother
 * to mix initiating sends and receives. Rather, all receives are
 * initiated before any send is initiated.
 *
 * <li> irecv_isend_packed Separates the network transfer from the
 * iteration of the MPI datatype by creating an MPI_PACKED buffer
 * filled via MPI_Pack internally.
 *
 * <li> irecv_isend_ddt_packed Builds on irecv_isend_packed but uses
 * OpenACC kernels instead of MPI_Pack. This usually results in
 * significant performance gains on GPUs.
 *
 * <li> mix_irecv_isend
 *
 * <li> neigh_alltoall Available when a robust implementation of MPI 3
 * neighbor collectives is provided. Uses MPI_Neighbor_alltoallw
 * instead of point-to-point communication. This may benefit from
 * pre-created data paths but creates additional MPI communicators
 * which may be too costly in highly dynamic use cases.
 *
 * </ul>
 *
 * \section default_autoconvert Automatic conversions of index vectors
 * to stripes
 *
 * When constructing an Xmap via \ref xt_xmap_dist_dir_new or
 * \ref xt_xmap_all2all_new, index lists passed in as \ref xt_idxvec will be
 * automatically converted to \ref xt_idxstripes if the size is above a
 * limit which defaults to 128 indices.
 *
 * The XT_CONFIG_DEFAULT_IDXVEC_AUTOCONVERT_SIZE environment variable
 * can be used to override this value.
 *
 * \section mthread_mode Internal multi-threading
 *
 * For an MPI that supports multiple threads calling into it (see
 * MPI_Init_thread and MPI_THREAD_MULTIPLE), it can be more efficient
 * to use all multiple send and/or receive operations in parallel.
 *
 * Set the XT_CONFIG_DEFAULT_MULTI_THREAD_MODE environment variable to
 * "XT_MT_OPENMP" to enable this globally.
 * This is currently equivalent to adding an \ref Xt_config parameter which had
 * \ref xt_config_set_redist_mthread_mode with parameter \ref XT_MT_OPENMP
 * called on it.
 *
 * \section sort_algorithm Sort algorithm choice
 *
 * The default sort algorithm in YAXT is currenlty quicksort because
 * it does not normally need more than @f$O(\log N)@f$ memory. It is used
 * to speed up the computation of index list intersections when the
 * inputs are index vectors. Since almost sorted inputs lead to
 * exaggerated run-time resuls, it can be beneficial to use the
 * Mergesort algorithm on such index lists. Call \ref
 * xt_config_set_sort_algorithm_by_id on an object of class \ref Xt_config
 * and use it in calls to xmap constructors to change the sort algorithm.
 *
 * By setting the environment variable
 * XT_CONFIG_DEFAULT_SORT_ALGORITHM to "mergesort" (case-insensitive)
 * this can be made the global default.
 *
 * \section stripe_align Stripe alignment
 *
 * When computing the positions corresponding to an intersection
 * during the construction of an Xmap, the usual strategy for index
 * lists above a certain size (the limit described at \ref
 * default_autoconvert ) is to map multiple positions at the same time
 * by describing them as position extents (\ref Xt_pos_ext).
 *
 * In case only very few or even only one adjacent positions are ever
 * mapped together this is very inefficient. In such situations the
 * alternative strategy of mapping each index individually is faster.
 * Mapping individually can be set with \ref
 * xt_config_set_xmap_stripe_align and an object of type
 * Xt_config. Also this can be requested globally be setting the
 * environment variable
 * XT_CONFIG_DEFAULT_XMAP_STRIPE_ALIGN to "one_by_one" or the value 0.
 * The default is "auto" or 2 and stripe alignment can be enforced
 * with the values "always" or 1.
 *
 * \section mem_saving Memory saving optimizations
 *
 * In some situations, the default to use extra data structures to
 * speed up e.g. searching can be counter-productive. In those cases
 * the default can be changed to prefer expensive but
 * memory-conserving algorithms.
 *
 * This can be changed globally by setting the environment variable
 * XT_CONFIG_DEFAULT_MEM_SAVING to 1.
 */
void
xt_config_defaults_init(void)
{
  const char *config_env = getenv("XT_CONFIG_DEFAULT_EXCHANGE_METHOD");
  if (config_env) {
    int exchanger_id = xt_exchanger_id_by_name(config_env);
    if (exchanger_id != -1)
      xt_config_set_exchange_method(&xt_default_config, exchanger_id);
    else
      fprintf(stderr, "warning: Unexpected value "
              "for XT_CONFIG_DEFAULT_EXCHANGE_METHOD=%s\n", config_env);
  }
  config_env = getenv("XT_CONFIG_DEFAULT_IDXVEC_AUTOCONVERT_SIZE");
  if (config_env) {
    const char *endptr;
    long v;
    errno = 0;
    if (!strcmp(config_env, "INT_MAX")) {
      v = INT_MAX;
      endptr = config_env + 7;
    } else {
      char *endptr_;
      v = strtol(config_env, &endptr_, 0);
      endptr = endptr_;
    }
    if ((errno == ERANGE && (v == LONG_MAX || v == LONG_MIN))
        || (errno != 0 && v == 0)) {
      perror("failed to parse value of "
             "XT_CONFIG_DEFAULT_IDXVEC_AUTOCONVERT_SIZE environment variable");
    } else if (endptr == config_env) {
      fputs("malformed value of XT_CONFIG_DEFAULT_IDXVEC_AUTOCONVERT_SIZE"
            " environment variable, no digits were found\n",
            stderr);
    } else if (v < 1 || v > INT_MAX) {
      fprintf(stderr, "value of XT_CONFIG_DEFAULT_IDXVEC_AUTOCONVERT_SIZE"
              " environment variable (%ld) out of range [1,%d]\n",
              v, INT_MAX);
    } else
      xt_config_set_idxvec_autoconvert_size(&xt_default_config, (int)v);
  }
  config_env = getenv("XT_CONFIG_DEFAULT_MULTI_THREAD_MODE");
  if (config_env) {
    char *endptr;
    long v = strtol(config_env, &endptr, 0);
    if (endptr != config_env) {
      if ((errno == ERANGE && (v == LONG_MAX || v == LONG_MIN))
          || (errno != 0 && v == 0)) {
        perror("failed to parse value of "
               "XT_CONFIG_DEFAULT_MULTI_THREAD_MODE environment variable");
        goto dont_set_mt_mode;
      } else if (v < XT_MT_NONE || v > XT_MT_OPENMP) {
        fprintf(stderr, "numeric value of XT_CONFIG_DEFAULT_MULTI_THREAD_MODE"
                " environment variable (%ld) out of range [0,%d]\n",
                v, XT_MT_OPENMP);
        goto dont_set_mt_mode;
      } else if (*endptr) {
        fprintf(stderr, "trailing text '%s' found after numeric value (%*s) in "
                "XT_CONFIG_DEFAULT_MULTI_THREAD_MODE environment variable\n",
                endptr, (int)(endptr-config_env), config_env);
        goto dont_set_mt_mode;
      }
    } else {
      if (!strcasecmp(config_env, "XT_MT_OPENMP")) {
#ifndef _OPENMP
        fputs("multi-threaded operation requested via "
              "XT_CONFIG_DEFAULT_MULTI_THREAD_MODE, but OpenMP support is not"
              " compiled in!\n", stderr);
        goto dont_set_mt_mode;
#else
        v = XT_MT_OPENMP;
#endif
      } else if (!strcasecmp(config_env, "XT_MT_NONE")) {
        v = XT_MT_NONE;
      } else {
        fputs("unexpected value of XT_CONFIG_DEFAULT_MULTI_THREAD_MODE"
              " environment variable, unrecognized text or numeral\n",
              stderr);
        goto dont_set_mt_mode;
      }
    }
    xt_config_set_redist_mthread_mode(&xt_default_config, (int)v);
  }
dont_set_mt_mode:;

  config_env = getenv("XT_CONFIG_DEFAULT_SORT_ALGORITHM");
  if (config_env) {
    char *endptr;
    long v = strtol(config_env, &endptr, 0);
    if (endptr != config_env) {
      if ((errno == ERANGE && (v == LONG_MAX || v == LONG_MIN))
          || (errno != 0 && v == 0)) {
        perror("failed to parse value of "
               "XT_CONFIG_DEFAULT_SORT_ALGORITHM environment variable");
        goto dont_set_sort_algorithm;
      } else if (v < 0 || v > XT_MERGESORT) {
        fprintf(stderr, "numeric value of XT_CONFIG_DEFAULT_SORT_ALGORITHM"
                " environment variable (%ld) out of range [0,%d]\n",
                v, XT_MERGESORT);
        goto dont_set_sort_algorithm;
      } else if (*endptr) {
        fprintf(stderr, "trailing text '%s' found after numeric value (%*s) in "
                "XT_CONFIG_DEFAULT_SORT_ALGORITHM environment variable\n",
                endptr, (int)(endptr-config_env), config_env);
        goto dont_set_sort_algorithm;
      }
    } else {
      if (!strcasecmp(config_env, "QUICKSORT")) {
        v = XT_QUICKSORT;
      } else if (!strcasecmp(config_env, "MERGESORT")) {
        v = XT_MERGESORT;
      } else {
        fputs("unexpected value of XT_CONFIG_DEFAULT_SORT_ALGORITHM"
              " environment variable, unrecognized text or numeral\n",
              stderr);
        goto dont_set_sort_algorithm;
      }
    }
    xt_config_set_sort_algorithm_by_id(&xt_default_config, (int)v);
  }
dont_set_sort_algorithm:;
  {
    const char *evn = "XT_CONFIG_DEFAULT_XMAP_STRIPE_ALIGN";
    config_env = getenv(evn);
    if (!config_env)
      config_env = getenv((evn = "XT_CONFIG_DEFAULT_DIST_DIR_STRIPE_ALIGNMENT"));
    if (config_env) {
      const char *endptr;
      long v;
      errno = 0;
      if (!strcmp(config_env, "auto")) {
        v = 2;
        endptr = config_env + 4;
      } else if (!strcmp(config_env, "one_by_one")) {
        v = 0;
        endptr = config_env + 10;
      } else if (!strcmp(config_env, "always")) {
        v = 1;
        endptr = config_env + 6;
      } else {
        char *endptr_;
        v = strtol(config_env, &endptr_, 0);
        endptr = endptr_;
      }
      if ((errno == ERANGE && (v == LONG_MAX || v == LONG_MIN))
          || (errno != 0 && v == 0)) {
        fprintf(stderr, "warning: failed to parse value of "
                "environment variable %s: %s", evn, strerror(errno));
      } else if (endptr == config_env) {
        fprintf(stderr, "warning: malformed value of environment variable "
                "%s, no digits or symbolic constant found\n", evn);
      } else if (v < 0 || v > 2) {
        fprintf(stderr, "value of environment variable %s (%ld) out of range [0,2]\n", evn, v);
      } else
        xt_config_set_xmap_stripe_align(&xt_default_config, (int)v);
    }
  }

  {
    static const char evn[] = "XT_CONFIG_DEFAULT_MEM_SAVING";
    config_env = getenv(evn);
    if (config_env) {
      const char *endptr;
      long v;
      if (!strcasecmp(config_env, "FORCE_NOSORT")) {
        v = 1;
        endptr = config_env + 12;
      } else {
        char *endptr_;
        v = strtol(config_env, &endptr_, 0);
        endptr = endptr_;
      }
      if ((errno == ERANGE && (v == LONG_MAX || v == LONG_MIN))
          || (errno != 0 && v == 0)) {
        perror("failed to parse value of "
               "XT_CONFIG_DEFAULT_MEM_SAVING environment variable");
      } else if (endptr == config_env) {
        fprintf(stderr, "warning: malformed value of environment variable "
                "%s, no digits or symbolic constant found\n", evn);
      } else if (v < 0 || v > XT_MAX_MEM_SAVING) {
        fprintf(stderr, "numeric value of XT_CONFIG_DEFAULT_MEM_SAVING"
                " environment variable (%ld) out of range [0,%d]\n",
                v, XT_MAX_MEM_SAVING);
      } else if (*endptr) {
        fprintf(stderr, "trailing text '%s' found after value (%*s) in "
                "%s environment variable\n",
                endptr, (int)(endptr-config_env), evn, config_env);
      } else
        xt_config_set_mem_saving(&xt_default_config, (int)v);
    }
  }
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
