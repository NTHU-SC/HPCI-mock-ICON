/**
 * @file xt_config.h
 * @brief opaque configuration object for settings where the default
 * needs to be overridden
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

#ifndef XT_CONFIG_H
#define XT_CONFIG_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <xt/xt_xmap_dist_dir_bucket_gen.h>

typedef struct Xt_config_ *Xt_config;

/**
 * constructor for configuration object
 *
 * @return returns a configuration object where every setting is set
 * to the corresponding default.
 */
Xt_config xt_config_new(void);

/**
 * destructor of configuration objects
 *
 * @param[in,out] config configuration object to destroy
 */
void xt_config_delete(Xt_config config);

enum Xt_exchangers {
  xt_exchanger_irecv_send,
  xt_exchanger_irecv_isend,
  xt_exchanger_irecv_isend_packed,
  xt_exchanger_mix_isend_irecv,
  xt_exchanger_neigh_alltoall,
  xt_exchanger_irecv_isend_ddt_packed,
};

/**
 * set exchanger to use when the \a config object is passed to constructors
 * @param[in,out] config configuration object to modify
 * @param method an entry from enum Xt_exchangers to signify the
 * desired exchanger for data transfers
 */
void xt_config_set_exchange_method(Xt_config config, int method);

/**
 * get exchanger used when the \a config object is passed to constructors
 * @param[in] config configuration object to query
 * @return an entry from \a Xt_exchangers representing the method of
 * data transfer used
 */
int xt_config_get_exchange_method(Xt_config config);

/**
 * map exchanger name string to method id from \a Xt_exchangers
 * @param[in] name string that is supposed to match the part of the
 * corresponding enum after xt_exchanger_
 * @return for the string "irecv_send", the value of
 * xt_exchanger_irecv_send will be returned, for strings matching no
 * known exchanger, -1 will be returned
 */
int
xt_exchanger_id_by_name(const char *name);

/**
 * query size above which index lists of vector type will be converted
 *
 * For many operations it makes sense to first compress large index
 * vectors to stripes before continuing further computations.
 *
 * @param[in] config   configuration object to query
 * @return             size of vectors at which conversion happens
 */
int
xt_config_get_idxvec_autoconvert_size(Xt_config config);

/**
 * set size above which index lists of vector type will be converted
 *
 * For many operations it makes sense to first compress large index
 * vectors to stripes before continuing further computations. This
 * function sets the size of vectors at which this conversion
 * happens for operations that are called with the configuration
 * object as parameter.
 *
 * @param[in,out] config   configuration object to modify
 * @param[in]     cnvsize  size of vectors at which conversion happens
 */
void
xt_config_set_idxvec_autoconvert_size(Xt_config config, int cnvsize);

enum Xt_mthread_mode {
  /* xt_redist_[as]_exchange calls will be single-threaded */
  XT_MT_NONE = 0,
  /* xt_redist_[as]_exchange calls will open an OpenMP parallel region */
  XT_MT_OPENMP = 1,
};

/**
 * query multi-thread mode of message passing
 *
 * @param[in] config   configuration object to query
 * @return a value matching one of the enum Xt_mthread_mode members above
 */
int
xt_config_get_redist_mthread_mode(Xt_config config);

/**
 * set multi-thread mode of message passing
 *
 * @param[in,out] config   configuration object to modify
 * @param[in]     mode one of the enum Xt_mthread_mode members above
 */
void
xt_config_set_redist_mthread_mode(Xt_config config, int mode);

enum Xt_sort_algorithm {
  /** use default, in-place algorithm */
  XT_QUICKSORT,
  /** use merge sort with allocation to trade space for better
   * worst-case behaviour */
  XT_MERGESORT,
};

/**
 * query sorting algorithm suite
 *
 * @param[in] config   configuration object to query
 * @return a value matching one of the enum Xt_sort_algorithm members above
 */
int
xt_config_get_sort_algorithm_id(Xt_config config);

/**
 * set sorting algorithm suite
 *
 * @param[in,out] config  configuration object to modify
 * @param[in]     algo    one of the enum Xt_sort_algorithm members above
 */
void
xt_config_set_sort_algorithm_by_id(Xt_config config, int algo);

/**
 * Set memory conservation parameter.
 *
 * @param[in,out]   config configuration object to modify
 * @param       memconserve switch on methods of memory saving
 *              currently supported values: 0 and 1
 *
 * When set to 1, operations parameterized by config will
 * aggressively try to trade memory savings for higher computational
 * cost.
 *
 * Currently, by not sorting index lists used in xmap construction.
 *
 * When set to 0, use default trade-offs for space vs. time complexity.
 *
 */
void
xt_config_set_mem_saving(Xt_config config, int memconserve);

/**
 * Get memory conservation parameter.
 * @param[in]   config configuration object to query
 */
int
xt_config_get_mem_saving(Xt_config config);

/**
 * Query configured bucket generator
 */
Xt_xmdd_bucket_gen
xt_config_get_xmdd_bucket_gen(Xt_config config);

/**
 * Change generator for buckets in distributed directory
 * The default is a generator that tiles the range of indices of all
 * index lists involved in xmap creation in 1D fashion with stripes.
 */
void
xt_config_set_xmdd_bucket_gen(Xt_config config,
                              Xt_xmdd_bucket_gen bucket_gen_iface);

/**
 * Set xmap stripe alignment parameter.
 *
 * @param[in,out]   config configuration object to modify
 * @param           use_stripe_align
 *                        0: detect element positions one by one,
 *                        1: align full stripes to form position extents
 *                        2: automatically choose element of stripe alignment
 *
 * When set to 1, the computations in the xmap
 * constructor can be much more expensive if stripes are very short.
 * When set to 0, use very simple algorithm that positions each
 * element individually but will use memory proportional to the local
 * index list sizes.
 *
 * By default the xmap constructor switches automatically to
 * stripe alignment if it appears profitable, which might be untrue if
 * the number of resulting stripes/extents is rather high.
 *
 */
void
xt_config_set_xmap_stripe_align(Xt_config config, int use_stripe_align);

/**
 * Get xmap stripe alignment parameter.
 * @param[in]   config configuration object to query
 */
int
xt_config_get_xmap_stripe_align(Xt_config config);


#endif

/*
 * Local Variables:
 * c-basic-offset: 2
 * coding: utf-8
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
