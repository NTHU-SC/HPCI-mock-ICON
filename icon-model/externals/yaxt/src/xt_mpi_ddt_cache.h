/**
 * @file xt_mpi_ddt_cache.h
 *
 * @copyright Copyright  (C)  2023 Jörg Behrens <behrens@dkrz.de>
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
#ifndef XT_MPI_DDT_CACHE_H
#define XT_MPI_DDT_CACHE_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <mpi.h>

#include "core/ppm_visibility.h"

struct Xt_mpiddt_list_entry;

struct Xt_mpiddt_list {
  struct Xt_mpiddt_list_entry *entries;
  MPI_Datatype *struct_dt;      /* points to concatenation of arrays of the datatypes
                                 * used for each struct */
  size_t num_entries, size_entries;
};

#define Xt_mpiddt_empty_list ((struct Xt_mpiddt_list){NULL,NULL,0,0})
/**
 * lookup MPI contiguous datatype in cache or create if not yet present
 * @param[in,out] ddt_list cache of already created derived datatypes
 * @param[in] count number of elements in contiguous datatype
 * @param[in] oldtype element type to create contiguous sequence
 * datatype of
 * @param[in] comm communicator to use for coordination
 */
PPM_DSO_INTERNAL MPI_Datatype
Xt_mpi_ddt_cache_acquire_contiguous(
  struct Xt_mpiddt_list *ddt_list,
  int count, MPI_Datatype oldtype,
  MPI_Comm comm);

/**
 * lookup MPI vector datatype in cache or create if not yet present
 * @param[in,out] ddt_list cache of already created derived datatypes
 * @param[in] count number of repeated blocks in vector datatype
 * @param[in] blocklength number of contiguous elements per repeated block
 * @param[in] stride vector stride in unit of elements from start of
 *                   one block to next block
 * @param[in] oldtype element type to create derived datatype of
 * @param[in] comm communicator to use for coordination
 */
PPM_DSO_INTERNAL MPI_Datatype
Xt_mpi_ddt_cache_acquire_vector(
  struct Xt_mpiddt_list *ddt_list,
  int count, int blocklength, int stride, MPI_Datatype oldtype,
  MPI_Comm comm);

/**
 * lookup MPI hvector datatype in cache or create if not yet present
 * @param[in,out] ddt_list cache of already created derived datatypes
 * @param[in] count number of repeated blocks in hvector datatype
 * @param[in] blocklength number of contiguous elements per repeated block
 * @param[in] stride hvector stride in bytes from start of
 *                   one block to next block
 * @param[in] oldtype element type to create derived datatype of
 * @param[in] comm communicator to use for coordination
 */
PPM_DSO_INTERNAL MPI_Datatype
Xt_mpi_ddt_cache_acquire_hvector(
  struct Xt_mpiddt_list *ddt_list,
  int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
  MPI_Comm comm);

/**
 * lookup MPI indexed block datatype in cache or create if not yet present
 * @param[in,out] ddt_list cache of already created derived datatypes
 * @param[in] count number of repeated blocks in indexed block datatype
 * @param[in] blocklength number of contiguous elements per repeated block
 * @param[in] disp displacement of each block in unit of elements from start
 * @param[in] oldtype element type to create derived datatype of
 * @param[in] comm communicator to use for coordination
 */
PPM_DSO_INTERNAL MPI_Datatype
Xt_mpi_ddt_cache_acquire_indexed_block(
  struct Xt_mpiddt_list *ddt_list,
  int count, int blocklength, const int disp[count], MPI_Datatype oldtype,
  MPI_Comm comm);

/**
 * lookup MPI hindexed block datatype in cache or create if not yet present
 * @param[in,out] ddt_list cache of already created derived datatypes
 * @param[in] count number of repeated blocks in hindexed block datatype
 * @param[in] blocklength number of contiguous elements per repeated block
 * @param[in] disp displacement of each block in bytes
 * @param[in] oldtype element type to create derived datatype of
 * @param[in] comm communicator to use for coordination
 */
PPM_DSO_INTERNAL MPI_Datatype
Xt_mpi_ddt_cache_acquire_hindexed_block(
  struct Xt_mpiddt_list *ddt_list,
  int count, int blocklength, const MPI_Aint disp[count], MPI_Datatype oldtype,
  MPI_Comm comm);

/**
 * lookup MPI indexed datatype in cache or create if not yet present
 * @param[in,out] ddt_list cache of already created derived datatypes
 * @param[in] count number of blocks in indexed datatype
 * @param[in] blocklength number of contiguous elements for each block
 * @param[in] disp displacement of each block in unit of elements
 * @param[in] oldtype element type to create derived datatype of
 * @param[in] comm communicator to use for coordination
 */
PPM_DSO_INTERNAL MPI_Datatype
Xt_mpi_ddt_cache_acquire_indexed(
  struct Xt_mpiddt_list *ddt_list,
  int count, const int blocklength[count], const int disp[count],
  MPI_Datatype oldtype, MPI_Comm comm);

/**
 * lookup MPI hindexed datatype in cache or create if not yet present
 * @param[in,out] ddt_list cache of already created derived datatypes
 * @param[in] count number of blocks in hindexed datatype
 * @param[in] blocklength number of contiguous elements for each block
 * @param[in] disp displacement of each block in bytes
 * @param[in] oldtype element type to create derived datatype of
 * @param[in] comm communicator to use for coordination
 */
PPM_DSO_INTERNAL MPI_Datatype
Xt_mpi_ddt_cache_acquire_hindexed(
  struct Xt_mpiddt_list *ddt_list,
  int count, const int blocklength[count], const MPI_Aint disp[count],
  MPI_Datatype oldtype, MPI_Comm comm);

/**
 * lookup MPI struct datatype in cache or create if not yet present
 * @param[in,out] ddt_list cache of already created derived datatypes
 * @param[in] count number of blocks in struct datatype
 * @param[in] blocklength number of contiguous elements for each block
 * @param[in] disp displacement of each block in bytes
 * @param[in] oldtype element type for each block
 * @param[in] comm communicator to use for coordination
 */
PPM_DSO_INTERNAL MPI_Datatype
Xt_mpi_ddt_cache_acquire_struct(
  struct Xt_mpiddt_list *ddt_list,
  int count, const int blocklength[count], const MPI_Aint disp[count],
  const MPI_Datatype oldtype[count], MPI_Comm comm);


/**
 * reduce reference counter for MPI datatype in cache
 * @param[in,out] ddt_list cache of already created derived datatypes
 * @param[in] dt derived datatype to mark as used one instance less
 * @param[in] comm communicator to use for coordination and error reporting
 */
PPM_DSO_INTERNAL void
Xt_mpi_ddt_cache_entry_release(struct Xt_mpiddt_list *ddt_list,
                               MPI_Datatype *dt, MPI_Comm comm);

/**
 * Remove ddt cache data structure
 *
 * This function will also call MPI_Type_free for all datatypes in
 * cache that have 0 reference count
 * @param[in,out] ddt_list cache of already created derived datatypes
 * @param[in] comm communicator to use for coordination and error reporting
  */
PPM_DSO_INTERNAL void
Xt_mpi_ddt_cache_free(struct Xt_mpiddt_list *ddt_list,
                      MPI_Comm comm);

#include "xt/xt_redist.h"

PPM_DSO_INTERNAL void
Xt_mpi_ddt_cache_check_retention(struct Xt_mpiddt_list *ddt_list,
                                 size_t nmsg,
                                 struct Xt_redist_msg msgs[nmsg]);


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
