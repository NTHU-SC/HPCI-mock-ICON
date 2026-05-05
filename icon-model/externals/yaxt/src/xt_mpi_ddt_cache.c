/**
 * @file xt_mpi_ddt_cache.c
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
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <assert.h>
#include <stdint.h>
#include <string.h>

#include <mpi.h>

#include "xt/xt_mpi.h"
#include "xt_arithmetic_util.h"
#include "core/cksum.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "xt_mpi_ddt_wrap.h"
#include "xt_mpi_ddt_cache.h"

#if ! HAVE_DECL___BUILTIN_CLZL       \
  && (HAVE_DECL___LZCNT && SIZEOF_LONG == SIZEOF_INT               \
      || HAVE_DECL___LZCNT64 && SIZEOF_LONG == 8 && CHAR_BIT == 8)
#include <intrin.h>
#endif


struct Xt_mpi_contiguous_arg_desc {
  int count;
  MPI_Datatype oldtype;
};

struct Xt_mpi_vector_arg_desc {
  int count, blocklength, stride;
  MPI_Datatype oldtype;
};

struct Xt_mpi_hvector_arg_desc {
  int count, blocklength;
  MPI_Aint stride;
  MPI_Datatype oldtype;
};

struct Xt_mpi_indexed_block_arg_desc {
  int count, blocklength;
  uint32_t disp_hash;
  MPI_Datatype oldtype;
};

struct Xt_mpi_indexed_arg_desc {
  int count;
  uint32_t blocklength_hash, disp_hash;
  MPI_Datatype oldtype;
};

struct Xt_mpi_struct_arg_desc {
  int count;
  uint32_t blocklength_hash, disp_hash, oldtype_hash;
};

struct Xt_mpiddt_list_entry {
  union {
    struct Xt_mpi_contiguous_arg_desc contiguous;
    struct Xt_mpi_vector_arg_desc vector;
    struct Xt_mpi_hvector_arg_desc hvector;
    struct Xt_mpi_indexed_block_arg_desc indexed_block;
    struct Xt_mpi_indexed_arg_desc indexed;
    struct Xt_mpi_struct_arg_desc struct_dt;
  } args;
  MPI_Datatype cached_dt;
  int combiner, use_count;
};

static struct Xt_mpiddt_list_entry *
grow_ddt_list(struct Xt_mpiddt_list *ddt_list)
{
  size_t size_entries = ddt_list->size_entries;
  ddt_list->size_entries = size_entries = size_entries ? size_entries * 2 : 8;
  return ddt_list->entries
    = xrealloc(ddt_list->entries, size_entries * sizeof (*ddt_list->entries));
}

#define GROW_DDT_LIST(ddt_list)                                 \
  do {                                                          \
    if (ddt_list->num_entries == ddt_list->size_entries)        \
      entries = grow_ddt_list(ddt_list);                        \
  } while (0)

static inline void
free_dt_unless_named(MPI_Datatype *dt, MPI_Comm comm)
{
  int num_integers, num_addresses, num_datatypes, combiner;
  xt_mpi_call(MPI_Type_get_envelope(*dt, &num_integers,
                                    &num_addresses, &num_datatypes, &combiner), comm);
  if (combiner != MPI_COMBINER_NAMED)
    xt_mpi_call(MPI_Type_free(dt), comm);
}


MPI_Datatype
Xt_mpi_ddt_cache_acquire_contiguous(
  struct Xt_mpiddt_list *ddt_list,
  int count, MPI_Datatype oldtype,
  MPI_Comm comm)
{
  MPI_Datatype dt;
  if (ddt_list) {
    struct Xt_mpiddt_list_entry *restrict entries = ddt_list->entries;
    size_t num_entries = ddt_list->num_entries;
    for (size_t i = 0; i < num_entries; ++i)
      if (entries[i].combiner == MPI_COMBINER_CONTIGUOUS) {
        struct Xt_mpi_contiguous_arg_desc *args
          = &entries[i].args.contiguous;
        if (args->count == count && args->oldtype == oldtype) {
          entries[i].use_count += 1;
          dt = entries[i].cached_dt;
          goto dt_is_set;
        }
      }
    GROW_DDT_LIST(ddt_list);
    xt_mpi_call(MPI_Type_contiguous(count, oldtype, &dt), comm);
    entries[num_entries] = (struct Xt_mpiddt_list_entry){
      .args.contiguous = (struct Xt_mpi_contiguous_arg_desc){
        .count = count, .oldtype = oldtype,
      },
      .cached_dt = dt,
      .use_count = 1,
      .combiner = MPI_COMBINER_CONTIGUOUS,
    };
    ddt_list->num_entries = num_entries + 1;
  } else {
    xt_mpi_call(MPI_Type_contiguous(count, oldtype, &dt), comm);
  }
dt_is_set:
  return dt;
}

MPI_Datatype
Xt_mpi_ddt_cache_acquire_vector(
  struct Xt_mpiddt_list *ddt_list,
  int count, int blocklength, int stride, MPI_Datatype oldtype,
  MPI_Comm comm)
{
  MPI_Datatype dt;
  if (ddt_list) {
    struct Xt_mpiddt_list_entry *restrict entries = ddt_list->entries;
    size_t num_entries = ddt_list->num_entries;
    for (size_t i = 0; i < num_entries; ++i)
      if (entries[i].combiner == MPI_COMBINER_VECTOR) {
        struct Xt_mpi_vector_arg_desc *args = &entries[i].args.vector;
        if (args->count == count && args->blocklength == blocklength
            && args->oldtype == oldtype && args->stride == stride) {
          entries[i].use_count += 1;
          dt = entries[i].cached_dt;
          goto dt_is_set;
        }
      }
    GROW_DDT_LIST(ddt_list);
    xt_mpi_call(MPI_Type_vector(count, blocklength, stride, oldtype, &dt),
                comm);
    entries[num_entries] = (struct Xt_mpiddt_list_entry){
      .args.vector = (struct Xt_mpi_vector_arg_desc){
        .count = count, .blocklength = blocklength,
        .stride = stride, .oldtype = oldtype,
      },
      .cached_dt = dt,
      .use_count = 1,
      .combiner = MPI_COMBINER_VECTOR,
    };
    ddt_list->num_entries = num_entries + 1;
  } else {
    xt_mpi_call(MPI_Type_vector(count, blocklength, stride, oldtype, &dt),
                comm);
  }
dt_is_set:
  return dt;
}

MPI_Datatype
Xt_mpi_ddt_cache_acquire_hvector(
  struct Xt_mpiddt_list *ddt_list,
  int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
  MPI_Comm comm)
{
  MPI_Datatype dt;
  if (ddt_list) {
    struct Xt_mpiddt_list_entry *restrict entries = ddt_list->entries;
    size_t num_entries = ddt_list->num_entries;
    for (size_t i = 0; i < num_entries; ++i)
      if (entries[i].combiner == MPI_COMBINER_HVECTOR) {
        struct Xt_mpi_hvector_arg_desc *args = &entries[i].args.hvector;
        if (args->count == count && args->blocklength == blocklength
            && args->oldtype == oldtype && args->stride == stride) {
          entries[i].use_count += 1;
          dt = entries[i].cached_dt;
          goto dt_is_set;
        }
      }
    GROW_DDT_LIST(ddt_list);
    xt_mpi_call(MPI_Type_create_hvector(count, blocklength, stride,
                                        oldtype, &dt), comm);
    entries[num_entries] = (struct Xt_mpiddt_list_entry){
      .args.hvector = (struct Xt_mpi_hvector_arg_desc){
        .count = count, .blocklength = blocklength,
        .stride = stride, .oldtype = oldtype,
      },
      .cached_dt = dt,
      .use_count = 1,
      .combiner = MPI_COMBINER_HVECTOR,
    };
    ddt_list->num_entries = num_entries + 1;
  } else {
    xt_mpi_call(MPI_Type_create_hvector(count, blocklength, stride,
                                        oldtype, &dt), comm);
  }
dt_is_set:
  return dt;
}

MPI_Datatype
Xt_mpi_ddt_cache_acquire_indexed_block(
  struct Xt_mpiddt_list *ddt_list,
  int count, int blocklength, const int disp[count], MPI_Datatype oldtype,
  MPI_Comm comm)
{
  MPI_Datatype dt;
  if (ddt_list) {
    size_t disp_size = (count > 0 ? (size_t)count : (size_t)0) * sizeof (*disp);
    uint32_t disp_hash = Xt_memcrc((const void *)disp, disp_size);
    struct Xt_mpiddt_list_entry *restrict entries = ddt_list->entries;
    size_t num_entries = ddt_list->num_entries;
    int *disp_cmp = NULL;
    for (size_t i = 0; i < num_entries; ++i)
      if (entries[i].combiner == MPI_COMBINER_INDEXED_BLOCK) {
        struct Xt_mpi_indexed_block_arg_desc *args
          = &entries[i].args.indexed_block;
        if (args->count == count && args->blocklength == blocklength
            && args->disp_hash == disp_hash && args->oldtype == oldtype) {
          if (!disp_cmp)
            disp_cmp = xmalloc(disp_size + 2 * sizeof (int));
          MPI_Datatype cached_dt = entries[i].cached_dt, oldtype_;
          xt_mpi_call(MPI_Type_get_contents(cached_dt, count + 2, 0, 1,
                                            disp_cmp, NULL, &oldtype_), comm);
          free_dt_unless_named(&oldtype_, comm);
          if (memcmp(disp, disp_cmp+2, disp_size))
            continue;
          entries[i].use_count += 1;
          dt = cached_dt;
          goto dt_is_set;
        }
      }
    GROW_DDT_LIST(ddt_list);
    Xt_Type_create_indexed_block(count, blocklength,
                                 disp, oldtype, &dt, comm);
    entries[num_entries] = (struct Xt_mpiddt_list_entry){
      .args.indexed_block = (struct Xt_mpi_indexed_block_arg_desc){
        .count = count, .blocklength = blocklength,
        .disp_hash = disp_hash, .oldtype = oldtype,
      },
      .cached_dt = dt,
      .use_count = 1,
      .combiner = MPI_COMBINER_INDEXED_BLOCK,
    };
    ddt_list->num_entries = num_entries + 1;
  dt_is_set:
    free(disp_cmp);
  } else {
    Xt_Type_create_indexed_block(count, blocklength,
                                 disp, oldtype, &dt, comm);
  }
  return dt;
}

MPI_Datatype
Xt_mpi_ddt_cache_acquire_hindexed_block(
  struct Xt_mpiddt_list *ddt_list,
  int count, int blocklength, const MPI_Aint disp[count], MPI_Datatype oldtype,
  MPI_Comm comm)
{
  MPI_Datatype dt;
  if (ddt_list) {
    size_t count_ = count > 0 ? (size_t)count : (size_t)0,
      disp_size = count_ * sizeof (*disp);
    uint32_t disp_hash = Xt_memcrc((const void *)disp, disp_size);
    struct Xt_mpiddt_list_entry *restrict entries = ddt_list->entries;
    size_t num_entries = ddt_list->num_entries;
    MPI_Aint *disp_cmp = NULL;
    for (size_t i = 0; i < num_entries; ++i)
      if (entries[i].combiner == MPI_COMBINER_HINDEXED_BLOCK) {
        struct Xt_mpi_indexed_block_arg_desc *args
          = &entries[i].args.indexed_block;
        if (args->count == count && args->blocklength == blocklength
            && args->disp_hash == disp_hash && args->oldtype == oldtype) {
#if MPI_VERSION < 3
#define disp_size (disp_size + (count_ + 2) * sizeof (int))
#endif
          if (!disp_cmp)
            disp_cmp = xmalloc(disp_size);
          MPI_Datatype cached_dt = entries[i].cached_dt, oldtype_;
#if MPI_VERSION >= 3
          int icmp[2];
#else
#undef disp_size
          int *icmp = (void *)(disp_cmp + count_);
#endif
          xt_mpi_call(MPI_Type_get_contents(cached_dt, 2 + count, count, 1,
                                            icmp, disp_cmp, &oldtype_), comm);
          free_dt_unless_named(&oldtype_, comm);
          if (memcmp(disp, disp_cmp, disp_size))
            continue;
          entries[i].use_count += 1;
          dt = cached_dt;
          goto dt_is_set;
        }
      }
    GROW_DDT_LIST(ddt_list);
    Xt_Type_create_hindexed_block(count, blocklength,
                                  disp, oldtype, &dt, comm);
    entries[num_entries] = (struct Xt_mpiddt_list_entry){
      .args.indexed_block = (struct Xt_mpi_indexed_block_arg_desc){
        .count = count, .blocklength = blocklength,
        .disp_hash = disp_hash, .oldtype = oldtype,
      },
      .cached_dt = dt,
      .use_count = 1,
      .combiner = MPI_COMBINER_HINDEXED_BLOCK,
    };
    ddt_list->num_entries = num_entries + 1;
  dt_is_set:
    free(disp_cmp);
  } else {
    Xt_Type_create_hindexed_block(count, blocklength,
                                  disp, oldtype, &dt, comm);
  }
  return dt;
}


MPI_Datatype
Xt_mpi_ddt_cache_acquire_indexed(
  struct Xt_mpiddt_list *ddt_list,
  int count, const int blocklength[count], const int disp[count],
  MPI_Datatype oldtype, MPI_Comm comm)
{
  MPI_Datatype dt;
  if (ddt_list) {
    size_t asize = (count > 0 ? (size_t)count : (size_t)0) * sizeof (int);
    uint32_t disp_hash = Xt_memcrc((const void *)disp, asize),
      blocklength_hash = Xt_memcrc((const void *)blocklength, asize);
    struct Xt_mpiddt_list_entry *restrict entries = ddt_list->entries;
    size_t num_entries = ddt_list->num_entries;
    int *acmp = NULL;
    for (size_t i = 0; i < num_entries; ++i)
      if (entries[i].combiner == MPI_COMBINER_INDEXED) {
        struct Xt_mpi_indexed_arg_desc *args
          = &entries[i].args.indexed;
        if (args->count == count && args->blocklength_hash == blocklength_hash
            && args->disp_hash == disp_hash && args->oldtype == oldtype) {
          if (!acmp)
            acmp = xmalloc(2 * asize + sizeof (int));
          MPI_Datatype cached_dt = entries[i].cached_dt, oldtype_;
          xt_mpi_call(MPI_Type_get_contents(cached_dt, 2 * count + 1, 0, 1,
                                            acmp, NULL, &oldtype_), comm);
          free_dt_unless_named(&oldtype_, comm);
          if (memcmp(blocklength, acmp+1, asize)
              || memcmp(disp, acmp+count+1, asize))
            continue;
          entries[i].use_count += 1;
          dt = cached_dt;
          goto dt_is_set;
        }
      }
    GROW_DDT_LIST(ddt_list);
    Xt_Type_indexed(count, blocklength,
                    disp, oldtype, &dt, comm);
    entries[num_entries] = (struct Xt_mpiddt_list_entry){
      .args.indexed = (struct Xt_mpi_indexed_arg_desc){
        .count = count, .blocklength_hash = blocklength_hash,
        .disp_hash = disp_hash, .oldtype = oldtype,
      },
      .cached_dt = dt,
      .use_count = 1,
      .combiner = MPI_COMBINER_INDEXED,
    };
    ddt_list->num_entries = num_entries + 1;
  dt_is_set:
    free(acmp);
  } else {
    Xt_Type_indexed(count, blocklength, disp, oldtype, &dt, comm);
  }
  return dt;
}

MPI_Datatype
Xt_mpi_ddt_cache_acquire_hindexed(
  struct Xt_mpiddt_list *ddt_list,
  int count, const int blocklength[count], const MPI_Aint disp[count],
  MPI_Datatype oldtype, MPI_Comm comm)
{
  MPI_Datatype dt;
  if (ddt_list) {
    size_t count_ = count > 0 ? (size_t)count : (size_t)0,
      disp_size = count_ * sizeof (*disp),
      blocklength_size = count_ * sizeof (*blocklength);
    uint32_t disp_hash = Xt_memcrc((const void *)disp, disp_size),
      blocklength_hash = Xt_memcrc((const void *)blocklength, blocklength_size);
    struct Xt_mpiddt_list_entry *restrict entries = ddt_list->entries;
    size_t num_entries = ddt_list->num_entries;
    MPI_Aint *disp_cmp = NULL;
    for (size_t i = 0; i < num_entries; ++i)
      if (entries[i].combiner == MPI_COMBINER_HINDEXED) {
        struct Xt_mpi_indexed_arg_desc *args
          = &entries[i].args.indexed;
        if (args->count == count && args->blocklength_hash == blocklength_hash
            && args->disp_hash == disp_hash && args->oldtype == oldtype) {
          if (!disp_cmp)
            disp_cmp = xmalloc(sizeof (int) + blocklength_size + disp_size);
          MPI_Datatype cached_dt = entries[i].cached_dt, oldtype_;
          xt_mpi_call(MPI_Type_get_contents(
                        cached_dt, count + 1, count, 1,
                        (void *)(disp_cmp+count_), disp_cmp, &oldtype_), comm);
          free_dt_unless_named(&oldtype_, comm);
          if (memcmp(blocklength, disp_cmp+count_+1, blocklength_size)
              || memcmp(disp, disp_cmp, disp_size))
            continue;
          entries[i].use_count += 1;
          dt = cached_dt;
          goto dt_is_set;
        }
      }
    GROW_DDT_LIST(ddt_list);
    Xt_Type_create_hindexed(count, blocklength,
                            disp, oldtype, &dt, comm);
    entries[num_entries] = (struct Xt_mpiddt_list_entry){
      .args.indexed = (struct Xt_mpi_indexed_arg_desc){
        .count = count, .blocklength_hash = blocklength_hash,
        .disp_hash = disp_hash, .oldtype = oldtype,
      },
      .cached_dt = dt,
      .use_count = 1,
      .combiner = MPI_COMBINER_HINDEXED,
    };
    ddt_list->num_entries = num_entries + 1;
  dt_is_set:
    free(disp_cmp);
  } else {
    Xt_Type_create_hindexed(count, blocklength,
                            disp, oldtype, &dt, comm);
  }
  return dt;
}


MPI_Datatype
Xt_mpi_ddt_cache_acquire_struct(
  struct Xt_mpiddt_list *ddt_list,
  int count, const int blocklength[count],
  const MPI_Aint disp[count],
  const MPI_Datatype oldtype[count], MPI_Comm comm)
{
  MPI_Datatype dt;
  if (ddt_list) {
    size_t count_ = (count > 0 ? (size_t)count : (size_t)0),
      disp_size = count_ * sizeof (disp[0]),
      blocklength_size = count_ * sizeof (blocklength[0]),
      oldtype_size = count_ * sizeof (oldtype[0]);
    uint32_t disp_hash = Xt_memcrc((const void *)disp, disp_size),
      blocklength_hash = Xt_memcrc((const void *)blocklength, blocklength_size),
      oldtype_hash = Xt_memcrc((const void *)oldtype, oldtype_size);
    struct Xt_mpiddt_list_entry *restrict entries = ddt_list->entries;
    size_t num_entries = ddt_list->num_entries;
    MPI_Aint *disp_cmp = NULL;
    MPI_Datatype *oldtype_contents = NULL,
      *oldtype_cmp = ddt_list->struct_dt;
    int *icmp = NULL;
    for (size_t i = 0; i < num_entries; ++i)
      if (entries[i].combiner == MPI_COMBINER_STRUCT) {
        struct Xt_mpi_struct_arg_desc *args
          = &entries[i].args.struct_dt;
        if (args->count == count && args->blocklength_hash == blocklength_hash
            && args->disp_hash == disp_hash
            && args->oldtype_hash == oldtype_hash) {
          if (!disp_cmp) {
            disp_cmp = xmalloc(disp_size + oldtype_size
                               + sizeof (int) + blocklength_size);
            oldtype_contents = (void *)(disp_cmp + count_);
            icmp = (void *)(oldtype_contents + count_);
          }
          MPI_Datatype cached_dt = entries[i].cached_dt;
          xt_mpi_call(MPI_Type_get_contents(cached_dt, count+1, count, count,
                                            icmp, disp_cmp, oldtype_contents), comm);
          int oldtypes_mismatch = memcmp(oldtype, oldtype_cmp, oldtype_size);
          for (size_t j = 0; j < count_; ++j)
            free_dt_unless_named(oldtype_contents+j, comm);
          assert(icmp[0] == count);
          if (!oldtypes_mismatch && !memcmp(blocklength, icmp+1, blocklength_size)
              && !memcmp(disp, disp_cmp, disp_size)) {
            entries[i].use_count += 1;
            dt = cached_dt;
            goto dt_is_set;
          }
        }
        oldtype_cmp += args->count;
      }
    GROW_DDT_LIST(ddt_list);
    size_t struct_dt_size = (size_t)(oldtype_cmp - ddt_list->struct_dt),
      struct_dt_size_p2 = next_2_pow(struct_dt_size),
      struct_dt_needed = struct_dt_size + count_ + (oldtype_cmp == ddt_list->struct_dt);
    if (struct_dt_needed > struct_dt_size_p2) {
      ddt_list->struct_dt
        = xrealloc(ddt_list->struct_dt,
                   next_2_pow(struct_dt_needed) * sizeof (*oldtype_cmp));
      oldtype_cmp = ddt_list->struct_dt + struct_dt_size;
    }
    memcpy(oldtype_cmp, oldtype, oldtype_size);
    Xt_Type_create_struct(count, blocklength, disp, oldtype, &dt, comm);
    entries[num_entries] = (struct Xt_mpiddt_list_entry){
      .args.struct_dt = (struct Xt_mpi_struct_arg_desc){
        .count = count, .blocklength_hash = blocklength_hash,
        .disp_hash = disp_hash, .oldtype_hash = oldtype_hash,
      },
      .cached_dt = dt,
      .use_count = 1,
      .combiner = MPI_COMBINER_STRUCT,
    };
    ddt_list->num_entries = num_entries + 1;
  dt_is_set:
    free(disp_cmp);
  } else {
    Xt_Type_create_struct(count, blocklength, disp, oldtype, &dt, comm);
  }
  return dt;
}

void
Xt_mpi_ddt_cache_entry_release(struct Xt_mpiddt_list *ddt_list,
                               MPI_Datatype *dt, MPI_Comm comm)
{
  if (ddt_list) {
    struct Xt_mpiddt_list_entry *restrict entries = ddt_list->entries;
    size_t num_entries = ddt_list->num_entries;
    MPI_Datatype dt_ = *dt;
    for (size_t i = 0; i < num_entries; ++i)
      if (entries[i].cached_dt == dt_) {
#ifndef NDEBUG
        int new_use_count =
#endif
          --entries[i].use_count;
        assert(new_use_count >= 0);
        /**
         * @todo: implement heuristic to free datatypes going unused
         * in a while
         */
        *dt = MPI_DATATYPE_NULL;
        return;
      }
  }
  xt_mpi_call(MPI_Type_free(dt), comm);
}

void
Xt_mpi_ddt_cache_free(struct Xt_mpiddt_list *ddt_list,
                      MPI_Comm comm)
{
  if (ddt_list) ; else return;
  struct Xt_mpiddt_list_entry *restrict entries = ddt_list->entries;
  size_t num_entries = ddt_list->num_entries;
  for (size_t i = 0; i < num_entries; ++i)
    if (!entries[i].use_count)
      xt_mpi_call(MPI_Type_free(&entries[i].cached_dt), comm);
  free(entries);
  free(ddt_list->struct_dt);
  ddt_list->struct_dt = NULL;
  ddt_list->entries = NULL;
  ddt_list->num_entries = 0;
  ddt_list->size_entries = 0;
}

void
Xt_mpi_ddt_cache_check_retention(struct Xt_mpiddt_list *ddt_list,
                                 size_t nmsg,
                                 struct Xt_redist_msg msgs[nmsg])
{
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  for (size_t i = 0, n = ddt_list->num_entries; i < n; ++i)
    if (ddt_list->entries[i].use_count) {
      MPI_Datatype cached_dt = ddt_list->entries[i].cached_dt;
      for (size_t j = 0; j < nmsg; ++j)
        if (msgs[j].datatype == cached_dt)
          goto use_count_is_fine;
      char buf[256];
      sprintf(buf, "%d: cache inconsistency: In-use marked datatype "
              "encountered that is not in any mesage!\n", world_rank);
      Xt_abort(Xt_default_comm, buf, "xt_mpi_ddt_cache.c", __LINE__);
    use_count_is_fine:
      ;
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
