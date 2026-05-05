/**
 * @file xt_idxlist_collection.c
 *
 * @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
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
#include <config.h>
#endif

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "mpi.h"

#include "xt/xt_core.h"
#include "xt/xt_idxlist.h"
#include "xt_idxlist_internal.h"
#include "xt/xt_idxlist_collection.h"
#include "xt_idxlist_collection_internal.h"
#include "xt/xt_idxempty.h"
#include "xt/xt_config.h"
#include "xt_config_internal.h"
#include "xt_idxvec_internal.h"
#include "xt/xt_idxstripes.h"
#include "xt_idxstripes_internal.h"
#include "xt/xt_mpi.h"
#include "xt_idxlist_unpack.h"
#include "xt_stripe_util.h"
#include "core/core.h"
#include "core/ppm_xfuncs.h"
#include "ensure_array_size.h"

static const char filename[] = "xt_idxlist_collection.c";

static void
idxlist_collection_delete(Xt_idxlist data);

static size_t
idxlist_collection_get_pack_size(Xt_idxlist data, MPI_Comm comm);

static void
idxlist_collection_pack(Xt_idxlist data, void *buffer, int buffer_size,
                        int *position, MPI_Comm comm);

static Xt_idxlist
idxlist_collection_copy(Xt_idxlist idxlist);

static Xt_idxlist
idxlist_collection_sorted_copy(Xt_idxlist idxlist, Xt_config config);

static void
idxlist_collection_get_indices(Xt_idxlist idxlist, Xt_int *indices);

static const Xt_int *
idxlist_collection_get_indices_const(Xt_idxlist idxlist);

static int
idxlist_collection_get_num_index_stripes(Xt_idxlist idxlist);

static void
idxlist_collection_get_index_stripes(Xt_idxlist idxlist,
                                     struct Xt_stripe *restrict stripes,
                                     size_t num_stripes_alloc);

static int
idxlist_collection_get_index_at_position(Xt_idxlist idxlist, int position,
                                         Xt_int * index);

static int
idxlist_collection_get_position_of_index(Xt_idxlist idxlist, Xt_int index,
                                         int * position);

static int
idxlist_collection_get_position_of_index_off(Xt_idxlist idxlist, Xt_int index,
                                             int * position, int offset);

static Xt_int
idxlist_collection_get_min_index(Xt_idxlist idxlist);

static Xt_int
idxlist_collection_get_max_index(Xt_idxlist idxlist);

static int
idxlist_collection_get_sorting(Xt_idxlist idxlist);

static const struct xt_idxlist_vtable idxlist_collection_vtable = {
  .delete                      = idxlist_collection_delete,
  .get_pack_size               = idxlist_collection_get_pack_size,
  .pack                        = idxlist_collection_pack,
  .copy                        = idxlist_collection_copy,
  .sorted_copy			= idxlist_collection_sorted_copy,
  .get_indices                 = idxlist_collection_get_indices,
  .get_indices_const           = idxlist_collection_get_indices_const,
  .get_num_index_stripes	= idxlist_collection_get_num_index_stripes,
  .get_index_stripes           = idxlist_collection_get_index_stripes,
  .get_index_at_position       = idxlist_collection_get_index_at_position,
  .get_indices_at_positions    = NULL,
  .get_position_of_index       = idxlist_collection_get_position_of_index,
  .get_positions_of_indices    = NULL,
  .get_position_of_index_off   = idxlist_collection_get_position_of_index_off,
  .get_positions_of_indices_off = NULL,
  .get_min_index               = idxlist_collection_get_min_index,
  .get_max_index               = idxlist_collection_get_max_index,
  .get_sorting			= idxlist_collection_get_sorting,
  .get_bounding_box            = NULL,
  .idxlist_pack_code           = COLLECTION,
};

typedef struct Xt_idxlist_collection_ *Xt_idxlist_collection;

struct Xt_idxlist_collection_ {

  struct Xt_idxlist_ parent;

  int num_idxlists;
  unsigned flags;
  Xt_int min, max;
  Xt_int *index_array_cache;
  Xt_idxlist idxlists[];
};

Xt_idxlist
xt_idxlist_collection_new(Xt_idxlist *idxlists, int num_idxlists) {
  // ensure that yaxt is initialized
  assert(xt_initialized());

  /* todo: find collections in idxlists and prevent hierarchical
   * copy, flatten list instead  */
  Xt_idxlist result;
  size_t num_non_empty_idxlists = 0, first_non_empty_idxlist = (size_t)-1;
  for (int i = 0; i < num_idxlists; ++i) {
    int num_indices_of_list = xt_idxlist_get_num_indices(idxlists[i]);
    if (first_non_empty_idxlist == (size_t)-1 && num_indices_of_list > 0)
      first_non_empty_idxlist = (size_t)i;
    num_non_empty_idxlists +=
      (idxlists[i]->vtable == &idxlist_collection_vtable)
      ? (size_t)((Xt_idxlist_collection)idxlists[i])->num_idxlists
      : (size_t)(num_indices_of_list > 0);
  }
  if (num_non_empty_idxlists  > 1)
  {
    long long num_indices
      = xt_idxlist_get_num_indices(idxlists[first_non_empty_idxlist]);
    Xt_int min = xt_idxlist_get_min_index(idxlists[first_non_empty_idxlist]),
      max = xt_idxlist_get_max_index(idxlists[first_non_empty_idxlist]),
      prev_max = max;
    unsigned ntrans_dn = 0, ntrans_up = 0,
      nsort_asc = 0, nsort_dsc = 0;
    switch (xt_idxlist_get_sorting(idxlists[first_non_empty_idxlist])) {
    case 1:
      ++nsort_asc;
      break;
    case -1:
      ++nsort_dsc;
      break;
    case 0:
      break;
#  if HAVE_DECL___BUILTIN_UNREACHABLE
    default:
      __builtin_unreachable();
#  endif
    }
    for (int i = (int)first_non_empty_idxlist+1; i < num_idxlists; ++i) {
      int num_indices_of_list = xt_idxlist_get_num_indices(idxlists[i]);
      if (num_indices_of_list > 0) {
        Xt_int tmp_min, tmp_max;
        if ((tmp_min = xt_idxlist_get_min_index(idxlists[i])) < min)
          min = tmp_min;
        if ((tmp_max = xt_idxlist_get_max_index(idxlists[i])) > max)
          max = tmp_max;
        ntrans_dn += (tmp_min < prev_max);
        ntrans_up += (tmp_min > prev_max);
        prev_max = tmp_max;
        int sort = xt_idxlist_get_sorting(idxlists[i]);
        switch (sort) {
        case 1:
          ++nsort_asc;
          break;
        case -1:
          ++nsort_dsc;
          break;
        case 0:
          break;
#  if HAVE_DECL___BUILTIN_UNREACHABLE
        default:
          __builtin_unreachable();
#  endif
        }
        num_indices += num_indices_of_list;
      }
    }
    assert(num_indices <= INT_MAX);
    if (min != max) {
      Xt_idxlist_collection collectionlist
        = xmalloc(sizeof (*collectionlist)
                  + num_non_empty_idxlists
                  * sizeof (collectionlist->idxlists[0]));

      collectionlist->num_idxlists = (int)num_non_empty_idxlists;
      collectionlist->index_array_cache = NULL;

      for (size_t i = first_non_empty_idxlist, j = 0;
           i < (size_t)num_idxlists; ++i) {
        if (idxlists[i]->vtable == &idxlist_collection_vtable) {
          Xt_idxlist_collection comp_collection
            = (Xt_idxlist_collection)idxlists[i];
          size_t num_comp_idxlists
            = (size_t)comp_collection->num_idxlists;
          Xt_idxlist *comp_lists = comp_collection->idxlists;
          for (size_t k = 0; k < num_comp_idxlists; ++k)
            collectionlist->idxlists[j+k] = xt_idxlist_copy(comp_lists[k]);
          j += num_comp_idxlists;
        } else if (xt_idxlist_get_num_indices(idxlists[i]) > 0)
          collectionlist->idxlists[j++] = xt_idxlist_copy(idxlists[i]);
      }
      Xt_idxlist_init(&collectionlist->parent, &idxlist_collection_vtable,
                      (int)num_indices);
      collectionlist->flags = XT_SORT_FLAGS(ntrans_up + nsort_asc,
                                            ntrans_dn + nsort_dsc) & sort_mask;
      result = (Xt_idxlist)collectionlist;
    } else /* min == max => all values identical */ {
      result =
        xt_idxstripes_new(&(struct Xt_stripe){
            .start = min, .stride = 0, .nstrides = (int)num_indices }, 1);
    }
  }
  else if (num_non_empty_idxlists == 1)
    result = xt_idxlist_copy(idxlists[first_non_empty_idxlist]);
  else /* num_idxlists == 0 */
    result = xt_idxempty_new();
  return result;
}

static void
idxlist_collection_delete(Xt_idxlist data) {

   Xt_idxlist_collection collectionlist = (Xt_idxlist_collection)data;

   int num_lists = collectionlist->num_idxlists;
   for (int i = 0; i < num_lists; ++i)
      xt_idxlist_delete(collectionlist->idxlists[i]);

   free(collectionlist->index_array_cache);
   free(collectionlist);
}

static size_t
idxlist_collection_get_pack_size(Xt_idxlist data, MPI_Comm comm) {

   Xt_idxlist_collection collectionlist = (Xt_idxlist_collection)data;

   int size_header, num_lists = collectionlist->num_idxlists;
   size_t size_idxlists = 0;

   xt_mpi_call(MPI_Pack_size(2, MPI_INT, comm, &size_header), comm);

   for (int i = 0; i < num_lists; ++i)
      size_idxlists
        += xt_idxlist_get_pack_size(collectionlist->idxlists[i], comm);

   return (size_t)size_header + size_idxlists;
}

static void
idxlist_collection_pack(Xt_idxlist data, void *buffer, int buffer_size,
                        int *position, MPI_Comm comm) {

   Xt_idxlist_collection collectionlist = (Xt_idxlist_collection)data;
   int num_lists = collectionlist->num_idxlists;
   int header[2] = { COLLECTION, num_lists };

   xt_mpi_call(MPI_Pack(header, 2, MPI_INT, buffer,
                        buffer_size, position, comm), comm);

   for (int i = 0; i < num_lists; ++i)
      xt_idxlist_pack(collectionlist->idxlists[i], buffer, buffer_size,
                      position, comm);
}

Xt_idxlist
xt_idxlist_collection_unpack(void *buffer, int buffer_size, int *position,
                             MPI_Comm comm) {

  int num_lists;
  xt_mpi_call(MPI_Unpack(buffer, buffer_size, position,
                         &num_lists, 1, MPI_INT, comm), comm);

  Xt_idxlist_collection collectionlist
    = xmalloc(sizeof (*collectionlist)
              + (size_t)num_lists * sizeof (collectionlist->idxlists[0]));

  collectionlist->index_array_cache = NULL;
  collectionlist->num_idxlists = num_lists;

  long long num_indices = 0;
  for (int i = 0; i < num_lists; ++i) {
    collectionlist->idxlists[i] = xt_idxlist_unpack(buffer, buffer_size,
                                                    position, comm);
    num_indices += xt_idxlist_get_num_indices(collectionlist->idxlists[i]);
  }

  assert(num_indices <= INT_MAX);
  Xt_idxlist_init(&collectionlist->parent, &idxlist_collection_vtable,
                  (int)num_indices);
  return (Xt_idxlist)collectionlist;
}

Xt_idxlist
xt_idxlist_collection_get_intersection(Xt_idxlist XT_UNUSED(idxlist_src),
                                       Xt_idxlist XT_UNUSED(idxlist_dst),
                                       Xt_config XT_UNUSED(config)) {

  return NULL;
}

static Xt_idxlist
idxlist_collection_copy(Xt_idxlist idxlist) {

   Xt_idxlist_collection collectionlist = (Xt_idxlist_collection)idxlist;

   return xt_idxlist_collection_new(collectionlist->idxlists,
                                    collectionlist->num_idxlists);
}

struct range_list_sort
{
  Xt_int min, max;
  int pos;
};

static int
range_list_cmp(const void *a, const void *b)
{
  const struct range_list_sort *ra = a, *rb = b;
  int ret;
  if (ra->min != rb->min) {
    ret = (ra->min > rb->min) - (ra->min < rb->min);
  } else {
    ret = (ra->max > rb->max) - (ra->max < rb->max);
  }
  return ret;
}

struct Xt_stripes {
  int num_stripes;
  struct Xt_stripe *stripes;
};

static Xt_idxlist
coll_get_sorted_stripes(Xt_idxlist_collection collectionlist,
                        struct range_list_sort *list_sorter,
                        Xt_config config)
{
  /* accumulate stripes into this buffer */
  struct Xt_stripe *stripes_accum = NULL;
  size_t stripes_accum_array_size = 0;
  size_t num_accum_stripes;
  const Xt_idxlist *unsorted_idxlists = collectionlist->idxlists;
  {
    int num_accum_stripes_;
    xt_idxlist_get_index_stripes(unsorted_idxlists[list_sorter[0].pos], &stripes_accum,
                                 &num_accum_stripes_);
    num_accum_stripes = (size_t)num_accum_stripes_;
  }
  size_t num_idxlists = (size_t)collectionlist->num_idxlists;
  struct Xt_stripe *stripe_buf = NULL; size_t stripe_buf_size = 0;
  /* obtain stripes, locally mostly sorted if reasonably available */
  for (size_t i = 1; i < num_idxlists; ++i) {
    Xt_idxlist component_list, unsorted_list
      = unsorted_idxlists[list_sorter[i].pos];
    Xt_int *sorted_indices = NULL;
    switch (unsorted_list->vtable->idxlist_pack_code) {
    case EMPTY:                 /* xt_idxempty should not usually be part
                                 * of a collection, but better safe
                                 * than sorry, I guess */
      Xt_abort(Xt_default_comm, "internal error", filename, __LINE__);
      /* these are either too expensive or pointless to sort, just
       * use their ranges as is */
    case COLLECTION:
    case STRIPES:
      component_list = unsorted_list;
      break;
    case VECTOR:
      {
        struct Xt_config_ derived_config = *config;
        derived_config.flags |= xt_force_nosort;
        size_t num_indices
          = (size_t)xt_idxlist_get_num_indices(unsorted_list);
        const Xt_int *indices
          = xt_idxvec_get_sorted_vector(unsorted_list, &derived_config);
        if (!indices) {
          const Xt_int *vector
            = xt_idxlist_get_indices_const(unsorted_list);
           indices = sorted_indices
             = xmalloc(num_indices * sizeof (*sorted_indices));
           for (size_t j = 0; j < num_indices; ++j)
             sorted_indices[j] = vector[j];
          config->sort_funcs->sort_xt_int(sorted_indices, num_indices);
        }
        component_list
          = xt_idxvec_prealloc_new(indices, (int)num_indices);
      }
      break;
    default:
      component_list
        = xt_idxlist_sorted_copy_custom(unsorted_list, config);
    }
    size_t num_stripes
      = (size_t)(xt_idxlist_get_num_index_stripes(component_list));
    if (num_stripes > stripe_buf_size) {
      stripe_buf = xrealloc(stripe_buf, num_stripes * sizeof (stripe_buf[0]));
      stripe_buf_size = num_stripes;
    }
    xt_idxlist_get_index_stripes_keep_buf(component_list,
                                          stripe_buf, stripe_buf_size);
    if (component_list != unsorted_list) {
      xt_idxlist_delete(component_list);
      free(sorted_indices);
    }
    ENSURE_ARRAY_SIZE(stripes_accum, stripes_accum_array_size,
                      num_accum_stripes + (size_t)stripe_buf_size);

    num_accum_stripes
      += xt_stripes_merge_copy((size_t)stripe_buf_size,
                               stripes_accum + num_accum_stripes,
                               stripe_buf,
                               num_accum_stripes > 0);

  }
  /* sort/shuffle the entirety of stripes */
  free(stripe_buf);
  Xt_idxlist sorted_copy
    = xt_idxstripes_sort_new((size_t)num_accum_stripes, stripes_accum, config);
  free(stripes_accum);
  return sorted_copy;
}


static Xt_idxlist
idxlist_collection_sorted_copy(Xt_idxlist idxlist, Xt_config config) {

   Xt_idxlist_collection collectionlist = (Xt_idxlist_collection)idxlist;
   size_t num_idxlists = (size_t)collectionlist->num_idxlists;
   Xt_int prev_max = XT_INT_MIN;
   Xt_idxlist *idxlists = collectionlist->idxlists;
   struct range_list_sort *list_sorter;
   Xt_idxlist_collection sorted_collection
     = xmalloc(sizeof (*sorted_collection)
               + num_idxlists * sizeof(sorted_collection->idxlists[0])
               + num_idxlists * sizeof (struct range_list_sort));
   sorted_collection->num_idxlists = (int)num_idxlists;
   sorted_collection->index_array_cache = NULL;
   list_sorter
     = (struct range_list_sort *)(sorted_collection->idxlists+num_idxlists);
   Xt_idxlist result;
   bool component_overlap = false;
   for (size_t i = 0; i < num_idxlists; ++i) {
     Xt_int min = list_sorter[i].min = xt_idxlist_get_min_index(idxlists[i]);
     component_overlap |= (min < prev_max);
     prev_max = list_sorter[i].max = xt_idxlist_get_max_index(idxlists[i]);
     list_sorter[i].pos = (int)i;
   }
   if (component_overlap) {
     /* individual lists not sorted or even overlapping each other */
     /* 1. sort the lists according to minimal element */
     qsort(list_sorter, num_idxlists, sizeof (*list_sorter),
           range_list_cmp);
     /* 2. re-check overlap after sorting */
     prev_max = XT_INT_MIN;
     component_overlap = false;
     for (size_t i = 0; i < num_idxlists; ++i) {
       Xt_int min = list_sorter[i].min;
       component_overlap |= (min < prev_max);
       prev_max = list_sorter[i].max;
     }
   }
   if (!component_overlap) {
     /* the component lists are in a sequence either naturally or
      * after sorting without problematic overlaps, a concatenation of
      * their sorted copies is all that's needed. */
     Xt_idxlist *sorted_lists = sorted_collection->idxlists;
     for (size_t i = 0; i < num_idxlists; ++i)
       sorted_lists[i]
         /* list_sorter[i].pos will be an identity map if the lists
          * were already sorted */
         = xt_idxlist_sorted_copy_custom(idxlists[list_sorter[i].pos],
                                         config);
     sorted_collection
       = xrealloc(sorted_collection, sizeof (*sorted_collection)
                  + num_idxlists * sizeof(sorted_collection->idxlists[0]));
     Xt_idxlist_init(&sorted_collection->parent, &idxlist_collection_vtable,
                     xt_idxlist_get_num_indices(idxlist));
     result = (Xt_idxlist)sorted_collection;
   } else {
     /* components overlap in a way that cannot be resolved by sorting
      * only the lists, need to create stripes and sort those */
     result = coll_get_sorted_stripes(collectionlist, list_sorter, config);
     free(sorted_collection);
   }
   return result;
}

static void
idxlist_collection_get_indices(Xt_idxlist idxlist, Xt_int *indices) {

   Xt_idxlist_collection collectionlist = (Xt_idxlist_collection)idxlist;
   /// \todo use memcpy with index_array_cache if available
   int offlist = 0, num_lists = collectionlist->num_idxlists;

   for (int i = 0; i < num_lists; ++i) {

      xt_idxlist_get_indices(collectionlist->idxlists[i], indices+offlist);
      offlist += xt_idxlist_get_num_indices(collectionlist->idxlists[i]);
   }
}

static const Xt_int *
idxlist_collection_get_indices_const(Xt_idxlist idxlist) {

  Xt_idxlist_collection collection = (Xt_idxlist_collection)idxlist;

  if (collection->index_array_cache) return collection->index_array_cache;

  unsigned num_indices = (unsigned)xt_idxlist_get_num_indices(idxlist);

  Xt_int *tmp_index_array
    = xmalloc(num_indices * sizeof (collection->index_array_cache[0]));

  idxlist_collection_get_indices(idxlist, tmp_index_array);

  collection->index_array_cache = tmp_index_array;

  return collection->index_array_cache;
}

static int
idxlist_collection_get_num_index_stripes(Xt_idxlist idxlist)
{
  Xt_idxlist_collection collectionlist = (Xt_idxlist_collection)idxlist;
  int num_lists = collectionlist->num_idxlists;
  size_t num_stripes = 0;
  for (int i = 0; i < num_lists; ++i)
    num_stripes
      += (size_t)xt_idxlist_get_num_index_stripes(collectionlist->idxlists[i]);
  assert(num_stripes <= INT_MAX);
  return (int)num_stripes;
}

static void
idxlist_collection_get_index_stripes(Xt_idxlist idxlist,
                                     struct Xt_stripe *restrict stripes,
                                     size_t num_stripes_alloc) {

  Xt_idxlist_collection collectionlist = (Xt_idxlist_collection)idxlist;

  size_t num_stripes = 0;

  size_t num_lists = (size_t)(collectionlist->num_idxlists);
  Xt_idxlist *idxlists = collectionlist->idxlists;
  for (size_t i = 0; i < num_lists; ++i) {
    size_t num_stripes_of_list
      = (size_t)(xt_idxlist_get_num_index_stripes(idxlists[i]));
    assert(num_stripes_alloc - num_stripes >= num_stripes_of_list);
    xt_idxlist_get_index_stripes_keep_buf(idxlists[i], stripes + num_stripes,
                                          num_stripes_alloc - num_stripes);
    num_stripes += num_stripes_of_list;
  }
}

static int
idxlist_collection_get_index_at_position(Xt_idxlist idxlist, int position,
                                         Xt_int * index) {

  Xt_idxlist_collection collectionlist = (Xt_idxlist_collection)idxlist;
  int num_lists = collectionlist->num_idxlists;

  for (int i = 0; i < num_lists; ++i) {
    int n = xt_idxlist_get_num_indices(collectionlist->idxlists[i]);
    if (position >= n)
      position -= n;
    else {
      return xt_idxlist_get_index_at_position(collectionlist->idxlists[i],
                                              position, index);
    }
  }
  return 1;

}

static int
idxlist_collection_get_position_of_index_off(Xt_idxlist idxlist, Xt_int index,
                                             int * position, int offset) {

  Xt_idxlist_collection collectionlist = (Xt_idxlist_collection)idxlist;

  int curr_num_indices = 0;

  int idxlist_offsets = 0;

  assert(offset >= 0);

  int i = 0, num_lists = collectionlist->num_idxlists;

  do {
    idxlist_offsets += curr_num_indices;
    curr_num_indices = xt_idxlist_get_num_indices(collectionlist->idxlists[i]);
  } while (idxlist_offsets + curr_num_indices <= offset && ++i < num_lists);

  offset -= idxlist_offsets;

  for (;i < num_lists; ++i)
    if (!xt_idxlist_get_position_of_index_off(collectionlist->idxlists[i],
                                              index, position, offset)) {
      *position += idxlist_offsets;
      return 0;
    } else {
      idxlist_offsets
        += xt_idxlist_get_num_indices(collectionlist->idxlists[i]);
      offset = 0;
    }

  return 1;
}

static int
idxlist_collection_get_position_of_index(Xt_idxlist idxlist, Xt_int index,
                                         int * position) {

  return idxlist_collection_get_position_of_index_off(idxlist, index,
                                                      position, 0);
}

static Xt_int
idxlist_collection_get_min_index(Xt_idxlist idxlist) {

  Xt_idxlist_collection collectionlist = (Xt_idxlist_collection)idxlist;

  size_t num_lists = (size_t)collectionlist->num_idxlists;
  assert(collectionlist->num_idxlists > 0);

  Xt_int tmp_min, min = XT_INT_MAX;

  for (size_t i = 0; i < num_lists; ++i)
    if ((tmp_min = xt_idxlist_get_min_index(collectionlist->idxlists[i])) < min)
      min = tmp_min;

  return min;
}

static Xt_int
idxlist_collection_get_max_index(Xt_idxlist idxlist) {

  Xt_idxlist_collection collectionlist = (Xt_idxlist_collection)idxlist;

  size_t num_lists = (size_t)collectionlist->num_idxlists;
  assert(collectionlist->num_idxlists > 0);

  Xt_int tmp_max, max = XT_INT_MIN;

  for (size_t i = 0; i < num_lists; ++i)
    if ((tmp_max = xt_idxlist_get_max_index(collectionlist->idxlists[i])) > max)
      max = tmp_max;

  return max;
}

static int
idxlist_collection_get_sorting(Xt_idxlist idxlist)
{
  unsigned sort_flags = (((Xt_idxlist_collection)idxlist)->flags) & sort_mask;
  return (int)sort_flags-(sort_flags < 3);
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
