/**
   @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

   SPDX-License-Identifier: BSD-3-Clause

   Please see the file LICENSE in the root of the source tree for this code.
   Where software is supplied by third parties, it is indicated in the
   headers of the routines. **/

#include "comin.h"
#include <stdio.h>

void comin_main() {
  /* Test local / global index mappings:

     We get the global cell index of this MPI rank's cell with
     local index 42... and then re-translate this index to the local
     index.
  */
  int local_idx = 42;
  int jg        = 1;
  const int* glb;
  int arr_size[1];
  comin_descrdata_get_domain_cells_glb_index(jg, &glb, arr_size);
  if (arr_size[0] > local_idx) {
    int glb_index = glb[local_idx];
    int lcl       = comin_descrdata_index_lookup_glb2loc_cell(1, glb_index) -
              1; // subtract one due to fortran indexing
    comin_print_info_f("global index: %d", glb_index);
    comin_print_info_f("local index: %d (reference)    %d (looked up)",
                       local_idx, lcl);
    if (lcl != local_idx)
      comin_plugin_finish("test_index_mapping", "Check failed.");
  } else {
    comin_plugin_finish("test_index_mapping", "insufficient local cells.");
  }
}
