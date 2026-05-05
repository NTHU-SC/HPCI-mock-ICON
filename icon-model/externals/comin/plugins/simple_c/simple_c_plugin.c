/*
  Example plugin for the ICON Community Interface (ComIn)
  with basic (not MPI-parallel) callbacks and accessing variables and
  descriptive data structures.

  Note that in order to demonstrate ComIn's language interoperability,
  a similary plugin has been implemented in FORTRAN, see the subdirectory
  "simple".

  @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

  SPDX-License-Identifier: BSD-3-Clause

  Please see the file LICENSE in the root of the source tree for this code.
  Where software is supplied by third parties, it is indicated in the
  headers of the routines.
*/

#include <comin.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// points to the pointer to the buffer (buffer swapping)
void *pres, *simple_c_var, *simple_c_tracer;

void simple_c_constructor() {
  unsigned int major, minor, patch;
  comin_setup_get_version(&major, &minor, &patch);
  comin_print_info_f("ComIn v%u.%u.%u simple_c_constructor called!", major,
                     minor, patch);

  t_comin_var_descriptor desc;
  for (void *it = comin_var_get_descr_list_head(); it != NULL;
       it       = comin_var_get_descr_list_next(it)) {
    comin_var_get_descr_list_var_desc(it, &desc);
    comin_print_info_f("found variable: %s on domain id %d", desc.name,
                       desc.id);
  }

  // test explicit free'ing of a var descriptor list iterator:
  {
    t_comin_var_descr_list_iterator *descr_it = comin_var_get_descr_list_head();
    descr_it = comin_var_get_descr_list_next(descr_it);
    comin_var_descr_list_iterator_delete(&descr_it);
  }

  t_comin_var_descriptor pres_d     = {.name = "pres", .id = 1};
  t_comin_entry_point before_output = EP_ATM_WRITE_OUTPUT_BEFORE;
  pres = comin_var_get(1, &before_output, pres_d, COMIN_FLAG_READ);
  if (pres == NULL) {
    comin_plugin_finish("simple_c_plugin", "Internal error!");
  }

  t_comin_var_descriptor var_d = {.name = "simple_c_var", .id = 1};
  simple_c_var = comin_var_get(1, &before_output, var_d, COMIN_FLAG_WRITE);
  if (simple_c_var == NULL) {
    comin_plugin_finish("simple_c_plugin", "Internal error!");
  }

  t_comin_var_descriptor tracer_d = {.name = "simple_c_tracer", .id = 1};
  simple_c_tracer =
      comin_var_get(1, &before_output, tracer_d, COMIN_FLAG_WRITE);
  if (simple_c_tracer == NULL) {
    comin_plugin_finish("simple_c_plugin", "Internal error!");
  }

  comin_print_info("Metadata:");
  // iterator over the metadata:
  void *end_it = comin_metadata_get_iterator_end(pres_d);
  void *it     = comin_metadata_get_iterator_begin(pres_d);
  for (; !comin_metadata_iterator_compare(it, end_it);
       comin_metadata_iterator_next(it)) {
    const char *key = comin_metadata_iterator_get_key(it);
    int type        = comin_metadata_get_typeid(pres_d, key);
    comin_print_info_f("%s type:%i", key, type);
    if (type == COMIN_METADATA_TYPEID_INTEGER) {
      int val = -999;
      comin_metadata_get_integer(pres_d, key, &val);
      comin_print_info_f(" value: %i", val);
    }
  }
  comin_metadata_iterator_delete(it);
  comin_metadata_iterator_delete(end_it);

  // test "dim_semantics" auxiliary function, used for
  // interpretation of array dimensions:
  int dim_semantics[5];
  comin_var_get_dim_semantics(pres, dim_semantics);
  if (dim_semantics[0] != COMIN_DIM_SEMANTICS_NPROMA) {
    comin_plugin_finish("simple_c_plugin", "Internal error!");
  }
  if (dim_semantics[1] != COMIN_DIM_SEMANTICS_LEVEL) {
    comin_plugin_finish("simple_c_plugin", "Internal error!");
  }
  if (dim_semantics[2] != COMIN_DIM_SEMANTICS_BLOCK) {
    comin_plugin_finish("simple_c_plugin", "Internal error!");
  }
}

void simple_c_diagfct() {
  int jg;
  comin_print_info("simple_c_diagfct called!");
  jg = comin_current_get_domain_id();
  comin_print_info_f("currently on domain %i", jg);
  int pres_shape[5], tracer_shape[5];
  comin_var_get_shape(pres, pres_shape);
  double *simple_c_var_data = comin_var_get_ptr_double(simple_c_var);
  double *pres_data         = comin_var_get_ptr_double(pres);
  for (int i = 0; i < pres_shape[0] * pres_shape[1] * pres_shape[2] *
                          pres_shape[3] * pres_shape[4];
       ++i) {
    simple_c_var_data[i] = pres_data[i] + 42.;
  }
  comin_var_get_shape(simple_c_tracer, tracer_shape);
  double *simple_c_tracer_data = comin_var_get_ptr_double(simple_c_tracer);
  for (int i = 0; i < tracer_shape[0] * tracer_shape[1] * tracer_shape[2];
       ++i) {
    simple_c_tracer_data[i] /= 1337.;
  }
}

void simple_c_destructor() { comin_print_info("simple_c_destructor called!"); }

void comin_main() {
  int ilen                = -1;
  const char *plugin_name = NULL;
  comin_current_get_plugin_name(&plugin_name, &ilen);

  int plugin_id = comin_current_get_plugin_id();
  comin_print_info_f("plugin %.*s has id %d", ilen, plugin_name, plugin_id);
  t_comin_var_descriptor simple_var_d = {.name = "simple_c_var", .id = 1};
  comin_var_request_add(simple_var_d, true);

  comin_metadata_set_integer(simple_var_d, "zaxis_id", COMIN_ZAXIS_3D);
  comin_metadata_set_logical(simple_var_d, "restart", false);
  comin_metadata_set_logical(simple_var_d, "tracer", false);
  comin_metadata_set_integer(simple_var_d, "tracer_vlimit", 0);
  comin_metadata_set_integer(simple_var_d, "tracer_hlimit", 0);

  t_comin_var_descriptor simple_tracer_d = {.name = "simple_c_tracer",
                                            .id   = -1};
  comin_var_request_add(simple_tracer_d, false);
  comin_metadata_set_integer(simple_tracer_d, "zaxis_id", COMIN_ZAXIS_3D);
  comin_metadata_set_logical(simple_tracer_d, "restart", false);
  comin_metadata_set_logical(simple_tracer_d, "tracer", true);
  comin_metadata_set_integer(simple_tracer_d, "tracer_vlimit", 0);
  comin_metadata_set_integer(simple_tracer_d, "tracer_hlimit", 0);

  comin_callback_register(EP_SECONDARY_CONSTRUCTOR, &simple_c_constructor);
  comin_callback_register(EP_ATM_WRITE_OUTPUT_BEFORE, &simple_c_diagfct);
  comin_callback_register(EP_DESTRUCTOR, &simple_c_destructor);

  char ep_name[COMIN_MAX_LEN_EP_NAME + 1];
  comin_callback_get_ep_name(EP_DESTRUCTOR, ep_name);
  if (strncmp(ep_name, "EP_DESTRUCTOR", (size_t)13)) {
    char output_text[255];
    sprintf(output_text, "Expected EP_DESTRUCTOR; got |%s|\n", ep_name);
    comin_plugin_finish("simple_c: comin_main", output_text);
  }
  // TODO: Access descriptive data structures

  /* access to comin descriptive global data (exemplary)*/
  int n_dom = comin_descrdata_get_global_n_dom();
  comin_print_info_f("n_dom: %d", n_dom);
  int max_dom = comin_descrdata_get_global_max_dom();
  comin_print_info_f("max_dom: %d", max_dom);
  int nproma = comin_descrdata_get_global_nproma();
  comin_print_info_f("nproma: %d", nproma);
  int min_rlcell_int = comin_descrdata_get_global_min_rlcell_int();
  comin_print_info_f("min_rlcell_int: %d", min_rlcell_int);
}
