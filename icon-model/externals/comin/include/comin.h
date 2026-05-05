/**
   @file comin.h
   @brief C interface for the ICON Community Interface

   @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

   SPDX-License-Identifier: BSD-3-Clause

   Please see the file LICENSE in the root of the source tree for this code.
   Where software is supplied by third parties, it is indicated in the
   headers of the routines. **/

#ifndef COMIN_H
#define COMIN_H

#include "comin_global.inc"
#include "comin_version.inc"
#include <stdarg.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

#ifdef DOXYGEN
// This is a dirty workaround to make C symbols referable if there
// is a Fortran symbol with the same name
namespace comin_h {
#endif

#ifdef __cplusplus
extern "C" {
#endif

  /// @defgroup c_interface C Interface
  /// @{

  typedef struct t_comin_var_descriptor {
    char name[COMIN_MAX_LEN_VAR_NAME + 1];
    int id;
  } t_comin_var_descriptor;

  //< Opaque data type for a ComIn variable handle.
  typedef struct t_comin_var_handle t_comin_var_handle;

  //< Opaque data type for a variable list iterator.
  typedef struct t_comin_var_descr_list_iterator
      t_comin_var_descr_list_iterator;

  //< Opaque data type for a metadata iterator.
  typedef struct t_comin_metadata_iterator t_comin_metadata_iterator;

  typedef enum t_comin_entry_point {
    EP_SECONDARY_CONSTRUCTOR = 1,
    EP_ATM_YAC_DEFCOMP_BEFORE,
    EP_ATM_YAC_DEFCOMP_AFTER,
    EP_ATM_YAC_SYNCDEF_BEFORE,
    EP_ATM_YAC_SYNCDEF_AFTER,
    EP_ATM_YAC_ENDDEF_BEFORE,
    EP_ATM_YAC_ENDDEF_AFTER,
    EP_ATM_INIT_FINALIZE,
    EP_ATM_TIMELOOP_BEFORE,
    EP_ATM_TIMELOOP_START,
    EP_ATM_TIMELOOP_END,
    EP_ATM_TIMELOOP_AFTER,
    EP_ATM_INTEGRATE_BEFORE,
    EP_ATM_INTEGRATE_START,
    EP_ATM_INTEGRATE_END,
    EP_ATM_INTEGRATE_AFTER,
    EP_ATM_WRITE_OUTPUT_BEFORE,
    EP_ATM_WRITE_OUTPUT_AFTER,
    EP_ATM_CHECKPOINT_BEFORE,
    EP_ATM_CHECKPOINT_AFTER,
    EP_ATM_ADVECTION_BEFORE,
    EP_ATM_ADVECTION_AFTER,
    EP_ATM_PHYSICS_BEFORE,
    EP_ATM_PHYSICS_AFTER,
    EP_ATM_NUDGING_BEFORE,
    EP_ATM_NUDGING_AFTER,
    EP_ATM_SURFACE_BEFORE,
    EP_ATM_SURFACE_AFTER,
    EP_ATM_TURBULENCE_BEFORE,
    EP_ATM_TURBULENCE_AFTER,
    EP_ATM_MICROPHYSICS_BEFORE,
    EP_ATM_MICROPHYSICS_AFTER,
    EP_ATM_CONVECTION_BEFORE,
    EP_ATM_CONVECTION_AFTER,
    EP_ATM_RADIATION_BEFORE,
    EP_ATM_RADIATION_AFTER,
    EP_ATM_RADHEAT_BEFORE,
    EP_ATM_RADHEAT_AFTER,
    EP_ATM_GWDRAG_BEFORE,
    EP_ATM_GWDRAG_AFTER,
    EP_FINISH,
    EP_DESTRUCTOR
  } t_comin_entry_point;

  typedef enum t_comin_zaxis {
    COMIN_ZAXIS_UNDEF   = -1,
    COMIN_ZAXIS_NONE    = 0,
    COMIN_ZAXIS_2D      = 1,
    COMIN_ZAXIS_3D      = 2,
    COMIN_ZAXIS_3D_HALF = 3,
  } t_comin_zaxis;

  typedef enum t_comin_var_access_flag {
    COMIN_FLAG_NONE      = 0,
    COMIN_FLAG_READ      = 1 << 1,
    COMIN_FLAG_WRITE     = 1 << 2,
    COMIN_FLAG_SYNC_HALO = 1 << 3,
    COMIN_FLAG_DEVICE    = 1 << 4,
  } t_comin_var_access_flag;

  typedef enum t_comin_error_code {
    COMIN_SUCCESS = 0,
    COMIN_INFO,
    COMIN_WARNING,
    COMIN_ERROR_STATUS,
    COMIN_ERROR_CALLBACK_REGISTER_OUTSIDE_PRIMARYCONSTRUCTOR,
    COMIN_ERROR_CALLBACK_COMPLETE,
    COMIN_ERROR_CALLBACK_EP_ID_UNKNOWN,
    COMIN_ERROR_DESCRDATA_SET_FCT_GLB2LOC,
    COMIN_ERROR_DESCRDATA_FINALIZE,
    COMIN_ERROR_METADATA_SET_OUTSIDE_PRIMARYCONSTRUCTOR,
    COMIN_ERROR_METADATA_KEY_NOT_FOUND,
    COMIN_ERROR_METADATA_GET_INSIDE_PRIMARYCONSTRUCTOR,
    COMIN_ERROR_SETUP_FINALIZE,
    COMIN_ERROR_SETUP_COMIN_ALREADY_INITIALIZED,
    COMIN_ERROR_PLUGIN_INIT_COMIN_VERSION,
    COMIN_ERROR_PLUGIN_INIT_PRECISION,
    COMIN_ERROR_PLUGIN_INIT_STATE_INITIALIZED,
    COMIN_ERROR_SETUP_ERRHANDLER_NOT_ASSOCIATED,
    COMIN_ERROR_SETUP_ERRHANDLER_NOT_SET,
    COMIN_ERROR_SETUP_PRECISION_TEST_FAILED,
    COMIN_ERROR_VAR_REQUEST_AFTER_PRIMARYCONSTRUCTOR,
    COMIN_ERROR_VAR_REQUEST_EXISTS_IS_LMODEXCLUSIVE,
    COMIN_ERROR_VAR_REQUEST_EXISTS_REQUEST_LMODEXCLUSIVE,
    COMIN_ERROR_VAR_DESCRIPTOR_NOT_FOUND,
    COMIN_ERROR_VAR_ITEM_NOT_ASSOCIATED,
    COMIN_ERROR_FIELD_NOT_ALLOCATED,
    COMIN_ERROR_POINTER_NOT_ASSOCIATED,
    COMIN_ERROR_TRACER_REQUEST_NOT_FOR_ALL_DOMAINS,
    COMIN_ERROR_VAR_SYNC_DEVICE_MEM_NOT_ASSOCIATED,
    COMIN_ERROR_VAR_SYNC_HALO_NOT_ASSOCIATED,
    COMIN_ERROR_VAR_GET_OUTSIDE_SECONDARY_CONSTRUCTOR,
    COMIN_ERROR_VAR_GET_NO_DEVICE,
    COMIN_ERROR_VAR_GET_VARIABLE_NOT_FOUND,
    COMIN_ERROR_VAR_GET_CONTAINER_CAN_NOT_HALO_SYNCHRONIZED,
    COMIN_ERROR_VAR_GET_IRREGULAR_VAR_CAN_NOT_HALO_SYNCHRONIZED,
    COMIN_ERROR_VAR_SYNC_HALO_NOT_SUPPORTED_ZAXIS,
    COMIN_ERROR_VAR_METADATA_INCONSISTENT_TYPE,
    COMIN_ERROR_FATAL,
  } t_comin_error_code;

  typedef enum t_comin_hgrid_id {
    COMIN_HGRID_UNSTRUCTURED_CELL   = 1,
    COMIN_HGRID_UNSTRUCTURED_EDGE   = 2,
    COMIN_HGRID_UNSTRUCTURED_VERTEX = 3
  } t_comin_hgrid_id;

  typedef enum t_comin_metadata_typeid {
    COMIN_METADATA_TYPEID_UNDEFINED = -1,
    COMIN_METADATA_TYPEID_INTEGER,
    COMIN_METADATA_TYPEID_REAL,
    COMIN_METADATA_TYPEID_CHARACTER,
    COMIN_METADATA_TYPEID_LOGICAL
  } t_comin_metadata_typeid;

  typedef enum t_comin_var_datatype {
    COMIN_VAR_DATATYPE_DOUBLE = 1,
    COMIN_VAR_DATATYPE_FLOAT,
    COMIN_VAR_DATATYPE_INT,
  } t_comin_var_datatype;

  typedef enum t_comin_dim_semantics {
    COMIN_DIM_SEMANTICS_UNDEF     = 1, // not explicitly defined
    COMIN_DIM_SEMANTICS_NPROMA    = 2, // nproma
    COMIN_DIM_SEMANTICS_BLOCK     = 3, // blocked data layout
    COMIN_DIM_SEMANTICS_UNBLOCK   = 4, // unblocked data layout
    COMIN_DIM_SEMANTICS_LEVEL     = 5, // vertical axis dimension
    COMIN_DIM_SEMANTICS_CONTAINER = 6, // container
    COMIN_DIM_SEMANTICS_OTHER     = 7, // special cases (tracers, tiles, etc)
    COMIN_DIM_SEMANTICS_UNUSED    = 8  // unused dimension
  } t_comin_dim_semantics;

  static const int COMIN_DOMAIN_OUTSIDE_LOOP = -1;

  t_comin_entry_point comin_current_get_ep();
  int comin_current_get_domain_id();
  int comin_current_get_plugin_id();
  void comin_current_get_plugin_name(char const** val, int* len);
  void comin_current_get_plugin_options(char const** val, int* len);
  void comin_current_get_plugin_comm(char const** val, int* len);
  void comin_current_get_datetime(char const** val, int* len);

  int comin_parallel_get_plugin_mpi_comm();
  int comin_parallel_get_host_mpi_comm();
  int comin_parallel_get_host_mpi_rank();

  void comin_plugin_finish(const char* routine, const char* text);
  void comin_print_debug(const char* msg);
  void comin_print_info(const char* msg);
  void comin_print_warning(const char* msg);

  void comin_print_debug_f(const char* fmt, ...);
  void comin_print_info_f(const char* fmt, ...);
  void comin_print_warning_f(const char* fmt, ...);

  void comin_error_get_message(t_comin_error_code error_code, char category[11],
                               char message[COMIN_MAX_LEN_ERR_MESSAGE]);
  void comin_error_check(t_comin_error_code error_code, const char* scope);
  void comin_error_set_errors_return(bool errors_return);
  t_comin_error_code comin_error_get();
  void comin_error_reset();

  void comin_var_request_add(t_comin_var_descriptor var_descriptor,
                             bool lmodexclusive);
  t_comin_var_handle* comin_var_get(int context_len,
                                    t_comin_entry_point* context,
                                    t_comin_var_descriptor var_descriptor,
                                    int flag);
  void* comin_var_get_ptr(t_comin_var_handle* handle);
  double* comin_var_get_ptr_double(t_comin_var_handle* handle);
  float* comin_var_get_ptr_float(t_comin_var_handle* handle);
  int* comin_var_get_ptr_int(t_comin_var_handle* handle);
  void* comin_var_get_device_ptr(t_comin_var_handle* handle);
  double* comin_var_get_device_ptr_double(t_comin_var_handle* handle);
  float* comin_var_get_device_ptr_float(t_comin_var_handle* handle);
  int* comin_var_get_device_ptr_int(t_comin_var_handle* handle);
  void comin_var_get_shape(t_comin_var_handle* handle, int shape[5]);
  void comin_var_get_dim_semantics(t_comin_var_handle* handle,
                                   int dim_semantics[5]);
  void comin_var_get_ncontained(t_comin_var_handle* handle, int* ncontained);
  void comin_var_get_descriptor(t_comin_var_handle* handle,
                                t_comin_var_descriptor* descr);

  t_comin_var_descr_list_iterator* comin_var_get_descr_list_head();
  t_comin_var_descr_list_iterator*
  comin_var_get_descr_list_next(t_comin_var_descr_list_iterator* current);
  void
  comin_var_get_descr_list_var_desc(t_comin_var_descr_list_iterator* current,
                                    t_comin_var_descriptor* var_desc);
  void
  comin_var_descr_list_iterator_delete(t_comin_var_descr_list_iterator** it);

  typedef void (*t_comin_callback_function)();
  void comin_callback_register(t_comin_entry_point entry_point_id,
                               t_comin_callback_function fct_ptr);
  void comin_callback_get_ep_name(t_comin_entry_point iep,
                                  char out_ep_name[COMIN_MAX_LEN_EP_NAME + 1]);

  t_comin_metadata_typeid
  comin_metadata_get_typeid(t_comin_var_descriptor var_descriptor,
                            const char* key);
  void comin_metadata_set_integer(t_comin_var_descriptor var_descriptor,
                                  const char* key, int val);
  void comin_metadata_set_logical(t_comin_var_descriptor var_descriptor,
                                  const char* key, bool val);
  void comin_metadata_set_real(t_comin_var_descriptor var_descriptor,
                               const char* key, double val);
  void comin_metadata_set_character(t_comin_var_descriptor var_descriptor,
                                    const char* key, char const* val);
  void comin_metadata_get_integer(t_comin_var_descriptor var_descriptor,
                                  const char* key, int* val);
  void comin_metadata_get_logical(t_comin_var_descriptor var_descriptor,
                                  const char* key, bool* val);
  void comin_metadata_get_real(t_comin_var_descriptor var_descriptor,
                               const char* key, double* val);
  void comin_metadata_get_character(t_comin_var_descriptor var_descriptor,
                                    const char* key, char const** val,
                                    int* len);

  t_comin_metadata_iterator*
  comin_metadata_get_iterator_begin(t_comin_var_descriptor var_descriptor);
  t_comin_metadata_iterator*
  comin_metadata_get_iterator_end(t_comin_var_descriptor var_descriptor);
  const char* comin_metadata_iterator_get_key(t_comin_metadata_iterator* it);
  bool comin_metadata_iterator_compare(t_comin_metadata_iterator* it1,
                                       t_comin_metadata_iterator* it2);
  void comin_metadata_iterator_next(t_comin_metadata_iterator* it);
  void comin_metadata_iterator_delete(t_comin_metadata_iterator* it);

  double comin_descrdata_get_timesteplength(int jg);
  int comin_descrdata_get_index(int j);
  int comin_descrdata_get_block(int j);
  void comin_descrdata_get_cell_indices(int jg, int i_blk, int i_startblk,
                                        int i_endblk, int* i_startidx,
                                        int* i_endidx, int irl_start,
                                        int irl_end);
  void comin_descrdata_get_edge_indices(int jg, int i_blk, int i_startblk,
                                        int i_endblk, int* i_startidx,
                                        int* i_endidx, int irl_start,
                                        int irl_end);
  void comin_descrdata_get_vert_indices(int jg, int i_blk, int i_startblk,
                                        int i_endblk, int* i_startidx,
                                        int* i_endidx, int irl_start,
                                        int irl_end);
  int comin_descrdata_get_cell_npromz(int jg);
  int comin_descrdata_get_edge_npromz(int jg);
  int comin_descrdata_get_vert_npromz(int jg);
  int comin_descrdata_index_lookup_glb2loc_cell(int jg, int global_idx);
  void comin_descrdata_get_simulation_interval_exp_start(char const** val,
                                                         int* len);
  void comin_descrdata_get_simulation_interval_exp_stop(char const** val,
                                                        int* len);
  void comin_descrdata_get_simulation_interval_run_start(char const** val,
                                                         int* len);
  void comin_descrdata_get_simulation_interval_run_stop(char const** val,
                                                        int* len);

  /// returns version info.
  static inline void comin_setup_get_version(unsigned int* major,
                                             unsigned int* minor,
                                             unsigned int* patch) {
    (*major) = COMIN_VERSION_MAJOR;
    (*minor) = COMIN_VERSION_MINOR;
    (*patch) = COMIN_VERSION_PATCH;
  }

  // ---------------------------------------------------------

  // The following **internal** struct is used to describe the
  // descrdata structure dynamically, mainly for use in the
  // python_adapter.
  struct comin_descrdata_property_t {
    const char* name;
    void* get_function;
    const char* datatype; // c datatype
    int ndims;
    bool has_jg;
    const struct comin_descrdata_property_t* subtypes;
  };

#ifdef __cplusplus
} // extern C
#endif

/* Header extension for get of grid data and domain routines generated by python
 * script (comin_header_c_ext_descrdata_get_domain.h.py) in utils/.  */
#include "comin_header_c_ext_descrdata_query_domain.h"

/* Header extension for get of global data routines generated by python script
 * (comin_header_c_ext_descrdata_get_global.h.py) in utils/.  */
#include "comin_header_c_ext_descrdata_query_global.h"

/// @}

#ifdef DOXYGEN
}; // namespace comin_h
#endif

#endif
