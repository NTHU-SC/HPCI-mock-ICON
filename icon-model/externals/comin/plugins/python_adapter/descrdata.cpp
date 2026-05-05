/* @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

   SPDX-License-Identifier: BSD-3-Clause

   Please see the file LICENSE in the root of the source tree for this code.
   Where software is supplied by third parties, it is indicated in the
   headers of the routines. */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "util.h"

#include "comin.h"
#include "descrdata_properties.h"
#include "exception.h"

static PyObject* py_comin_current_get_plugin_info(PyObject* /*self*/,
                                                  PyObject* /*args*/) {
  int plugin_id = comin_current_get_plugin_id();
  py_comin_check_error();
  int ilen_name           = -1;
  const char* plugin_name = NULL;
  comin_current_get_plugin_name(&plugin_name, &ilen_name);
  py_comin_check_error();

  int ilen_options           = -1;
  const char* plugin_options = NULL;
  comin_current_get_plugin_options(&plugin_options, &ilen_options);
  py_comin_check_error();

  int ilen_comm           = -1;
  const char* plugin_comm = NULL;
  comin_current_get_plugin_comm(&plugin_comm, &ilen_comm);
  py_comin_check_error();

  return Py_BuildValue("{siss#ss#ss#}", "id", plugin_id, "name", plugin_name,
                       (Py_ssize_t)ilen_name, "options", plugin_options,
                       (Py_ssize_t)ilen_options, "comm", plugin_comm,
                       (Py_ssize_t)ilen_comm);
}

static PyObject* py_comin_descrdata_get_domain(PyObject* /*self*/,
                                               PyObject* /*args*/) {
  return py_comin_descrdata_properties_to_dict(
      comin_descrdata_domain_properties);
}

static PyObject* py_comin_descrdata_get_global(PyObject* /*self*/,
                                               PyObject* /*args*/) {
  return py_comin_descrdata_properties_to_dict(
      comin_descrdata_global_properties);
}

static PyObject*
py_comin_descrdata_get_simulation_interval(PyObject* /*self*/,
                                           PyObject* /*args*/) {
  int ilen_exp_start    = -1;
  const char* exp_start = NULL;
  comin_descrdata_get_simulation_interval_exp_start(&exp_start,
                                                    &ilen_exp_start);
  py_comin_check_error();
  int ilen_exp_stop    = -1;
  const char* exp_stop = NULL;
  comin_descrdata_get_simulation_interval_exp_stop(&exp_stop, &ilen_exp_stop);
  py_comin_check_error();
  int ilen_run_start    = -1;
  const char* run_start = NULL;
  comin_descrdata_get_simulation_interval_run_start(&run_start,
                                                    &ilen_run_start);
  py_comin_check_error();
  int ilen_run_stop    = -1;
  const char* run_stop = NULL;
  comin_descrdata_get_simulation_interval_run_stop(&run_stop, &ilen_run_stop);
  py_comin_check_error();

  return Py_BuildValue("{ss#ss#ss#ss#}", "exp_start", exp_start,
                       (Py_ssize_t)ilen_exp_start, "exp_stop", exp_stop,
                       (Py_ssize_t)ilen_exp_stop, "run_start", run_start,
                       (Py_ssize_t)ilen_run_start, "run_stop", run_stop,
                       (Py_ssize_t)ilen_run_stop);
}

static PyObject* py_comin_descrdata_get_timesteplength(PyObject* /*self*/,
                                                       PyObject* args,
                                                       PyObject* kwargs) {
  int j;
  static char const* kwlist[] = {"j", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "i", (char**)&kwlist, &j))
    return NULL;
  double t = comin_descrdata_get_timesteplength(j);
  py_comin_check_error();
  return PyFloat_FromDouble(t);
}

static PyObject* py_comin_descrdata_get_index(PyObject* /*self*/,
                                              PyObject* args,
                                              PyObject* kwargs) {
  int j;
  static char const* kwlist[] = {"j", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "i", (char**)&kwlist, &j))
    return NULL;
  int i = comin_descrdata_get_index(j);
  py_comin_check_error();
  return PyLong_FromLong(i);
}

static PyObject* py_comin_descrdata_get_block(PyObject* /*self*/,
                                              PyObject* args,
                                              PyObject* kwargs) {
  int j;
  static char const* kwlist[] = {"j", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "i", (char**)&kwlist, &j))
    return NULL;
  int b = comin_descrdata_get_block(j);
  py_comin_check_error();
  return PyLong_FromLong(b);
}

static PyObject* py_comin_descrdata_get_cell_indices(PyObject* /*self*/,
                                                     PyObject* args,
                                                     PyObject* kwargs) {
  int jg, i_blk, i_startblk, i_endblk, irl_start, irl_end;
  static char const* kwlist[] = {
      "jg", "i_blk", "i_startblk", "i_endblk", "irl_start", "irl_end", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "iiiiii", (char**)&kwlist, &jg,
                                   &i_blk, &i_startblk, &i_endblk, &irl_start,
                                   &irl_end))
    return NULL;
  int i_startidx, i_endidx;
  comin_descrdata_get_cell_indices(jg, i_blk, i_startblk, i_endblk, &i_startidx,
                                   &i_endidx, irl_start, irl_end);
  py_comin_check_error();
  return Py_BuildValue("ii", i_startidx, i_endidx);
}

static PyObject* py_comin_descrdata_get_cell_npromz(PyObject* /*self*/,
                                                    PyObject* args,
                                                    PyObject* kwargs) {
  int jg;
  static char const* kwlist[] = {"jg", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "i", (char**)&kwlist, &jg))
    return NULL;
  int npz = comin_descrdata_get_cell_npromz(jg);
  py_comin_check_error();
  return PyLong_FromLong(npz);
}

static PyObject* py_comin_descrdata_get_edge_npromz(PyObject* /*self*/,
                                                    PyObject* args,
                                                    PyObject* kwargs) {
  int jg;
  static char const* kwlist[] = {"jg", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "i", (char**)&kwlist, &jg))
    return NULL;
  int npz = comin_descrdata_get_edge_npromz(jg);
  py_comin_check_error();
  return PyLong_FromLong(npz);
}

static PyObject* py_comin_descrdata_get_vert_npromz(PyObject* /*self*/,
                                                    PyObject* args,
                                                    PyObject* kwargs) {
  int jg;
  static char const* kwlist[] = {"jg", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "i", (char**)&kwlist, &jg))
    return NULL;
  int npz = comin_descrdata_get_vert_npromz(jg);
  py_comin_check_error();
  return PyLong_FromLong(npz);
}

static PyObject* py_comin_current_get_datetime(PyObject* /*self*/,
                                               PyObject* /*args*/) {
  int ilen_datetime_str    = -1;
  const char* datetime_str = NULL;
  comin_current_get_datetime(&datetime_str, &ilen_datetime_str);
  py_comin_check_error();
  return PyUnicode_FromStringAndSize(datetime_str, ilen_datetime_str);
}

static PyObject* py_comin_current_get_domain_id(PyObject* /*self*/,
                                                PyObject* /*args*/) {
  int jg = comin_current_get_domain_id();
  py_comin_check_error();
  if (jg < 0)
    return PyErr_Format(PyExc_RuntimeError,
                        "comin_current_get_domain_id failed. (jg=%d)", jg);
  return PyLong_FromLong(jg);
}

static PyObject*
py_comin_descrdata_index_lookup_glb2loc_cell(PyObject* /*self*/, PyObject* args,
                                             PyObject* kwargs) {
  int jg, global_idx;
  static char const* kwlist[] = {"jg", "global_idx", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ii", (char**)&kwlist, &jg,
                                   &global_idx))
    return NULL;
  int loc = comin_descrdata_index_lookup_glb2loc_cell(jg, global_idx);
  py_comin_check_error();
  return PyLong_FromLong(loc);
}

static PyObject* py_comin_setup_get_version(PyObject* /*self*/,
                                            PyObject* /*args*/) {
  unsigned int major, minor, patch;
  comin_setup_get_version(&major, &minor, &patch);
  py_comin_check_error();
  return Py_BuildValue("iii", major, minor, patch);
}

static PyObject* py_comin_parallel_get_host_mpi_comm(PyObject* /*self*/,
                                                     PyObject* /*args*/) {
  int comm = comin_parallel_get_host_mpi_comm();
  py_comin_check_error();
  return PyLong_FromLong((long)comm);
}
static PyObject* py_comin_parallel_get_plugin_mpi_comm(PyObject* /*self*/,
                                                       PyObject* /*args*/) {
  int comm = comin_parallel_get_plugin_mpi_comm();
  py_comin_check_error();
  return PyLong_FromLong((long)comm);
}
static PyObject* py_comin_parallel_get_host_mpi_rank(PyObject* /*self*/,
                                                     PyObject* /*args*/) {
  int rank = comin_parallel_get_host_mpi_rank();
  py_comin_check_error();
  return PyLong_FromLong((long)rank);
}

static PyObject* py_comin_print_debug(PyObject* /*self*/, PyObject* args,
                                      PyObject* kwargs) {
  PyObject* pymsg;
  static char const* kwlist[] = {"msg", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", (char**)&kwlist, &pymsg))
    return NULL;
  PyObject* msg_str = PyObject_Str(pymsg);
  if (!PyUnicode_Check(msg_str)) {
    return PyErr_Format(PyExc_ValueError, "Conversion to string failed");
  }
  const char* val_str = PyUnicode_AsUTF8(msg_str);
  if (!val_str)
    return NULL;
  comin_print_debug(val_str);
  py_comin_check_error();
  Py_DECREF(msg_str);
  Py_RETURN_NONE;
}

static PyObject* py_comin_print_info(PyObject* /*self*/, PyObject* args,
                                     PyObject* kwargs) {
  PyObject* pymsg;
  static char const* kwlist[] = {"msg", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", (char**)&kwlist, &pymsg))
    return NULL;
  PyObject* msg_str = PyObject_Str(pymsg);
  if (!PyUnicode_Check(msg_str)) {
    return PyErr_Format(PyExc_ValueError, "Conversion to string failed");
  }
  const char* val_str = PyUnicode_AsUTF8(msg_str);
  if (!val_str)
    return NULL;
  comin_print_info(val_str);
  py_comin_check_error();
  Py_DECREF(msg_str);
  Py_RETURN_NONE;
}

static PyObject* py_comin_print_warning(PyObject* /*self*/, PyObject* args,
                                        PyObject* kwargs) {
  PyObject* pymsg;
  static char const* kwlist[] = {"msg", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", (char**)&kwlist, &pymsg))
    return NULL;
  PyObject* msg_str = PyObject_Str(pymsg);
  if (!PyUnicode_Check(msg_str)) {
    return PyErr_Format(PyExc_ValueError, "Conversion to string failed");
  }
  const char* val_str = PyUnicode_AsUTF8(msg_str);
  if (!val_str)
    return NULL;
  comin_print_warning(val_str);
  py_comin_check_error();
  Py_DECREF(msg_str);
  Py_RETURN_NONE;
}

const std::vector<PyMethodDef> py_comin_descrdata_methods() {
  return {
      {"_current_get_plugin_info",
       py_comin_func_wrapper<py_comin_current_get_plugin_info>, METH_NOARGS,
       ""},
      {"_descrdata_eval_property",
       (PyCFunction)py_comin_func_wrapper<py_comin_descrdata_eval_property>,
       METH_VARARGS | METH_KEYWORDS, ""},
      {"_descrdata_get_domain",
       py_comin_func_wrapper<py_comin_descrdata_get_domain>, METH_NOARGS, ""},
      {"_descrdata_get_global",
       py_comin_func_wrapper<py_comin_descrdata_get_global>, METH_NOARGS, ""},
      {"_descrdata_get_simulation_interval",
       py_comin_func_wrapper<py_comin_descrdata_get_simulation_interval>,
       METH_NOARGS, ""},
      {"descrdata_get_timesteplength",
       (PyCFunction)
           py_comin_func_wrapper<py_comin_descrdata_get_timesteplength>,
       METH_VARARGS | METH_KEYWORDS,
       "C function signature: `void double "
       "comin_descrdata_get_timesteplength(int jg)`"},
      {"descrdata_get_index",
       (PyCFunction)py_comin_func_wrapper<py_comin_descrdata_get_index>,
       METH_VARARGS | METH_KEYWORDS,
       "C function signature: `int comin_descrdata_get_index(int j)`"},
      {"descrdata_get_block",
       (PyCFunction)py_comin_func_wrapper<py_comin_descrdata_get_block>,
       METH_VARARGS | METH_KEYWORDS,
       "C function signature: `int comin_descrdata_get_block(int j)`"},
      {"descrdata_get_cell_indices",
       (PyCFunction)py_comin_func_wrapper<py_comin_descrdata_get_cell_indices>,
       METH_VARARGS | METH_KEYWORDS,
       "C function signature: `void comin_descrdata_get_cell_indices(int jg, "
       "int i_blk, int i_startblk, int i_endblk, int* i_startidx, int* "
       "i_endidx, int irl_start, int irl_end)`"},
      {"descrdata_get_cell_npromz",
       (PyCFunction)py_comin_func_wrapper<py_comin_descrdata_get_cell_npromz>,
       METH_VARARGS | METH_KEYWORDS,
       "C function signature: `int comin_descrdata_get_cell_npromz(int jg)`"},
      {"descrdata_get_edge_npromz",
       (PyCFunction)py_comin_func_wrapper<py_comin_descrdata_get_edge_npromz>,
       METH_VARARGS | METH_KEYWORDS,
       "C function signature: `int comin_descrdata_get_edge_npromz(int jg)`"},
      {"descrdata_get_vert_npromz",
       (PyCFunction)py_comin_func_wrapper<py_comin_descrdata_get_vert_npromz>,
       METH_VARARGS | METH_KEYWORDS,
       "C function signature: `int comin_descrdata_get_vert_npromz(int jg)`"},
      {"current_get_datetime",
       py_comin_func_wrapper<py_comin_current_get_datetime>, METH_NOARGS,
       "C function signature: `void comin_current_get_datetime(char "
       "const**,int*,int*);`"},
      {"current_get_domain_id",
       py_comin_func_wrapper<py_comin_current_get_domain_id>, METH_NOARGS,
       "C function signature: `int comin_current_get_domain_id()`"},
      {"descrdata_index_lookup_glb2loc_cell",
       (PyCFunction)
           py_comin_func_wrapper<py_comin_descrdata_index_lookup_glb2loc_cell>,
       METH_VARARGS | METH_KEYWORDS,
       "C function signature: `int "
       "comin_descrdata_index_lookup_glb2loc_cell(int jg, int global_idx)`"},
      {"setup_get_version", py_comin_func_wrapper<py_comin_setup_get_version>,
       METH_NOARGS, "returns (major, minor, patch) version info"},
      {"parallel_get_host_mpi_comm",
       py_comin_func_wrapper<py_comin_parallel_get_host_mpi_comm>, METH_NOARGS,
       "C function signature: `int comin_parallel_get_host_mpi_comm()`"},
      {"parallel_get_plugin_mpi_comm",
       py_comin_func_wrapper<py_comin_parallel_get_plugin_mpi_comm>,
       METH_NOARGS,
       "C function signature: `int comin_parallel_get_host_mpi_comm()`"},
      {"parallel_get_host_mpi_rank",
       py_comin_func_wrapper<py_comin_parallel_get_host_mpi_rank>, METH_NOARGS,
       "C function signature: `int comin_parallel_get_host_mpi_rank()`"},
      {"print_debug", (PyCFunction)py_comin_func_wrapper<py_comin_print_debug>,
       METH_VARARGS | METH_KEYWORDS,
       "C function signature: `void comin_print_debug(const char* msg)`"},
      {"print_info", (PyCFunction)py_comin_func_wrapper<py_comin_print_info>,
       METH_VARARGS | METH_KEYWORDS,
       "C function signature: `void comin_print_info(const char* msg)`"},
      {"print_warning",
       (PyCFunction)py_comin_func_wrapper<py_comin_print_warning>,
       METH_VARARGS | METH_KEYWORDS,
       "C function signature: `void comin_print_warning(const char* msg)`"},
  };
}
