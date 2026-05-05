/* @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

   SPDX-License-Identifier: BSD-3-Clause

   Please see the file LICENSE in the root of the source tree for this code.
   Where software is supplied by third parties, it is indicated in the
   headers of the routines. */

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "comin.h"

#include "exception.h"
#include "util.h"
#include "variables.h"

using namespace std::string_literals;

static PyObject* py_comin_var_request_add(PyObject* /*self*/, PyObject* args,
                                          PyObject* kwargs) {
  t_comin_var_descriptor var_descr;
  static char const* kwlist[] = {"var_descriptor", "lmodexclusive", NULL};
  char* name;
  Py_ssize_t len;
  int lmodexclusive;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "(s#i)p", (char**)&kwlist,
                                   &name, &len, &(var_descr.id),
                                   &lmodexclusive))
    return NULL;

  if (len > COMIN_MAX_LEN_VAR_NAME)
    return PyErr_Format(PyExc_ValueError, "Variable name to long!");

  strncpy(var_descr.name, name, COMIN_MAX_LEN_VAR_NAME);

  comin_var_request_add(var_descr, lmodexclusive);
  py_comin_check_error("var_descr: "s + std::string(name) + ", "s +
                       std::to_string(var_descr.id));
  Py_RETURN_NONE;
}

static PyObject* py_comin_metadata_set(PyObject* /*self*/, PyObject* args,
                                       PyObject* kwargs) {
  t_comin_var_descriptor var_descr;
  char *key, *name;
  Py_ssize_t len;
  PyObject* val;
  static char const* kwlist[] = {"var_descriptor", "key", "value", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "(s#i)sO", (char**)&kwlist,
                                   &name, &len, &(var_descr.id), &key, &val))
    return NULL;

  if (len > COMIN_MAX_LEN_VAR_NAME)
    return PyErr_Format(PyExc_ValueError, "Variable name to long!");

  strncpy(var_descr.name, name, COMIN_MAX_LEN_VAR_NAME);

  if (PyBool_Check(val)) {
    comin_metadata_set_logical(var_descr, key, val == Py_True);
    py_comin_check_error("var_descr: "s + std::string(name) + ", "s +
                         std::to_string(var_descr.id) +
                         " key: " + std::string(key));
    Py_RETURN_NONE;
  }
  if (PyLong_Check(val)) {
    int overflow = 0;
    long val_l   = PyLong_AsLongAndOverflow(val, &overflow);
    if (val_l > std::numeric_limits<int>::max() ||
        val_l < std::numeric_limits<int>::min() || overflow != 0) {
      return PyErr_Format(PyExc_ValueError, "value is out of the range of int");
    }
    comin_metadata_set_integer(var_descr, key, val_l);
    py_comin_check_error("var_descr: "s + std::string(name) + ", "s +
                         std::to_string(var_descr.id) +
                         " key: " + std::string(key));
    Py_RETURN_NONE;
  }
  if (PyFloat_Check(val)) {
    double val_d = PyFloat_AsDouble(val);
    if (val_d == -1. && PyErr_Occurred())
      return PyErr_Format(PyExc_ValueError, "value is not a float");
    comin_metadata_set_real(var_descr, key, val_d);
    py_comin_check_error("var_descr: "s + std::string(name) + ", "s +
                         std::to_string(var_descr.id) +
                         " key: " + std::string(key));
    Py_RETURN_NONE;
  }
  if (PyUnicode_Check(val)) {
    const char* val_str = PyUnicode_AsUTF8(val);
    if (val_str == NULL)
      return PyErr_Format(PyExc_ValueError,
                          "value cannot be converted to char* (UTF-8)");
    comin_metadata_set_character(var_descr, key, val_str);
    py_comin_check_error("var_descr: "s + std::string(name) + ", "s +
                         std::to_string(var_descr.id) +
                         " key: " + std::string(key));
    Py_RETURN_NONE;
  }
  PyTypeObject* type = val->ob_type;
  return PyErr_Format(PyExc_ValueError,
                      "comin_metadata_set is not implemented for the provided "
                      "object type (type=%s, name=%s, id=%d, key=%s)",
                      type->tp_name, var_descr.name, var_descr.id, key);
}

static PyObject* py_comin_metadata_get(PyObject* /*self*/, PyObject* args,
                                       PyObject* kwargs) {
  t_comin_var_descriptor var_desc;
  static char const* kwlist[] = {"var_descriptor", "key", NULL};
  char* name;
  Py_ssize_t len;
  char* key;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "(s#i)s", (char**)&kwlist,
                                   &name, &len, &(var_desc.id), &key))
    return NULL;
  strncpy(var_desc.name, name, len + 1);
  int type = comin_metadata_get_typeid(var_desc, key);
  py_comin_check_error("var_descr: "s + std::string(name) + ", "s +
                       std::to_string(var_desc.id) +
                       " key: " + std::string(key));
  if (type == COMIN_METADATA_TYPEID_LOGICAL) {
    bool val;
    comin_metadata_get_logical(var_desc, key, &val);
    py_comin_check_error("var_descr: "s + std::string(name) + ", "s +
                         std::to_string(var_desc.id) +
                         " key: " + std::string(key));
    if (val)
      Py_RETURN_TRUE;
    else
      Py_RETURN_FALSE;
  }
  if (type == COMIN_METADATA_TYPEID_INTEGER) {
    int val = 0;
    comin_metadata_get_integer(var_desc, key, &val);
    py_comin_check_error("var_descr: "s + std::string(name) + ", "s +
                         std::to_string(var_desc.id) +
                         " key: " + std::string(key));
    return PyLong_FromLong(val);
  }
  if (type == COMIN_METADATA_TYPEID_REAL) {
    double val = 0.;
    comin_metadata_get_real(var_desc, key, &val);
    py_comin_check_error("var_descr: "s + std::string(name) + ", "s +
                         std::to_string(var_desc.id) +
                         " key: " + std::string(key));
    return PyFloat_FromDouble(val);
  }
  if (type == COMIN_METADATA_TYPEID_CHARACTER) {
    const char* val = NULL;
    int len         = -1;
    comin_metadata_get_character(var_desc, key, &val, &len);
    py_comin_check_error("var_descr: "s + std::string(name) + ", "s +
                         std::to_string(var_desc.id) +
                         " key: " + std::string(key));
    return PyUnicode_FromStringAndSize(val, len);
  }
  return PyErr_Format(PyExc_KeyError, "Key not found: %s", key);
}

static PyObject* py_comin_var_get(PyObject* /*self*/, PyObject* args,
                                  PyObject* kwargs) {
  PyObject* context;
  t_comin_var_descriptor var_desc;
  int flag = -1;

  static char const* kwlist[] = {"context", "var_descriptor", "flag", NULL};
  char* name;
  Py_ssize_t len;
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O(s#i)|i", (char**)&kwlist,
                                   &context, &name, &len, &(var_desc.id),
                                   &flag))
    return NULL;
  strncpy(var_desc.name, name, len + 1);

  if (!PyList_Check(context))
    return PyErr_Format(PyExc_ValueError,
                        "PyCominVar requires list of integers as contexts");
  std::vector<int> icontext(PyList_Size(context));
  for (size_t i = 0; i < icontext.size(); ++i)
    icontext[i] = (int)PyLong_AsLong(PyList_GetItem(context, i));

  t_comin_var_handle* var = comin_var_get(
      icontext.size(), (t_comin_entry_point*)icontext.data(), var_desc, flag);
  py_comin_check_error("var_descr: "s + std::string(name) + ", "s +
                       std::to_string(var_desc.id));
  return PyCapsule_New(var, "var", NULL);
}

static PyObject* py_comin_var_get_buffer(PyObject* self, PyObject* args,
                                         PyObject* kwargs) {
  PyObject* handle_cap;
  static char const* kwlist[] = {"handle", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", (char**)&kwlist,
                                   &handle_cap))
    return NULL;
  t_comin_var_handle* handle =
      (t_comin_var_handle*)PyCapsule_GetPointer(handle_cap, "var");
  t_comin_var_descriptor descr;
  comin_var_get_descriptor(handle, &descr);

  void* ptr = comin_var_get_ptr(handle);
  py_comin_check_error();
  if (ptr == NULL)
    return PyErr_Format(PyExc_ValueError, "comin_var_get_ptr failed");
  int shape[5];
  comin_var_get_shape(handle, shape);
  py_comin_check_error();
  Py_buffer buffer;
  int type = -1;
  comin_metadata_get_integer(descr, "datatype", &type);
  py_comin_check_error("'datatype' not set for "s + std::string(descr.name) +
                       ", id="s + std::to_string(descr.id));
  switch (type) {
  case COMIN_VAR_DATATYPE_DOUBLE:
    fill_buffer<double>(&buffer, ptr, shape, 5, 0);
    break;
  case COMIN_VAR_DATATYPE_FLOAT:
    fill_buffer<float>(&buffer, ptr, shape, 5, 0);
    break;
  case COMIN_VAR_DATATYPE_INT:
    fill_buffer<int>(&buffer, ptr, shape, 5, 0);
    break;
  default:
    throw std::runtime_error("Unknown datatype for "s +
                             std::string(descr.name) + ", id="s +
                             std::to_string(descr.id));
  };
  buffer.obj = self;
  return PyMemoryView_FromBuffer(&buffer);
}

static PyObject* py_comin_var_get_device_ptr(PyObject* /*self*/, PyObject* args,
                                             PyObject* kwargs) {
  PyObject* handle_cap;
  static char const* kwlist[] = {"handle", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", (char**)&kwlist,
                                   &handle_cap))
    return NULL;
  t_comin_var_handle* handle =
      (t_comin_var_handle*)PyCapsule_GetPointer(handle_cap, "var");

  void* ptr = comin_var_get_device_ptr(handle);
  py_comin_check_error();
  return PyLong_FromVoidPtr(ptr);
}

static PyObject* py_comin_var_get_dim_semantics(PyObject* /*self*/,
                                                PyObject* args,
                                                PyObject* kwargs) {
  PyObject* handle_cap;
  static char const* kwlist[] = {"handle", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", (char**)&kwlist,
                                   &handle_cap))
    return NULL;
  t_comin_var_handle* handle =
      (t_comin_var_handle*)PyCapsule_GetPointer(handle_cap, "var");
  int dim_semantics[5];
  comin_var_get_dim_semantics(handle, dim_semantics);
  py_comin_check_error();
  return Py_BuildValue("iiiii", dim_semantics[0], dim_semantics[1],
                       dim_semantics[2], dim_semantics[3], dim_semantics[4]);
}

static PyObject* py_comin_var_get_ncontained(PyObject* /*self*/, PyObject* args,
                                             PyObject* kwargs) {
  PyObject* handle_cap;
  static char const* kwlist[] = {"handle", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", (char**)&kwlist,
                                   &handle_cap))
    return NULL;
  t_comin_var_handle* handle =
      (t_comin_var_handle*)PyCapsule_GetPointer(handle_cap, "var");
  int ncontained;
  comin_var_get_ncontained(handle, &ncontained);
  py_comin_check_error();
  return PyLong_FromLong(ncontained);
}

static PyObject* py_comin_var_get_descriptor(PyObject* /*self*/, PyObject* args,
                                             PyObject* kwargs) {
  PyObject* handle_cap;
  static char const* kwlist[] = {"handle", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", (char**)&kwlist,
                                   &handle_cap))
    return NULL;
  t_comin_var_handle* handle =
      (t_comin_var_handle*)PyCapsule_GetPointer(handle_cap, "var");
  t_comin_var_descriptor descr;
  comin_var_get_descriptor(handle, &descr);
  py_comin_check_error();
  return Py_BuildValue("si", descr.name, descr.id);
}

static PyObject* py_comin_var_get_descr_list_head(PyObject* /*self*/,
                                                  PyObject* /*args*/) {
  void* head = comin_var_get_descr_list_head();
  py_comin_check_error();
  return PyCapsule_New(head, "var_descr_list", NULL);
}

static PyObject* py_comin_var_get_descr_list_next(PyObject* /*self*/,
                                                  PyObject* args,
                                                  PyObject* kwargs) {
  PyObject* py_current;
  static char const* kwlist[] = {"current", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", (char**)&kwlist,
                                   &py_current))
    return NULL;
  t_comin_var_descr_list_iterator* it =
      (t_comin_var_descr_list_iterator*)PyCapsule_GetPointer(py_current,
                                                             "var_descr_list");
  void* next = comin_var_get_descr_list_next(it);
  py_comin_check_error();
  if (next == NULL)
    Py_RETURN_NONE;
  else
    return PyCapsule_New(next, "var_descr_list", NULL);
}

static PyObject* py_comin_var_get_descr_list_var_desc(PyObject* /*self*/,
                                                      PyObject* args,
                                                      PyObject* kwargs) {
  PyObject* py_current;
  static char const* kwlist[] = {"current", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", (char**)&kwlist,
                                   &py_current))
    return NULL;
  t_comin_var_descr_list_iterator* current =
      (t_comin_var_descr_list_iterator*)PyCapsule_GetPointer(py_current,
                                                             "var_descr_list");
  t_comin_var_descriptor var_desc;
  comin_var_get_descr_list_var_desc(current, &var_desc);
  py_comin_check_error();
  return Py_BuildValue("si", var_desc.name, var_desc.id);
}

static PyObject* py_comin_var_descr_list_iterator_delete(PyObject* /*self*/,
                                                         PyObject* args,
                                                         PyObject* kwargs) {
  PyObject* py_current;
  static char const* kwlist[] = {"current", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", (char**)&kwlist,
                                   &py_current))
    return NULL;
  t_comin_var_descr_list_iterator* it =
      (t_comin_var_descr_list_iterator*)PyCapsule_GetPointer(py_current,
                                                             "var_descr_list");
  comin_var_descr_list_iterator_delete(&it);
  py_comin_check_error();
  Py_RETURN_NONE;
}

static PyObject* py_comin_metadata_get_iterator_begin(PyObject* /*self*/,
                                                      PyObject* args,
                                                      PyObject* kwargs) {
  t_comin_var_descriptor var_desc;
  char* name;
  Py_ssize_t len;
  static char const* kwlist[] = {"var_descriptor", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "(s#i)", (char**)&kwlist,
                                   &name, &len, &(var_desc.id)))
    return NULL;
  strncpy(var_desc.name, name, COMIN_MAX_LEN_VAR_NAME);
  void* it = comin_metadata_get_iterator_begin(var_desc);
  py_comin_check_error("var_descr: "s + std::string(var_desc.name) + ", "s +
                       std::to_string(var_desc.id));
  return PyCapsule_New(it, "metadata_iterator", NULL);
}

static PyObject* py_comin_metadata_get_iterator_end(PyObject* /*self*/,
                                                    PyObject* args,
                                                    PyObject* kwargs) {
  t_comin_var_descriptor var_desc;
  char* name;
  Py_ssize_t len;
  static char const* kwlist[] = {"var_descriptor", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "(s#i)", (char**)&kwlist,
                                   &name, &len, &(var_desc.id)))
    return NULL;
  strncpy(var_desc.name, name, COMIN_MAX_LEN_VAR_NAME);
  void* it = comin_metadata_get_iterator_end(var_desc);
  py_comin_check_error("var_descr: "s + std::string(var_desc.name) + ", "s +
                       std::to_string(var_desc.id));
  return PyCapsule_New(it, "metadata_iterator", NULL);
}

static PyObject* py_comin_metadata_iterator_get_key(PyObject* /*self*/,
                                                    PyObject* args,
                                                    PyObject* kwargs) {
  PyObject* py_it;
  static char const* kwlist[] = {"it", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", (char**)&kwlist, &py_it))
    return NULL;
  t_comin_metadata_iterator* it =
      (t_comin_metadata_iterator*)PyCapsule_GetPointer(py_it,
                                                       "metadata_iterator");
  const char* key = comin_metadata_iterator_get_key(it);
  py_comin_check_error();
  return Py_BuildValue("s", key);
}

static PyObject* py_comin_metadata_iterator_compare(PyObject* /*self*/,
                                                    PyObject* args,
                                                    PyObject* kwargs) {
  PyObject *py_it1, *py_it2;
  static char const* kwlist[] = {"it1", "it2", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OO", (char**)&kwlist, &py_it1,
                                   &py_it2))
    return NULL;
  t_comin_metadata_iterator* it1 =
      (t_comin_metadata_iterator*)PyCapsule_GetPointer(py_it1,
                                                       "metadata_iterator");
  t_comin_metadata_iterator* it2 =
      (t_comin_metadata_iterator*)PyCapsule_GetPointer(py_it2,
                                                       "metadata_iterator");
  bool result = comin_metadata_iterator_compare(it1, it2);
  py_comin_check_error();
  if (result)
    Py_RETURN_TRUE;
  else
    Py_RETURN_FALSE;
}

static PyObject* py_comin_metadata_iterator_next(PyObject* /*self*/,
                                                 PyObject* args,
                                                 PyObject* kwargs) {
  PyObject* py_it;
  static char const* kwlist[] = {"it", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", (char**)&kwlist, &py_it))
    return NULL;
  t_comin_metadata_iterator* it =
      (t_comin_metadata_iterator*)PyCapsule_GetPointer(py_it,
                                                       "metadata_iterator");
  comin_metadata_iterator_next(it);
  py_comin_check_error();
  Py_RETURN_NONE;
}

static PyObject* py_comin_metadata_iterator_delete(PyObject* /*self*/,
                                                   PyObject* args,
                                                   PyObject* kwargs) {
  PyObject* py_it;
  static char const* kwlist[] = {"it", NULL};
  if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O", (char**)&kwlist, &py_it))
    return NULL;
  t_comin_metadata_iterator* it =
      (t_comin_metadata_iterator*)PyCapsule_GetPointer(py_it,
                                                       "metadata_iterator");
  comin_metadata_iterator_delete(it);
  py_comin_check_error();
  Py_RETURN_NONE;
}

std::vector<PyMethodDef> py_comin_variables_methods() {
  return {
      {"var_request_add",
       (PyCFunction)py_comin_func_wrapper<py_comin_var_request_add>,
       METH_VARARGS | METH_KEYWORDS,
       "Request the host model to add a variable, arguments:"
       " (name string, domain id), lmodexclusive"},
      {"metadata_get",
       (PyCFunction)py_comin_func_wrapper<py_comin_metadata_get>,
       METH_VARARGS | METH_KEYWORDS,
       "retrieve metadata, arguments: (name string, domain id) , metadata key"},
      {"_metadata_set",
       (PyCFunction)py_comin_func_wrapper<py_comin_metadata_set>,
       METH_VARARGS | METH_KEYWORDS, ""},
      {"_var_get", (PyCFunction)py_comin_func_wrapper<py_comin_var_get>,
       METH_VARARGS | METH_KEYWORDS, ""},
      {"_var_get_device_ptr",
       (PyCFunction)py_comin_func_wrapper<py_comin_var_get_device_ptr>,
       METH_VARARGS | METH_KEYWORDS, ""},
      {"_var_get_buffer",
       (PyCFunction)py_comin_func_wrapper<py_comin_var_get_buffer>,
       METH_VARARGS | METH_KEYWORDS, ""},
      {"_var_get_dim_semantics",
       (PyCFunction)py_comin_func_wrapper<py_comin_var_get_dim_semantics>,
       METH_VARARGS | METH_KEYWORDS, ""},
      {"_var_get_ncontained",
       (PyCFunction)py_comin_func_wrapper<py_comin_var_get_ncontained>,
       METH_VARARGS | METH_KEYWORDS, ""},
      {"_var_get_descriptor",
       (PyCFunction)py_comin_func_wrapper<py_comin_var_get_descriptor>,
       METH_VARARGS | METH_KEYWORDS, ""},
      {"_var_get_descr_list_head",
       py_comin_func_wrapper<py_comin_var_get_descr_list_head>, METH_NOARGS,
       ""},
      {"_var_get_descr_list_next",
       (PyCFunction)py_comin_func_wrapper<py_comin_var_get_descr_list_next>,
       METH_VARARGS | METH_KEYWORDS, ""},
      {"_var_get_descr_list_var_desc",
       (PyCFunction)py_comin_func_wrapper<py_comin_var_get_descr_list_var_desc>,
       METH_VARARGS | METH_KEYWORDS, ""},
      {"_var_descr_list_iterator_delete",
       (PyCFunction)
           py_comin_func_wrapper<py_comin_var_descr_list_iterator_delete>,
       METH_VARARGS | METH_KEYWORDS, ""},
      {"_metadata_get_iterator_begin",
       (PyCFunction)py_comin_func_wrapper<py_comin_metadata_get_iterator_begin>,
       METH_VARARGS | METH_KEYWORDS, ""},
      {"_metadata_get_iterator_end",
       (PyCFunction)py_comin_func_wrapper<py_comin_metadata_get_iterator_end>,
       METH_VARARGS | METH_KEYWORDS, ""},
      {"_metadata_iterator_get_key",
       (PyCFunction)py_comin_func_wrapper<py_comin_metadata_iterator_get_key>,
       METH_VARARGS | METH_KEYWORDS, ""},
      {"_metadata_iterator_compare",
       (PyCFunction)py_comin_func_wrapper<py_comin_metadata_iterator_compare>,
       METH_VARARGS | METH_KEYWORDS, ""},
      {"_metadata_iterator_next",
       (PyCFunction)py_comin_func_wrapper<py_comin_metadata_iterator_next>,
       METH_VARARGS | METH_KEYWORDS, ""},
      {"_metadata_iterator_delete",
       (PyCFunction)py_comin_func_wrapper<py_comin_metadata_iterator_delete>,
       METH_VARARGS | METH_KEYWORDS, ""},
  };
}
