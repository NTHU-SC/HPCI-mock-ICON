#ifndef EXCPEPTION_H
#define EXCPEPTION_H

#include <memory>
#include <stdexcept>
#include <string>

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <comin.h>

extern PyObject* PyExc_ComInError;

class ComInError : public std::runtime_error {
  std::string _what;

  static std::string format_message(t_comin_error_code error_code,
                                    std::string add_info) {
    char category[11];
    char message[COMIN_MAX_LEN_ERR_MESSAGE];
    comin_error_get_message(error_code, category, message);
    return std::string(category) + ": " + std::string(message) + "\n" +
           add_info;
  }

public:
  ComInError(t_comin_error_code error_code, std::string add_info)
      : std::runtime_error(format_message(error_code, add_info)) {}
};

void py_comin_check_error(std::string add_info = "");

template <PyCFunction Fun>
PyObject* py_comin_func_wrapper(PyObject* self, PyObject* args) {
  try {
    return Fun(self, args);
  } catch (ComInError& e) {
    return PyErr_Format(PyExc_ComInError, e.what());
  }
}

template <PyCFunctionWithKeywords Fun>
PyObject* py_comin_func_wrapper(PyObject* self, PyObject* args,
                                PyObject* kwargs) {
  try {
    return Fun(self, args, kwargs);
  } catch (ComInError& e) {
    return PyErr_Format(PyExc_ComInError, e.what());
  }
}

#endif
