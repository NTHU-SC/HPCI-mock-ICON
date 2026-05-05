/* @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

   SPDX-License-Identifier: BSD-3-Clause

   Please see the file LICENSE in the root of the source tree for this code.
   Where software is supplied by third parties, it is indicated in the
   headers of the routines. */

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <array>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <dlfcn.h>

#include "comin.h"
#include "util.h"

#include "comin.py.h"
#include "config.h"

#include "callbacks.h"
#include "descrdata.h"
#include "exception.h"
#include "variables.h"

static struct PyModuleDef py_comin_module = {
    PyModuleDef_HEAD_INIT, "_comin", /* name of module */
    NULL,                            /* module documentation, may be NULL */
    -1,  /* size of per-interpreter state of the module,
            or -1 if the module keeps state in global variables. */
    NULL /*m_methods - will be set in main before calling PyModule_Create*/
};

static std::vector<PyMethodDef> pyCominMethods;

PyMODINIT_FUNC PyInit_comin(void) {
  // collect methods from other cpp files:
  for (const auto& m : py_comin_callbacks_methods())
    pyCominMethods.push_back(m);
  for (const auto& m : py_comin_variables_methods())
    pyCominMethods.push_back(m);
  for (const auto& m : py_comin_descrdata_methods())
    pyCominMethods.push_back(m);
  pyCominMethods.push_back({NULL, NULL, 0, NULL}); /* Sentinel */

  py_comin_module.m_methods = pyCominMethods.data();
  PyObject* pSelf           = PyModule_Create(&py_comin_module);
  if (pSelf == NULL)
    return NULL;

  PyExc_ComInError =
      PyErr_NewException("comin.ComInError", PyExc_RuntimeError, NULL);
  Py_INCREF(PyExc_ComInError);
  PyModule_AddObject(pSelf, "ComInError", PyExc_ComInError);

  return pSelf;
}

static char* extract_filename(const char* str, int size) {
  PyObject* pShlex = PyImport_ImportModule("shlex");
  PyObject* pDict  = PyModule_GetDict(pShlex); // returns a borrowed reference
  PyObject* pSplit = PyDict_GetItemString(
      pDict, (char*)"split"); // returns a borrowed reference
  PyObject* pList = PyObject_CallFunction(pSplit, "s#", str, (Py_ssize_t)size);
  PyObject* pFilename =
      PyList_GetItem(pList, 0); // returns a borrowed reference
  char* filename = strdup(PyUnicode_AsUTF8(pFilename));
  Py_DECREF(pList);
  Py_DECREF(pShlex);
  return filename;
}

static std::string exec(std::string cmd) {
  std::array<char, 128> buffer;
  std::string result;
  std::unique_ptr<FILE, void (*)(FILE*)> pipe(
      popen(cmd.c_str(), "r"),
      [](FILE* f) -> void { std::ignore = pclose(f); });
  if (!pipe) {
    comin_plugin_finish("python_adapter",
                        "Cannot execute given python executable.");
  }
  while (fgets(buffer.data(), static_cast<int>(buffer.size()), pipe.get()) !=
         nullptr) {
    result += buffer.data();
  }
  return result;
}

static void python_version_check(std::string python_exe) {
  std::string py_version_str = exec(python_exe + " --version");
  // assuming that the output is something like "Python 3.13.0"
  py_version_str =
      py_version_str.substr(py_version_str.find(" ") + 1); // truncate "Python "
  int major =
      std::atoi(py_version_str.substr(0, py_version_str.find(".")).c_str());
  py_version_str =
      py_version_str.substr(py_version_str.find(".") + 1); // truncate "3."
  int minor =
      std::atoi(py_version_str.substr(0, py_version_str.find(".") + 1).c_str());
  if (major != PY_MAJOR_VERSION || minor != PY_MINOR_VERSION) {
    comin_plugin_finish("python_adapter::python_version_check",
                        ("Versions of the embedded interpreter (" +
                         std::to_string(PY_MAJOR_VERSION) + "." +
                         std::to_string(PY_MINOR_VERSION) + ")" +
                         " and the python executable (" +
                         std::to_string(major) + "." + std::to_string(minor) +
                         ") are not compatible.")
                            .c_str());
  }
}

static void initialize_python() {
  PyStatus status;
  PyConfig config;
  PyConfig_InitPythonConfig(&config);
  std::string python3_exe = Python3_EXECUTABLE;
  if (const char* python3_exe_env = std::getenv("COMIN_PYTHON_EXECUTABLE"))
    python3_exe = python3_exe_env;

  python_version_check(python3_exe);

  status =
      PyConfig_SetBytesString(&config, &config.executable, python3_exe.c_str());
  if (PyStatus_Exception(status)) {
    comin_plugin_finish(
        "python_adapter::comin_main",
        ("Connot set config.executable for python initialization." +
         std::string(status.err_msg))
            .c_str());
    return;
  }

  status = Py_InitializeFromConfig(&config);
  if (PyStatus_Exception(status)) {
    comin_plugin_finish(
        "python_adapter::comin_main",
        ("Python initialization failed. (" + std::string(status.err_msg) + ")")
            .c_str());
    return;
  }
  PyConfig_Clear(&config);
}

static int py_comin_instance_counter = 0;

extern "C" {
  void comin_main() {

    using namespace std::string_literals;

    int mpi_rank = comin_parallel_get_host_mpi_rank();
    if (mpi_rank == 0)
      std::cerr << "setup_python_adapter" << std::endl;

    if (!Py_IsInitialized()) {
      // this is a workaround for a python problem that occurs if python
      // is embedded in a library that is loaded with dlopen.
      // See for example:
      // https://bugs.python.org/issue4434
      // https://stackoverflow.com/questions/8302810/undefined-symbol-in-c-when-loading-a-python-shared-library
      // https://stackoverflow.com/questions/64295279/so-type-plugin-with-embedded-python-interpreter-how-to-link-load-libpython
      char libname[17];
      sprintf(libname, "libpython3.%d.so", PY_MINOR_VERSION);
      void* handle = dlopen(libname, RTLD_LAZY | RTLD_GLOBAL);
      if (handle == nullptr) {
        if (mpi_rank == 0) {
          std::cerr << "Cannot load " << libname << ": " << dlerror()
                    << std::endl;
          comin_plugin_finish(__func__, "Plugin cannot be loaded.");
        } else {
          if (mpi_rank == 0)
            std::cerr << "Python adapter: " << libname << " loaded!"
                      << std::endl;
        }
      }
      dlclose(handle);

      PyImport_AppendInittab("_comin", &PyInit_comin);

      initialize_python();

      std::string comin_py_c((char*)comin_py, comin_py_len);
      PyObject* compiled_comin =
          Py_CompileString(comin_py_c.c_str(), "comin.py", Py_file_input);
      PyImport_ExecCodeModule("comin", compiled_comin);
    }

    py_comin_instance_counter++;
    comin_callback_register(EP_DESTRUCTOR, []() {
      py_comin_generic_callback();
      py_comin_instance_counter--;
      if (py_comin_instance_counter == 0)
        Py_Finalize();
    });

    // Dictionary to store the plugins global variables
    // independently from other python plugins
    PyObject* globals = PyDict_New();

    int ilen                     = -1;
    const char* plugin_options_c = NULL;
    comin_current_get_plugin_options(&plugin_options_c, &ilen);
    char* filename = extract_filename(plugin_options_c, ilen);

    // add directory of filename to sys.path
    std::filesystem::path fpath = filename;
    auto parent_path = std::filesystem::canonical(fpath.parent_path());
    PyRun_SimpleString(
        ("import sys; sys.path.insert(0,\""s + parent_path.c_str() + "\")"s)
            .c_str());

    // errors need to be handled by py_comin_check_error
    comin_error_set_errors_return(true);
    try {
      if (mpi_rank == 0)
        std::cerr << "Running python script " << filename << std::endl;
      FILE* file = fopen(filename, "r");
      if (file == NULL)
        throw std::runtime_error("Cannot read "s + std::string(filename));
      else {
        PyObject* pyResult =
            PyRun_File(file, filename, Py_file_input, globals, globals);
        Py_XDECREF(pyResult);
        fclose(file);
        if (PyErr_Occurred()) {
          PyErr_Print();
          throw std::runtime_error("Error while executing "s +
                                   std::string(filename));
        }
        PyRun_SimpleString(
            ("sys.path.remove(\""s + parent_path.c_str() + "\")"s).c_str());
      }
    } catch (std::exception& err) {
      std::cerr << err.what() << std::endl;
      comin_plugin_finish("python_adapter::comin_main",
                          "Error while executing script");
      return;
    }
  }
}
