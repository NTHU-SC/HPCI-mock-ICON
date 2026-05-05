#include "exception.h"

PyObject* PyExc_ComInError;

void py_comin_check_error(std::string add_info) {
  t_comin_error_code error_code = comin_error_get();
  if (error_code == COMIN_SUCCESS)
    return;
  comin_error_reset();
  throw ComInError(error_code, add_info);
}
