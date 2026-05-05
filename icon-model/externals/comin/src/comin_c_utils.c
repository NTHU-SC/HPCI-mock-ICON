/**
   @file comin_c_utils.c
   @brief C interface for the ICON Community Interface

   @authors 12/2024 :: ICON Community Interface  <comin@icon-model.org>

   SPDX-License-Identifier: BSD-3-Clause

   Please see the file LICENSE in the root of the source tree for this code.
   Where software is supplied by third parties, it is indicated in the
   headers of the routines. **/

#include <stdarg.h>
#include <stdio.h>

#include <comin.h>

void comin_print_debug_f(const char* fmt, ...) {
  va_list argp;
  va_start(argp, fmt);
  va_list args2;
  va_copy(args2, argp);
  char buf[1 + vsnprintf(NULL, 0, fmt, argp)];
  va_end(argp);
  vsnprintf(buf, sizeof(buf), fmt, args2);
  va_end(args2);
  comin_print_debug(buf);
}

void comin_print_info_f(const char* fmt, ...) {
  va_list argp;
  va_start(argp, fmt);
  va_list args2;
  va_copy(args2, argp);
  char buf[1 + vsnprintf(NULL, 0, fmt, argp)];
  va_end(argp);
  vsnprintf(buf, sizeof(buf), fmt, args2);
  va_end(args2);
  comin_print_info(buf);
}

void comin_print_warning_f(const char* fmt, ...) {
  va_list argp;
  va_start(argp, fmt);
  va_list args2;
  va_copy(args2, argp);
  char buf[1 + vsnprintf(NULL, 0, fmt, argp)];
  va_end(argp);
  vsnprintf(buf, sizeof(buf), fmt, args2);
  va_end(args2);
  comin_print_warning(buf);
}
