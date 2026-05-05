// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERP_METHOD_FILE_H
#define INTERP_METHOD_FILE_H

#include "interp_method.h"

// YAC PUBLIC HEADER START

enum yac_interp_file_on_missing_file {
  YAC_INTERP_FILE_MISSING_ERROR = 0, //!< abort on missing file
  YAC_INTERP_FILE_MISSING_CONT  = 1, //!< continue on missing file
};

enum yac_interp_file_on_success {
  YAC_INTERP_FILE_SUCCESS_STOP = 0, //!< prevents following interpolation method
                                    //!< from computating further weights
  YAC_INTERP_FILE_SUCCESS_CONT = 1, //!< continue weight computation with
                                    //!< following interpolation methods
};

#define YAC_WEIGHT_FILE_VERSION_STRING "yac weight file 1.0"

#define YAC_INTERP_FILE_WEIGHT_FILE_NAME_DEFAULT (NULL)
#define YAC_INTERP_FILE_ON_MISSING_FILE_DEFAULT (YAC_INTERP_FILE_MISSING_ERROR)
#define YAC_INTERP_FILE_ON_SUCCESS_DEFAULT (YAC_INTERP_FILE_SUCCESS_CONT)

/**
 * Construtor for an interpolation method that reads in a weight file
 * @param[in] weight_file_name name of the weight file
 * @param[in] on_missing_file  specifies how YAC should behave if no file was
 *                             found
 * @param[in] on_success       specifies how YAC should behave in case a weight
 *                             file was successfully read
 */
struct interp_method * yac_interp_method_file_new(
  char const * weight_file_name,
  enum yac_interp_file_on_missing_file on_missing_file,
  enum yac_interp_file_on_success on_success);

// YAC PUBLIC HEADER STOP

#endif // INTERP_METHOD_FILE_H
