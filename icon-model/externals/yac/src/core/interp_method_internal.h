// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef INTERP_METHOD_INTERNAL_H
#define INTERP_METHOD_INTERNAL_H

#include "interp_method.h"
#include "dist_grid_internal.h"
#include "interp_grid_internal.h"
#include "interp_weights_internal.h"

struct interp_method_vtable {

  size_t (*do_search)( // returns number of target points interpolated by
                       // this method
    struct interp_method * method,       // pointer to interpolation method
    struct yac_interp_grid * grid,       // interpolation grid
    size_t * tgt_points,                 // local indices for target points
                                         // that are to be interpolated
                                         // (do_search has to reorder entries,
                                         //  such that all target points
                                         //  interpolated by this method are at
                                         //  the front of the array
    size_t count,                        // number of target points to be
                                         // interpolated
    struct yac_interp_weights * weights, // interpolation weights
    int * interpolation_complete);       // if != 0, do not do any
                                         // interpolation
                                         // set to 1 if remaining interpolation
                                         // stack is not supposed to compute any
                                         // more weights
  void (*delete)(struct interp_method * method);
};

struct interp_method {

  struct interp_method_vtable *vtable;
};

#endif // INTERP_METHOD_INTERNAL_H
