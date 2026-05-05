// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef YAC_XMAP_H
#define YAC_XMAP_H

#include "yaxt.h"
#include <remote_point.h>

struct yac_xmap_;
typedef struct yac_xmap_ * yac_xmap;

yac_xmap yac_xmap_from_point_infos(
  struct remote_point_infos * point_infos, size_t count, MPI_Comm comm);

Xt_redist yac_xmap_generate_redist(yac_xmap xmap, MPI_Datatype base_type);

void yac_xmap_delete(yac_xmap xmap);

#endif // YAC_XMAPL_H
