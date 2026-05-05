// ICON
//
// ---------------------------------------------------------------
// Copyright (C) 2004-2026, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
// Contact information: icon-model.org
//
// See AUTHORS.TXT for a list of authors
// See LICENSES/ for license information
// SPDX-License-Identifier: BSD-3-Clause
// ---------------------------------------------------------------

///
/// @file
/// @brief This header file contains the definitions required
/// for the data/procedure exchange with ICON Fortran code
///
//----------------------------

#ifndef RAGNAROK_SUPPORT_ICON_BRIDGE_RAGNAROK_H_
#define RAGNAROK_SUPPORT_ICON_BRIDGE_RAGNAROK_H_

namespace ragnarok {

bool is_initialized();
void init();
const char* retrieve_kokkos_version(int&);

}  // namespace ragnarok

extern "C" {
void init_ragnarok();
char* retrieve_kokkos_version_c(int&);
}

#endif  // RAGNAROK_SUPPORT_ICON_BRIDGE_RAGNAROK_H_
