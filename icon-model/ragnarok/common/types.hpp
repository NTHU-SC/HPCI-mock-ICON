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
#ifndef RAGNAROK_COMMON_TYPES_H_
#define RAGNAROK_COMMON_TYPES_H_

#include <Kokkos_Core.hpp>

// unmanaged mutable views
template <typename T>
using View1D = Kokkos::View<T*, Kokkos::DefaultExecutionSpace, Kokkos::MemoryUnmanaged>;

template <typename T>
using View2D = Kokkos::View<T**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace, Kokkos::MemoryUnmanaged>;

template <typename T>
using View3D = Kokkos::View<T***, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace, Kokkos::MemoryUnmanaged>;

//  unmanaged immutable views
template <typename T>
using ConstView1D = Kokkos::View<const T*, Kokkos::DefaultExecutionSpace, Kokkos::MemoryUnmanaged>;

template <typename T>
using ConstView2D =
    Kokkos::View<const T**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace, Kokkos::MemoryUnmanaged>;

template <typename T>
using ConstView3D =
    Kokkos::View<const T***, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace, Kokkos::MemoryUnmanaged>;

//  unmanaged mutable host views
template <typename T>
using HostView1D = Kokkos::View<T*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;

template <typename T>
using HostView2D = Kokkos::View<T**, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;

template <typename T>
using HostView3D = Kokkos::View<T***, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;

//  managed mutable views
template <typename T>
using ManagedView1D = Kokkos::View<T*, Kokkos::DefaultExecutionSpace, Kokkos::MemoryManaged>;

template <typename T>
using ManagedView2D = Kokkos::View<T**, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace, Kokkos::MemoryManaged>;

template <typename T>
using ManagedView3D = Kokkos::View<T***, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace, Kokkos::MemoryManaged>;

// Kokkos default memory space
using MemorySpace   = Kokkos::DefaultExecutionSpace::memory_space;

using RangePolicy2D = Kokkos::MDRangePolicy<Kokkos::Rank<2, Kokkos::Iterate::Left, Kokkos::Iterate::Left>>;

#endif  // RAGNAROK_COMMON_TYPES_H_
