#!/bin/sh

# SPDX-FileCopyrightText: 2025 DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
# SPDX-License-Identifier: BSD-3-Clause

set -e

script_dir=`echo "$0" | sed 's@[^/]*$@@'`
(unset CDPATH) >/dev/null 2>&1 && unset CDPATH
cd "$script_dir"

exec autoreconf -fvi
