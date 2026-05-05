"""
Test plugin for the ICON Community Interface (ComIn)

This simple test plugin shows how to use the basic features of
ComIn analogous to simple_c_plugin and simple_fortran_plugin.

@authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>

SPDX-License-Identifier: BSD-3-Clause

Please see the file LICENSE in the root of the source tree for this code.
Where software is supplied by third parties, it is indicated in the
headers of the routines.
"""

import numpy as np
import numpy.ma as npma
import comin
from itertools import product

rank = comin.parallel_get_host_mpi_rank()
# first domain
glob = comin.descrdata_get_global()
domains = [comin.descrdata_get_domain(jg) for jg in range(1, glob.n_dom + 1)]

# create owner mask
owner_masks = [
    np.asarray(domains[jg - 1].cells.decomp_domain) == 0
    for jg in range(1, glob.n_dom + 1)
]

promz_c_masks = [
    np.full((glob.nproma, domains[jg - 1].cells.nblks), False)
    for jg in range(1, glob.n_dom + 1)
]
for jg in range(1, glob.n_dom + 1):
    promz_c_masks[jg - 1][comin.descrdata_get_cell_npromz(jg) :, -1] = True

dim2comin_zaxis = {2: comin.COMIN_ZAXIS_2D, 3: comin.COMIN_ZAXIS_3D}

comin_dtypes = [
    comin.COMIN_VAR_DATATYPE_DOUBLE,
    comin.COMIN_VAR_DATATYPE_FLOAT,
    comin.COMIN_VAR_DATATYPE_INT,
]

var_descriptors = [
    (
        (f"test_var_{dim}D_{dtype}_{sync_mode}", jg),
        {"zaxis_id": dim2comin_zaxis[dim], "datatype": dtype},
        sync_mode,
    )
    for dim, dtype, jg, sync_mode in product(
        [2, 3], comin_dtypes, range(1, glob.n_dom + 1), ["after_write", "before_read"]
    )
]

for var_descr, metadata, _ in var_descriptors:
    comin.var_request_add(var_descr, False)
    comin.metadata_set(var_descr, **metadata)


@comin.register_callback(comin.EP_SECONDARY_CONSTRUCTOR)
def simple_python_constructor():
    global test_vars_before_sync, test_vars_after_sync

    test_vars_before_sync = [
        comin.var_get(
            [comin.EP_ATM_TIMELOOP_START],
            var_descr,
            comin.COMIN_FLAG_WRITE
            | (comin.COMIN_FLAG_SYNC_HALO if sm == "after_write" else 0),
        )
        for var_descr, _, sm in var_descriptors
    ]

    test_vars_after_sync = [
        comin.var_get(
            [comin.EP_ATM_TIMELOOP_END],
            var_descr,
            comin.COMIN_FLAG_READ
            | (comin.COMIN_FLAG_SYNC_HALO if sm == "before_read" else 0),
        )
        for var_descr, _, sm in var_descriptors
    ]


@comin.register_callback(comin.EP_ATM_TIMELOOP_START)
def fill_test_var():
    # fill test_var
    for test_var in test_vars_before_sync:
        comin.print_info(f"Filling {test_var.descriptor}")
        name, jg = test_var.descriptor
        if comin.COMIN_DIM_SEMANTICS_LEVEL in test_var.dim_semantics:  # 3d
            mask = np.repeat(
                owner_masks[jg - 1][:, None, :], domains[jg - 1].nlev, axis=1
            )
        else:  # 2d
            mask = owner_masks[jg - 1]
        test_var_np = np.asarray(test_var)
        test_var_np[mask] = rank + 1
        test_var_np[~mask] = -1.0


@comin.register_callback(comin.EP_ATM_TIMELOOP_END)
def check_halo_sync():
    for test_var in test_vars_after_sync:
        comin.print_info(f"Checking {test_var.descriptor}")
        name, jg = test_var.descriptor
        if comin.COMIN_DIM_SEMANTICS_LEVEL in test_var.dim_semantics:  # 3d
            test_var_np = npma.array(
                test_var,
                mask=np.repeat(
                    promz_c_masks[jg - 1][:, None, :], domains[jg - 1].nlev, axis=1
                ),
            )
            mask = np.repeat(
                owner_masks[jg - 1][:, None, :], domains[jg - 1].nlev, axis=1
            )
        else:  # 2d
            test_var_np = npma.array(test_var, mask=promz_c_masks[jg - 1])
            mask = owner_masks[jg - 1]
        assert npma.all(
            test_var_np >= 1
        ), "Invalid value. (Non-owner cell was probably not overwritten."
        assert npma.all(
            test_var_np[mask] == (rank + 1)
        ), "Invaid value. Owner cell was overwritten."
