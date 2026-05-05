"""
Test plugin for the ICON Community Interface (ComIn)

Creates variables with different datatypes and checks the corresponding buffer

@authors 12/2024 :: ICON Community Interface  <comin@icon-model.org>

SPDX-License-Identifier: BSD-3-Clause

Please see the file LICENSE in the root of the source tree for this code.
Where software is supplied by third parties, it is indicated in the
headers of the routines.
"""

import comin
import numpy as np

# types to test
comin_datatypes = {
    "dp": comin.COMIN_VAR_DATATYPE_DOUBLE,
    "sp": comin.COMIN_VAR_DATATYPE_FLOAT,
    "i": comin.COMIN_VAR_DATATYPE_INT,
}
numpy_datatypes = {
    "dp": np.double,
    "sp": np.single,
    "i": np.intc,
}

for name, comin_dt in comin_datatypes.items():
    descr = ("test_var_" + name, 1)
    comin.var_request_add(descr, False)
    comin.metadata_set(descr, datatype=comin_dt)


var_handles = {}


@comin.EP_SECONDARY_CONSTRUCTOR
def simple_python_constructor():
    global var_handles
    flag = comin.COMIN_FLAG_READ | comin.COMIN_FLAG_WRITE
    var_handles = {
        name: comin.var_get(
            [comin.EP_ATM_TIMELOOP_START], ("test_var_" + name, 1), flag
        )
        for name, comin_dt in comin_datatypes.items()
    }


@comin.EP_ATM_TIMELOOP_START
def check():
    for name, var_handle in var_handles.items():
        assert (
            np.asarray(var_handle).dtype == numpy_datatypes[name]
        ), f"Check for {name} unsuccessful ({np.asarray(var_handle).dtype} != {numpy_datatypes[name]})"
        comin.print_info(f"Check for {name} successful ({numpy_datatypes[name]})")
