import comin

import numpy as np
from itertools import product

glob = comin.descrdata_get_global()
domains = [comin.descrdata_get_domain(jg) for jg in range(1, glob.n_dom + 1)]

test_space = product(
    range(1, glob.n_dom + 1),
    [comin.COMIN_ZAXIS_2D, comin.COMIN_ZAXIS_3D, comin.COMIN_ZAXIS_3D_HALF],
    [
        comin.COMIN_HGRID_UNSTRUCTURED_CELL,
        comin.COMIN_HGRID_UNSTRUCTURED_EDGE,
        comin.COMIN_HGRID_UNSTRUCTURED_VERTEX,
    ],
)

for jg, zaxis, hgrid in test_space:
    descr = (f"{zaxis}_{hgrid}", jg)
    comin.var_request_add(descr, False)
    comin.metadata_set(descr, zaxis_id=zaxis, hgrid_id=hgrid)


@comin.EP_SECONDARY_CONSTRUCTOR
def secondary_constructor():
    global var_handles
    var_handles = {
        (jg, zaxis, hgrid): comin.var_get(
            [comin.EP_ATM_TIMELOOP_BEFORE],
            (f"{zaxis}_{hgrid}", jg),
            comin.COMIN_FLAG_READ,
        )
        for jg, zaxis, hgrid in test_space
    }


@comin.EP_ATM_TIMELOOP_BEFORE
def check_shapes():
    global var_handles
    for jg, zaxis, hgrid in test_space:
        var = var_handles[(jg, zaxis, hgrid)]
        shape = np.shape(var)
        dim_sema = var.dim_semantics

        assert shape[dim_sema.index(comin.COMIN_DIM_SEMANTICS_NPROMA)] == glob.nproma

        assert (
            shape[dim_sema.index(comin.COMIN_DIM_SEMANTICS_BLOCK)]
            == {
                comin.COMIN_HGRID_UNSTRUCTURED_CELL: domains[jg].cells.nblks,
                comin.COMIN_HGRID_UNSTRUCTURED_EDGE: domains[jg].edges.nblks,
                comin.COMIN_HGRID_UNSTRUCTURED_VERTEX: domains[jg].verts.nblks,
            }[hgrid]
        )

        dim = len(np.squeeze(var).shape)
        if zaxis == comin.COMIN_ZAXIS_2D:
            assert comin.COMIN_DIM_SEMANTICS_LEVEL not in dim_sema
            assert dim == 2
        else:
            assert (
                shape[dim_sema.index(comin.COMIN_DIM_SEMANTICS_LEVEL)]
                == {
                    comin.COMIN_ZAXIS_3D: glob.nlev,
                    comin.COMIN_ZAXIS_3D: glob.nlev + 1,
                }[zaxis]
            )
            assert dim == 3
