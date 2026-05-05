#!/usr/bin/env python3

# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

##
# @file test_raw_exchange.py
# @test
# A test for the raw exchange with python

import numpy as np
from scipy.sparse import csr_matrix

from yac import (
    YAC,
    def_calendar,
    Calendar,
    Reg2dGrid,
    Location,
    TimeUnit,
    Field,
    InterpolationStack,
    ConservNormalizationType,
    Reduction,
)

def_calendar(Calendar.PROLEPTIC_GREGORIAN)

yac = YAC()

yac.def_datetime("2020-01-01T00:00:00", "2020-01-02T00:00:00")

comp1, comp2 = yac.def_comps(["comp1", "comp2"])

grid1 = Reg2dGrid("grid1", [-1, 0, 1], [-1, 0, 1])
grid1.set_core_mask([1, 1, 0, 1], Location.CELL)
points1 = grid1.def_points(Location.CELL, [-0.5, 0.5], [-0.5, 0.5])

grid2 = Reg2dGrid("grid2", [-1, -1/3, 1/3, 1], [-1, -1/3, 1/3, 1])
points2 = grid2.def_points(Location.CELL, [-2/3, 0.0, 2/3], [-2/3, 0.0, 2/3])

field1 = Field.create("field1", comp1, points1, 1, "1", TimeUnit.HOUR)
field2 = Field.create("field2", comp2, points2, 1, "1", TimeUnit.HOUR)
field3 = Field.create("field3", comp2, points2, 1, "1", TimeUnit.HOUR)
field4 = Field.create("field4", comp2, points2, 1, "1", TimeUnit.HOUR)

interp = InterpolationStack()
interp.add_conservative(1, 0, 1, ConservNormalizationType.DESTAREA)
interp.add_fixed(42.0)

couple_kwargs = {"src_comp": "comp1",
                 "src_grid": "grid1",
                 "src_field": "field1",
                 "tgt_comp": "comp2",
                 "tgt_grid": "grid2",
                 "coupling_timestep": "60",
                 "timeunit": TimeUnit.MINUTE,
                 "time_reduction": Reduction.TIME_NONE,
                 "interp_stack": interp}

yac.def_couple(
    **couple_kwargs,
    tgt_field="field2",
    use_raw_exchange=True
)

yac.def_couple(
    **couple_kwargs,
    tgt_field="field3"
)

yac.def_couple(
    **couple_kwargs,
    tgt_field="field4",
    use_raw_exchange=True
)

yac.enddef()

raw_data = field2.get_raw_interp_weights_data()
indptr = np.insert(np.cumsum(raw_data.num_src_per_tgt), 0, 0)

raw_data_csr = field4.get_raw_interp_weights_data_csr()
W = csr_matrix((raw_data_csr.weights, raw_data_csr.src_idx, raw_data_csr.src_indptr))

buf = np.empty(shape=(1, raw_data.src_field_buffer_sizes[0]), dtype=np.float64)

i = 0.
while field2.datetime < yac.end_datetime:
    field1.put(np.arange(i, i + 4., dtype=np.float64))

    buf[:] = 0.0
    field2.get_raw(buf)
    f2_buf = np.empty(shape=(1, field2.size), dtype=np.float64)
    f2_buf[0, raw_data.tgt_idx_fixed] = raw_data.fixed_values[0]
    for tgt_idx, start_idx, end_idx in zip(raw_data.wgt_tgt_idx, indptr[:-1], indptr[1:]):
        f2_buf[0, tgt_idx] = raw_data.weights[start_idx:end_idx]@buf[0, raw_data.src_idx[start_idx: end_idx]]

    f3_buf, info = field3.get()
    assert np.allclose(f2_buf, f3_buf)

    field4.get_raw(buf)
    f4_buf = W@buf[0, :]
    f4_buf[raw_data_csr.tgt_idx_fixed] = raw_data_csr.fixed_values[0]
    assert np.allclose(f2_buf[0, :], f4_buf)

    i += 1
