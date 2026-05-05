#!/usr/bin/env python3

# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

from yac import (YAC, Reg2dGrid, Field, Location, TimeUnit)
import numpy as np


class NoiseGenerator:
    def __init__(self, timestep, yac=None, collection_size=1):
        self.yac = yac or YAC.default_instance
        self.comp = self.yac.predef_comp("noisegenerator")
        self.timestep = timestep
        self.collection_size = collection_size

    def setup(self):
        global grid, noise_field
        x = np.linspace(0, 2*np.pi, 360)
        y = np.linspace(-0.5*np.pi, 0.5*np.pi, 180)
        grid = Reg2dGrid("noise_grid", x, y)
        points = grid.def_points(Location.CORNER, x, y)

        noise_field = Field.create("noise", self.comp, points, self.collection_size,
                                   self.timestep, TimeUnit.ISO_FORMAT)

    def def_couples(self):
        pass

    def step(self):
        print("NoiseGenerator: Lets make some noise!!!")
        noise_field.put(np.random.rand(self.collection_size, grid.nbr_corners))
        return noise_field.datetime
