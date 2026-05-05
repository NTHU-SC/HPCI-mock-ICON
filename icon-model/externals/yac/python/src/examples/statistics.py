#!/usr/bin/env python3

# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

from yac import *
import numpy as np

class Statistics:

    def __init__(self, variables, comp_name = "statistics",
                 grid_name = "stat_grid",
                 bounds = [0, 2*np.pi, -0.5*np.pi, 0.5*np.pi],
                 resolution = (360,180),
                 yac = None):

        self.yac = yac or YAC.default_instance
        self.comp_name = comp_name
        self.grid_name = grid_name
        self.comp = self.yac.predef_comp(comp_name)
        self.variables = variables
        self.bounds = bounds
        self.resolution = resolution

    def setup(self):
        self.x = np.linspace(self.bounds[0],self.bounds[1],self.resolution[0])
        self.y = np.linspace(self.bounds[2],self.bounds[3],self.resolution[1])
        grid = Reg2dGrid(self.grid_name, self.x, self.y)
        self.points = grid.def_points(Location.CORNER, self.x, self.y)

    def def_couples(self):
        nnn = InterpolationStack()
        nnn.add_nnn(NNNReductionType.AVG, 1, 0., 1.)

        self.fields = []
        for var in self.variables:
            timestep = self.yac.get_field_timestep(*var)
            collection_size = self.yac.get_field_collection_size(*var)
            self.fields.append(Field.create(var[2], self.comp, self.points, 1,
                                            timestep, TimeUnit.ISO_FORMAT))
            self.yac.def_couple(*var,
                                self.comp_name, self.grid_name, var[2],
                                timestep, TimeUnit.ISO_FORMAT, 0, nnn)

    def step(self):
        for field in self.fields:
            time = field.datetime
            buf, info = field.get()
            print(f"{self.comp_name}:{self.grid_name}:{field.name} at {time}: "
                  f"min {np.min(buf)} max {np.max(buf)} mean {np.mean(buf)}")
        return field.datetime
