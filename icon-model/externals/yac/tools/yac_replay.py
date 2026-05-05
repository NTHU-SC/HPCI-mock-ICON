#!/usr/bin/env python3
# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
import uxarray as ux
import xarray as xr
import yac
import argparse
import logging
from mpi4py import MPI # < some MPI impl. need this to be imported explicitly


parser = argparse.ArgumentParser("yac_replay",
                                 description="""Replay simulation from a dataset.
The given files are loaded with `uxarray.open_dataset`.
It is currently not supported to run this component in parallel.""")
parser.add_argument("gridfile", type=str,
                    default="path to the gridfile")
parser.add_argument("datafile", type=str, nargs='+',
                    default="path to the dataset")
parser.add_argument("--compname", type=str, default="replay",
                    help="Name for the yac component (default: replay)")
parser.add_argument("--gridname", type=str, default="replay_grid",
                    help="Name for the yac grid (default: replay_grid)")
parser.add_argument("--log-level", default=logging.WARNING, type=lambda x: getattr(logging, x.upper()),
                    help="Configure the logging level.")
parser.add_argument("--coupling-config", type=str,
                    help="If given the yac config is read with yac.read_config_yaml")
parser.add_argument("--cell-mask", type=str,
                    help="Expression for the cell mask used when defining the fields. All variables in the grid file can be used.")
parser.add_argument("--edge-mask", type=str,
                    help="Expression for the edge mask used when defining the fields. All variables in the grid file can be used.")
parser.add_argument("--vertex-mask", type=str,
                    help="Expression for the vertex mask used when defining the fields. All variables in the grid file can be used.")

args = parser.parse_args()

log_handler = logging.StreamHandler()
log_handler.setFormatter(logging.Formatter("%(asctime)-10s %(name)-10s %(levelname)-8s %(message)s"))
logging.basicConfig(level=args.log_level, handlers=[log_handler])

logging.info("Instantiate yac")
y = yac.YAC()

if args.coupling_config:
    y.read_config_yaml(args.coupling_config)

logging.info(f"open dataset: {args.gridfile=}, {args.datafile=}")
ds = ux.open_mfdataset(args.gridfile, args.datafile)

logging.info(f"define component: {args.compname}")
comp = y.def_comp(args.compname)

fields = []

cell_to_edge = np.asarray(ds.uxgrid.face_edge_connectivity).reshape((-1,), order="C")
cell_to_edge = cell_to_edge[cell_to_edge != ux.INT_FILL_VALUE]

grid = yac.UnstructuredGridEdge(args.gridname,
                                ds.uxgrid.n_nodes_per_face,
                                np.deg2rad(ds.uxgrid.node_lon),
                                np.deg2rad(ds.uxgrid.node_lat),
                                cell_to_edge, ds.uxgrid.edge_node_connectivity)

cell_point_id = grid.def_points(yac.Location.CELL,
                                np.deg2rad(ds.uxgrid.face_lon),
                                np.deg2rad(ds.uxgrid.face_lat))
edge_point_id = grid.def_points(yac.Location.EDGE,
                                np.deg2rad(ds.uxgrid.edge_lon),
                                np.deg2rad(ds.uxgrid.edge_lat))
vertex_point_id = grid.def_points(yac.Location.CORNER,
                                  np.deg2rad(ds.uxgrid.node_lon),
                                  np.deg2rad(ds.uxgrid.node_lat))

# define the masks
grid_ds = xr.open_dataset(args.gridfile)
if args.cell_mask is not None:
    mask_values = eval(args.cell_mask, globals=None, locals=grid_ds.variables)
    cell_mask = grid.def_mask(yac.Location.CELL,
                              mask_values, "replay_cell_mask")
else:
    cell_mask = None
if args.edge_mask is not None:
    mask_values = eval(args.edge_mask, globals=None, locals=grid_ds.variables)
    edge_mask = grid.def_mask(yac.Location.EDGE,
                              mask_values, "replay_edge_mask")
else:
    edge_mask = None
if args.vertex_mask is not None:
    mask_values = eval(args.vertex_mask, globals=None, locals=grid_ds.variables)
    vertex_mask = grid.def_mask(yac.Location.CORNER,
                                mask_values, "replay_vertex_mask")
else:
    vertex_mask = None

dt = np.diff(ds.coords["time"])
assert np.all(dt[0] == dt), "Time coordinates are not equidistant"
dt = dt[0]

for varname, var in ds.variables.items():
    logging.info(f"checking variable {varname}")
    if "time" not in var.dims:
        logging.info(f"{varname} was skipped: lacking time coordinate")
        continue
    spatial_axis = {"n_face", "n_edge", "n_node"} & set(var.dims)
    if len(spatial_axis) != 1:
        logging.info(
            f"{varname} was skipped: exactly one of n_face, n_edge or n_node must be a coordinate.")
        continue
    vaxis = set(var.dims) - {"time", "n_face", "n_edge", "n_node"}
    if len(vaxis) > 1:
        logging.info(f"{varname} was skipped: too many coordinates")
        continue
    if len(vaxis) > 0:
        zcoord, = vaxis
        collection_size = np.asarray(ds[zcoord]).shape[0]
    else:
        collection_size = 1

    spatial_axis, = spatial_axis
    point_id, mask = {"n_face": (cell_point_id, cell_mask),
                      "n_edge": (edge_point_id, edge_mask),
                      "n_node": (vertex_point_id, vertex_mask)}[spatial_axis]
    fields.append((
        yac.Field.create(varname,
                         comp,
                         point_id,
                         collection_size,
                         str(int(dt / np.timedelta64(1, 'ms'))), yac.TimeUnit.MILLISECOND, mask),
        var.transpose("time", *list(vaxis), spatial_axis)
    ))

logging.info("calling enddef")
y.enddef()

fields = [(field, var) for field, var in fields if field.role == yac.ExchangeType.SOURCE]
logging.info(f"{len(fields)} fields are coupled")

t0 = np.searchsorted(ds.coords["time"], np.datetime64(y.start_datetime))
assert ds.coords["time"][t0] == np.datetime64(y.start_datetime), "starttime not in time axis"

while True:
    for field, var in fields:
        logging.info(f"processing {field.name} at {field.datetime}")
        t = np.datetime64(field.datetime)
        t_idx = np.searchsorted(ds["time"], t)
        assert ds["time"][t_idx] == t
        info = field.put(var.isel(time=t_idx).data)
    if info == yac.Action.PUT_FOR_RESTART:
        break

logging.info("done")
