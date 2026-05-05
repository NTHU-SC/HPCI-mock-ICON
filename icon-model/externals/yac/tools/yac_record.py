#!/usr/bin/env python3
# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
import uxarray as ux
import yac
import argparse
import logging
from isodate import parse_duration
from mpi4py import MPI # < some MPI impl. need this to be imported explicitly

parser = argparse.ArgumentParser("yac_replay",
                                 description="""Replay simulation from a dataset.
The given files are loaded with `uxarray.open_dataset`.
It is currently not supported to run this component in parallel.""")
parser.add_argument("gridfile", type=str,
                    default="path to the gridfile")
parser.add_argument("datafile", type=str,
                    default="path to the dataset")
parser.add_argument("variables", type=str, nargs="+",
                    help="Variable names to be recorded")
parser.add_argument("--compname", type=str, default="record",
                    help="Name for the yac component (default: record)")
parser.add_argument("--gridname", type=str, default="record_grid",
                    help="Name for the yac grid (default: record_grid)")
parser.add_argument("--log-level", default=logging.WARNING, type=lambda x: getattr(logging, x.upper()),
                    help="Configure the logging level.")
parser.add_argument("--coupling-config", type=str, default=None,
                    help="If given the yac config is read with yac.read_config_yaml")
parser.add_argument("--points", type=str, nargs="*", choices=("cell", "edge", "vertex"),
                    default=[],
                    help="can be used to specify the entities on which the data should be recorded.")

args = parser.parse_args()
log_handler = logging.StreamHandler()
log_handler.setFormatter(logging.Formatter("%(asctime)-10s %(name)-10s %(levelname)-8s %(message)s"))
logging.basicConfig(level=args.log_level, handlers=[log_handler])

logging.info("Instantiate yac")
y = yac.YAC()

logging.info(f"reading config file {args.coupling_config!r}")

if args.coupling_config:
    y.read_config_yaml(args.coupling_config)

logging.info(f"open dataset: {args.gridfile!r}, {args.datafile!r}")

uxgrid = ux.open_grid(args.gridfile)

logging.info(f"define component: {args.compname}")
comp = y.def_comp(args.compname)

fields = []


cell_to_edge = np.asarray(uxgrid.face_edge_connectivity).reshape((-1,), order="C")
cell_to_edge = cell_to_edge[cell_to_edge != ux.INT_FILL_VALUE]

logging.info("defining grid")
grid = yac.UnstructuredGridEdge(args.gridname,
                                uxgrid.n_nodes_per_face,
                                np.deg2rad(uxgrid.node_lon),
                                np.deg2rad(uxgrid.node_lat),
                                cell_to_edge, uxgrid.edge_node_connectivity)

cell_point_id = grid.def_points(yac.Location.CELL,
                                np.deg2rad(uxgrid.face_lon),
                                np.deg2rad(uxgrid.face_lat))
edge_point_id = grid.def_points(yac.Location.EDGE,
                                np.deg2rad(uxgrid.edge_lon),
                                np.deg2rad(uxgrid.edge_lat))
vertex_point_id = grid.def_points(yac.Location.CORNER,
                                  np.deg2rad(uxgrid.node_lon),
                                  np.deg2rad(uxgrid.node_lat))

if len(args.points) == 1:
    args.points = len(args.variables)*[args.points[0]]
args.points = args.points or len(args.variables)*["cell"]
assert len(args.points) == len(args.variables), "If points are speified it must be given once or the same number as variables"

logging.info("calling sync_def")

y.sync_def()

for varname in args.variables:
    assert y.get_field_role(args.compname, args.gridname, varname) == yac.ExchangeType.TARGET, \
        f"{varname!r} is not coupled to a target. Please define a coupling or remove it from the variables."

start = np.datetime64(y.start_datetime)
end = np.datetime64(y.end_datetime)


first_source = y.get_field_source(args.compname, args.gridname, args.variables[0])
dt = np.timedelta64(parse_duration(y.get_field_timestep(*first_source)))

logging.info(f"Timestep is {dt}")

# create the dataset
ds = ux.UxDataset(uxgrid=uxgrid)
ds["time"] = np.arange(start, end+dt, dt)

for varname, points in zip(args.variables, args.points):
    logging.info(f"configuring {varname!r}")
    source = y.get_field_source(args.compname, args.gridname, varname)
    logging.info(f"coupled to {source!r}")
    assert dt == np.timedelta64(parse_duration(y.get_field_timestep(*source))), \
        "Source fields of variables have incosistent timesteps"
    collection_size = y.get_field_collection_size(*source)
    logging.info(f"collection size: {collection_size}")
    if collection_size > 1:
        zaxis = [f"zaxis_{varname}"]
    else:
        zaxis = []
    logging.info(f"zaxis: {zaxis}")
    spatial_dim, spatial_size, point_id = \
        {"cell": ("n_face", uxgrid.n_face, cell_point_id),
         "edge": ("n_edge", uxgrid.n_edge, edge_point_id),
         "vertex": ("n_node", uxgrid.n_node, vertex_point_id)}[points]
    data = np.nan*np.ones((len(ds["time"]), *(len(zaxis)*[collection_size]), spatial_size))
    ds = ds.assign({varname:
                    (("time", *zaxis, spatial_dim), data)})
    fields.append((
        yac.Field.create(varname,
                         comp,
                         point_id,
                         collection_size,
                         str(int(dt / np.timedelta64(1, 'ms'))), yac.TimeUnit.MILLISECOND),
        ds[varname]
    ))

logging.info("calling enddef")
y.enddef()

while True:
    for field, var in fields:
        logging.info(f"processing {field.name} at {field.datetime}")
        t = np.datetime64(field.datetime)
        t_idx = np.searchsorted(ds["time"], t)
        assert ds["time"][t_idx] == t, f"{t} not found in timeaxis"
        data, info = field.get()
        var.isel(time=t_idx).data[:] = data
        logging.info(f"info: {info!r}")
    if info == yac.Action.GET_FOR_RESTART:
        break

ds.to_netcdf(args.datafile)
