# @authors 11/2023 :: ICON Community Interface  <comin@icon-model.org>
#
# SPDX-License-Identifier: BSD-3-Clause
#
# Please see the file LICENSE in the root of the source tree for this code.
# Where software is supplied by third parties, it is indicated in the
# headers of the routines.

from _comin import *  # noqa: F403

import _comin
from dataclasses import dataclass
import shlex as _shlex
import sys as _sys
import collections.abc as _abc


def register_callback(ep):
    def __callback(fun):
        _comin._callback_register(ep, fun)
        return fun

    return __callback


COMIN_ZAXIS_UNDEF = -1
COMIN_ZAXIS_NONE = 0
COMIN_ZAXIS_2D = 1
COMIN_ZAXIS_3D = 2
COMIN_ZAXIS_3D_HALF = 3


class _variable:
    def __init__(self, handle):
        self._handle = handle
        try:
            import numpy as _np

            self.np = _np
        except ImportError:
            print("Warning: cant import numpy", file=_sys.stderr)

    def __array__(self):
        return self.np.asarray(_comin._var_get_buffer(self._handle))

    @property
    def __cuda_array_interface__(self):
        host_buf = _comin._var_get_buffer(self._handle)
        return {
            "shape": host_buf.shape,
            "typestr": self.np.dtype(host_buf.format).str,
            "data": (_comin._var_get_device_ptr(self._handle), False),
            "version": 3,
            "strides": host_buf.strides,
        }

    @property
    def dim_semantics(self):
        return _comin._var_get_dim_semantics(self._handle)

    def to_xarray(self):
        import xarray as xr

        unused_dims = tuple(
            i
            for i, s in enumerate(self.dim_semantics)
            if s == COMIN_DIM_SEMANTICS_UNUSED
        )
        data = self.np.squeeze(self.np.asarray(self), axis=unused_dims)
        return xr.DataArray(
            data,
            dims=tuple(
                _comin_dim_semantics_get_name(s)
                for s in self.dim_semantics
                if s != COMIN_DIM_SEMANTICS_UNUSED
            ),
            name=self.descriptor[0],
            attrs=dict(metadata(self.descriptor)),
        )

    @property
    def ncontained(self):
        return _comin._var_get_ncontained(self._handle)

    @property
    def to_3d(self):
        ds = self.dim_semantics
        if COMIN_DIM_SEMANTICS_LEVEL in ds:
            pos = (
                ds.index(COMIN_DIM_SEMANTICS_NPROMA),
                ds.index(COMIN_DIM_SEMANTICS_LEVEL),
                ds.index(COMIN_DIM_SEMANTICS_BLOCK),
            )
        else:
            pos = (
                ds.index(COMIN_DIM_SEMANTICS_NPROMA),
                ds.index(COMIN_DIM_SEMANTICS_BLOCK),
            )
        missing_dims = {*range(5)}.difference({*pos})
        return self.np.asarray(self).transpose(*pos, *missing_dims)[..., 0, 0]

    @property
    def descriptor(self):
        return _comin._var_get_descriptor(self._handle)


def var_get(context, var_descriptor, flag):
    """get variable object, arguments: [entry point], (name string, domain id), access flag)"""
    return _variable(_comin._var_get(context, var_descriptor, flag))


class _EntryPoint:
    def __init__(self, ep_int):
        self.ep = ep_int

    def __call__(self, fun):
        _comin._callback_register(self.ep, fun)
        return fun

    def __index__(self):
        return self.ep


for ep in range(1, _comin._EP_DESTRUCTOR() + 1):
    name = _comin.callback_get_ep_name(ep)
    vars()[name] = _EntryPoint(ep)


@dataclass
class plugin_info:
    id: int
    name: str
    options: str
    comm: str

    @property
    def args(self):
        """
        Extract the argument from the options string like as the
        script was called from the command line. This is supposed to
        be passed to `argparse.ArgumentParser.parse_args`.
        """
        return _shlex.split(self.options)[1:]


def current_get_plugin_info():
    """returns object describing the current plugin"""
    return plugin_info(**_comin._current_get_plugin_info())


class _descrdata:
    def __init__(self, properties, jg=0):
        self.properties = properties
        self.jg = jg

    def __dir__(self):
        return self.properties.keys()

    def __getattr__(self, key):
        val = _comin._descrdata_eval_property(self.properties[key], jg=self.jg)
        if isinstance(val, dict):
            return _descrdata(val, jg=self.jg)
        else:
            return val


def descrdata_get_domain(jg):
    """returns descriptive data for a given domain, arguments: jg"""
    return _descrdata(_comin._descrdata_get_domain(), jg=jg)


def descrdata_get_global():
    """returns global descriptive data object"""
    return _descrdata(_comin._descrdata_get_global())


def var_descr_list():
    """List of exposed variables (descriptors)"""
    current = _comin._var_get_descr_list_head()
    try:
        while current is not None:
            yield _comin._var_get_descr_list_var_desc(current)
            current = _comin._var_get_descr_list_next(current)
    finally:
        if current is not None:
            _comin._var_descr_list_iterator_delete(current)


def metadata_set(var_descriptor, **kwargs):
    """sets metadata for a requested field, arguments: name string, domain id, metadata key, metadata value"""
    for n, v in kwargs.items():
        _comin._metadata_set(var_descriptor, n, v)


class metadata(_abc.Mapping):
    def __init__(self, var_descr):
        self.descr = var_descr

    def __getitem__(self, key):
        return _comin.metadata_get(self.descr, key)

    def __iter__(self):
        it = _comin._metadata_get_iterator_begin(self.descr)
        end = _comin._metadata_get_iterator_end(self.descr)
        while not _comin._metadata_iterator_compare(it, end):
            yield _comin._metadata_iterator_get_key(it)
            _comin._metadata_iterator_next(it)
        _comin._metadata_iterator_delete(it)
        _comin._metadata_iterator_delete(end)

    def __len__(self):
        # impl. could be improved
        return len(self.__iter__)


@dataclass
class simulation_interval:
    exp_start: str
    exp_stop: str
    run_start: str
    run_stop: str


def descrdata_get_simulation_interval():
    """returns simulation intervals: exp_start, exp_stop, run_start, run_stop"""
    return simulation_interval(**_comin._descrdata_get_simulation_interval())


COMIN_FLAG_NONE = 0
COMIN_FLAG_READ = 1 << 1
COMIN_FLAG_WRITE = 1 << 2
COMIN_FLAG_SYNC_HALO = 1 << 3
COMIN_FLAG_DEVICE = 1 << 4

COMIN_HGRID_UNSTRUCTURED_CELL = 1
COMIN_HGRID_UNSTRUCTURED_EDGE = 2
COMIN_HGRID_UNSTRUCTURED_VERTEX = 3

COMIN_VAR_DATATYPE_DOUBLE = 1
COMIN_VAR_DATATYPE_FLOAT = 2
COMIN_VAR_DATATYPE_INT = 3

COMIN_DIM_SEMANTICS_UNDEF = 1
COMIN_DIM_SEMANTICS_NPROMA = 2
COMIN_DIM_SEMANTICS_BLOCK = 3
COMIN_DIM_SEMANTICS_UNBLOCK = 4
COMIN_DIM_SEMANTICS_LEVEL = 5
COMIN_DIM_SEMANTICS_CONTAINER = 6
COMIN_DIM_SEMANTICS_OTHER = 7
COMIN_DIM_SEMANTICS_UNUSED = 8


def _comin_dim_semantics_get_name(ds):
    return {
        COMIN_DIM_SEMANTICS_UNDEF: "undef",
        COMIN_DIM_SEMANTICS_NPROMA: "nproma",
        COMIN_DIM_SEMANTICS_BLOCK: "block",
        COMIN_DIM_SEMANTICS_UNBLOCK: "unblock",
        COMIN_DIM_SEMANTICS_LEVEL: "level",
        COMIN_DIM_SEMANTICS_CONTAINER: "container",
        COMIN_DIM_SEMANTICS_OTHER: "other",
        COMIN_DIM_SEMANTICS_UNUSED: "unused",
    }[ds]
