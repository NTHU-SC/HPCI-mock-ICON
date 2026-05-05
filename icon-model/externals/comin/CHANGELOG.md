# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 0.4.0 (2025-11-24)

### Added

- Child/parent relationship between domains (parent_id, child_id) added. !457
- Python plugins are now able to import modules in the same directory. !451
- Add a CHANGELOG. !454
- Documentation: fix pointer update for multi-timelevel fields
- Documentation: specify an external microphysics scheme by a ComIn plugin.
- Documentation: specify gas concentrations for ecRad radiation in a ComIn plugin.
- User guide: Add section for gathering publications on and with ComIn
- Update authors file
- Documentation: update user guide documentation on GPUs. !474
- Add utility function to convert a variable to an xarray.DataArray (!476)
- Documentation: Update ComIn standalone setup documentation for Levante (!478)
- Documentation: ask new users to "star" the Gitlab repo. (!477)
- Documentation: add WindWaker plugin to applist. (!481)

### Changed

- Variable lookups for metadata, etc., use hash tables. !462
- Testing: Reference are written out by the process itself instead of relying on MPI !482

### Fixed

- add missing include in comin_keyval.cpp !468
- Error handling if an error occured outside of an entry point in comin. !453
- Add YAC entrypoints to plugin interface
- workaround for NVIDIA deep copy problem. !472
- Fix `number_of_grid_used` descriptive data. !473

### Removed

- Remove `comin.plugin_finish` from python API. !450

## 0.3.0 (2025-04-04)

No changelog before v0.3.0
