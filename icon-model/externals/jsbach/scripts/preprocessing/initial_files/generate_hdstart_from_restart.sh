#!/bin/bash

# ICON-Land
#
# ---------------------------------------
# Copyright (C) 2013-2026, MPI-M, MPI-BGC
#
# Contact: icon-model.org
# Authors: AUTHORS.md
# See LICENSES/ for license information
# SPDX-License-Identifier: BSD-3-Clause
# ---------------------------------------

#_____________________________________________________________________________
# Script to construct an initial condition file for the HD model from a
# restart file of an existing experiment.
#
# The restart can be either a single NetCDF file or a multi-file restart directory.
# No remapping between different grids is performed.
#
# Requirements: cdo and nco
#_____________________________________________________________________________

set -eu

# Either restart file or directory of multi-file restart files
src_restart=/work/mh1421/m220053/merge-to-icon-mpim/feature-new-r2b4-land-data-rev/experiments/jsbalone-intel/restart/jsbalone-intel_restart_jsbalone_19920101.nc

# In case of a single restart file: corresponding grid file
src_grid=/pool/data/ICON/grids/public/mpim/0049/icon_grid_0049_R02B04_G.nc

# Definition of the output file name
hd_label=r2b4_0049_sc_hfrac_s_v1  # label corresponding to hdpara file used in experiment
output_hdstart_file=hdstart_${hd_label}.nc

CDO="cdo -s -P 8"

function finish {
  rm temp?_$$.nc *ncap* >& /dev/null
}
trap finish EXIT

# Select HD reservoir variables
if [[ -d $src_restart ]]; then
  $CDO collgrid,hd_overlflow_res_box,hd_riverflow_res_box,hd_baseflow_res_box $src_restart/patch1_*.nc temp0_$$.nc
else
  $CDO setgrid,${src_grid} -selvar,hd_overlflow_res_box,hd_riverflow_res_box,hd_baseflow_res_box $src_restart temp0_$$.nc
fi

# Remove all global namelist attributes
ncatted -O -a nml_\*,global,d,c, temp0_$$.nc temp1_$$.nc

# Renaming dimensions
if [[ $(ncdump -h temp1_$$.nc | grep 'cells = ') != "" ]]; then
  ncrename -d cells,cell -d layers_1,bresnum -d layers_5,rresnum temp1_$$.nc temp2_$$.nc
else
  ncrename               -d layers_1,bresnum -d layers_5,rresnum temp1_$$.nc temp2_$$.nc
fi

# Remove time dimension
ncwa -O -a y,time temp2_$$.nc temp3_$$.nc
ncatted -O -a cell_methods,,d,, temp3_$$.nc

# Renaming variables
${CDO} chname,hd_riverflow_res_box,FRFMEM,hd_baseflow_res_box,FGMEM,hd_overlflow_res_box,FLFMEM temp3_$$.nc temp4_$$.nc
# In case if dimension cell has been renamed to x rename it again to cell
if [[ $(ncdump -h temp4_$$.nc | grep 'cell = ') == "" ]]; then
  ncrename -O -d x,cell temp4_$$.nc
fi

# Define dimension oresnum and copy FLFMEM to FLFMEM_new, with dimensions (oresnum,cell) and good long_name attribute
ncap2 -O -s 'defdim("oresnum",$bresnum.size); FLFMEM_new[$oresnum,$cell]=0.0;FLFMEM_new(:,:)=FLFMEM(:,:);FLFMEM_new@long_name="content of the overflow reservoir"' \
            temp4_$$.nc temp5_$$.nc
ncks -x -v FLFMEM temp5_$$.nc temp6_$$.nc
ncrename -v FLFMEM_new,FLFMEM temp6_$$.nc ${output_hdstart_file}

echo "---"
echo " Generation of ${output_hdstart_file} done"
echo "---"
