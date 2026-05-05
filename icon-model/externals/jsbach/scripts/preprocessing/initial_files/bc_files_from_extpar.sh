#!/bin/ksh

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

#-----------------------------------------------------------------------------
# Generate ICON-Land boundary condition (bc) files based on extpar data
# and other sources.
#
# The script is part of a script series to generate jsbach initial files
# started by master script "create_jsbach_ini_files.sh".
#
#-----------------------------------------------------------------------------
set -e

# File names and paths exported from master script create_icon-land_ini_files.sh
path_extpar=${extpar_dir}        # Directory name of the extpar file
extpar_file=${extpar_file}       # Extpar file name

icon_grid=${icon_grid}
atmGridID=${atmGridID}
refinement=${refinement}

path_bc=${bc_file_dir}
year_list=${year_list}
work_dir=${work_dir}

export pool_prepare=/pool/data/JSBACH/prepare/

prog=$(basename $0)
scripts_dir=${scripts_dir}
clean_up=true  # Remove intermediate files

# Tools and commands
cdo="$cdo -f nc4 -b F64"
export rm=/usr/bin/rm
export cp=/usr/bin/cp

# Which of the initial files should be generated?
generate_bc_land_frac=true
generate_bc_land_soil=true
generate_bc_land_phys=true
generate_bc_land_sso=true
generate_ic_land_soil=true

# Information for history attribute
#vg git_repo=$(git remote -v | head -1 | cut -f2 | cut -f1 -d' ')
#vg git_rev=$(git rev-parse --short HEAD)
#vg git_branch=$(git log --pretty='format:%h %D' --first-parent | grep HEAD | cut -f4 -d' ' | sed 's/,//')
#vg history_att="$(date): Generated with $prog (${git_repo}:${git_branch} rev. ${git_rev}) by $(whoami)"

min_fract=${min_fract} # minimum land grid cell fraction other than 0.

start_year=$(echo ${year_list} | cut -f1 -d" ")
year=${start_year}

[[ -d ${path_bc} ]]   || mkdir -p ${path_bc}
[[ -d ${work_dir} ]]  || mkdir -p ${work_dir}
cd ${work_dir}
echo '------------------------------------------------------------------'
echo "${prog}: Working directory: $(pwd)"
echo "${prog}: Output  directory: ${path_bc}"
echo "${prog}: Extpar file used: ${path_extpar}/${extpar_file}"
echo "${prog}: Generation of boundary condition files: start_year=${start_year}"
echo '------------------------------------------------------------------'

ln -fs ${path_extpar}/${extpar_file}   ${extpar_file}

#----------------------------------------------------------------------------------------
# Preparations
#----------------------------------------------------------------------------------------
#
# -- Calculate masks --
#
# notsea: land + lake mask
# -------------------------
if [[ ${coupled} == "true" ]]; then
  # Get fractional mask from the fractional mask file depending on the ocean grid
  ${cdo} chname,cell_sea_land_mask,notsea -setgrid,${icon_grid} \
    -setrtoc,${max_fract},1.1,1 -setrtoc,-1,${min_fract},0 ${fractional_mask} notsea.nc
else
  # Get fractional mask from extpar
  if [[ $(cdo -s showvar ${extpar_file} | grep 'FR_LAND') != "" ]]; then
    # Extpar file contains FR_LAND and FR_LAKE
    ${cdo} expr,'notsea=FR_LAND+FR_LAKE' ${extpar_file} NOTSEA_extpar.nc
    if [[ $(${cdo} outputf,%1.8g -fldmax NOTSEA_extpar.nc) > 1. ]]; then
      echo " ERROR: Extpar file ${initial_extpar_file} "
      echo "     is outdated: FR_LAND and FR_LAKE do not match."
      echo "     Please switch to a newer version or generate the file anew."
      exit 1
    fi
    ${cdo} -setrtoc,-1,${min_fract},0 -setrtoc,${max_fract},1.1,1 -setmisstoc,0 \
        -setgrid,${icon_grid} NOTSEA_extpar.nc notsea.nc
  elif [[ $(cdo -s showvar ${extpar_file} | grep 'FR_LAND_TOPO') != "" ]]; then
    # Extpar file contains FR_LAND_TOPO including lake fractions (typically in older MPIM extpar file)
    ${cdo} -setname,notsea -setrtoc,-1,${min_fract},0 -setrtoc,${max_fract},1.1,1 \
        -setgrid,${icon_grid} -selvar,FR_LAND_TOPO ${extpar_file} notsea.nc
  else
    echo "Land sea mask in ${extpar_file} not found. It is neither FR_LAND nor FR_LAND_TOPO"
    exit 1
  fi
fi
ncatted -a long_name,notsea,m,c,'Fraction of land+lake' notsea.nc

# glac: integer glacier mask
# ---------------------------
#
# The glacier mask in extpar is fractional. We need a 1/0 glacier mask: 1 if the glacier fraction
# relative to the land fraction is greater than 0.5; 0 otherwise.
# We have to consider, that the land sea mask in extpar is generally different from the land sea
# mask generated here for the specific model setup.
${cdo} -setgrid,${icon_grid} -selvar,ICE ${extpar_file} ICE_extpar.nc
${cdo} -setgrid,${icon_grid} -expr,'notsea=FR_LAND+FR_LAKE' ${extpar_file} NOTSEA_extpar.nc
# Glacier fraction rel. to the original extpar data land fraction
${cdo} div ICE_extpar.nc NOTSEA_extpar.nc ICE_extpar_rel-land.nc
if [[ ${coupled} == true ]]; then
  # Extrapolation of the glacier mask to all non-land cells
  ${cdo} -gtc,0.5 -setmisstodis ICE_extpar_rel-land.nc ICE_extpar.filled.nc
   # Adaptation of the glacier mask to the land sea mask of the coupled setup
  ${cdo} setmisstoc,0 -ifthen notsea.nc ICE_extpar.filled.nc glac.tmp
else
  ${cdo} setmisstoc,0 -ifthen notsea.nc -gtc,0.5 ICE_extpar_rel-land.nc glac.tmp
fi
${cdo} -setvar,glac -setmisstoc,0 -gec,0.5 glac.tmp glac.nc

# notglac: not glacier mask
# ------------------
${cdo} -setvar,notglac -addc,1. -mulc,-1. glac.nc notglac.nc

# non-glac-land-frac: fractional mask of non-glacier land
# ------------------
${cdo} -expr,"nonglacland=max(0,notsea)" -sub notsea.nc glac.nc non-glac-land-frac.nc

# non-glac-land-mask: 1: cells with non-glacier land fraction, 0 otherwise
# ------------------
${cdo} setvar,nonglacland -gtc,0 non-glac-land-frac.nc non-glac-land-mask.nc

# land: 1: cell with land fraction, 0 otherwise
# ------------------
${cdo} setvar,land -gtc,0 -add non-glac-land-mask.nc glac.nc land.nc

# ocean-or-glac: 1: complete ocean or glacier cells, 0: non-glac land
# ------------------
${cdo} setvar,ocean-or-glac -addc,1. -mulc,-1. non-glac-land-mask.nc ocean-or-glac.nc

# Helper array: zero everywhere (cdo const does not work with irregular grids.)
${cdo} -setname,constant -sub notsea.nc notsea.nc zero.nc

#
# -- Calculate remapping weights --
#
# Note:
# - All 05 deg grids are expected to be identical (lon: 0.25 -> 359.75; lat: -89.75 -> 89.75)
# - The non-glacier-land mask is expected to be identical for all variables (additional
#   missing values do not matter)
# - Depression storage would be available on resolutions up to 15 arcsec if needed
#
# Remapping weights: regular 0.5 degree grid -> Icon grid
# These weights are used with
#  - ${pool_prepare}/05/vegmax_6_05_0-360.nc
#  - ${pool_prepare}/05/soil_parameters_05.nc
#  - ${pool_prepare}/05/soil_parameter_05_apr2024.nc          -> invertlat
#  - ${pool_prepare}/05/wise_som-fract_05deg_jsbsoillayers.nc -> sellonlatbox...
#  - ${pool_prepare}/05/ECHAM6/05_jan_surf.nc
#  - ${pool_prepare}/05/albedo_05.nc                          -> sellonlatbox...
#  - ${pool_prepare}/05/depr_stor_cop_05.nc                   -> invertlat

bisect=${refinement#*B}
if [[ $bisect -ge 8  ]]; then
  ulimit -s unlimited
  ulimit -v unlimited
fi
if [[ $bisect -le 5 ]]; then
  export remap_scheme=ycon
else
  export remap_scheme=dis
fi

# 0.5 degree grid non-glacier-land mask (used for extrapolation)
${cdo} setname,nonglacland -invertlat -setmisstoc,0 -ifthenc,1 \
    -selvar,heat_capacity ${pool_prepare}/05/soil_parameter_05_apr2024.nc \
    05_non-glac-land.nc

# Remapping weights
${cdo} gen${remap_scheme},${icon_grid} 05_non-glac-land.nc wgt05.nc

# Clean up
if [[ ${clean_up} == true ]]; then
  ${rm} *.tmp*
  ${rm} ICE_extpar*.nc NOTSEA_extpar.nc
fi

# ----------------------------------------------------------------------------
#  1. bc fract files
# -----------------------------
if [[ ${generate_bc_land_frac} == true ]]; then
  ${scripts_dir}/generate_bc_land_fractions.sh
  cp ${scripts_dir}/generate_bc_land_fractions.sh ${bc_file_dir}/scripts
fi

# ----------------------------------------------------------------------------
#  2. bc soil file
# -----------------------------
if [[ ${generate_bc_land_soil} == true ]]; then
  ${scripts_dir}/generate_bc_land_soil.sh
  cp ${scripts_dir}/generate_bc_land_soil.sh ${bc_file_dir}/scripts
fi

# ----------------------------------------------------------------------------
#  3. bc physics file
# -----------------------------
if [[ ${generate_bc_land_phys} == true ]]; then
  ${scripts_dir}/generate_bc_land_phys.sh
  cp ${scripts_dir}/generate_bc_land_phys.sh ${bc_file_dir}/scripts
fi

# ----------------------------------------------------------------------------
#  4. bc surface orography file
# -----------------------------
if [[ ${generate_bc_land_sso} == true ]]; then
  ${scripts_dir}/generate_bc_land_sso.sh
  cp ${scripts_dir}/generate_bc_land_sso.sh ${bc_file_dir}/scripts
fi

# ----------------------------------------------------------------------------
#  5. ic soil file
# -----------------------------
if [[ ${generate_ic_land_soil} == true ]]; then
  ${scripts_dir}/generate_ic_land_soil.sh
  cp ${scripts_dir}/generate_ic_land_soil.sh ${bc_file_dir}/scripts
fi

# ----------------------------------------------------------------------------
# Clean up

if [[ ${clean_up} == true ]]; then
 ${rm} ${extpar_file}
 ${rm} -f ${atmGridID}_*.nc 05_*.nc
 cd ..
 #vg rm -rf ${work_dir}
fi

echo "${prog}:  done"
echo "------------------------------------------------------------------"

exit 0
