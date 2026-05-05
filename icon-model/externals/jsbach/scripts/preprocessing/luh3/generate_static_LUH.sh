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

#--------------------------------------------------------------------
### Batch Queuing System is SLURM
#SBATCH --output=generate_static_LUH.o%j
#SBATCH --error=generate_static_LUH.o%j
#SBATCH --account=mj0143
#SBATCH --partition=interactive
#SBATCH --mem-per-cpu=3000
#SBATCH --exclusive    #-- exclusive required when working with the transitions file
#SBATCH --ntasks=1

# script for input preprocessing
# set -x
module unload cdo
module load cdo/2.5.0-gcc-11.2.0

scriptName="generate_static_LUH.sh"
echo "--> ${scriptName}"

# change into the directory where this script is located; thus, you can start this script from anywhere (no need to cd into a specific directory)
cd "$( dirname "$0" )"
path_to_postproc_maps_dir=$(pwd)
# set cdo command to be silent, create no history and use double precision
cdo="cdo -s --no_history -b 64"

# path and name of the target output folder
output_bc_static_dir=YOUROUTPUTFOLDER
grid_choice=GRID_CHOICE

# path to input files
input_bc_static_dir=/work/mj0143/icon_land_data_pool/preprocessing/land_use_data/rev001
input_bc_static_crops_file=${input_bc_static_dir}/IQ_croptypes_r360x720.nc
input_bc_static_slash_file=${input_bc_static_dir}/slash_fractions_map_BEF_medium_gridcell.nc
r2b4=/pool/data/ICON/grids/public/mpim/0049/icon_grid_0049_R02B04_G.nc
r2b5=/pool/data/ICON/grids/public/edzw/icon_grid_0030_R02B05_G.nc

# selection of the target grid file based on grid_choice
if [[ "$grid_choice" == "r2b4" ]]; then
  grid_file=${r2b4}
  echo "Using grid R2B4"
elif [[ "$grid_choice" == "r2b5" ]]; then
  grid_file=${r2b5}
  echo "Using grid R2B5"
else
  echo "You have not specified grid choice"
fi

# creation of the output subdirectory if not existing
mkdir -p "${output_bc_static_dir}/static_file"
out="${output_bc_static_dir}/static_file"

# conversion of different variables into specified grid
${cdo} remaplaf,"${grid_file}" "${input_bc_static_crops_file}" "${out}/crops.${grid_choice}.nc"

${cdo} -P 2 \
    -setmisstoc,0.2 \
    -remapycon,${grid_file} \
    -shiftx,720,cyclic \
    -invertlatdata \
    -setgrid,r1440x720 \
    -setattribute,slash@long_name='forest harvest slash' \
    "${input_bc_static_slash_file}" \
    "${out}/slash_fractions_map_BEF_medium_gridcell_${grid_choice}.nc"


# merging all the files into one input file
${cdo} merge "${out}/*.nc" "${out}/bc_land_quincy_land_use_static_${grid_choice}.nc"
# removal of the intermediate files
rm "${out}/crops.${grid_choice}.nc" "${out}/slash_fractions_map_BEF_medium_gridcell_${grid_choice}.nc"
