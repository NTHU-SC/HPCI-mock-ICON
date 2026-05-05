#!/bin/bash
#
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
#
# Script to generate the HD parameter file for simulations with ICON-Land and
# internal HD.  It is based on the instructions given in
# https://gitlab.dkrz.de/jsbach/jsbach/-/wikis/Documentation/ICON-Land-input-data/Icon-Land-initial-file-generation/low-res-hdpara-file-generation
#
#------------------------------------------------------------------------------
set -e
prog=$(basename $0)

tmpdir=${work_dir}/hd_tmp
mkdir -p ${tmpdir}; rm -f ${tmpdir}/*.nc
cd ${tmpdir}

cdo="$cdo"

echo "=============================================================="
echo "===  Generation of the HD parameter file"
echo "=============================================================="

# Find out, if this is a high resolution setup (R02B07 or higher)
hd_high_res=false
if [[   $(echo ${atmRes} | cut -c3 ) -gt 2  \
     || $((10#$(echo ${atmRes} | cut -c5-6 ))) -gt 6 ]]; then
  hd_high_res=true
  echo "$0: HD parameter file generation for high resolution setup: "${atmRes}
  echo " "
fi

# 1. Setup the repository with icon hdpara generation tools
[[ ! -d icon_hdpara_generation_tools ]] || rm -rf icon_hdpara_generation_tools
git clone https://github.com/ThomasRiddick/DynamicHD.git icon_hdpara_generation_tools
cd icon_hdpara_generation_tools
git checkout icon_hd_tools_version_1.5.1
cp -r ${scripts_dir}/dynamic_hd_tools/*  Dynamic_HD_bash_scripts/parameter_generation_scripts/

# 2. Generate the necessary mamba environment (can be slow; needs to be done only once)
if [[ $(mamba info -e | grep "dyhdenv_mamba") == "" ]]; then
  ./Dynamic_HD_bash_scripts/regenerate_mamba_environment.sh
fi
if [[ $(julia -E 'import Pkg
                  Pkg.status()' | grep -c "NetCDF") -eq 0 ]]; then
  julia -E 'import Pkg
            Pkg.add("NetCDF")'
fi
if [[ $(julia -E 'import Pkg
                  Pkg.status()' | grep -c "ArgParse") -eq 0 ]]; then
  julia -E 'import Pkg
            Pkg.add("ArgParse")'
fi

# 3. Activate the mamba environment
source activate dyhdenv_mamba

# 4. Compile code and generate example run script configuration files
make

# 5. Further preparations

# Use the fractional land sea mask from the bc_land_frac file.
# Note: Due to the min- and max_fract definition, the fractional mask in the bc_land_frac
# file ('notsea') is generally not identical with the fractional mask in the fractional_mask
# file. The expected variable name for the HD parameter generation tools is 'cell_sea_land_mask'.
if [[ $(cdo showvar ${hd_fractional_mask} | grep notsea ) != "" ]]; then
  $cdo setvar,cell_sea_land_mask -selvar,notsea ${hd_fractional_mask} ${tmpdir}/fractional_mask.nc
  fractional_lsmask_filepath=${tmpdir}/fractional_mask.nc
else
  if [[ $(cdo showvar ${hd_fractional_mask} | grep cell_sea_land_mask ) != "" ]]; then
    fractional_lsmask_filepath=${hd_fractional_mask}
  else
    echo "$0: ERROR: ICON-Land bc file ${hd_fractional_mask} "\
              "does not exist or does not contain variable 'notsea' nor 'cell_sea_land_mask'."
    echo ""
    exit 1
  fi
fi

# Change format of AtmRes: e.g. R02B04 -> r2b4
icon_atmo_grid_res=$(echo ${atmRes} | tr -s 'RB' 'rb' | tr -d 0)

# 6. Create and edit a run script configuration file
cd run
sed     "s/icon_atmo_grid_id=.*/icon_atmo_grid_id=${atmGridID}/"          examples/r2b3_example.cfg \
  | sed "s/icon_ocean_grid_id=.*/icon_ocean_grid_id=${oceGridID}/"                                  \
  | sed "s/icon_atmo_grid_res=.*/icon_atmo_grid_res=${icon_atmo_grid_res}/"                         \
  | sed "s:icon_grid_filepath=.*:icon_grid_filepath=${icon_grid}:"                                  \
  | sed "s:fractional_lsmask_filepath=.*:fractional_lsmask_filepath=${fractional_lsmask_filepath}:" \
    > ${icon_atmo_grid_res}.cfg
if [[ ${hd_high_res} == true ]]; then
  sed   "s:orography_filepath=.*:orography_filepath=${hd_orography_file}:" ${icon_atmo_grid_res}.cfg \
   > ${icon_atmo_grid_res}.tmp
  mv ${icon_atmo_grid_res}.tmp  ${icon_atmo_grid_res}.cfg
fi

~/.conda/envs/dyhdenv_mamba/bin/python ../utils/run_utilities/mkproject.py ${icon_atmo_grid_res}.cfg

# 7. Run the script (use an appropriate interactive node)
./${icon_atmo_grid_res}.run

echo "***"
echo "*** Please ignore HDF5 error messages !!"
echo "***"

# 8. Save the new hdpara file and clean up
hdpara_file=$(ls ../projects/${icon_atmo_grid_res}/output/hdpara_${icon_atmo_grid_res}_${atmGridID}_${oceGridID}_*v?.nc)
hdpara_file_new=$(echo ${hdpara_file##*/} | tr -s _ )  # remove path and extra '_' in case of empty oceGridID
cp ${hdpara_file} ${hd_file_dir}/${hdpara_file_new}

echo "----------------------------------------------------------------"
echo "$0: Generated HD parameter file:"
echo "     ${hd_file_dir}/${hdpara_file_new}"

if [[ ${keep_hd_output_dir} ]]; then
  mkdir -p ${hd_file_dir}/hd_output_dir
  mv ../projects/${icon_atmo_grid_res}/output/* ${hd_file_dir}/hd_output_dir
  echo "Output from hdpara generation, temporarily kept for debugging etc." > ${hd_file_dir}/hd_output_dir/README
  echo "$0: Additional HD output:"
  echo "     ${hd_file_dir}/hd_output_dir"
fi
echo "----------------------------------------------------------------"

cd ..
rm -rf ${tmpdir}

