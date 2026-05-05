#!/usr/bin/bash

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
# This script processes nitrogen deposition data from CMIP7 input datasets, combining wet and dry components of
# NHx and NOy fluxes into a single NetCDF file on the ICON grid.
# In addition, phosphorus depositions were added, based on
# Brahney et al. (2015, Global Biogeochem. Cycles,doi:10.1002/2015GB005137).

# load software module
module load nco

#----------------------------------------------------------------------


# Grid selection
grid_res_id=r2b5 # r2b4 or r2b5

grid_def_r2b4=/pool/data/ICON/grids/public/mpim/0049/icon_grid_0049_R02B04_G.nc
grid_def_r2b5=/pool/data/ICON/grids/public/edzw/icon_grid_0030_R02B05_G.nc

if [[ "$grid_res_id" == "r2b4" ]]; then
  grid_def_file=${grid_def_r2b4}
  institute_tag=mpim
  grid_id=0049
  echo "grid R2B4"
elif [[ "$grid_res_id" == "r2b5" ]]; then
  grid_def_file=${grid_def_r2b5}
  institute_tag=edzw
  grid_id=0030
  echo "grid R2B5"
else
  echo "ERROR! you have not specified supported grid id (r2b4/r2b5)"
  exit 1
fi

# Paths and filenames
input_path_n="/work/mj0143/icon_land_data_pool/preprocessing/np_depositions/cmip7/"
input_path_p="/work/mj0143/icon_land_data_pool/preprocessing/np_depositions/p_dep_files/"
output_path="./"
mkdir -p ${output_path}/temp
temp_path=${output_path}"/temp"

# Define years to process depositions
start_year=1850
n_dep_end_year=2022
p_dep_end_year=2019

# Define units
dep_unit="kg/m2/s"

# Enable verbose mode for debugging
# set -vx
# ---- Process NHx deposition ----
# Combine wet and dry NHx deposition
cdo --no_history -b F64 -add ${input_path_n}/wetnhx_input4MIPs_surfaceFluxes_CMIP_FZJ-CMIP-nitrogen-1-2_gn_185001-202212.nc \
                ${input_path_n}/drynhx_input4MIPs_surfaceFluxes_CMIP_FZJ-CMIP-nitrogen-1-2_gn_185001-202212.nc \
                ${temp_path}/foo.nc
# Fix missing values
cdo --no_history -b F64 -setmissval,-9.e33 -setmisstoc,0 -setmissval,Infinity -setmisstoc,0 -setmissval,nan ${temp_path}/foo.nc ${temp_path}/foo1.nc
# Select only the NHx variable, rename it, and set attributes
ncks -O -h -v wetnhx ${temp_path}/foo1.nc ${temp_path}/foo2.nc
ncrename -O -h -v wetnhx,NHx_deposition ${temp_path}/foo2.nc
ncatted -h -a long_name,NHx_deposition,o,c,"NHx_deposition (dry+wet)" \
        -h -a standard_name,,d,c,"" \
        -h -a original_name,,d,c,"" \
        -h -a bounds,,d,c,"" \
        -h -a units,NHx_deposition,o,c,"${dep_unit}" ${temp_path}/foo2.nc
# Remove boundary variables
ncks -O -h -x -v lat_bnds,lon_bnds,time_bnds ${temp_path}/foo2.nc ${temp_path}/foo3.nc
# Remap to target model grid
cdo --no_history -remapcon,${grid_def_file} ${temp_path}/foo3.nc ${temp_path}/nhx_deposition.nc
# Clean temporary files
rm ${temp_path}/foo*

# ---- Process NOy deposition ----
# Combine wet and dry NOy deposition
cdo --no_history -b F64 -add  ${input_path_n}/wetnoy_input4MIPs_surfaceFluxes_CMIP_FZJ-CMIP-nitrogen-1-2_gn_185001-202212.nc \
                 ${input_path_n}/drynoy_input4MIPs_surfaceFluxes_CMIP_FZJ-CMIP-nitrogen-1-2_gn_185001-202212.nc \
                ${temp_path}/foo.nc
# Fix missing values
cdo --no_history -b F64 -setmissval,-9.e33 -setmisstoc,0 -setmissval,Infinity -setmisstoc,0 -setmissval,nan ${temp_path}/foo.nc ${temp_path}/foo1.nc
# Select NOy variable, rename it, and set attributes
ncks -O -h -v wetnoy ${temp_path}/foo1.nc ${temp_path}/foo2.nc
ncrename -O -h -v wetnoy,NOy_deposition ${temp_path}/foo2.nc
ncatted -h -a long_name,NOy_deposition,o,c,"NOy_deposition (dry+wet)" \
        -h -a units,NOy_deposition,o,c,"${dep_unit}" \
        -h -a standard_name,,d,c,"" \
        -h -a original_name,,d,c,"" \
        -h -a bounds,,d,c,"" ${temp_path}/foo2.nc
# Remove boundary variables
ncks -O -h -x -v lat_bnds,lon_bnds,time_bnds ${temp_path}/foo2.nc ${temp_path}/foo3.nc
# Remap to target grid
cdo --no_history -remapcon,${grid_def_file} ${temp_path}/foo3.nc ${temp_path}/noy_deposition.nc

# Clean temporary files
rm ${temp_path}/foo*

# ---- Merge NHx and NOy into a single N-deposition file ----
cdo --no_history -merge ${temp_path}/nhx_deposition.nc ${temp_path}/noy_deposition.nc ${temp_path}/ndep.nc
# Split merged file into yearly files
cdo --no_history -splityear ${temp_path}/ndep.nc ${temp_path}/ndep_

# Loop through each year from start_year to end_year
year=${start_year}
while [ ${year} -le ${n_dep_end_year} ]; do
  # Define input and output paths
  ndep_in="${temp_path}/ndep_${year}.nc"
  npdep_out_pre="${temp_path}/bc_land_npdep_quincy_${year}_pre.nc"
  npdep_out="${output_path}/bc_land_npdep_quincy_${year}.nc"
  if [ ${year} -gt ${p_dep_end_year} ]; then
    pdep_in="${input_path_p}/bc_land_npdep_quincy_2019.nc"
  else
    pdep_in="${input_path_p}/bc_land_npdep_quincy_${year}.nc"
  fi
  # Select the 'pdep' variable from the files containing phosphorus depositions
  cdo --no_history -selvar,pdep "${pdep_in}" ${temp_path}/pdep_${year}_temp.nc
  # Remap the selected Pdep file to the target grid defined by grid_def_file
  cdo --no_history -remapcon,${grid_def_file} ${temp_path}/pdep_${year}_temp.nc ${temp_path}/pdep_${year}.nc
  # Merge yearly Ndep and Pdep
  cdo --no_history -merge "${ndep_in}" "${temp_path}/pdep_${year}.nc" "${npdep_out_pre}"
  ncatted -O -h \
  -a references,global,a,c,"; Brahney, J., N. Mahowald, D. S. Ward, A. P. Ballantyne, and J. C. Neff (2015), Is atmospheric phosphorus pollution altering global alpine Lake stoichiometry?, Global Biogeochem. Cycles, 29, 1369–1383, doi:10.1002/2015GB005137." \
  -a variable_id,global,d,, \
  "${npdep_out_pre}" "${npdep_out}"
  # Clean up temporary files
  rm -f "${temp_path}/pdep_${year}_temp.nc" "${temp_path}/pdep_${year}.nc"
  year=$((year + 1))
done
echo "Processing complete. Yearly NP-deposition files are located in ${output_path}."

# Remove temporary directory
rm -r ${temp_path}

exit 0

