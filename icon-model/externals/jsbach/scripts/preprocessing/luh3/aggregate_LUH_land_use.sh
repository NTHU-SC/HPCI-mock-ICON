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
#SBATCH --output=aggregate_LUH_land_use.o%j
#SBATCH --error=aggregate_LUH_land_use.o%j
#SBATCH --account=mj0143
#SBATCH --partition=interactive
#SBATCH --mem-per-cpu=3000
#SBATCH --exclusive    #-- exclusive required when working with the transitions file
#SBATCH --ntasks=1

module load nco
module unload cdo
module load cdo/2.5.0-gcc-11.2.0
# set ncatted command to create no history
ncatted="ncatted -h"
# set cdo command to be silent, create no history and use double precision
cdo="cdo -s --no_history -b 64"

scriptName="aggregate_LUH_land_use.sh"
#----------------------------------------------------------------------
# aggregate and map LUH area harvest variables and LUH fertilisation variables for use with IQ
echo "--> ${scriptName}"

# path and name of the target output folder
outfolder=./YOUROUTPUTFOLDER

# path to the LUH files
tag="CMIP7_forcings"
version="3-1-1"
transitionsFileWithPath=/work/bb1469/b383102/${tag}/RAW_DATA/land/multiple-transitions_input4MIPs_landState_CMIP_UofMD-landState-${version}_gn_0850-2023.nc
ifile_management=/work/bb1469/b383102/${tag}/RAW_DATA/land/multiple-management_input4MIPs_landState_CMIP_UofMD-landState-${version}_gn_0850-2024.nc
ifile_states=/work/bb1469/b383102/${tag}/RAW_DATA/land/multiple-states_input4MIPs_landState_CMIP_UofMD-landState-${version}_gn_0850-2024.nc

# grid ...
targetGrid=/pool/data/ICON/grids/public/mpim/0049/icon_grid_0049_R02B04_G.nc
#targetGrid=/pool/data/ICON/grids/public/edzw/icon_grid_0030_R02B05_G.nc   # dwd R2B5 grid

start_year=1850 # start year from which on to use data from the LUH file
end_year=2023   # end year until which on to use data from the LUH file

#-------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------
# Some global attributes to be set
history="Aggregate LUH land-use data (${scriptName}) and ${tag} LUH data version ${version}"
comment="Harvest: sum all harvest area; Fertilisation: sum c3 and c4 fertilisation, respectively"
contact="jnabel@bgc-jena.mpg.de"

#-------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------
varListFertiliser="fertl_c3ann,fertl_c3per,fertl_c3nfx,fertl_c4ann,fertl_c4per"
varListCrops="c3ann,c3per,c3nfx,c4ann,c4per"
cdoExprFertilC3='fertl_c3crops=(c3ann*fertl_c3ann+c3per*fertl_c3per+c3nfx*fertl_c3nfx)/(c3ann+c3per+c3nfx)'
cdoExprFertilC4='fertl_c4crops=(c4ann*fertl_c4ann+c4per*fertl_c4per)/(c4ann+c4per)'

# Note: previously we included all harvested areas, while we now only include harvest from (LUH) forest
#-- previously
# varListHarvest=primf_harv,primn_harv,secmf_harv,secyf_harv,secnf_harv,pltns_harv
# cdoExprHarvest="harvest_fract=primf_harv+primn_harv+secmf_harv+secyf_harv+secnf_harv+pltns_harv"
#-- now
varListHarvest=primf_harv,secmf_harv,secyf_harv,pltns_harv
cdoExprHarvest="harvest_fract=primf_harv+secmf_harv+secyf_harv+pltns_harv"
# Note: for now we include pltns_harv (but this is zero in the historical forcing anyway)
#       - have to further be considered when processing scenario data

if [ ! -d ${outfolder} ]; then
    mkdir ${outfolder}
fi

cd ${outfolder}

# aggregate harvest for selected years
${cdo} -selyear,${start_year}/${end_year} -selvar,${varListHarvest} ${transitionsFileWithPath} tmp_trans_${start_year}-${end_year}.nc
${cdo} -expr,${cdoExprHarvest} tmp_trans_${start_year}-${end_year}.nc tmp_harvest_${start_year}-${end_year}.nc

# restrict to a maximum of 0.2 - i.e. do not harvest more than 0.2 of a grid-cell in one year
${cdo} -gtc,0.2 tmp_harvest_${start_year}-${end_year}.nc tmp_mask_harvest_${start_year}-${end_year}.nc
${cdo} -mulc,0.2 tmp_mask_harvest_${start_year}-${end_year}.nc tmp_max_harvest_${start_year}-${end_year}.nc
${cdo} -ifthenelse tmp_mask_harvest_${start_year}-${end_year}.nc \
  tmp_max_harvest_${start_year}-${end_year}.nc tmp_harvest_${start_year}-${end_year}.nc tmp_res_harvest_${start_year}-${end_year}.nc

# get required crop related variables from the management and state files
${cdo} -selyear,${start_year}/${end_year} -selvar,${varListFertiliser} ${ifile_management} tmp_luh_fertiliser_${start_year}-${end_year}.nc
${cdo} -selyear,${start_year}/${end_year} -selvar,${varListCrops} ${ifile_states} tmp_luh_frac_${start_year}-${end_year}.nc
${cdo} -merge tmp_luh_fertiliser_${start_year}-${end_year}.nc tmp_luh_frac_${start_year}-${end_year}.nc tmp_merged_${start_year}-${end_year}.nc
# calc C3 crop fertliser
${cdo} -expr,${cdoExprFertilC3} tmp_merged_${start_year}-${end_year}.nc fertl_C3crops_${start_year}-${end_year}.nc
# calc C4 crop fertliser
${cdo} -expr,${cdoExprFertilC4} tmp_merged_${start_year}-${end_year}.nc fertl_C4crops_${start_year}-${end_year}.nc

# merge into one file for future work
${cdo} merge fertl_C3crops_${start_year}-${end_year}.nc fertl_C4crops_${start_year}-${end_year}.nc \
  tmp_res_harvest_${start_year}-${end_year}.nc tmp_land_use_${start_year}-${end_year}.nc

# remap and set miss to zero
${cdo} -setmisstoc,0 -remapycon,${targetGrid} tmp_land_use_${start_year}-${end_year}.nc tmp_land_use_mapped_${start_year}-${end_year}.nc

# remove attributes which are not correct (TODO: consider replacing instead of deleting)
${ncatted} -a long_name,harvest_fract,d,,  -a standard_name,harvest_fract,d,, tmp_land_use_mapped_${start_year}-${end_year}.nc
${ncatted} -a _FillValue,harvest_fract,d,, -a missing_value,harvest_fract,d,, tmp_land_use_mapped_${start_year}-${end_year}.nc
${ncatted} -a long_name,fertl_c4crops,d,, tmp_land_use_mapped_${start_year}-${end_year}.nc
${ncatted} -a _FillValue,fertl_c4crops,d,, -a missing_value,fertl_c4crops,d,, tmp_land_use_mapped_${start_year}-${end_year}.nc
${ncatted} -a long_name,fertl_c3crops,d,, tmp_land_use_mapped_${start_year}-${end_year}.nc
${ncatted} -a _FillValue,fertl_c3crops,d,, -a missing_value,fertl_c3crops,d,, tmp_land_use_mapped_${start_year}-${end_year}.nc

# set further global attributes
${ncatted} -O -a comment,global,o,c,"${comment}" tmp_land_use_mapped_${start_year}-${end_year}.nc
${ncatted} -O -a contact,global,o,c,"${contact}" tmp_land_use_mapped_${start_year}-${end_year}.nc

# split file to get annual files
${cdo} -splityear tmp_land_use_mapped_${start_year}-${end_year}.nc bc_land_use_quincy_

# reduce time dimension
files=$( ls bc_land_use_quincy_????.nc )

for file in ${files}; do
  ${cdo} --reduce_dim copy ${file} tmp.nc

  # set clean history
  ${ncatted} -O -a history,global,o,c,"${history}" tmp.nc

  mv tmp.nc ${file}
done

# cleanup folder
rm tmp*_${start_year}-${end_year}.nc
rm fertl_C3crops_${start_year}-${end_year}.nc fertl_C4crops_${start_year}-${end_year}.nc
