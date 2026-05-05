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
#SBATCH --output=calculate_mean_LUH_harvest.o%j
#SBATCH --error=calculate_mean_LUH_harvest.o%j
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

scriptName="calculate_mean_LUH_harvest.sh"
#----------------------------------------------------------------------
# calculate aggregated and mapped mean LUH area harvest for use with IQ

start_year=1850 # start year from which on to use data from the LUH transitions file
end_year=1869   # end year from which on to use data from the LUH transitions file

# path and name of the target output folder
outfolder=YOUROUTPUTFOLDER

# path to the file with original LUH harvest data
tag="CMIP7_forcings"
version="3-1-1"
transitionsFileWithPath=/work/bb1469/b383102/${tag}/RAW_DATA/land/multiple-transitions_input4MIPs_landState_CMIP_UofMD-landState-${version}_gn_0850-2023.nc

# grid ...
targetGrid=/pool/data/ICON/grids/public/mpim/0049/icon_grid_0049_R02B04_G.nc
#targetGrid=/pool/data/ICON/grids/public/edzw/icon_grid_0030_R02B05_G.nc   # dwd R2B5 grid

gridID=0049

#-------------------------------------------------------------------------------------------------------------
# Note: previously we included all harvested areas, while we now only include harvest from (LUH) forest
#-- previously
# varList=primf_harv,primn_harv,secmf_harv,secyf_harv,secnf_harv,pltns_harv
# cdoExprString="harvest_fract=primf_harv+primn_harv+secmf_harv+secyf_harv+secnf_harv+pltns_harv"
#-- now
varList=primf_harv,secmf_harv,secyf_harv,pltns_harv
cdoExprString="harvest_fract=primf_harv+secmf_harv+secyf_harv+pltns_harv"
# Note: for now we include pltns_harv (but this is zero in the historical forcing anyway)
#       - have to further be considered when processing scenario data

#----------------------------------------------------------------------
echo "--> ${scriptName}"

# Some global attributes to be set
history="Derived from LUH data using ${scriptName} and ${tag} LUH data version ${version}"
comment="${start_year} to ${end_year} mean of summed harvest area"
contact="jnabel@bgc-jena.mpg.de"

#----------------------------------------------------------------------
if [ ! -d ${outfolder} ]; then
    mkdir ${outfolder}
fi

cd ${outfolder}

# aggregate for selected years
${cdo} -selyear,${start_year}/${end_year} -selvar,${varList} ${transitionsFileWithPath} tmp_sel_${start_year}-${end_year}.nc
${cdo} -expr,${cdoExprString} tmp_sel_${start_year}-${end_year}.nc tmp_${start_year}-${end_year}.nc

# restrict to a maximum of 0.2 - i.e. do not harvest more than 0.2 of a grid-cell in one year
${cdo} -gtc,0.2 tmp_${start_year}-${end_year}.nc tmp_mask_harvest_${start_year}-${end_year}.nc
${cdo} -mulc,0.2 tmp_mask_harvest_${start_year}-${end_year}.nc tmp_max_harvest_${start_year}-${end_year}.nc
${cdo} -ifthenelse tmp_mask_harvest_${start_year}-${end_year}.nc \
    tmp_max_harvest_${start_year}-${end_year}.nc tmp_${start_year}-${end_year}.nc tmp_res_${start_year}-${end_year}.nc

# remap and set miss to zero
${cdo} -setmisstoc,0 -remapycon,${targetGrid} tmp_res_${start_year}-${end_year}.nc tmp_mapped_${start_year}-${end_year}.nc

# remove wrong attributes
${ncatted} -a long_name,harvest_fract,d,, -a standard_name,harvest_fract,d,, tmp_mapped_${start_year}-${end_year}.nc
${ncatted} -a _FillValue,harvest_fract,d,, -a missing_value,harvest_fract,d,, tmp_mapped_${start_year}-${end_year}.nc

# set some attributes
${ncatted} -O -a comment,global,o,c,"${comment}" tmp_mapped_${start_year}-${end_year}.nc
${ncatted} -O -a contact,global,o,c,"${contact}" tmp_mapped_${start_year}-${end_year}.nc

# calculate time mean
${cdo} --reduce_dim -timmean tmp_mapped_${start_year}-${end_year}.nc LUH_harvest_mean_${start_year}-${end_year}.nc

# set clean history
${ncatted} -O -a history,global,o,c,"${history}" LUH_harvest_mean_${start_year}-${end_year}.nc

# cleanup folder
rm tmp*_${start_year}-${end_year}.nc
