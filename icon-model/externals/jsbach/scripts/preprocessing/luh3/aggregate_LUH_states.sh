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
#SBATCH --output=aggregate_LUH_states.o%j
#SBATCH --error=aggregate_LUH_states.o%j
#SBATCH --account=mj0143
#SBATCH --partition=interactive
#SBATCH --mem-per-cpu=3000
#SBATCH --exclusive
#SBATCH --ntasks=1

#----------------------------------------------------------------------
# JN 05.03.25, script to aggregate LUH type states
# - with rangelands treated depending on the fstnf value from the static data (multiple-static ...)
# (if 1, i.e. forests in LUH, there is a conversion to pasture; if 0: no conversion i.e. stays natural)

module load nco
module unload cdo
module load cdo/2.5.0-gcc-11.2.0
# set ncatted command to create no history
ncatted="ncatted -h"
# set cdo command to be silent, create no history and use double precision
cdo="cdo -s --no_history -b 64"

scriptName=aggregate_LUH_states.sh

# usually the original grid (for further use in the icon-land init file scripts)
targetGrid="global_0.25"
# -> note: lats in LUH data have another orientation as compared to the one in global_0.25
#          therefore a remapping here makes sense to ensure that all 0.25 grids used have the same orientation!
# targetGrid=/pool/data/ICON/grids/public/mpim/0049/icon_grid_0049_R02B04_G.nc   # mpim R2B4 grid
# targetGrid=/pool/data/ICON/grids/public/edzw/icon_grid_0030_R02B05_G.nc   # dwd R2B5 grid

start_year=1850
end_year=2024

tag="CMIP7_forcings"
version="3-1-1"
mainStatesFile=/work/bb1469/b383102/${tag}/RAW_DATA/land/multiple-states_input4MIPs_landState_CMIP_UofMD-landState-${version}_gn_0850-2024.nc
fixedInputFile=/work/bb1469/b383102/${tag}/RAW_DATA/land/multiple-static_input4MIPs_landState_CMIP_UofMD-landState-${version}_gn.nc

outfolder=./YOUROUTPUTFOLDER

#-------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------
echo "--> ${scriptName}"

# Some global attributes to be set
luh3_preprocessing="Aggregated LUH states (using ${scriptName}) and ${tag} LUH data version ${version}"
luh3_preprocessing="${luh3_preprocessing}; contact: jnabel@bgc-jena.mpg.de"
luh3_comment="Aggregated to natural, c3&c4crops, urban and pasture."
luh3_comment="${luh3_comment} Rangelands on forests are treated as pasture (LUH static data: fstnf = 1), else: natural vegetation"

# Note: the "pltns" state is zero until 2024
#       prop. needs to be handled > 2024
# https://gitlab.dkrz.de/quincy-community/quincy-iq-scripts/iq-scripts/-/issues/66#note_308611 : merge with natural forest

if [ ! -d ${outfolder} ]; then
    mkdir ${outfolder}
fi

cd ${outfolder}

# get fstnf file
fstnfFile=${outfolder}/multiple-fixed_input4MIPs_fstnf.nc
${cdo} -selvar,fstnf ${fixedInputFile} ${fstnfFile}

# get selected years from multiple states file
statesFile=./multiple-fixed_input4MIPs_landState_${start_year}-${end_year}.nc
${cdo} -selyear,${start_year}/${end_year} ${mainStatesFile} ${statesFile}

#------------------------------------------------------------------------------
#-- get pasture and natural vegetation share from rangelands depending on fstnf
${cdo} -mul -selvar,range ${statesFile} ${fstnfFile} tmp_rangeland_to_pasture.nc
${cdo} -mul -selvar,range ${statesFile} -mulc,-1 -subc,1 ${fstnfFile} tmp_rangeland_to_nat_veg.nc

#------------------------------------------------------------------------------
#-- aggregate states

#-- aggregate natural vegetation: primf, primn, secdf, secdn and range* as well as pltns
cdoExprString="nat=primf+primn+secdf+secdn+pltns"
${cdo} -expr,${cdoExprString} ${statesFile} tmp_sum_nat_veg.nc
${cdo} -add tmp_sum_nat_veg.nc tmp_rangeland_to_nat_veg.nc tmp_aggregated_nat_${start_year}-${end_year}.nc
${ncatted} -a long_name,nat,d,, tmp_aggregated_nat_${start_year}-${end_year}.nc

#--- aggregate pasture: pastr range*
${cdo} -add -selvar,pastr ${statesFile} tmp_rangeland_to_pasture.nc tmp_aggregated_pastr_${start_year}-${end_year}.nc
${ncatted} -a long_name,pastr,d,, tmp_aggregated_pastr_${start_year}-${end_year}.nc

#-- aggregate c3 crops: c3ann c3per and c3nfx
cdoExprString="c3crops=c3ann+c3per+c3nfx"
${cdo} -expr,${cdoExprString} ${statesFile} tmp_aggregated_c3crops_${start_year}-${end_year}.nc
${ncatted} -a long_name,c3crops,d,, tmp_aggregated_c3crops_${start_year}-${end_year}.nc

#-- aggregate c4 crops: c4ann and c4per
cdoExprString="c4crops=c4ann+c4per"
${cdo} -expr,${cdoExprString} ${statesFile} tmp_aggregated_c4crops_${start_year}-${end_year}.nc
${ncatted} -a long_name,c4crops,d,, tmp_aggregated_c4crops_${start_year}-${end_year}.nc

#--- select urban
${cdo} -selvar,urban ${statesFile} tmp_aggregated_urban_${start_year}-${end_year}.nc
${ncatted} -a long_name,urban,d,, tmp_aggregated_urban_${start_year}-${end_year}.nc

# join the files
${cdo} merge tmp_aggregated_*_${start_year}-${end_year}.nc tmp_org_LUH_states_${start_year}-${end_year}.nc
rm tmp*_aggregated_*_${start_year}-${end_year}.nc tmp_sum_nat_veg.nc tmp_rangeland_to_nat_veg.nc tmp_rangeland_to_pasture.nc

#------------------------------------------------------------------------------
#-- remap to target grid
${cdo} -remapcon,${targetGrid} tmp_org_LUH_states_${start_year}-${end_year}.nc tmp_finalGrid_LUH_states_${start_year}-${end_year}.nc

#------------------------------------------------------------------------------
#-- scale to one (divide by tmp_fact_no_zeros to avoid division by zero)
cdoExprString="fact=c3crops+c4crops+pastr+nat+urban"
${cdo} -expr,${cdoExprString} tmp_finalGrid_LUH_states_${start_year}-${end_year}.nc tmp_fact_${start_year}-${end_year}.nc
${cdo} -div tmp_finalGrid_LUH_states_${start_year}-${end_year}.nc tmp_fact_${start_year}-${end_year}.nc LUH_states_${start_year}-${end_year}.nc

#------------------------------------------------------------------------------
# care for attributes
# set global attributes
${ncatted} -O -a ,global,d,, LUH_states_${start_year}-${end_year}.nc
${ncatted} -O -a luh3_comment,global,o,c,"${luh3_comment}" LUH_states_${start_year}-${end_year}.nc
${ncatted} -O -a luh3_preprocessing,global,o,c,"${luh3_preprocessing}" LUH_states_${start_year}-${end_year}.nc

#------------------------------------------------------------------------------
#-- split file to get annual files
${cdo} -splityear LUH_states_${start_year}-${end_year}.nc LUH_states_

# reduce time dimension
files=$( ls LUH_states_????.nc )

for file in ${files}; do
  ${cdo} --reduce_dim copy ${file} tmp.nc
  mv tmp.nc ${file}
done

rm LUH_states_${start_year}-${end_year}.nc ${statesFile}
rm tmp_org_LUH_states_${start_year}-${end_year}.nc tmp_finalGrid_LUH_states_${start_year}-${end_year}.nc tmp_fact_${start_year}-${end_year}.nc


