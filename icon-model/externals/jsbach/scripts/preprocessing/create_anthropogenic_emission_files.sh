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

#-----------------------------------------------------------------------------
# Scripts to combine aviation emission file into other anthropogenic emission files.
###############################################################################
### Batch Queuing System is SLURM
#SBATCH --job-name=create_anthropogenic_emission_files
#SBATCH --output=create_anthropogenic_emission_files.o%j
#SBATCH --error=create_anthropogenic_emission_files.o%j
#SBATCH --partition=compute
#SBATCH --mem-per-cpu=5120 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=07:30:00
#SBATCH --account mj0143
#============================================================================
#============================================================================
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
module load cdo nco

# 'User' definitions

# directory where the anthropogenic emisson data from the CMIP7 dataset is located
DataPath='/work/mj0143/icon_land_data_pool/preprocessing/CO2-em-anthro/'

# directory where the scratch data should be located
tmpPath='YOURWORKPLACE'

# directory where the new data should be located
workingPath='YOURWORKPLACE'

ATMO_GRID_ID=0049
ATMO_GRID_TYPE=R02B04
INSTITUTE_TAG=mpim

CO2_em_anthro_reference="Hoesly et al. (2025) CEDS v_2025_04_18 0.5 degree gridded data"
CO2_em_anthro_history="Combine anthropogenic emission data with aircraft emissions by summing over 25 vertical levels, then regrid the result to the ICON ${ATMO_GRID_TYPE} grid."
CO2_em_anthro_institution="Max Planck Institute for Biogeochemistry"
CO2_em_anthro_contact="YOUR CONTACT"
CO2_em_anthro_title="1750-2023 anthropogenic emission data"

# The raw data is provided in 50-year blocks as monthly data
# Data starting from the year 2000 covers the period 2000–2023
first_year=1750
last_year=2023

currentYear=${first_year}

while [[ ${currentYear} -le ${last_year} ]]; do

    if (( currentYear == 1750 )); then
        Data_start_year=1750
        Data_end_year=1799
    elif (( currentYear == 1800 )); then
        Data_start_year=1800
        Data_end_year=1849
    elif (( currentYear == 1850 )); then
        Data_start_year=1850
        Data_end_year=1899
    elif (( currentYear == 1900 )); then
        Data_start_year=1900
        Data_end_year=1949
    elif (( currentYear == 1950 )); then
        Data_start_year=1950
        Data_end_year=1999
    elif (( currentYear == 2000 )); then
        Data_start_year=2000
        Data_end_year=2023
    fi

    echo "currentYear is ${currentYear}"
    echo "Processing ${currentYear} with file ${Data_start_year}01-${Data_end_year}12"

    #Sum up the values across 25 layers
    cdo -s vertsum, -selyear,${currentYear} \
        ${DataPath}/CO2-em-AIR-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_${Data_start_year}01-${Data_end_year}12.nc \
        ${tmpPath}/aircraft_emission_sum_${currentYear}.nc
    cdo selyear,${currentYear} ${DataPath}/CO2-em-anthro_input4MIPs_emissions_CMIP_CEDS-CMIP-2025-04-18_gn_${Data_start_year}01-${Data_end_year}12.nc ${tmpPath}/tmp_anthro.nc

    #define new dimension "sector" to merge into the anthropogenic files
    ncap2 -O -s 'defdim("sector",1); CO2_em_anthro[$time,$sector,$lat,$lon]=float(CO2_em_AIR_anthro(:,:,:));' \
    ${tmpPath}/aircraft_emission_sum_${currentYear}.nc ${tmpPath}/air_tmp_${currentYear}.nc
    ncks -O -x -v CO2_em_AIR_anthro ${tmpPath}/air_tmp_${currentYear}.nc ${tmpPath}/air_tmp2_${currentYear}.nc
    ncap2 -O -s 'defdim("sector",9);' \
        ${tmpPath}/air_tmp2_${currentYear}.nc ${tmpPath}/air_tmp3_${currentYear}.nc

    #merge aircraft emission file into anthropogenic emission file
    cdo merge ${tmpPath}/tmp_anthro.nc ${tmpPath}/air_tmp3_${currentYear}.nc ${tmpPath}/combined_anthropogenic_with_aircraft_${currentYear}.nc
    rm -f ${tmpPath}/aircraft_emission_sum_${currentYear}.nc
    rm -f ${tmpPath}/air_tmp_${currentYear}.nc
    rm -f ${tmpPath}/air_tmp2_${currentYear}.nc
    rm -f ${tmpPath}/air_tmp3_${currentYear}.nc
    rm -f ${tmpPath}/tmp_anthro.nc

    (( currentYear += 1 ))

done

cdo mergetime ${tmpPath}/combined_anthropogenic_with_aircraft* ${workingPath}/combined_anthropogenic_with_aircraft_${first_year}_${last_year}.nc

#add "8: Aviation" into the sector ids
ncap2 -O -s 'defdim("sector",9); sector[$sector]={0,1,2,3,4,5,6,7,8}' -o ${tmpPath}/sector.nc
ncks -A ${tmpPath}/sector.nc ${workingPath}/combined_anthropogenic_with_aircraft_${first_year}_${last_year}.nc
ncatted -O -a ids,sector,o,c,"0: Agriculture; 1: Energy; 2: Industrial; 3: Transportation; 4: Residential, Commercial, Other; 5: Solvents production and application; 6: Waste; 7: International Shipping; 8: Aviation" ${workingPath}/combined_anthropogenic_with_aircraft_${first_year}_${last_year}.nc
rm -f ${tmpPath}/sector.nc
rm -f ${tmpPath}/combined_anthropogenic_with_aircraft*

#add regrid into the ICON grids
cdo -remapcon,/pool/data/ICON/grids/public/${INSTITUTE_TAG}/${ATMO_GRID_ID}/icon_grid_${ATMO_GRID_ID}_${ATMO_GRID_TYPE}_G.nc ${workingPath}/combined_anthropogenic_with_aircraft_${first_year}_${last_year}.nc ${workingPath}/combined_anthropogenic_with_aircraft_icongrids_${first_year}_${last_year}.nc

#sum up the values across sectors
ncwa -a sector -y sum ${workingPath}/combined_anthropogenic_with_aircraft_icongrids_${first_year}_${last_year}.nc ${workingPath}/bc_anthropogenic_co2_${first_year}-${last_year}_${ATMO_GRID_TYPE}.nc

# change global attributes on prepared anthropogenic emission file
ncatted -h -O -a ,global,d,, ${workingPath}/bc_anthropogenic_co2_${first_year}-${last_year}_${ATMO_GRID_TYPE}.nc
ncatted -h -O -a references,global,o,c,"${CO2_em_anthro_reference}" "${workingPath}/bc_anthropogenic_co2_${first_year}-${last_year}_${ATMO_GRID_TYPE}.nc"
ncatted -h -O -a history,global,o,c,"${CO2_em_anthro_history}" "${workingPath}/bc_anthropogenic_co2_${first_year}-${last_year}_${ATMO_GRID_TYPE}.nc"
ncatted -h -O -a institution,global,o,c,"${CO2_em_anthro_institution}" "${workingPath}/bc_anthropogenic_co2_${first_year}-${last_year}_${ATMO_GRID_TYPE}.nc"
ncatted -h -O -a contact,global,o,c,"${CO2_em_anthro_contact}" "${workingPath}/bc_anthropogenic_co2_${first_year}-${last_year}_${ATMO_GRID_TYPE}.nc"
ncatted -h -O -a title,global,o,c,"${CO2_em_anthro_title}" "${workingPath}/bc_anthropogenic_co2_${first_year}-${last_year}_${ATMO_GRID_TYPE}.nc"

rm -f ${workingPath}/combined_anthropogenic_with_aircraft_${first_year}_${last_year}.nc
rm -f ${workingPath}/combined_anthropogenic_with_aircraft_icongrids_${first_year}_${last_year}.nc