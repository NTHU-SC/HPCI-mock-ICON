#!/bin/bash
#
# adding pft info to netcdf file attributes
#----------------------------------------------------------------

# getting arguments
model=$1
output_file=$2

module load nco

# set ncatted command to create no history
ncatted="ncatted -h"

if [ ${model} == "jsbach" ]; then
  ${ncatted} -O -a long_name,TE,o,c,"Tropical evergreen trees" "${output_file}"
  ${ncatted} -O -a short_name,TE,o,c,"TE" "${output_file}"
  ${ncatted} -O -a long_name,TD,o,c,"Tropical deciduous trees" "${output_file}"
  ${ncatted} -O -a short_name,TD,o,c,"TD" "${output_file}"
  ${ncatted} -O -a long_name,ETE,o,c,"Extra-tropical evergreen trees" "${output_file}"
  ${ncatted} -O -a short_name,ETE,o,c,"ETE" "${output_file}"
  ${ncatted} -O -a long_name,ETD,o,c,"Extra-tropical deciduous trees" "${output_file}"
  ${ncatted} -O -a short_name,ETD,o,c,"ETD" "${output_file}"
  ${ncatted} -O -a long_name,RShrubs,o,c,"Raingreen shrubs" "${output_file}"
  ${ncatted} -O -a short_name,RShrubs,o,c,"RShrubs" "${output_file}"
  ${ncatted} -O -a long_name,DShrubs,o,c,"Deciduous shrubs" "${output_file}"
  ${ncatted} -O -a short_name,DShrubs,o,c,"DShrubs" "${output_file}"
  ${ncatted} -O -a long_name,H,o,c,"C3 grass" "${output_file}"
  ${ncatted} -O -a short_name,H,o,c,"TeH" "${output_file}"
  ${ncatted} -O -a long_name,HC4,o,c,"C4 grass" "${output_file}"
  ${ncatted} -O -a short_name,HC4,o,c,"TrH" "${output_file}"
  ${ncatted} -O -a long_name,bare,o,c,"bare" "${output_file}"
  ${ncatted} -O -a short_name,bare,o,c,"bare" "${output_file}"
elif [ ${model} == "quincy" ]; then
  ${ncatted} -O -a long_name,TrBE,o,c,"Moist broadleaved evergreen" "${output_file}"
  ${ncatted} -O -a short_name,TrBE,o,c,"BEM" "${output_file}"
  ${ncatted} -O -a long_name,TeBE,o,c,"Dry broadleaved evergreen" "${output_file}"
  ${ncatted} -O -a short_name,TeBE,o,c,"BED" "${output_file}"
  ${ncatted} -O -a long_name,TrBR,o,c,"Rain green broadleaved deciduous" "${output_file}"
  ${ncatted} -O -a short_name,TrBR,o,c,"BDR" "${output_file}"
  ${ncatted} -O -a long_name,TeBS,o,c,"Summer green broadleaved deciduous" "${output_file}"
  ${ncatted} -O -a short_name,TeBS,o,c,"BDS" "${output_file}"
  ${ncatted} -O -a long_name,NEEV,o,c,"Needle-leaved evergreen" "${output_file}"
  ${ncatted} -O -a short_name,NEEV,o,c,"NE" "${output_file}"
  ${ncatted} -O -a long_name,NEDE,o,c,"Summer green needle-leaved" "${output_file}"
  ${ncatted} -O -a short_name,NEDE,o,c,"NS" "${output_file}"
  ${ncatted} -O -a long_name,H,o,c,"C3 grass" "${output_file}"
  ${ncatted} -O -a short_name,H,o,c,"TeH" "${output_file}"
  ${ncatted} -O -a long_name,HC4,o,c,"C4 grass" "${output_file}"
  ${ncatted} -O -a short_name,HC4,o,c,"TrH" "${output_file}"
  ${ncatted} -O -a long_name,bare,o,c,"bare" "${output_file}"
  ${ncatted} -O -a short_name,bare,o,c,"BSO" "${output_file}"
fi
