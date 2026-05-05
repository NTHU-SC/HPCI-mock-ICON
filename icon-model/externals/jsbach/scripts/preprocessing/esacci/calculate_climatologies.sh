#!/bin/bash
#
# calculate climatologies for annual precip, tavg (0.56*tmin+0.44*tmax) and coldest month mean temperature
#----------------------------------------------------------------

# Extract arguments
climate_input_path=$1
climate_file=$2
clim_start_year=$3
clim_end_year=$4
tmin_tag=$5
tmin_name=$6
precip_target_unit=$7
output_dir_name=$8
forcing_type=$9
temporary_dir=${10}

# set cdo command to be silent and create no history
cdo="cdo -s --no_history -b 64"

# get required vars from files for all years over which to calculate the climatology and merge to one file
years_span=${clim_start_year}-${clim_end_year}
merge_files=""
year=${clim_start_year}
while [ ${year} -le ${clim_end_year} ]; do
  echo -n "."
  ${cdo} -selvar,tmin,tmax,precip ${climate_input_path}${climate_file}${year}.nc ${temporary_dir}/tmp_${year}.nc
  merge_files="${merge_files} ${temporary_dir}/tmp_${year}.nc"
  (( year = year + 1 ))
done
${cdo} -mergetime ${merge_files} ${temporary_dir}/tmp_${years_span}.nc
echo "."

# calculate tavg
${cdo} -setname,tavg -add -mulc,0.44 -selvar,tmax ${temporary_dir}/tmp_${years_span}.nc -mulc,0.56 -selvar,tmin ${temporary_dir}/tmp_${years_span}.nc ${temporary_dir}/tmp_tavg_${years_span}.nc
${cdo} -timmean ${temporary_dir}/tmp_tavg_${years_span}.nc ${output_dir_name}/tavg_${years_span}.nc
echo ">> - calculate average temperature"
# and the average of the temperature of the coldest month
${cdo} -setname,${tmin_name} -timmean -yearmin -monmean ${temporary_dir}/tmp_tavg_${years_span}.nc ${output_dir_name}/${tmin_tag}${years_span}.nc
echo ">> - calculate coldest month mean temperature"
# and mean annual precip
echo "Forcing is: ${forcing_type}"

if [ ${forcing_type} == "gswp3" ]; then
  # Convert from kg/m²/s to mm/day  (60x60x24)
  conversion_factor=86400
  ${cdo} -mulc,${conversion_factor} -setunit,${precip_target_unit} -timmean -yearsum -selvar,precip ${temporary_dir}/tmp_${years_span}.nc ${output_dir_name}/precip_${years_span}.nc
  echo ">> - calculated mean annual precipitation"
elif [ ${forcing_type} == "crujra" ]; then
  ${cdo} -setunit,${precip_target_unit} -timmean -yearsum -selvar,precip ${temporary_dir}/tmp_${years_span}.nc ${output_dir_name}/precip_${years_span}.nc
  echo ">> - calculate mean annual precipitation"
else 
  echo ">> no proper forcing provided"
fi