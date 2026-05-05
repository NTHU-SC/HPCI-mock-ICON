#!/bin/bash
#
# split broadleaved pfts from given map according to temperature and precip criteria
# [criteria taken from /Net/Groups/BSI/people/szaehle/Projects/GCP/LandUse/Synmap/gather_synmap_PFT.R]
#----------------------------------------------------------------

# Extract arguments
prepared_esacci_file=$1 # in and output file!
clim_start_year=$2
clim_end_year=$3
output_dir_name=$4
precip_threshold=$5
min_mean_temperature_threshold=$6
tmp_dir=$7

# set cdo command to be silent and create no history
cdo="cdo -s --no_history -b 64"

#=============================================================================================================
# --------------- main
#=============================================================================================================
BE_source="BREV"  # var name of broadleaved evergreen pft
E_rain_pft="TrBE" # broad leaved evergreen (rainforest)
E_xeric_pft="TeBE" # broad leaved evergreen (xeric forest)

BD_source="BRDE"  # var name of broadleaved deciduous pft
D_rain_pft="TrBR" # broad leaved deciduous (rain green)
D_summer_pft="TeBS" # broad leaved deciduous (summer green)

# Note: temperature is assumed to be in degC, precip in mm
mean_annual_precip_file="${output_dir_name}/precip_${clim_start_year}-${clim_end_year}.nc"
avg_min_mean_temp_file="${output_dir_name}/avg_of_coldest_month_mean_temperature_${clim_start_year}-${clim_end_year}.nc"

# get a help file with zeros only
${cdo} -mulc,0 -selvar,${BE_source} ${prepared_esacci_file} ${tmp_dir}/tmp_zeros.nc

# change source names to rainforest / rain green pft
${cdo} -chname,${BE_source},${E_rain_pft} -chname,${BD_source},${D_rain_pft} ${prepared_esacci_file} ${tmp_dir}/tmp_in_file_renamed_rain_pfts.nc

# rainforest vs xeric forest
${cdo} -mul -gtc,${precip_threshold} ${mean_annual_precip_file} \
        -gtc,${min_mean_temperature_threshold} ${avg_min_mean_temp_file} ${tmp_dir}/tmp_rainforest_mask.nc
${cdo} -ifthenelse ${tmp_dir}/tmp_rainforest_mask.nc \
        -selvar,${E_rain_pft} ${tmp_dir}/tmp_in_file_renamed_rain_pfts.nc ${tmp_dir}/tmp_zeros.nc ${tmp_dir}/tmp_TrBE.nc
${cdo} -mulc,-1 -addc,-1 ${tmp_dir}/tmp_rainforest_mask.nc ${tmp_dir}/tmp_rainforest_mask_negated.nc
${cdo} -setname,${E_xeric_pft} -ifthenelse ${tmp_dir}/tmp_rainforest_mask_negated.nc -selvar,${E_rain_pft} \
          ${tmp_dir}/tmp_in_file_renamed_rain_pfts.nc ${tmp_dir}/tmp_zeros.nc ${tmp_dir}/tmp_TeBE.nc

# get rain vs summergreen mask
${cdo} -gtc,${min_mean_temperature_threshold} ${avg_min_mean_temp_file} ${tmp_dir}/tmp_raingreen_mask.nc
${cdo} -ifthenelse ${tmp_dir}/tmp_raingreen_mask.nc \
         -selvar,${D_rain_pft} ${tmp_dir}/tmp_in_file_renamed_rain_pfts.nc ${tmp_dir}/tmp_zeros.nc ${tmp_dir}/tmp_TrBR.nc
${cdo} -mulc,-1 -addc,-1 ${tmp_dir}/tmp_raingreen_mask.nc ${tmp_dir}/tmp_raingreen_mask_negated.nc
${cdo} -setname,${D_summer_pft} -ifthenelse ${tmp_dir}/tmp_raingreen_mask_negated.nc \
          -selvar,${D_rain_pft} ${tmp_dir}/tmp_in_file_renamed_rain_pfts.nc ${tmp_dir}/tmp_zeros.nc ${tmp_dir}/tmp_TeBS.nc

# replace rainforest and rain green pft
${cdo} -merge ${tmp_dir}/tmp_TrBE.nc ${tmp_dir}/tmp_TrBR.nc ${tmp_dir}/tmp_new_rain_pft_values.nc
${cdo} -replace ${tmp_dir}/tmp_in_file_renamed_rain_pfts.nc ${tmp_dir}/tmp_new_rain_pft_values.nc ${tmp_dir}/tmp_in_file_updated_rain_pft_values.nc

# add xeric and summer green pft
mv ${prepared_esacci_file} ${tmp_dir}/tmp_prep_esacci_file_pre_broadleave_split
${cdo} -merge ${tmp_dir}/tmp_in_file_updated_rain_pft_values.nc ${tmp_dir}/tmp_TeBE.nc ${tmp_dir}/tmp_TeBS.nc ${prepared_esacci_file}

# assert new = old
# ${cdo} -sub -selvar,${BE_source} ${in_file} -add ${tmp_dir}/tmp_TeBE.nc ${tmp_dir}/tmp_TrBE.nc ${tmp_dir}/tmp_assert_zero_BE.nc
# ${cdo} -sub -selvar,${BD_source} ${in_file} -add ${tmp_dir}/tmp_TeBS.nc ${tmp_dir}/tmp_TrBR.nc ${tmp_dir}/tmp_assert_zero_BD.nc
echo ">> Broadleaved pfts splitted"


