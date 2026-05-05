#!/bin/bash
#
# split given map into c3 and c4 according to avg temperature criterion
#----------------------------------------------------------------

# Extract arguments
prepared_esacci_file=$1 # in and output file!
clim_start_year=$2
clim_end_year=$3
output_dir_name=$4
qpft=$5
temporary_dir=$6


# set cdo command to be silent and create no history
cdo="cdo -s --no_history -b 64"

# Note: temperature is assumed to be in degC!
mean_temp_file="${output_dir_name}/tavg_${clim_start_year}-${clim_end_year}.nc"

#=============================================================================================================
# --------------- helper functions
#=============================================================================================================
echo ">> split ${qpft} into c3 and c4"

# SZ: frac_c4 <- 0.8 * ( 1-exp(-((max(0,var2D-5.5))/10.0)**2))
# temperature threshold: 5.5degC
if [ ! -f ${temporary_dir}/tmp_fract_${qpft}_c4.nc ]; then
  # calculate if not already available
  ${cdo} -expr,'c4=0.8*(1-exp(-((max(0,tavg-5.5))/10.0)^2))' ${mean_temp_file} ${temporary_dir}/tmp_fract_${qpft}_c4.nc
  ${cdo} -setmissval,NaN ${temporary_dir}/tmp_fract_${qpft}_c4.nc ${temporary_dir}/tmp_fract_${qpft}_c4_nan.nc
  mv ${temporary_dir}/tmp_fract_${qpft}_c4_nan.nc ${temporary_dir}/tmp_fract_${qpft}_c4.nc
fi

# copy in temp dir in case that its requried for checking (keep_tmp_dir=Y)
cp ${prepared_esacci_file} ${temporary_dir}/tmp_prep_esacci_file_pre_c4_1

# calculate C3 and C4 PFTs
${cdo} -setname,${qpft}C4 -mul -selvar,${qpft} ${prepared_esacci_file} ${temporary_dir}/tmp_fract_${qpft}_c4.nc ${temporary_dir}/tmp_${qpft}_C4.nc
${cdo} -mul -selvar,${qpft} ${prepared_esacci_file} -mulc,-1 -addc,-1 ${temporary_dir}/tmp_fract_${qpft}_c4.nc ${temporary_dir}/tmp_${qpft}_C3.nc

# add C4 type and replace org by C3
if [ -f ${temporary_dir}/tmp_prep_esacci_file_pre_${qpft}_c4_2 ]; then
  rm ${temporary_dir}/tmp_prep_esacci_file_pre_${qpft}_c4_2
fi
${cdo} -merge ${prepared_esacci_file} ${temporary_dir}/tmp_${qpft}_C4.nc ${temporary_dir}/tmp_prep_esacci_file_pre_${qpft}_c4_2
rm ${prepared_esacci_file}
${cdo} -replace ${temporary_dir}/tmp_prep_esacci_file_pre_${qpft}_c4_2 ${temporary_dir}/tmp_${qpft}_C3.nc ${prepared_esacci_file}

# assert new C3 + C4 = old C3...
# ${cdo} -sub -selvar,${sourcePFT} ${in_file} -add tmp_C4.nc tmp_C3.nc tmp_assert_zero.nc

