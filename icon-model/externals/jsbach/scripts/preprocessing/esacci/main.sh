#!/usr/bin/bash
#
# script for esacci preprocessing, created files will be located in the folder specified in the config.txt
# please refer to the README (STEP 2) for further information
#----------------------------------------------------------------

#--------------------------------------------------------------------
### Batch Queuing System is SLURM
#SBATCH --output=esacci.o%j
#SBATCH --error=esacci.o%j
#SBATCH --account=mj0143
#SBATCH --partition=interactive
#SBATCH --ntasks=1

iq_scripts_repro=YOURSCRIPTSPATH
scriptDir=${iq_scripts_repro}/preprocessing/esacci/

# get settings from config file
ln -sf ${scriptDir}/config.txt ./config.txt
source ./config.txt
# add info to history about the scripts and the commit
gitCommitNum=$(git -C ${iq_scripts_repro} rev-parse --short=8 HEAD)
gitShortStatusFirstLetters=$(git -C ${iq_scripts_repro} status -s | head -n 1 | cut -c 1-2)
esacci_comment="${esacci_comment} -- using preprocessing/esacci/main.sh - commit ${gitCommitNum} (${gitShortStatusFirstLetters})"

module load nco
module unload cdo
module load cdo/2.5.0-gcc-11.2.0
# set ncatted command to create no history
ncatted="ncatted -h"

#=============================================================================================================
# --------------- set-up io
#=============================================================================================================
if [ ! -d ${scriptDir} ]; then
  echo ">> ERROR: ${scriptDir} does not exists"
  exit
fi

# make output and temporary directories, if not already there
if [ ! -d ${output_dir_name} ]; then
  mkdir -p ${output_dir_name}
  mkdir ${output_dir_name}/${model}
  echo ">> created: ${output_dir_name}"
else
  echo ">> WARNING: ${output_dir_name} already exists"
fi
if [ ! -d ${temporary_dir} ]; then
  mkdir -p ${temporary_dir}
else
  echo ">> WARNING: ${temporary_dir} already exists"
fi
# add the config to the output folder
cp config.txt ${output_dir_name}/${model}/config.txt

# ---- input files
# check for esacci file
if [[ ! -f ${input_file_esacci} ]]; then
  echo ">> ERROR: passed esacci file '"${input_file_esacci}"' does not exist! Please check."
  exit
fi

# ---- output file
prepared_esacci_file=${output_dir_name}/${model}/${model}_prep_ESACCI${esacci_year}.nc

#=============================================================================================================
# --------------- start of calculations
#=============================================================================================================
echo ">> preprocess esacci ${esacci_year} for ${model} with grid ${grid_name}"

# ---- calculate climatologies for annual precip, tavg (0.56*tmin+0.44*tmax) and coldest month mean temperature
if [ ! -f ${output_dir_name}/precip_${clim_start_year}-${clim_end_year}.nc ]; then
  echo ">> calculate climatologies"
  ${scriptDir}/calculate_climatologies.sh ${climate_input_path} ${climate_file} ${clim_start_year} ${clim_end_year} ${tmin_tag} ${tmin_name}  ${precip_target_unit} ${output_dir_name} ${forcing_type} ${temporary_dir}
else
  echo ">> Note: climatologies have already been calculated before and are available in ${output_dir_name}"
fi

# ---- prepare esacci file
# -- select pft variables, remap, select final pfts and rescale
echo ">> call select_esacci_pfts_and_remap.sh"
${scriptDir}/select_esacci_pfts_and_remap.sh ${input_file_esacci} ${prepared_esacci_file} \
                                               ${grid_name} ${esacci_out_vars} ${esacci_final_vars} ${temporary_dir}

# -- split ESACCI grasses into c3 and c4 according to avg temperature criterion (prepared_esacci_file is in and output)
echo ">> call split_into_c3_and_c4.sh"
qpft="H"
${scriptDir}/split_into_c3_and_c4.sh ${prepared_esacci_file} ${clim_start_year} ${clim_end_year} ${output_dir_name} ${qpft} ${temporary_dir}

# -- tree pft splits (prepared_esacci_file is in and output)
if [ ${model} == "jsbach" ]; then
  # -- for jsbach: split evergreen and deciduous in tropical and extratropical and shrubs in raingreen and deciudous, following the 1991-2020 Koppen-Geiger map
  echo ">> call split_jsbach_woody_pfts.sh (jsbach)"
  ${scriptDir}/split_jsbach_woody_pfts.sh ${prepared_esacci_file} ${jsb_KG_map} ${output_dir_name} ${temporary_dir}
elif [ ${model} == "quincy" ]; then
  # -- split broadleaved pfts from given map according to temperature and precip criteria
  echo ">> call split_broadleaved_pfts.sh (quincy)"
  ${scriptDir}/split_broadleaved_pfts.sh ${prepared_esacci_file} ${clim_start_year} ${clim_end_year} ${output_dir_name} ${precip_threshold} ${min_mean_temperature_threshold} ${temporary_dir}
fi

# set nan mask
${cdo} ifthen -expr,"mask=${expStringEsacciMask}" ${prepared_esacci_file} ${prepared_esacci_file} ${prepared_esacci_file}_nans
mv ${prepared_esacci_file}_nans ${prepared_esacci_file}

# change global attributes on prepared esacci file
${ncatted} -O -a ,global,d,, ${prepared_esacci_file}
${ncatted} -O -a esacci_reference,global,o,c,"${esacci_reference}" "${prepared_esacci_file}"
${ncatted} -O -a esacci_comment,global,o,c,"${esacci_comment}" "${prepared_esacci_file}"
${ncatted} -O -a esacci_preprocessing,global,o,c,"${esacci_preprocessing}" "${prepared_esacci_file}"
# add pft info
echo ">> call add_pft_info_to_attributes.sh"
${scriptDir}/add_pft_info_to_attributes.sh ${model} ${prepared_esacci_file}
echo ">> attributes of ${prepared_esacci_file} file have been changed, file ready!"

# ---- if demanded (in the config) remove temporary directory
if [ "${keep_tmp_dir}" == "N" ]; then
  rm -r ${temporary_dir}
fi
