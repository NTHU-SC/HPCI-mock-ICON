#!/usr/bin/bash
#
# select pfts required from the esacci tool output and remap to given grid
#----------------------------------------------------------------

# Extract arguments
input_file_esacci=$1
prepared_esacci_file=$2
grid_name=$3
esacci_out_vars=$4
esacci_final_vars=$5
temporary_dir=$6

# set cdo command to be silent and create no history
cdo="cdo -s --no_history -b 64"

###### selection of variables #######
temporary_file=${temporary_dir}/tmp_selected_esacci_pfts.nc
${cdo} -selvar,${esacci_out_vars} ${input_file_esacci} ${temporary_file}

###### map to selected grid ######
temporary_file_2=${temporary_dir}/tmp_target_grid_selected_esacci_pfts.nc
${cdo} -remapcon,${grid_name} ${temporary_file} ${temporary_file_2}
echo ">> - remapping with ${grid_name}, output file saved as tmp file ${temporary_file_2}"

###### selection of required pfts from ESACCi file #######
temporary_file_3=${temporary_dir}/tmp_select_from_mapped_selected_esacci_pfts.nc
${cdo} -selvar,${esacci_final_vars} ${temporary_file_2} ${temporary_file_3}
echo ">> - selected pfts from esacci output"

###### scale ######
cdoExprString=$(echo "${esacci_final_vars}" | tr , +)
cdoExprString="fact=${cdoExprString}"
${cdo} -expr,${cdoExprString} ${temporary_file_3} ${temporary_dir}/tmp_fract.nc
${cdo} -setvals,0,1 ${temporary_dir}/tmp_fract.nc ${temporary_dir}/tmp_fract_no_zeros.nc
${cdo} -div ${temporary_file_3} ${temporary_dir}/tmp_fract_no_zeros.nc ${prepared_esacci_file}
echo ">> - scale selected esacci pfts to 1.0, output file saved as "${prepared_esacci_file}