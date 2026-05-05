#!/bin/ksh

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
set -e
#-----------------------------------------------------------------------------
#   Generate
#     ICON-Land soil boundary conditions file
#-----------------------------------------------------------------------------
#
# The script is part of a script series to generate jsbach initial files
#   started by master script "create_jsbach_ini_files.sh".
#
#   bc_land_soil:
#   - Soil textures from extpar
#   - Soil properties from FAO
#   - Soil organic carbon fraction based on WISE data
#   - Deprecated variable maxmoist calculated.

# Variables exported from master script create_icon-land_ini_files.sh
extpar_file=${extpar_file}

icon_grid=${icon_grid}
atmGridID=${atmGridID}

path_bc=${bc_file_dir}
work_dir=${work_dir}
pool_prepare=${pool_prepare}

cdo=${cdo}
rm=$rm
cp=$cp

# Local variables
prog=$(basename $0)
file_bc_soil=bc_land_soil
clean_up=true  # Remove intermediate files

cd ${work_dir}
# ----------------------------------------------------------------------------
#  Soil data file
# -----------------------------
echo "${prog}: Generating ${file_bc_soil}:"

# The order of variables should not be changed
var_list="FR_SAND FR_SILT FR_CLAY SUB_FR_SAND SUB_FR_SILT SUB_FR_CLAY \
 root_depth maxmoist soil_depth fract_org_sl bclapp bclapp_mineral soil_field_cap \
 soil_field_cap_mineral heat_capacity heat_capacity_mineral heat_conductivity \
 heat_conductivity_mineral hyd_cond_sat hyd_cond_sat_mineral moisture_pot \
 moisture_pot_mineral pore_size_index pore_size_index_mineral soil_porosity \
 soil_porosity_mineral wilting_point wilting_point_mineral residual_water \
 residual_water_mineral fao"

# Preparations that need to be done outside the parallel loop
${cdo} splitvar ${pool_prepare}/05/soil_parameter_05_apr2024.nc   05_

# Process root depth, as it is needed for soil depth calculation
if [[ ${var_list} == *"soil_depth"* || ${var_list} == *"root_depth"* ]]; then
  # ----------------------
  var=root_depth
  # ----------------------
  echo "   ${var} ..."
  minimum_root_depth=0.15  # Note: root depth of GLCC data used in extpar ranges from 0.3 to 2.0 m.
  # Get root depth from extpar file and extrapolate land values to ocean area
  ${cdo} -setmisstodis -setmissval,0 -selvar,ROOTDP ${extpar_file} extpar_ROOTDP.nc
  # Set minimum root depth - also use this value for complete ocean cells. Root depth of glaciers is zero.
  ${cdo} -setattribute,root_depth@long_name="Root depth" -setattribute,root_depth@units="m" \
      -setname,root_depth -mul notglac.nc -maxc,${minimum_root_depth} \
      -mul non-glac-land-mask.nc extpar_ROOTDP.nc  ${atmGridID}_root_depth.nc
  # Copy for later usage in generate_ic_land_soil
  $cp ${atmGridID}_root_depth.nc ${atmGridID}_root_depth.keep
fi


for var in ${var_list}; do
  {
  echo "   ${var} ..."
  case ${var} in
  # ----------------------
    FR_SAND | FR_SILT | FR_CLAY )
  # ----------------------
      # Soil texture fractions in the upper 30 cm of the soil column
      # Convert from percent to fraction and remove 'institution' attribute (i.e. DWD)
      ${cdo} mulc,0.01 -selvar,${var} ${extpar_file} ${atmGridID}_${var}.tmp
      ncatted -a institution,global,d,c,sng       ${atmGridID}_${var}.tmp

      # Generate mask for valid data, used below to mask out grid cells
      # where any of the texture variables are missing
      ${cdo} gec,0. ${atmGridID}_${var}.tmp ${atmGridID}_${var}_MASK.nc
      ;;
  # ----------------------
    SUB_FR_SAND | SUB_FR_SILT | SUB_FR_CLAY )
  # ----------------------
      # Soil texture fractions in deeper soil (below 30 cm)
      # Convert from percent to fraction and remove institution attribute (i.e. DWD)
      ${cdo} mulc,0.01 -selvar,${var} ${extpar_file} ${atmGridID}_${var}.tmp
      ncatted -a institution,global,d,c,sng       ${atmGridID}_${var}.tmp

      # Generate mask for valid data, used below to mask out grid cells
      # where any of the texture variables are missing
      ${cdo} gec,0. ${atmGridID}_${var}.tmp ${atmGridID}_${var}_MASK.nc
      ;;
  # ----------------------
    root_depth )
  # ----------------------
      # Root depth already calculated above ...
      ;;
  # ----------------------
    maxmoist )
  # ----------------------
      # Maxmoist calculated below ...
      ;;
  # ----------------------
    soil_depth )
  # ----------------------
      # In case soil depth is smaller than root depth it needs to be adapted: As the root depth
      # data set in extpar is more recent than the soil data from echam we rather adapt soil
      # depth to root depth - and not vice versa as we did in previous versions.
      ${cdo} selvar,soildepth ${pool_prepare}/05/soil_parameters_05.nc 05_soil_depth.nc
      ${cdo} remap,${icon_grid},wgt05.nc -setmisstodis 05_soil_depth.nc \
          ${atmGridID}_soil_depth.tmp
      ${cdo} -setattribute,soil_depth@long_name="Soil depth until bedrock" -setattribute,soil_depth@units="m" \
          -setvar,soil_depth -setmisstoc,0. -ifthen non-glac-land-mask.nc \
          -max ${atmGridID}_soil_depth.tmp ${atmGridID}_root_depth.nc ${atmGridID}_soil_depth.nc
      ncatted -a code,${var},d,c, ${atmGridID}_${var}.nc
      # Copy for later usage in generate_ic_land_soil
      $cp ${atmGridID}_soil_depth.nc ${atmGridID}_soil_depth.keep
      ;;
  # ----------------------
    fract_org_sl )
  # ----------------------
      wise_data=${pool_prepare}/05/wise_som-fract_05deg_jsbsoillayers.nc
      ${cdo} -selvar,${var} -sellonlatbox,0,360,-90,90 $wise_data 05_${var}.tmp
      # Extrapolation
      ${cdo} -setmisstodis -ifthen 05_non-glac-land.nc 05_${var}.tmp 05_${var}.nc
      ${cdo} remap,${icon_grid},wgt05.nc 05_$var.nc  ${atmGridID}_${var}.tmp
      ${cdo} -setmisstoc,0 -ifthen non-glac-land-mask.nc ${atmGridID}_${var}.tmp \
          ${atmGridID}_${var}.nc
      ;;
  # ----------------------
    bclapp | bclapp_mineral | soil_field_cap | soil_field_cap_mineral | heat_capacity | heat_capacity_mineral \
    | heat_conductivity | heat_conductivity_mineral | hyd_cond_sat | hyd_cond_sat_mineral | moisture_pot \
    | moisture_pot_mineral | pore_size_index | pore_size_index_mineral | soil_porosity | soil_porosity_mineral \
    | wilting_point | wilting_point_mineral | residual_water | residual_water_mineral )
  # ----------------------
      ${cdo} -setvar,${var} -mul -remap,${icon_grid},wgt05.nc -invertlat -setmisstodis 05_${var}.nc \
                                 notglac.nc    ${atmGridID}_${var}.tmp
      case ${var} in
        bclapp* )
          ${cdo} -setmisstoc,4.5    -ifthen land.nc ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
          ;;
        soil_field_cap* )
          ${cdo} -setmisstoc,0.229  -ifthen land.nc ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
          # Copy for later usage in ic_soil script
          $cp ${atmGridID}_${var}.nc ${atmGridID}_${var}.keep
          ;;
        heat_capacity* )
          ${cdo} -setmisstoc,2.e+6  -ifthen land.nc ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
          ;;
        heat_conductivity* )
          ${cdo} -setmisstoc,7.     -ifthen land.nc ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
          ;;
        hyd_cond_sat* )
          ${cdo} -setmisstoc,5.e-6  -ifthen land.nc ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
          ;;
        moisture_pot* )
          ${cdo} -setmisstoc,-0.15  -ifthen land.nc ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
          ;;
        pore_size_index* )
          ${cdo} -setmisstoc,0.2    -ifthen land.nc ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
          ;;
        soil_porosity* )
          ${cdo} -setmisstoc,0.45   -ifthen land.nc ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
          ;;
        wilting_point* )
          ${cdo} -setmisstoc,0.15   -ifthen land.nc ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
          ;;
        residual_water* )
          ${cdo} -setmisstoc,0.05   -ifthen land.nc ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
          ;;
        * )
          ${cdo} -setmissval,-9.e33 -ifthen land.nc ${atmGridID}_${var}.tmp ${atmGridID}_${var}.nc
          ;;
      esac
      ;;
  # ----------------------
    fao )
  # ----------------------
      # We use old ECHAM5 data here! TODO: replace
      ${cdo} -remaplaf,${icon_grid} -selvar,FAO ${pool_prepare}/05/ECHAM6/05_jan_surf.nc \
          ${atmGridID}_fao.tmp
      ${cdo} setname,fao -setmisstoc,0 -ifthen non-glac-land-mask.nc ${atmGridID}_fao.tmp \
          ${atmGridID}_fao.nc
      ;;
  # ----------------------
    * )
  # ----------------------
      echo "ERROR: Variable ${var} not known to be a bc_land_soil file variable."
      exit 1
      ;;
  esac
  } &  # comment out '&' for serial processing
done
wait

# Further process soil texture data (needs to be done outside the above parallel loop)
#                -------------------
if [[ ${var_list} == *"FR_SAND"* || ${var_list} == *"FR_SILT"* || ${var_list} == *"FR_CLAY"* ]]; then
  if [[ -f ${atmGridID}_SUB_FR_SAND_MASK.nc \
     && -f ${atmGridID}_SUB_FR_SILT_MASK.nc \
     && -f ${atmGridID}_SUB_FR_CLAY_MASK.nc  ]]; then
    # Generate mask with grid cells where any of the soil texture data is missing
    ${cdo} -O ensmin \
        ${atmGridID}_SUB_FR_SAND_MASK.nc ${atmGridID}_SUB_FR_SILT_MASK.nc ${atmGridID}_SUB_FR_CLAY_MASK.nc \
        ${atmGridID}_SUB_SAND_SILT_CLAY_MASK.nc

    # Apply this mask to the texture data
    ${cdo} mul ${atmGridID}_SUB_FR_SAND.tmp ${atmGridID}_SUB_SAND_SILT_CLAY_MASK.nc ${atmGridID}_SUB_FR_SAND_0.nc
    ${cdo} mul ${atmGridID}_SUB_FR_SILT.tmp ${atmGridID}_SUB_SAND_SILT_CLAY_MASK.nc ${atmGridID}_SUB_FR_SILT_0.nc
    ${cdo} mul ${atmGridID}_SUB_FR_CLAY.tmp ${atmGridID}_SUB_SAND_SILT_CLAY_MASK.nc ${atmGridID}_SUB_FR_CLAY_0.nc

    # Fill data gaps with 50 percent sand, 25 percent silt and 25 percent clay fraction
    ${cdo} add ${atmGridID}_SUB_FR_SAND_0.nc -mulc,-0.5  -addc,-1. ${atmGridID}_SUB_SAND_SILT_CLAY_MASK.nc \
               ${atmGridID}_SUB_FR_SAND.nc
    ${cdo} add ${atmGridID}_SUB_FR_SILT_0.nc -mulc,-0.25 -addc,-1. ${atmGridID}_SUB_SAND_SILT_CLAY_MASK.nc \
               ${atmGridID}_SUB_FR_SILT.nc
    ${cdo} add ${atmGridID}_SUB_FR_CLAY_0.nc -mulc,-0.25 -addc,-1. ${atmGridID}_SUB_SAND_SILT_CLAY_MASK.nc \
               ${atmGridID}_SUB_FR_CLAY.nc
  else
    echo "ERROR: Deep soil texture variables incomplete"
    exit 1
  fi
  if [[ -f ${atmGridID}_FR_SAND_MASK.nc \
     && -f ${atmGridID}_FR_SILT_MASK.nc \
     && -f ${atmGridID}_FR_CLAY_MASK.nc  ]]; then
    # Generate mask with grid cells where any of the soil textures are missing.
    ${cdo} -O ensmin \
        ${atmGridID}_FR_SAND_MASK.nc ${atmGridID}_FR_SILT_MASK.nc ${atmGridID}_FR_CLAY_MASK.nc \
        ${atmGridID}_SAND_SILT_CLAY_MASK.nc

    # Apply this mask to the texture data
    ${cdo} mul ${atmGridID}_FR_SAND.tmp ${atmGridID}_SAND_SILT_CLAY_MASK.nc ${atmGridID}_FR_SAND_0.nc
    ${cdo} mul ${atmGridID}_FR_SILT.tmp ${atmGridID}_SAND_SILT_CLAY_MASK.nc ${atmGridID}_FR_SILT_0.nc
    ${cdo} mul ${atmGridID}_FR_CLAY.tmp ${atmGridID}_SAND_SILT_CLAY_MASK.nc ${atmGridID}_FR_CLAY_0.nc

    # Fill data gaps with 50 percent sand, 25 percent silt and 25 percent clay fraction
    ${cdo} add ${atmGridID}_FR_SAND_0.nc -mulc,-0.5  -addc,-1. ${atmGridID}_SAND_SILT_CLAY_MASK.nc \
               ${atmGridID}_FR_SAND.nc
    ${cdo} add ${atmGridID}_FR_SILT_0.nc -mulc,-0.25 -addc,-1. ${atmGridID}_SAND_SILT_CLAY_MASK.nc \
               ${atmGridID}_FR_SILT.nc
    ${cdo} add ${atmGridID}_FR_CLAY_0.nc -mulc,-0.25 -addc,-1. ${atmGridID}_SAND_SILT_CLAY_MASK.nc \
               ${atmGridID}_FR_CLAY.nc
  else
    echo "ERROR: Upper soil texture variables incomplete"
    exit 1
  fi
fi

# Maxmoist: water content of the root zone at field capacity
# --------
# Note: maxmoist is calculated here - outside the variable loop - because it depends on
#       soil_field_cap and root_depth generated above.
# Note: maxmoist is no longer read from initial files in recent ICON-Land hydrology versions
# as it depends on soil texture and organic layer fractions. It is calculated at runtime
# from field capacity and root depth. For backward compatibility we re-calculate it here to
# be consistent with root depth.

${cdo} -setattribute,maxmoist@long_name="Maximum amount of soil moisture" -setattribute,maxmoist@units="m" \
    -setname,maxmoist -mul ${atmGridID}_root_depth.nc ${atmGridID}_soil_field_cap.nc \
    ${atmGridID}_maxmoist.nc


# Merge newly calculated variables to the bc soil file
# -------------
[[ -f ${file_bc_soil}.tmp ]] && ${rm} ${file_bc_soil}.tmp
for var in ${var_list}; do
  if [[ ! -f ${file_bc_soil}.tmp ]]; then
    # first variable
    ${cp} ${atmGridID}_${var}.nc  ${file_bc_soil}.tmp
  else
    # following variable
    mv ${file_bc_soil}.tmp ${file_bc_soil}.tmp2
    ${cdo} -O merge ${file_bc_soil}.tmp2  ${atmGridID}_${var}.nc ${file_bc_soil}.tmp
  fi
done

${cdo} --no_history setattribute,history="${history_att}" ${file_bc_soil}.tmp \
       ${path_bc}/${file_bc_soil}.nc

# Clean up
if [[ ${clean_up} == true ]]; then
  ${rm} *.tmp*
fi
echo "${prog}:     ${file_bc_soil}.nc         done"

exit 0
