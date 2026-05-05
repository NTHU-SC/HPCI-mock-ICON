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
#     ICON-Land physical boundary conditions file
#-----------------------------------------------------------------------------
#
# The script is part of a script series to generate jsbach initial files
#   started by master script "create_jsbach_ini_files.sh"
#
#   bc_land_phys:
#   - Calculate albedo, roughness length, forest fraction, etc. based
#     on extpar data.
#
# Variables exported from master script create_icon-land_ini_files.sh
extpar_file=${extpar_file}       # Extpar file name

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
file_bc_phys=bc_land_phys
clean_up=true  # Remove intermediate files

cd ${work_dir}
# ----------------------------------------------------------------------------
#  Physical boundary conditions
# -----------------------------
echo "${prog}: Generating bc_land_phys:"

# The order of variables should not be changed
var_list="lai_clim veg_fract roughness_length \
 albedo albedo_veg_vis albedo_veg_nir albedo_soil_vis albedo_soil_nir \
 forest_fract skin_conductivity"

#  Preparations needed outside the parallel loop below
# -------------
# - For albedo
if [[ ${var_list} == *"albedo_"* ]]; then
  ${cdo} splitvar ${pool_prepare}/05/albedo_05.nc   tmp_05_

  # Define albedo compensate (as in jsbach_init_file.f90)
  #
  # Surface albedo is based on MODIS white-sky albedo data. White-sky albedo (reflectivity of diffuse radiation),
  # assumes the same solar zenith angle for all latitudes. "Albedo_compensate" corrects too low surface albedo
  # in the high latitudes resulting from always large zenith angles.
  ${cdo} const,0.0,${pool_prepare}/05/albedo_05.nc  albedo_comp_0.nc
  ${cdo} setclonlatbox,0.04,0,360,70,90   \
        -setclonlatbox,0.03,0,360,60,70   \
        -setclonlatbox,0.02,0,360,50,60   \
        -setclonlatbox,0.01,0,360,40,50   \
        -setclonlatbox,0.01,0,360,-50,-40 \
        -setclonlatbox,0.02,0,360,-60,-50 \
        -setclonlatbox,0.03,0,360,-70,-60 \
        -setclonlatbox,0.04,0,360,-90,-70 albedo_comp_0.nc  albedo_compensate.nc
  ${rm} albedo_comp_0.nc
fi

# - for NDVI vegetation fractions and LAI climatology
# ---------------------------------
# TODO: The NDVI data from extpar seams to include ocean fractions - it is not relative
#       to the land fraction, as expected by jsbach. (Dividing by FR_LAND did not work properly)
# For that reason it also does not make sense to extrapolate the data to complete ocean cells.
if [[ ${var_list} == *"lai_clim"* ]] || [[ ${var_list} == *"veg_fract"* ]]; then
  ${cdo} selvar,NDVI ${extpar_file} extpar_NDVI.nc
fi


for var in ${var_list}; do
  echo "   ${var} ..."
  {
  case ${var} in
  # ----------------------
    roughness_length )
  # ----------------------
      # Select extpar field of roughness length
      ${cdo} selvar,Z0 ${extpar_file} extpar_Z0.nc
      # Change variable name and attributes
      ${cdo} -setattribute,${var}@long_name="Surface roughness length" -setattribute,${var}@units="m" \
          -chname,Z0,${var} extpar_Z0.nc ${atmGridID}_roughness_length.nc
      ;;
  # ----------------------
    albedo )
  # ----------------------
      value_albedo_max=60. # maximum value of background albedo accepted [%]
      value_albedo_min=6.  # minimum value of background albedo accepted [%]
      # Remove very high albedo of incorrect glacier points and very low albedo at the coast and in snow areas
      ${cdo} selvar,ALB ${extpar_file} extpar_ALB.nc
      ${cdo} expr,"mask=ALB>${value_albedo_min} && ALB<${value_albedo_max}" extpar_ALB.nc extpar_ALB_mask.nc
      ${cdo} ifthen extpar_ALB_mask.nc extpar_ALB.nc extpar_ALB.tmp
      ${cdo} -setmisstodis extpar_ALB.tmp extpar_ALB_filled.nc

      # Average in time (weighted by days per month) and divide by 100 ([% reflection] ==> [albedo])
      ${cdo} -setattribute,${var}@long_name="Surface albedo" -setattribute,${var}@units="1" \
          -setname,albedo -divc,100. -yearmonmean extpar_ALB_filled.nc ${atmGridID}_albedo.nc
      ;;
  # ----------------------
    albedo_veg_vis | albedo_veg_nir | albedo_soil_vis | albedo_soil_nir )
  # ----------------------
      # Soil and vegetation albedo in visible and NIR range
      ${cdo} sellonlatbox,0,360,-90,90 \
        -add tmp_05_$var.nc albedo_compensate.nc 05_$var.nc

      [[ ${var} == albedo_veg_vis ]]  && longname="Vegetation albedo in the visible range"
      [[ ${var} == albedo_veg_nir ]]  && longname="Vegetation albedo in the NIR"
      [[ ${var} == albedo_soil_vis ]] && longname="Soil albedo in the visible range"
      [[ ${var} == albedo_soil_nir ]] && longname="Soil albedo in the NIR"
      ${cdo} -setattribute,${var}@long_name="${longname}" -setattribute,${var}@units="1" \
          -remap,${icon_grid},wgt05.nc 05_$var.nc ${atmGridID}_${var}.nc
      ncatted -a code,${var},d,c, ${atmGridID}_${var}.nc
      ;;
  # ----------------------
    forest_fract )
  # ----------------------
      # Deciduous forest
      ${cdo} selvar,FOR_D ${extpar_file} extpar_FOR_D.nc
      # Evergreen forest
      ${cdo} selvar,FOR_E ${extpar_file} extpar_FOR_E.nc
      # All forest
      ${cdo} add extpar_FOR_D.nc extpar_FOR_E.nc extpar_FOR.nc
      ${cdo} -setattribute,${var}@long_name="Forest fraction" -setattribute,${var}@units="1" \
          -chname,FOR_D,forest_fract extpar_FOR.nc ${atmGridID}_forest_fract.nc
      ;;
  # ----------------------
    veg_fract )
  # ----------------------
      # NDVI ==> veg_fract
      upper_bound_NDVI_veg_fract=0.92
      lower_bound_NDVI_veg_fract=0.12
      upper_bound_veg_fract=0.98
      fact_NDVI_veg_fract=1.225 # should be: upper_bound_veg_fract / (upper_bound_NDVI_veg_fract - lower_bound_NDVI_veg_fract)
      ${cdo} -maxc,${lower_bound_NDVI_veg_fract} -minc,${upper_bound_NDVI_veg_fract} \
            extpar_NDVI.nc extpar_NDVI_${var}_bound.nc
      ${cdo} -setattribute,${var}@long_name="Vegetation fraction" -setattribute,${var}@units="1" \
            -expr,"veg_fract=${fact_NDVI_veg_fract}*(max(0.,NDVI-${lower_bound_NDVI_veg_fract}))" \
            extpar_NDVI_${var}_bound.nc ${atmGridID}_veg_fract.nc
      ;;
  # ----------------------
    lai_clim )
  # ----------------------
      # NDVI ==> lai_clim
      upper_bound_NDVI_lai_clim=0.8
      lower_bound_NDVI_lai_clim=0.12
      upper_bound_lai_clim=6.
      fact_NDVI_lai_clim=8.8235294 # should be: upper_bound_lai_clim / (upper_bound_NDVI_lai_clim - lower_bound_NDVI_lai_clim)
      ${cdo} maxc,${lower_bound_NDVI_lai_clim} -minc,${upper_bound_NDVI_lai_clim} \
            extpar_NDVI.nc extpar_NDVI_${var}_bound.nc
      ${cdo} -setattribute,${var}@long_name="Leaf area index" -setattribute,${var}@units="1" \
            -expr,"lai_clim=${fact_NDVI_lai_clim}*(max(0.,NDVI-${lower_bound_NDVI_lai_clim}))" \
            extpar_NDVI_${var}_bound.nc ${atmGridID}_lai_clim.nc
      ;;
  # ----------------------
    skin_conductivity )
  # ----------------------
      # Skin layer conductivity
      ${cdo} setvar,${var} -ifthen notsea.nc -setmisstonn -setmissval,0 \
          -selvar,SKC ${extpar_file}  ${atmGridID}_${var}.nc
      ;;
  # ----------------------
    * )
  # ----------------------
      echo "ERROR: Variable ${var} not known to be a bc_land_phys file variable."
      exit 1
      ;;
  esac
  } &  # comment out '&' for serial processing
done
wait

# Merge newly calculated variables to the bc physics file
# -------------
[[ -f ${file_bc_phys}.tmp ]] && ${rm} ${file_bc_phys}.tmp
for var in ${var_list}; do
  if [[ ! -f ${file_bc_phys}.tmp ]]; then
    # first variable
    ${cp} ${atmGridID}_${var}.nc  ${file_bc_phys}.tmp
  else
    # following variable
    mv ${file_bc_phys}.tmp ${file_bc_phys}.tmp2
    ${cdo} -O merge ${file_bc_phys}.tmp2  ${atmGridID}_${var}.nc ${file_bc_phys}.tmp
  fi
done

${cdo} --no_history setattribute,history="${history_att}" ${file_bc_phys}.tmp \
       ${path_bc}/${file_bc_phys}.nc

# Clean up
if [[ ${clean_up} == true ]]; then
  ${rm} *.tmp*
  ${rm} extpar_*
fi
echo "${prog}:     ${file_bc_phys}.nc         done"

exit 0
