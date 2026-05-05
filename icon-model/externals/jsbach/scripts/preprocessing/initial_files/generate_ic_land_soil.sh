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
#     ICON-Land soil initial conditions file
#-----------------------------------------------------------------------------
#
# The script is part of a script series to generate jsbach initial files
#   started by master script "create_jsbach_ini_files.sh".
#
#   ic_land_soil:
#   - Initial soil mositure and snow cover

# Variables exported from master script create_icon-land_ini_files.sh
extpar_file=${extpar_file}

icon_grid=${icon_grid}
atmGridID=${atmGridID}

path_bc=${bc_file_dir}
work_dir=${work_dir}

cdo=${cdo}
rm=$rm
cp=$cp

# Local variables
prog=$(basename $0)
file_ic_soil=ic_land_soil
clean_up=true  # Remove intermediate files

cd ${work_dir}
# ----------------------------------------------------------------------------
#  Soil initial conditions
# -----------------------------
echo "${prog}: Generating ${file_ic_soil}:"

# The order of variables should not be changed
var_list="init_moist surf_temp layer_moist snow"

# Preparations
ini_mon=01

# Assume relative soil moisture from veg_ratio_max (1-desert)
# We define masks for specific ranges of veg_ration_max, i.e. 0.0-0.5, 0.5-0.6, etc. ...
${cdo} -setrtoc,-1.,0.0,0.0 -setrtoc,0.0,0.5,1.0 -setrtoc,0.5,1.0,-1. ${atmGridID}_veg_ratio_max.keep veg_ratio_max.00-05.tmp
${cdo} -setrtoc,0.0,0.5,0.0 -setrtoc,0.5,0.6,1.0 -setrtoc,0.6,1.0,0.0 ${atmGridID}_veg_ratio_max.keep veg_ratio_max.05-06.tmp
${cdo} -setrtoc,0.0,0.6,0.0 -setrtoc,0.6,0.9,1.0 -setrtoc,0.9,1.0,0.0 ${atmGridID}_veg_ratio_max.keep veg_ratio_max.06-09.tmp
${cdo} -setrtoc,0.0,0.9,0.0 -setrtoc,0.9,.95,1.0 -setrtoc,.95,1.0,0.0 ${atmGridID}_veg_ratio_max.keep veg_ratio_max.09-95.tmp
${cdo} -setrtoc,0.0,.95,0.0 -setrtoc,.99,1.0,1.0                      ${atmGridID}_veg_ratio_max.keep veg_ratio_max.95-10.tmp
rm -f init_moist_rel.tmp

# ... and assign a specific relative soil moisture values for each of the masks.
${cdo} mul non-glac-land-mask.nc -enssum -mulc,0.35 veg_ratio_max.00-05.tmp \
                                         -mulc,0.4  veg_ratio_max.05-06.tmp \
                                         -mulc,0.5  veg_ratio_max.06-09.tmp \
                                         -mulc,0.8  veg_ratio_max.09-95.tmp \
                                                    veg_ratio_max.95-10.tmp init_moist_rel.tmp

for var in ${var_list}; do
  {
  echo "   ${var} ..."
  case ${var} in
  # ----------------------
    init_moist )
  # ----------------------
      # Relative rootzone soil moisture calculated above - outside parallel loop
      # Convert relative to absolute soil moisture: [m3/m3] -> [m water equivalent]
      ${cdo} -setattribute,init_moist@long_name="Soil water content" \
          -setattribute,init_moist@units="m water equivalent" \
          -setmissval,-9.e+33 -setname,init_moist \
          -mul init_moist_rel.tmp \
          -mul ${atmGridID}_soil_field_cap.keep  ${atmGridID}_root_depth.keep \
          ${atmGridID}_${var}.nc
      ;;
  # ----------------------
    layer_moist )
  # ----------------------
      levels="0.065 0.319 1.232 4.134 9.834"

      above=0.
      for lev in ${levels}; do
         ldepth=$(echo "${lev} - ${above}" | bc)                    # Current layer depth
         ${cdo} minc,${lev} -subc,${above} ${atmGridID}_soil_depth.keep sdepth.tmp  # Soil depth in layer
         ${cdo} ifthenelse -gtc,0 sdepth.tmp \
               -mul init_moist_rel.tmp -mul ${atmGridID}_soil_field_cap.keep sdepth.tmp \
               zero.nc   ${var}_${lev}.tmp1
         ncecat -u soillev ${var}_${lev}.tmp1 ${var}_${lev}.tmp2    # Add new record dimension "soillev"
         ncecat ${var}_${lev}.tmp2 ${var}_${lev}.tmp3               # Convert soilev to non-record dim.
         ncwa -a record ${var}_${lev}.tmp3 ${var}_${lev}.tmp4       # Remove degenerate record dim.
         ${cdo} setlevel,${lev} ${var}_${lev}.tmp4 ${var}_${lev}.tmp  # Set level
         above=${lev}
      done
      $rm -f ${atmGridID}_${var}.nc
      ${cdo} -setattribute,${var}@long_name="Soil water content" \
          -setattribute,${var}@units="m" -setname,${var} \
          -setattribute,soillev@long_name="Soil layer (lower boundary)" \
          -setattribute,soillev@units="m" \
          -setmissval,-9.e+33 \
          -ifthen land.nc \
          -merge ${var}_?????.tmp ${atmGridID}_${var}.nc
      ;;
  # ----------------------
    snow )
  # ----------------------
      # Limit to 20 cm (water equivalent)
      ${cdo} -setattribute,${var}@long_name="Snow depth" -setattribute,${var}@units="m water equivalent" \
          -setname,${var} -setmissval,-9.e+33 -ifthen land.nc -mul non-glac-land-mask.nc -minc,0.2 -mulc,0.001 \
          -seldate,1111-${ini_mon}-11 -selvar,W_SNOW ${extpar_file} \
          ${atmGridID}_${var}.nc
      ;;
  # ----------------------
    surf_temp )
  # ----------------------
      ${cdo} -setname,${var} -setmissval,-9.e+33 -selvar,T_2M_CLIM ${extpar_file} \
          ${atmGridID}_${var}.nc
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

# Merge newly calculated variables to the ic soil file
# -------------
[[ -f ${file_ic_soil}.tmp ]] && ${rm} ${file_ic_soil}.tmp
for var in ${var_list}; do
  if [[ ! -f ${file_ic_soil}.tmp ]]; then
    # first variable
    ${cp} ${atmGridID}_${var}.nc  ${file_ic_soil}.tmp
  else
    # following variable
    mv ${file_ic_soil}.tmp ${file_ic_soil}.tmp2
    ${cdo} -O merge ${file_ic_soil}.tmp2  ${atmGridID}_${var}.nc ${file_ic_soil}.tmp
  fi
done

${cdo} --no_history setattribute,history="${history_att}" ${file_ic_soil}.tmp \
       ${path_bc}/${file_ic_soil}.nc

# Clean up
if [[ ${clean_up} == true ]]; then
  ${rm} *.tmp*
fi
echo "${prog}:     ${file_ic_soil}.nc          done"

exit 0
