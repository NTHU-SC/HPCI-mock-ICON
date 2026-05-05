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
#     ICON-Land orographic boundary conditions file
#-----------------------------------------------------------------------------
#
# The script is part of a script series to generate jsbach initial files
#   started by master script "create_jsbach_ini_files.sh"
#
#   bc_land_oro.nc
#   - Surface depression data prepared by T. Stacke
#   - Remaining orographic data from extpar

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
file_bc_sso=bc_land_sso
clean_up=true  # Remove intermediate files

cd ${work_dir}
# ----------------------------------------------------------------------------
#  Surface orography file
# -----------------------------
echo "${prog}: Generating ${file_bc_sso}:"

# The order of variables should not be changed
var_list="elevation oromea orostd orosig orogam orothe surf_depr_depth surf_depr_fract"

typeset -A sso_dict
sso_dict[elevation]=topography_c
sso_dict[oromea]=topography_c
sso_dict[orostd]=SSO_STDH
sso_dict[orosig]=SSO_SIGMA
sso_dict[orogam]=SSO_GAMMA
sso_dict[orothe]=SSO_THETA

# Preparation needed outside the parallel loop
if [[ ${var_list} == *"surf_depr_"* ]]; then
  ${cdo} splitvar ${pool_prepare}/05/depr_stor_cop_05.nc  05_
fi

for var in ${var_list}; do
  {
  echo "   ${var} ..."
  case ${var} in
  # ----------------------
    surf_depr_depth | surf_depr_fract )
  # ----------------------
      # Extrapolation
      ${cdo} -setmisstodis -ifthen 05_non-glac-land.nc -invertlat 05_${var}.nc 05_${var}.tmp
      # Remapping to the ICON grid; set glacier values to zero and ocean values to missing
      ${cdo} -setvar,${var} -mul -remap,${icon_grid},wgt05.nc 05_${var}.tmp \
                                 non-glac-land-mask.nc   ${atmGridID}_${var}.nc
      ;;
  # ----------------------
    elevation | oromea | orostd | orosig | orogam | orothe )
  # ----------------------
      ${cdo} setvar,${var} \
           -selvar,${sso_dict[${var}]} ${extpar_file} \
           ${atmGridID}_${var}.nc
      ;;
  # ----------------------
    * )
  # ----------------------
      echo "ERROR: Variable ${var} not known to be a bc_land_sso file variable."
      exit 1
      ;;
  esac
  } &  # comment out '&' for serial processing
done
wait

# Merge newly calculated variables to the bc orography file
# -------------
[[ -f ${file_bc_sso}.tmp ]] && ${rm} ${file_bc_sso}.tmp
for var in ${var_list}; do
  if [[ ! -f ${file_bc_sso}.tmp ]]; then
    # first variable
    ${cp} ${atmGridID}_${var}.nc  ${file_bc_sso}.tmp
  else
    # following variable
    mv ${file_bc_sso}.tmp ${file_bc_sso}.tmp2
    ${cdo} -O merge ${file_bc_sso}.tmp2  ${atmGridID}_${var}.nc ${file_bc_sso}.tmp
  fi
done

${cdo} --no_history setattribute,history="${history_att}" ${file_bc_sso}.tmp \
       ${path_bc}/${file_bc_sso}.nc

# Clean up
if [[ ${clean_up} == true ]]; then
  ${rm} *.tmp*
fi
echo "${prog}:     ${file_bc_sso}.nc          done"

exit 0
