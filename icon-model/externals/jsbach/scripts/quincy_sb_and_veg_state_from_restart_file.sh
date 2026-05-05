#!/usr/bin/env bash

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

#----------------------------------------------------------------------
# Collects those QUINCY vegetation and soil biogeochemistry variables from the restart file that define the
# state of the vegetation and of the soil biogeochemistry and stores these in an independent netcdf file.
#
# This is particularly required to inform a coupled simulation with a spun-up equilibrium state
# as achieved from a standalone ICON-Land with QUINCY simulation.
#
# Notes: 
# - The collected variables should represent the "long term" state of the vegetation and soil biogeochemistry.
# - They should be in sync with the variables read in mo_veg_init and mo_veg_sb (for each pft), respectively.
#
# Inspired by cpool_file_from_restart_files.sh

module unload cdo
module unload nco
module load nco/5.0.6-gcc-11.2.0
module load cdo/2.5.0-gcc-11.2.0

cdo="cdo -s -P 8"
#----------------------------------------------------------------------
# paths and filenames

input_path="./"
input_file="restart-file.nc"
output_path="./"
output_file="sb_and_veg_state_file.nc"

# usecase dependent number of pfts
nr_pft=13

# if true: include according product pool bgcms (Note: if "true" then the restart file is expected to contain these variables)
use_agriculture_with_prod_pools=true
use_sylviculture_with_prod_pools=true

# the name of the extra dimension in the icon output
name_of_soildim=layers_5

# selected soil model (simple_1d or jsm)
soilModel=simple_1d

# Assert that the output file does not already exist
if [ -f  ${output_path}/${output_file} ]; then
  echo "File ${output_path}/${output_file} already exists. Please delete or rename to proceed"
  exit
fi

#----------------------------------------------------------------------
#----------------------------------------------------------------------
# To derive the variable names as used in the restart file different lists / variables are created
# - for the bgcms
#   - all required elements
#   - all required compartments (veg and sb, respectively)
#   - bgcm main prefix (veg and sb, respectively)
#   - if required: product pool bgcms
# - for 'normal' variables
#   - process prefix (currently veg or sb)
#   - variable names
# - for all variables
#   - tile name postfix (currently only pfts)

#------------------
#-- elements
elements="C N P C13 C14 N15"

#------------------
#-- vegetation
# bgcm
veg_bgcm_prefix="veg_veg_bgcm_pool"
veg_compartments="leaf  fine_root  coarse_root  sap_wood  heart_wood  labile  reserve  fruit"

# other state variables (according to https://gitlab.dkrz.de/jsbach/jsbach/-/issues/294#note_341494, date 24.07.25)
veg_process_prefix="veg"
veg_variable_list="dens_ind  root_fraction_sl  mean_leaf_age
                   leaf2root_troot_mavg  unit_npp_troot_mavg  unit_uptake_n_pot_troot_mavg  unit_uptake_p_pot_troot_mavg
                   net_growth_tvegdyn_mavg  lai_tvegdyn_mavg"

# the seed bed is now located in an own bgcm
veg_bgcm_seedbed_prefix="veg_bgcm_seed_bed_pool"

# product pool bgcms - if the restart file contains product pool variables and these should be used
if [ ${use_agriculture_with_prod_pools} = true ] || [ ${use_sylviculture_with_prod_pools} = true ]; then
  veg_product_pool_bgcms=""
  if [ ${use_agriculture_with_prod_pools} = true ]; then
    veg_product_pool_bgcms="${veg_product_pool_bgcms} veg_bgcm_pp_crop"
  fi
  if [ ${use_sylviculture_with_prod_pools} = true ]; then
    veg_product_pool_bgcms="${veg_product_pool_bgcms} veg_bgcm_pp_fuel veg_bgcm_pp_paper veg_bgcm_pp_fiberboard"
    veg_product_pool_bgcms="${veg_product_pool_bgcms} veg_bgcm_pp_oirw veg_bgcm_pp_pv veg_bgcm_pp_sawnwood"
  fi
fi

#------------------
#-- soil
# bgcm
soil_bgcm_prefix="sb_sb_bgcm_pool"
soil_compartments="dom dom_assoc  sol_litter  pol_litter  woo_litter  mycorrhiza  microbial  residue  residue_assoc"

# other state variables (according to https://gitlab.dkrz.de/jsbach/jsbach/-/issues/294#note_341494, date 24.07.25)
sb_process_prefix="sb"
sb_variable_list="microbial_cue_eff_tmic_mavg  microbial_nue_eff_tmic_mavg  microbial_pue_eff_tmic_mavg"
sb_variable_list="${sb_variable_list} nh4_assoc po4_primary po4_assoc_slow"
if [[ "${soilModel}" == "simple_1d" ]]; then
  # currently no additional variables required
  sb_variable_list="${sb_variable_list}"
elif [[ "${soilModel}" == "jsm" ]]; then
  sb_variable_list="${sb_variable_list} enzyme_frac_poly_c_mavg enzyme_frac_poly_n_mavg enzyme_frac_poly_p_mavg"
else
  echo "Unexpected soil model scheme: ${soilModel}, please check!"
  exit 1
fi

#------------------
#-- tiles: all pfts
tile_list=""
i=1
while [ ${i} -le ${nr_pft} ]; do
  if [ ${i} -lt 10 ]; then
    tile_list="${tile_list} pft0${i}"
  else
    tile_list="${tile_list} pft${i}"
  fi
  (( i = i + 1 ))
done

#----------------------------------------------------------------------
#----------------------------------------------------------------------
# Generate list of to be extracted variables
var_list=""

# veg bgcm
for element in ${elements}; do
  for compartment in ${veg_compartments} ; do
    for tile in ${tile_list} ; do
      var_list="${var_list},${veg_bgcm_prefix}_${compartment}_${element}_${tile}"
    done
  done
done

# veg seedbed pool bgcms
for element in ${elements}; do
  for tile in ${tile_list} ; do
    var_list="${var_list},${veg_process_prefix}_${veg_bgcm_seedbed_prefix}_${element}_${tile}"
  done
done

# veg product pool bgcms
if [ ${use_agriculture_with_prod_pools} = true ] || [ ${use_sylviculture_with_prod_pools} = true ]; then
  for product_pool_bgcm in ${veg_product_pool_bgcms}; do
    for element in ${elements}; do
      for tile in ${tile_list} ; do
        var_list="${var_list},${veg_process_prefix}_${product_pool_bgcm}_${element}_${tile}"
      done
    done
  done
fi

# veg variables
for var in ${veg_variable_list} ; do
  for tile in ${tile_list} ; do
    var_list="${var_list},${veg_process_prefix}_${var}_${tile}"
  done
done

# sb bgcm
for element in ${elements}; do
  for compartment in ${soil_compartments} ; do
    for tile in ${tile_list} ; do
      var_list="${var_list},${soil_bgcm_prefix}_${compartment}_${element}_${tile}"
    done
  done
done

# sb variables
for var in ${sb_variable_list} ; do
  for tile in ${tile_list} ; do
    var_list="${var_list},${sb_process_prefix}_${var}_${tile}"
  done
done

#----------------------------------------------------------------------
#----------------------------------------------------------------------
# extract the variables from restart file
${cdo} --no_history -b 64 --reduce_dim -selvar${var_list} ${input_path}/${input_file} ${output_path}/${output_file}
#----------------------------------------------------------------------
# ... rename dimension to expected name
which ncatted || {
. ${MODULESHOME}/init/bash
  module load nco
}
ncrename -h -d cells,cell ${output_path}/${output_file}
ncrename -d ${name_of_soildim},soillev ${output_path}/${output_file}

#------------------------------------------------------------------------------
# set global attributes
ncatted -h -O -a ,global,d,, ${output_path}/${output_file}
ncatted -h -O -a comment,global,o,c,"States extracted from restart file ${input_file}" ${output_path}/${output_file}