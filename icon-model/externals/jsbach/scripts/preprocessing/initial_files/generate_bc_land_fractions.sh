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
#  Generate
#     ICON-Land land cover fractions
#
# The script is part of a script series to generate ICON-Land initial files
# started by master script "create_jsbach_ini_files.sh".
#-----------------------------------------------------------------------------

# Variables exported from master script create_icon-land_ini_files.sh
extpar_file=${extpar_file}       # Extpar file name

icon_grid=${icon_grid}
atmGridID=${atmGridID}
remap_scheme=${remap_scheme}

path_bc=${bc_file_dir}
year_list=${year_list}
work_dir=${work_dir}
pool_prepare=${pool_prepare}

cdo=${cdo}
rm=$rm
cp=$cp

# Local variables
prog=$(basename $0)
file_bc_frac=bc_land_frac
clean_up=true  # Remove intermediate files

min_fract=${min_fract} # minimum land grid cell fraction other than 0.

# Definition of small_fract
#  With jsbach, all PFT tiles must have a minimum fraction of small_fract.
#  With quincy, PFT tiles fractions smaller than small_fract will be set to zero.
small_fract_jsb=1.e-10
small_fract_qcy=1.e-06

# Input files for cover fractions: natural (potential) vegetation and land use states
luh2_states_file=${pool_prepare}/T255/LUH2v2h_states_T255_all-oceans_no-dynveg.nc # without .gz
luh3_states_root=${pool_prepare}/025/LUH3/r0001/LUH_states   # without ${year}.nc
prep_esacci_jsb=${pool_prepare}/025/esacci/r0001/jsbach_prep_ESACCI2022.nc
prep_esacci_qcy=${pool_prepare}/025/esacci/r0001/quincy_prep_ESACCI2022.nc

prep_025_wgts_file=${prep_esacci_jsb}

start_year=$(echo ${year_list} | cut -f1 -d" ")
year=${start_year}

cd ${work_dir}
# ----------------------------------------------------------------------------
#  Land cover fractions
# -----------------------------
echo "${prog}: Generating land cover fraction files"

# ----------------------------------------------------------------------------
#  1. bc fract without PFT fractions
# -----------------------------
# The order of variables should not be changed
var_list="notsea fract_glac sea fract_lake fract_land fract_veg veg_ratio_max land lake glac"

# Preparations that need to be done outside the parallel variables loop.
#--------------

# Lake fraction
#---------------
# We use lake fractions from the FLake global lake database. The Caspian Sea is
# not counted as lake in this data set. Also other smaller water bodies in that region
# are missing. So lake fractions of Caspian Sea/Aral Lake area are replaced with
# Globcover data (extpar variable LU_CLASS_FRACTION, class 21).

# Select extpar field of lake fraction (from flake)
${cdo} selvar,FR_LAKE ${extpar_file} extpar_FR_LAKE.nc
# Select extpar field of lake fraction from Globcover
${cdo} -setlevel,0 -setgrid,${icon_grid} -sellevel,21 -selvar,LU_CLASS_FRACTION ${extpar_file} \
    extpar_LU21.nc
#   => lake fractions for Caspian Sea and Aral Lake area, missing value elsewhere
${cdo} masklonlatbox,46.,62.,36.,48. extpar_LU21.nc extpar_LU_CLASS_FRACTION_lev21_casp.nc
# mask: 1 for Caspian Sea and Aral Lake, missing value elsewhere
${cdo} gtc,-1. extpar_LU_CLASS_FRACTION_lev21_casp.nc extpar_LU_CLASS_FRACTION_mask_casp.nc
# mask: 0 for Caspian Sea and Aral Lake area, 1 elsewhere
${cdo} subc,1. -setmisstoc,2. extpar_LU_CLASS_FRACTION_mask_casp.nc \
    extpar_LU_CLASS_FRACTION_mask_mul_casp.nc
# lake fraction (from flake) excluding lakes in the Caspian Sea and Aral Lake area
${cdo} mul extpar_FR_LAKE.nc extpar_LU_CLASS_FRACTION_mask_mul_casp.nc \
    extpar_fr_lake_clean.nc
# add Caspian Sea and Aral Lake from Globcover to lake fractions from flake elsewhere
${cdo} add extpar_fr_lake_clean.nc -setmisstoc,0. extpar_LU_CLASS_FRACTION_lev21_casp.nc \
    extpar_fr_lake_casp.nc

# ICON-Land grid cells with ocean fraction cannot have a lake fraction. Thus lake
# fractions of coastal cells need to be set to zero.

# Define a mask of land grid boxes without ocean
${cdo} gec,1. notsea.nc extpar_fr_land_mask.nc
# Lake fraction just for grid boxes with no ocean ==> this is the new lake fraction
${cdo} -setattribute,${var}@long_name="Fraction of lake tile rel. to land tile" \
    -chname,FR_LAKE,fract_lake -mul extpar_fr_lake_casp.nc extpar_fr_land_mask.nc \
    ${atmGridID}_fract_lake.nc

for var in ${var_list}; do
  {
  echo "   ${var} ..."
  case ${var} in
  # ----------------------
    notsea )
  # ----------------------
      # Calculated above in preparations section
      ${cdo} -setattribute,${var}@long_name="Land fraction of the grid cell (incl. lakes)" \
          -setattribute,${var}@units="1" notsea.nc ${atmGridID}_notsea.nc
      ;;
  # ----------------------
    sea )
  # ----------------------
      ${cdo} -setattribute,${var}@long_name="Ocean fraction of the grid cell" \
          -setvar,sea -setrtoc,-1,0,0 -mulc,-1 -subc,1 notsea.nc ${atmGridID}_sea.nc
      ;;
  # ----------------------
    fract_glac )
  # ----------------------
      # Glacier fraction relative to the land tile
      ${cdo} -setattribute,${var}@long_name="Fraction of glacier tile rel. to land tile" \
          -setattribute,${var}@units="1" -setvar,fract_glac glac.nc ${atmGridID}_fract_glac.nc
      ;;
  # ----------------------
    glac )
  # ----------------------
      # The variables without "fract_" are considered relative to the grid cell. For the glacier fraction
      # this actually does not make a difference.
      ${cdo} -setattribute,${var}@long_name="Glacier fraction of the grid cell" \
          glac.nc ${atmGridID}_glac.nc
      ;;
  # ----------------------
    fract_veg )
  # ----------------------
      # Vegetation fraction relative to the land fraction, i.e. all non-glacier cells
      ${cdo} -setattribute,${var}@long_name="Fraction of vegetated tile rel. to land tile" \
          -setvar,fract_veg -ltc,0.5 glac.nc ${atmGridID}_fract_veg.nc
      ;;
  # ----------------------
    fract_lake )
  # ----------------------
      # Already generated above, outside this loop
      ;;
  # ----------------------
    lake )
  # ----------------------
      # The variables without "fract_" are considered relative to the grid box. For the lake fraction
      # this actually does not make a difference.
      ${cdo} -setattribute,${var}@long_name="Lake fraction of the grid cell" \
          -chname,fract_lake,lake ${atmGridID}_fract_lake.nc ${atmGridID}_lake.nc
      ;;
  # ----------------------
    fract_land )
  # ----------------------
      # Land fraction 'fract_land' - relative to the box tile, i.e. the part of the surface handled
      # by ICON-Land (i.e. 'notsea'). It is the fraction that is not lake and is 1 over the ocean.
      ${cdo} -setattribute,${var}@long_name="Fraction of land tile rel. to box tile" \
          -chname,fract_lake,fract_land -addc,1. -mulc,-1. ${atmGridID}_fract_lake.nc \
          ${atmGridID}_fract_land.nc
      ;;
  # ----------------------
    land )
  # ----------------------
      # The variables without "fract_" are considered relative to the grid cell. In contrast to
      # fract_land, land is zero for ocean grid cells.
      ${cdo} -setattribute,${var}@long_name="Land fraction to the grid cell" \
          -chname,notsea,land -sub notsea.nc ${atmGridID}_fract_lake.nc ${atmGridID}_land.nc
      ;;
  # ----------------------
     veg_ratio_max )
  # ----------------------
      # Note: veg_ratio_max is needed as input variable only with jsbach, not with quincy.
      # But initial soil moisture calculations are based on (the jsbach) veg_ratio_max also
      # with quincy (compare generate_ic_land_soil.sh). We thus need to calculate it here
      # for both setups.
      # The variable is based on esacci data (1 - bare; following the cross walking table
      # used for jsbach).

      # Remapping matrix for entire grid - without missing values
      if [[ ! -f wgt025.nc ]]; then
        ${cdo} gen${remap_scheme},${icon_grid} -setmisstoc,0. ${prep_025_wgts_file} wgt025.nc
      fi
      ${cdo} mul -remap,${icon_grid},wgt025.nc -setmisstodis ${prep_esacci_jsb} \
                 non-glac-land-mask.nc    ${atmGridID}_prep_esacci_12pfts.nc

      ${cdo} selvar,bare  ${atmGridID}_prep_esacci_12pfts.nc ${atmGridID}_esacci_jsb_bare.nc

      # Veg_ratio_max: 1 - bare; minimum vegetated fraction is "small_fract"
      ${cdo} -setvar,veg_ratio_max -mul non-glac-land-mask.nc -maxc,${small_fract_jsb} \
          -addc,1. -mulc,-1. ${atmGridID}_esacci_jsb_bare.nc  ${atmGridID}_veg_ratio_max.nc
      ${rm} ${atmGridID}_esacci_jsb_bare.nc

      # Keep veg_ratio_max for usage below in 12 PFT setup and for late use in
      # generate_ic_land_soil.nc
      cp ${atmGridID}_veg_ratio_max.nc ${atmGridID}_veg_ratio_max.keep
      ;;
  # ----------------------
    * )
  # ----------------------
      echo "ERROR: Variable ${var} not known to be a bc_land_fract file variable."
      exit 1
      ;;
  esac
  } &  # comment out '&' for serial processing
done
wait

# Merge newly calculated variables to the bc fraction file
# -------------
[[ -f ${file_bc_frac}.nc ]] && ${rm} ${file_bc_frac}.nc
for var in ${var_list}; do
  if [[ ! -f ${file_bc_frac}.nc ]]; then
    # first variable
    ${cp} ${atmGridID}_${var}.nc  ${file_bc_frac}.nc
  else
    # following variable
    mv ${file_bc_frac}.nc ${file_bc_frac}.tmp
    ${cdo} -O merge ${file_bc_frac}.tmp  ${atmGridID}_${var}.nc ${file_bc_frac}.nc
  fi
done

${cdo} --no_history setattribute,history="${history_att}" ${file_bc_frac}.nc \
       ${path_bc}/${file_bc_frac}.nc

# Consistency checks
${cdo} expr,"test=notsea+sea-1."            ${file_bc_frac}.nc  zero1.nc
${cdo} expr,"test=land+lake+sea-1."         ${file_bc_frac}.nc  zero2.nc
${cdo} expr,"test=fract_glac-glac"          ${file_bc_frac}.nc  zero3.nc
${cdo} expr,"test=fract_lake-lake"          ${file_bc_frac}.nc  zero4.nc
${cdo} expr,"test=fract_lake+fract_land-1." ${file_bc_frac}.nc  zero5.nc
${cdo} expr,"test=fract_veg+fract_glac-1."  ${file_bc_frac}.nc  zero6.nc
for test in zero1 zero2 zero3 zero4 zero5 zero6; do
  if   [[ $(cdo output -fldmin ${test}.nc | tr -d ' ') != 0 ]] \
    || [[ $(cdo output -fldmax ${test}.nc | tr -d ' ') != 0 ]]; then
    echo "ERROR: Consistency check failed for test ${test}.nc: "
    ${cdo} infon ${test}.nc
    exit 1
  fi
done

# Clean up
if [[ ${clean_up} == true ]]; then
  for var in ${var_list}; do; ${rm} -f ${atmGridID}_${var}.nc; done
  ${rm} -f extpar_*.nc
  ${rm} *.tmp*
fi

echo "${prog}:     ${file_bc_frac}.nc             done"

# ----------------------------------------------------------------------------
#  1b. bc fract with PFT fractions
# -----------------------------
npfts_list="12 13"

# Currently, only PFT fractions depend on the year
for year in ${year_list}; do
  for npfts in $npfts_list; do

    typeset -Z2 ilct ipft

    #
    # Generate PFT fractions - depending on the uscase
    #

    case $npfts in
    # ----------------------------------------------------------------------------
      12 | 13 )   # 12 PFTs - jsbach; 13 PFTs - quincy setup
    # ----------------------------------------------------------------------------
        [[ ${npfts} == 12 ]] && land_setup=jsbach
        [[ ${npfts} == 13 ]] && land_setup=quincy

        if [[ ${land_setup} == jsbach ]]; then
          pft_tag="${npfts}pfts"
          nat_pfts="pft01 pft02 pft03 pft04 pft05 pft06 pft07 pft08"
          ant_pfts="pft09 pft10 pft11 pft12"
        elif [[ ${land_setup} == quincy ]]; then
          pft_tag="iq_${npfts}pfts"
          # As pft13 has a natural (bare) and an anthropogenic (urban) fraction we add here an additional
          # preliminary PFT: pft14 for urban area. The urban will be added to the bare fraction below.
          nat_pfts="pft01 pft02 pft03 pft04 pft05 pft06 pft07 pft08 pft13"
          ant_pfts="pft09 pft10 pft11 pft12 pft14"
        fi
        pft_list="${nat_pfts} ${ant_pfts}"

        # Set variables depending on land setup
        [[ ${land_setup} == jsbach ]] && small_fract=${small_fract_jsb} && prep_esacci_file=${prep_esacci_jsb}
        [[ ${land_setup} == quincy ]] && small_fract=${small_fract_qcy} && prep_esacci_file=${prep_esacci_qcy}

        # Remapping matrix for entire grid - without missing values - for LUH states and natural vegetation
        if [[ ! -f wgt025.nc ]]; then
          ${cdo} gen${remap_scheme},${icon_grid} -setmisstoc,0. ${prep_025_wgts_file} wgt025.nc
        fi

        # Land use states
        #-----------------------

        # The sum of the LUH state variables (c3crops c4crops nat pastr urban) equals one for all
        # non-glacier land cells. Bare land is part of natural.
        if [[ ! -f ${atmGridID}_LUH3_states_${year}.nc ]]; then
          ${cdo} mul -remap,${icon_grid},wgt025.nc -setmisstodis ${luh3_states_root}_${year}.nc \
                     non-glac-land-mask.nc  ${atmGridID}_LUH3_states_${year}.nc

          ${cdo} -selvar,c3crops ${atmGridID}_LUH3_states_${year}.nc LUH3_c3crops_${year}.nc
          ${cdo} -selvar,c4crops ${atmGridID}_LUH3_states_${year}.nc LUH3_c4crops_${year}.nc
          ${cdo} -selvar,pastr   ${atmGridID}_LUH3_states_${year}.nc LUH3_pastr_${year}.nc
          ${cdo} -selvar,urban   ${atmGridID}_LUH3_states_${year}.nc LUH3_urban_${year}.nc
          ${cdo} -selvar,nat     ${atmGridID}_LUH3_states_${year}.nc LUH3_nat_${year}.nc

          # With jsbach we also need the maximum possible fractions of anthropogenic types, i.e. as
          # if there was no natural land cover.
          ${cdo} maxc,${small_fract} -expr,"anthro=1.-nat" ${atmGridID}_LUH3_states_${year}.nc LUH3_anthro_${year}.nc
          ${cdo} -div LUH3_c3crops_${year}.nc LUH3_anthro_${year}.nc LUH3_c3crops_max_${year}.nc
          ${cdo} -div LUH3_c4crops_${year}.nc LUH3_anthro_${year}.nc LUH3_c4crops_max_${year}.nc
          ${cdo} -div LUH3_pastr_${year}.nc   LUH3_anthro_${year}.nc LUH3_pastr_max_${year}.nc
          ${cdo} -div LUH3_urban_${year}.nc   LUH3_anthro_${year}.nc LUH3_urban_max_${year}.nc
        fi

        # Potential vegetation
        #-----------------------

        # Potential (natural) vegetation needs to be derived only once - in the first year.
        # Note: The sum of the natural vegetation types (including bare) equals one for all non-glacier
        #       land cells.
        if [[ ${year} == ${start_year} ]]; then

          if [[ ! -f ${atmGridID}_prep_esacci_${pft_tag}.nc ]]; then
            ${cdo} mul -remap,${icon_grid},wgt025.nc -setmisstodis ${prep_esacci_file} \
                   non-glac-land-mask.nc    ${atmGridID}_prep_esacci_${pft_tag}.nc
          fi

          ${cdo} splitvar  ${atmGridID}_prep_esacci_${pft_tag}.nc ${atmGridID}_esacci_${pft_tag}_
          ${rm} ${atmGridID}_prep_esacci_${pft_tag}.nc

          # C3/C4 ratios   (H: C3 grass; HC4: C4 grass)
          # Needed to define C3/C4 ratio of pastures corresponding to the C3/C4 ration of grasslands
          # Note: If neither C3 nor C4 grasses exist, we do an extrapolation.
          ${cdo} -setvar,c4_ratio -mul non-glac-land-mask.nc -setmisstodis \
              -div ${atmGridID}_esacci_${pft_tag}_HC4.nc \
                   -add ${atmGridID}_esacci_${pft_tag}_H.nc ${atmGridID}_esacci_${pft_tag}_HC4.nc  c4_ratio_${pft_tag}.nc
          ${cdo} -setvar,c3_ratio -addc,1. -mulc,-1. c4_ratio_${pft_tag}.nc  c3_ratio_${pft_tag}.nc

          #  Natural PFTs
          #------------------
          for pft in ${nat_pfts}; do
            {
            #---------------------------------------
            if [[ ${land_setup} == "jsbach" ]]; then
            #---------------------------------------
              case $pft in
                pft01 ) # Tropical evergreen trees
                  $cdo setname,fract_${pft} ${atmGridID}_esacci_${pft_tag}_TE.nc ${atmGridID}_${pft}.${npfts}.keep
                  ncatted -a long_name,fract_${pft},c,c,'Fraction of tropical evergreen tree tile rel. to veg tile' \
                      ${atmGridID}_${pft}.${npfts}.keep
                  ;;
                pft02 ) # Tropical deciduous trees
                  $cdo setname,fract_${pft} ${atmGridID}_esacci_${pft_tag}_TD.nc ${atmGridID}_${pft}.${npfts}.keep
                  ncatted -a long_name,fract_${pft},c,c,'Fraction of tropical deciduous tree tile rel. to veg tile' \
                      ${atmGridID}_${pft}.${npfts}.keep
                  ;;
                pft03 ) # Temperate broadl. evergreen + evergreen_conifer
                  $cdo setname,fract_${pft} ${atmGridID}_esacci_${pft_tag}_ETE.nc ${atmGridID}_${pft}.${npfts}.keep
                  ncatted -a long_name,fract_${pft},c,c,'Fraction of temperate evergreen tree tile rel. to veg tile' \
                      ${atmGridID}_${pft}.${npfts}.keep
                  ;;
                pft04 ) # Temperate broadl. deciduous + deciduous conifer
                  $cdo setname,fract_${pft} ${atmGridID}_esacci_${pft_tag}_ETD.nc ${atmGridID}_${pft}.${npfts}.keep
                  ncatted -a long_name,fract_${pft},c,c,'Fraction of temperate deciduous tree tile rel. to veg tile' \
                      ${atmGridID}_${pft}.${npfts}.keep
                  ;;
                pft05 ) # Raingreen shrub
                  $cdo setname,fract_${pft} ${atmGridID}_esacci_${pft_tag}_RShrubs.nc ${atmGridID}_${pft}.${npfts}.keep
                  ncatted -a long_name,fract_${pft},c,c,'Fraction of raingreen shrub tile rel. to veg tile' \
                      ${atmGridID}_${pft}.${npfts}.keep
                  ;;
                pft06 ) # Deciduous shrub
                  $cdo setname,fract_${pft} ${atmGridID}_esacci_${pft_tag}_DShrubs.nc ${atmGridID}_${pft}.${npfts}.keep
                  ncatted -a long_name,fract_${pft},c,c,'Fraction of deciduous shrub tile rel. to veg tile' \
                      ${atmGridID}_${pft}.${npfts}.keep
                  ;;
                pft07 ) # C3-grass
                  $cdo setname,fract_${pft} ${atmGridID}_esacci_${pft_tag}_H.nc ${atmGridID}_${pft}.${npfts}.keep
                  ncatted -a long_name,fract_${pft},c,c,'Fraction of C3-grass tile rel. to veg tile' \
                      ${atmGridID}_${pft}.${npfts}.keep
                  ;;
                pft08 ) # C4-grass
                  $cdo setname,fract_${pft} ${atmGridID}_esacci_${pft_tag}_HC4.nc ${atmGridID}_${pft}.${npfts}.keep
                  ncatted -a long_name,fract_${pft},c,c,'Fraction of C4-grass tile rel. to veg tile' \
                      ${atmGridID}_${pft}.${npfts}.keep
                  ;;
                * )
                  echo "ERROR: pft ${pft} not expected in ${land_setup} setup with ${npfts} PFTs."; exit 1
                  ;;
              esac
            #---------------------------------------
            elif [[ ${land_setup} == "quincy" ]]; then
            #---------------------------------------
              case $pft in
                pft01 ) # Tropical broadleaf evergreen trees
                  $cdo setname,fract_${pft} ${atmGridID}_esacci_${pft_tag}_TrBE.nc ${atmGridID}_${pft}.${npfts}.keep
                  ncatted -a long_name,fract_${pft},c,c,'Fraction of moist broadleaved evergreen rel. to veg tile' \
                          -a short_name,fract_${pft},c,c,'BEM' ${atmGridID}_${pft}.${npfts}.keep
                  ;;
                pft02 ) # Temperate broadleaf evergreen trees
                  $cdo setname,fract_${pft} ${atmGridID}_esacci_${pft_tag}_TeBE.nc ${atmGridID}_${pft}.${npfts}.keep
                  ncatted -a long_name,fract_${pft},c,c,'Fraction of dry broadleaved evergreen rel. to veg tile' \
                          -a short_name,fract_${pft},c,c,'BED' ${atmGridID}_${pft}.${npfts}.keep
                  ;;
                pft03 ) # Tropical broadleaf raingreen trees
                  $cdo setname,fract_${pft} ${atmGridID}_esacci_${pft_tag}_TrBR.nc ${atmGridID}_${pft}.${npfts}.keep
                  ncatted -a long_name,fract_${pft},c,c,'Fraction of rain green broadleaved deciduous rel. to veg tile' \
                          -a short_name,fract_${pft},c,c,'BDR'   ${atmGridID}_${pft}.${npfts}.keep
                  ;;
                pft04 ) # Temperate broadleaf deciduous trees
                  $cdo setname,fract_${pft} ${atmGridID}_esacci_${pft_tag}_TeBS.nc ${atmGridID}_${pft}.${npfts}.keep
                  ncatted -a long_name,fract_${pft},c,c,'Fraction of summer green broadleaved deciduous rel. to veg tile' \
                          -a short_name,fract_${pft},c,c,'BDS'   ${atmGridID}_${pft}.${npfts}.keep
                  ;;
                pft05 ) # Needleleaf evergreen trees
                  $cdo setname,fract_${pft} ${atmGridID}_esacci_${pft_tag}_NEEV.nc ${atmGridID}_${pft}.${npfts}.keep
                  ncatted -a long_name,fract_${pft},c,c,'Fraction of needle-leaved evergreen rel. to veg tile' \
                          -a short_name,fract_${pft},c,c,'NE'  ${atmGridID}_${pft}.${npfts}.keep
                  ;;
                pft06 ) # Needleleaf deciduous trees
                  $cdo setname,fract_${pft} ${atmGridID}_esacci_${pft_tag}_NEDE.nc ${atmGridID}_${pft}.${npfts}.keep
                  ncatted -a long_name,fract_${pft},c,c,'Fraction of summer green needle-leaved rel. to veg tile' \
                          -a short_name,fract_${pft},c,c,'NS' ${atmGridID}_${pft}.${npfts}.keep
                  ;;
                pft07 ) # C3-grass
                  $cdo setname,fract_${pft} ${atmGridID}_esacci_${pft_tag}_H.nc ${atmGridID}_${pft}.${npfts}.keep
                  ncatted -a long_name,fract_${pft},c,c,'Fraction of C3 grass tile rel. to veg tile' \
                          -a short_name,fract_${pft},c,c,'TeH' ${atmGridID}_${pft}.${npfts}.keep
                  ;;
                pft08 ) # C4-grass
                  $cdo setname,fract_${pft} ${atmGridID}_esacci_${pft_tag}_HC4.nc ${atmGridID}_${pft}.${npfts}.keep
                  ncatted -a long_name,fract_${pft},c,c,'Fraction of C4 grass tile rel. to veg tile' \
                          -a short_name,fract_${pft},c,c,'TrH' ${atmGridID}_${pft}.${npfts}.keep
                  ;;
                pft13 ) # bare  (preliminary; urban will be added below - after scaling)
                  if [[ ${land_setup} != quincy ]]; then
                    echo "ERROR: pft ${pft} not expected in ${land_setup} setup."; exit 1
                  fi
                  $cdo setname,fract_${pft} ${atmGridID}_esacci_${pft_tag}_bare.nc ${atmGridID}_${pft}.${npfts}.keep
                  ncatted -a long_name,fract_${pft},c,c,'Fraction of bare tile rel. to veg tile' \
                          -a short_name,fract_${pft},c,c,'BSO' ${atmGridID}_${pft}.${npfts}.keep
                  ;;
                * )
                  echo "ERROR: pft ${pft} not expected in ${land_setup} setup with ${npfts} PFTs."; exit 1
                  ;;
              esac
            else
              echo "ERROR: Unknown land_setup: ${land_setup}."; exit 1
            fi
            } &  # comment out '&' for serial processing
          done
          wait
        fi  # first year


        #  Anthropogenic PFTs
        #---------------------
        for pft in ${ant_pfts}; do
          {
          case $pft in
            pft09 ) # C3-pasture
              $cdo -setname,fract_${pft} -mul c3_ratio_${pft_tag}.nc LUH3_pastr_${year}.nc ${atmGridID}_${pft}.${npfts}.keep
              ncatted -a long_name,fract_${pft},c,c,'Fraction of C3 pasture tile rel. to veg tile' \
                      -a short_name,fract_${pft},c,c,'TeP' ${atmGridID}_${pft}.${npfts}.keep
              ;;
            pft10 ) # C4-pasture
              $cdo setname,fract_${pft} -mul c4_ratio_${pft_tag}.nc LUH3_pastr_${year}.nc ${atmGridID}_${pft}.${npfts}.keep
              ncatted -a long_name,fract_${pft},c,c,'Fraction of C4 pasture tile rel. to veg tile' \
                      -a short_name,fract_${pft},c,c,'TrP' ${atmGridID}_${pft}.${npfts}.keep
              ;;
            pft11 ) # C3-crop
              $cdo setname,fract_${pft} LUH3_c3crops_${year}.nc ${atmGridID}_${pft}.${npfts}.keep
              ncatted -a long_name,fract_${pft},c,c,'Fraction of C3-crop tile rel. to veg tile' \
                      -a short_name,fract_${pft},c,c,'TeC' ${atmGridID}_${pft}.${npfts}.keep
              ;;
            pft12 ) # C4-crop
              $cdo setname,fract_${pft} LUH3_c4crops_${year}.nc ${atmGridID}_${pft}.${npfts}.keep
              ncatted -a long_name,fract_${pft},c,c,'Fraction of C4-crop tile rel. to veg tile' \
                      -a short_name,fract_${pft},c,c,'TrC' ${atmGridID}_${pft}.${npfts}.keep
              ;;
            pft14 ) # urban (preliminary; will be added to bare below)
              if [[ ${land_setup} != quincy ]]; then
                echo "ERROR: pft ${pft} not expected in ${land_setup} setup."; exit 1
              fi
              $cdo setname,fract_${pft} LUH3_urban_${year}.nc ${atmGridID}_${pft}.${npfts}.keep
              ncatted -a long_name,fract_${pft},c,c,'Fraction of urban tile rel. to veg tile' \
                      -a short_name,fract_${pft},c,c,'UAR' ${atmGridID}_${pft}.${npfts}.keep
              ;;
            * )
              echo "ERROR: pft ${pft} not expected in ${land_setup} setup with ${npfts} PFTs."; exit 1
              ;;
          esac
          } &  # comment out '&' for serial processing
        done
        wait
        ;;
    # ----------------------------------------------------
      * )
    # ----------------------------------------------------
        echo "ERROR: No setup for ${npfts} PFTs available."; exit 1
        ;;
    esac

    # ----------------------------------------------------------------------------------------------
    #  Scaling: The sum of all PFT fractions must be one for all cells
    # ---------
    #
    # As the sum of all (natural) PFTs in the prep_esacci files is one, adding anthropogenic fractions
    # from the LUH data leads to total land cover fractions greater than one. Thus scaling is needed:
    #
    #  - Anthropogenic fractions are kept as defined in the LUH data file.
    #  - Natural PFT fractions are reduced accordingly, while the ratio between them is kept. We do not
    #    apply a pasture rule or such.
    #  - Bare fraction
    #     - In jsbach, PFT fractions are defined relative to the vegetated grid cell fraction of the
    #       veg tile ("1-veg_ratio_max") which is constant for all years. There is no explicit bare PFT.
    #     - In quincy, the veg tile also comprises a bare PFT fraction. Here, the bare PFT not only
    #       represents natural bare land but also urban land fractions.
    #  - Minimum fraction
    #     - With jsbach a minimum PFT fraction of small_fract is needed on all non-glacier land cells,
    #       otherwise NLCC (formerly dynveg) cannot be used.
    #     - With quincy, fractions can be zero, but very small fractions (between small_fract and zero)
    #       are not allowed.
    #  - In the current jsbach and quincy usecases the PFT fractions in bc_land_frac are defined relative
    #    to the 'veg' tile fraction, which is 1 everywhere, except for glacier cells - also for ocean
    #    cells. We thus need to assign the ocean area to one of the PFTs. We use pft07 (C3 grass) with
    #    jsbach and the bare tile with quincy, but this has no physical relevance.
    #
    # Due to the different handling of bare fractions and small_fract scaling works differently
    # for the different PFT setups.
    # ----------------------------------------------------------------------------------------------

    if [[ ${land_setup} == jsbach ]]; then

      # A vegetated fraction smaller than ntiles*small_fract cannot be accounted as vegetated, as all PFT tiles
      # need a minimum fraction of small_fract.
      sf=$(echo ${small_fract} | awk '{printf "%.13f\n", $1}') # Convert exponential number to float (needed with bc)
      min_veg=$(echo "${npfts} * ${sf}" | bc)
      ${cdo} -setname,has_veg -gtc,${min_veg} ${atmGridID}_veg_ratio_max.keep  has_veg.nc

      # Scaling with veg_ratio_max (i.e. 1 - bare)
      for pft in ${nat_pfts}; do
        {
        ${cdo} -div ${atmGridID}_${pft}.${npfts}.keep \
                    ${atmGridID}_veg_ratio_max.keep ${atmGridID}_${pft}.tmp1
        } &  # comment out '&' for serial processing
      done
      wait

      # We need to limit the anthropogenic fractions, otherwise fractions become huge in cells with little vegetation.
      ${cdo} mul LUH3_pastr_max_${year}.nc c3_ratio_${pft_tag}.nc LUH3_c3pastr_max_${year}.nc
      ${cdo} mul LUH3_pastr_max_${year}.nc c4_ratio_${pft_tag}.nc LUH3_c4pastr_max_${year}.nc
      ${cdo} -min -div ${atmGridID}_pft09.${npfts}.keep ${atmGridID}_veg_ratio_max.keep \
                  LUH3_c3pastr_max_${year}.nc   ${atmGridID}_pft09.tmp1
      ${cdo} -min -div ${atmGridID}_pft10.${npfts}.keep ${atmGridID}_veg_ratio_max.keep \
                  LUH3_c4pastr_max_${year}.nc   ${atmGridID}_pft10.tmp1
      ${cdo} -min -div ${atmGridID}_pft11.${npfts}.keep ${atmGridID}_veg_ratio_max.keep \
                  LUH3_c3crops_max_${year}.nc   ${atmGridID}_pft11.tmp1
      ${cdo} -min -div ${atmGridID}_pft12.${npfts}.keep ${atmGridID}_veg_ratio_max.keep \
                  LUH3_c4crops_max_${year}.nc   ${atmGridID}_pft12.tmp1

      # Scaling: The sum of the natural fractions from escci need to be scaled to the natural 'nat' fraction of LUH3
      ${rm} -f sum_ant_luh.nc sum_nat_luh.nc
      ${cdo} enssum ${atmGridID}_pft09.tmp1 ${atmGridID}_pft1[0-2].tmp1  sum_ant_luh.nc
      ${cdo} -maxc,${small_fract} -addc,1. -mulc,-1.                     sum_ant_luh.nc  sum_nat_luh.nc # 1-sum_ant_luh

      # Mask for grid cells with a relevant natural vegetation according to the LUH natural vegetation fraction
      # and according to veg_ratio_max (based on esacci)
      ${cdo} -setname,has_nat_luh -mul has_veg.nc -gtc,${min_veg} sum_nat_luh.nc has_nat_luh.nc

      ${cp} sum_nat_luh.nc scaling.nc

      # Apply the scaling to all natural PFTs
      for pft in ${nat_pfts}; do
        {
        ${cdo} -mul -mul ${atmGridID}_${pft}.tmp1 scaling.nc has_nat_luh.nc ${atmGridID}_${pft}.tmp2
        } &
      done
      for pft in ${ant_pfts}; do
        {
        cp ${atmGridID}_${pft}.tmp1  ${atmGridID}_${pft}.tmp2
        } &
      done
      wait
      # Now, the sum of all PFT fractions is between 0 and 1. We now need to handle small_fract.

      # Find out relevant and irrelevant PFTs
      ${rm} -f pft??_relevant.nc pft??_irrelevant.nc msk_relevant_pft??.nc
      for pft in ${pft_list}; do
        {
        ${cdo} -gtc,${small_fract} ${atmGridID}_${pft}.tmp2 msk_relevant_${pft}.nc
        ${cdo} setmisstoc,0. -ifthen    msk_relevant_${pft}.nc ${atmGridID}_${pft}.tmp2 ${pft}_relevant.nc
        ${cdo} setmisstoc,0. -ifnotthen msk_relevant_${pft}.nc ${atmGridID}_${pft}.tmp2 ${pft}_irrelevant.nc
        } &
      done
      wait

      # Sum up relevant natural / anthropogenic PFT fractions
      ${rm} -f sum_nat_relevant.nc sum_ant_relevant.nc num_irrelevant.nc
      ${cdo} enssum pft0[1-8]_relevant.nc                            sum_nat_relevant.nc
      ${cdo} enssum pft09_relevant.nc pft1?_relevant.nc              sum_ant_relevant.nc
      ${cdo} -mulc,-1. -subc,${npfts} -enssum msk_relevant_pft??.nc  num_irrelevant.nc

      # Scaling of relevant natural PFTs
      # With jsbach the total fraction of relevant natural PFTs should be
      #    sum_nat_relevant = 1. - sum_ant_relevant - num_irrelevant*small_fract
      ${cdo} -addc,1. -mulc,-1. -add sum_ant_relevant.nc \
                                     -mulc,${small_fract} num_irrelevant.nc  sum_nat_relevant_new.nc

      # Scaling factor for the relevant natural PFTs
      ${cdo} -setmisstoc,0. -div sum_nat_relevant_new.nc sum_nat_relevant.nc scaling_nat.nc

      # Apply the scaling to the relevant natural PFTs
      for pft in ${nat_pfts}; do
        {
        ${cdo} setmisstoc,${small_fract} -ifthen has_nat_luh.nc \
                  -mul ${atmGridID}_${pft}.tmp2 scaling_nat.nc \
               ${atmGridID}_${pft}.tmp3
        } &
      done
      wait

      # In case there is no relevant natural but anthropogenic vegetation, the anthropogenic PFTs are scaled.
      # Target: sum_ant_relevant = 1. - num_irrelevant*small_fract
      ${cdo} addc,1. -mulc,-1. -mulc,${small_fract} num_irrelevant.nc sum_ant_relevant_new.nc
      ${cdo} setname,only_ant -sub has_veg.nc has_nat_luh.nc only_ant.nc
      ${cdo} -ifthenelse only_ant.nc -div sum_ant_relevant_new.nc sum_ant_relevant.nc has_veg.nc scaling_ant.nc
      for pft in ${ant_pfts}; do
        {
        ${cdo} -mul ${atmGridID}_${pft}.tmp2 scaling_ant.nc ${atmGridID}_${pft}.tmp3
        } &
      done
      wait

      # In case there is no relevant vegetation at all, we define c3 or c4 grasses depending on c3/c4 ratios.
      #    fract_gras = 1. - num_irrelevant*small_fract
      grass_max=$(echo "1. - (${npfts} -1) * ${sf}" | bc)
      ${cdo} -ifthenelse has_veg.nc ${atmGridID}_pft07.tmp3 \
                                    -mul c3_ratio_${pft_tag}.nc -addc,${grass_max} zero.nc  ${atmGridID}_pft07.tmp4
      mv ${atmGridID}_pft07.tmp4 ${atmGridID}_pft07.tmp3
      ${cdo} -ifthenelse has_veg.nc ${atmGridID}_pft08.tmp3 \
                                    -mul c4_ratio_${pft_tag}.nc -addc,${grass_max} zero.nc  ${atmGridID}_pft08.tmp4
      mv ${atmGridID}_pft08.tmp4 ${atmGridID}_pft08.tmp3

      # Set minimum fraction for non-glacier land cells
      for pft in ${pft_list}; do
        {
        $cdo mul -maxc,${small_fract} ${atmGridID}_${pft}.tmp3 non-glac-land-mask.nc ${atmGridID}_${pft}.tmp4
        } &
      done
      wait

      # Final corrections in case number of irrelevant fractions changed
      ${cdo} -mul non-glac-land-mask.nc -subc,1 -enssum ${atmGridID}_pft??.tmp4 err.nc
      # Find out relevant PFTs
      for pft in ${pft_list}; do
        {
        ${cdo} -gtc,${min_veg} ${atmGridID}_${pft}.tmp4 msk_relevant_${pft}.nc
        } &
      done
      wait

      ${rm} -f num_relevant.nc
      ${cdo} enssum msk_relevant_pft??.nc num_relevant.nc
      for pft in ${pft_list}; do
        {
        ${cdo} ifthenelse msk_relevant_${pft}.nc \
               -sub ${atmGridID}_${pft}.tmp4 -div err.nc num_relevant.nc \
               ${atmGridID}_${pft}.tmp4   ${atmGridID}_${pft}.nc
        } &
      done
      wait

    elif [[ ${land_setup} == quincy ]]; then
      # With quincy we get the scaling factor from the LUH3 natural vegetation fraction 'nat'.
      ${cp} LUH3_nat_${year}.nc scaling.nc
      # Apply the scaling to all natural PFTs
      for pft in ${nat_pfts}; do
        {
        ${cdo} mul ${atmGridID}_${pft}.${npfts}.keep scaling.nc ${atmGridID}_${pft}.tmp
        } &
      done
      # No scaling of the anthropogenic PFTs
      for pft in ${ant_pfts}; do
        {
        ${cp} ${atmGridID}_${pft}.${npfts}.keep  ${atmGridID}_${pft}.tmp
        } &
      done
      wait

      # Set fractions to zero on non-glacier land
      for pft in ${pft_list}; do
        {
        $cdo mul ${atmGridID}_${pft}.tmp non-glac-land-mask.nc ${atmGridID}_${pft}.tmp1
        } &
      done
      wait

      # The sum of these PFTs (${atmGridID}_pft??.tmp1) is one on all non-glacier land cells.
      # Now, we have to remove all fractions below small_fract and rescale the relevant natural
      # vegetation accordingly.
      for pft in ${pft_list}; do
        {
        ${cdo} -gtc,${small_fract}  ${atmGridID}_${pft}.tmp1  msk_relevant_${pft}.nc
        ${cdo} setmisstoc,0. -ifthen    msk_relevant_${pft}.nc ${atmGridID}_${pft}.tmp1 ${pft}_relevant.nc
        ${cdo} setmisstoc,0. -ifnotthen msk_relevant_${pft}.nc ${atmGridID}_${pft}.tmp1 ${pft}_irrelevant.nc
        } &
      done
      wait

      # Scaling factor: (sum_nat_relevant + sum_irrelevant) / sum_nat_relevant
      ${rm} -f sum_nat_relevant.nc sum_ant_relevant.nc sum_irrelevant.nc
      ${cdo} enssum pft0[1-8]_relevant.nc pft13_relevant.nc  sum_nat_relevant.nc
      ${cdo} enssum pft09_relevant.nc pft1[0-2]_relevant.nc pft14_relevant.nc  sum_ant_relevant.nc
      ${cdo} enssum pft??_irrelevant.nc                      sum_irrelevant.nc
      ${cdo} -div -add sum_nat_relevant.nc sum_irrelevant.nc sum_nat_relevant.nc scaling_nat.nc
      ${cdo} -div -add sum_ant_relevant.nc sum_irrelevant.nc sum_ant_relevant.nc scaling_ant.nc

      # All PFTs: Set fractions below small_fract to zero. (no scaling)
      for pft in ${pft_list}; do
        {
        ${cdo} -setrtoc,0.,${small_fract},0. ${atmGridID}_${pft}.tmp1 ${atmGridID}_${pft}.tmp2
        } &
      done
      wait

      # Scaling of natural PFTs
      for pft in ${nat_pfts}; do
        {
        ${cdo} -ifthenelse -gtc,0. sum_nat_relevant.nc \
               -mul ${atmGridID}_${pft}.tmp2 scaling_nat.nc \
                ${atmGridID}_${pft}.tmp2     ${atmGridID}_${pft}.nc
        } &
      done
      wait

      # Anthropogenic PFTs: In rare cases, when no relevant natural vegetation exists, we need to scale
      #                     the relevant anthropogenic PFTs
      for pft in ${ant_pfts}; do
        {
        ${cdo} -ifthenelse -lec,0. sum_nat_relevant.nc \
               -mul ${atmGridID}_${pft}.tmp2 scaling_ant.nc \
                ${atmGridID}_${pft}.tmp2     ${atmGridID}_${pft}.nc
        } &
      done
      wait

      # With quincy add urban (pft14) to the bare (pft13) fraction.
      ${cdo} -add ${atmGridID}_pft13.nc ${atmGridID}_pft14.nc ${atmGridID}_pft13.tmp
      mv ${atmGridID}_pft13.tmp ${atmGridID}_pft13.nc
      ${rm} ${atmGridID}_pft14.nc  # otherwise test will fail
    fi

    # As the sum of the PFT fractions needs to be one, also entire glacier, lake and ocean cells
    # need to be assigned to a PFT.
    if [[ ${land_setup} == jsbach ]]; then
      # With jsbach all non-land grid cells are assigned to pft07.
      ${cdo} add ${atmGridID}_pft07.nc ocean-or-glac.nc pft07.tmp && mv pft07.tmp ${atmGridID}_pft07.nc
    elif [[ ${land_setup} == quincy ]]; then
      # With quincy all non-land grid cells are assigned to pft13 (bare).
      ${cdo} add ${atmGridID}_pft13.nc ocean-or-glac.nc pft13.tmp && mv pft13.tmp ${atmGridID}_pft13.nc
    fi

    # Test
    ${rm} -f sum_pft_${npfts}_${year}.nc
    ${cdo} enssum ${atmGridID}_pft??.nc sum_pft_${npfts}_${year}.nc
    if   [[ $(cdo output -fldmin sum_pft_${npfts}_${year}.nc | tr -d ' ') != 1 ]] \
      || [[ $(cdo output -fldmax sum_pft_${npfts}_${year}.nc | tr -d ' ') != 1 ]]; then
      echo "ERROR: Sum of PFT fraction differs from 1. (compare sum_pft_${npfts}_${year}.nc).:"
      exit 1
    fi

    # Generate bc fract file including PFT fractions
    # -------------
    # Merge PFT fractions
    [[ -f pft_fracts.nc ]] && ${rm} pft_fracts.nc
    ipft=01
    while [[ ${ipft} -le ${npfts} ]]; do
      pft=pft${ipft}
      if [[ ! -f pft_fracts.nc ]]; then
        # first variable
        ${cp} ${atmGridID}_${pft}.nc  pft_fracts.nc
      else
        # following variable
        mv pft_fracts.nc pft_fracts.tmp
        ${cdo} -O merge pft_fracts.tmp ${atmGridID}_${pft}.nc pft_fracts.nc
      fi
      (( ipft = ipft + 1 ))
    done

    if [[ ${npfts} == 12 ]]; then
      ${cp} ${file_bc_frac}.nc ${file_bc_frac}_${npfts}.nc
    elif [[ ${npfts} == 13 ]]; then
      # Veg_ratio_max is not needed with Quincy and is removed here to avoid confusion.
      ${cdo} delvar,veg_ratio_max ${file_bc_frac}.nc ${file_bc_frac}_${npfts}.nc
    fi

    # Merge PFT fractions into bc fractions file generated above
    [[ -f ${file_bc_frac}_${pft_tag}.nc ]] && ${rm} ${file_bc_frac}_${pft_tag}.nc
    ${cdo} -O merge ${file_bc_frac}_${npfts}.nc pft_fracts.nc ${file_bc_frac}_${pft_tag}.nc

    ${cdo} --no_history setattribute,history="${history_att}" ${file_bc_frac}_${pft_tag}.nc \
        ${path_bc}/${file_bc_frac}_${pft_tag}_${year}.nc

    # Add global attributes from LUH and esacci pre-preocessed files
    ncks -A -x -h ${prep_esacci_file}            ${path_bc}/${file_bc_frac}_${pft_tag}_${year}.nc
    ncks -A -x -h ${luh3_states_root}_${year}.nc ${path_bc}/${file_bc_frac}_${pft_tag}_${year}.nc

    echo "${prog}:     ${file_bc_frac}_${pft_tag}_${year}.nc done"

    # Clean up
    if [[ ${clean_up} == true ]]; then
      ${rm} *.tmp*
      ${rm} sum_ant_relevant.nc sum_nat_relevant.nc scaling.nc
      ${rm} -f num_irrelevant.nc sum_irrelevant.nc
      ${rm} ${atmGridID}_pft??.nc
      ${rm} -f ${atmGridID}_pft09.${npfts}.keep ${atmGridID}_pft1[0-2].${npfts}.keep ${atmGridID}_pft14.${npfts}.keep
      ${rm} -f msk_relevant_pft??.nc pft??_irrelevant.nc pft??_relevant.nc
      [[ ${year} == ${start_year} ]] && ${rm} -f lct0000??.nc
      ${rm} -f ${atmGridID}_LUH?_states_${year}.nc ${atmGridID}_esacci_${pft_tag}*.nc ${atmGridID}_bare_${npfts}.nc
    fi

  done   # npfts
done   # year

# Clean up
if [[ ${clean_up} == true ]]; then
  ${rm} -f ${atmGridID}_LUH2_states.nc
  ${rm} LUH?_*.nc ${atmGridID}_pft??.??.keep
  ${rm} -f c3_ratio_*.nc c4_ratio_*.nc ${atmGridID}_bare_??.keep \
            ${atmGridID}_veg_ratio_max.nc
fi

echo "${prog}:  done"

exit 0
