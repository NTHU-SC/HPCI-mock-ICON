#!/bin/bash
#
# only used for jsbach
# split evergreen and deciduous in tropical and extratropical and shrubs in raingreen and deciduous 
# following the 1991-2020 Koppen-Geiger map
#   the split is conducted as done for tropical and extratropical in Georgievski & Hagemann (2019) 
#   (https://link.springer.com/article/10.1007/s00704-018-2675-2)
#   "Tropical WHERE (kgclim==1 .OR. kgclim==2 .OR. kgclim==3 .OR.kgclim==4 .OR. kgclim==6) else extra-tropical"
#----------------------------------------------------------------

# Extract arguments
prepared_esacci_file=$1 # in and output file!
jsb_KG_map=$2
output_dir_name=$3
temporary_dir=$4

# set cdo command to be silent and create no history
cdo="cdo -s --no_history -b 64"

#=============================================================================================================
# --------------- main
#=============================================================================================================
E_source="E"  # var name of evergreen pft
TE_pft="TE" # tropical evergreen
ETE_pft="ETE" # extra-tropical evergreen

D_source="D"  # var name of deciduous pft
TD_pft="TD" # tropical deciduous
ETD_pft="ETD" # extra-tropical deciduous

S_source="S"  # var name of shrub pft
RS_pft="RShrubs" # raingreen shrubs
DS_pft="DShrubs" # deciduous shrubs

# get a help file with zeros only
${cdo} -mulc,0 -selvar,${E_source} ${prepared_esacci_file} ${temporary_dir}/tmp_zeros.nc

# change source names to tropical evergreen / deciduous / raingreen
${cdo} -chname,${E_source},${TE_pft} -chname,${D_source},${TD_pft} -chname,${S_source},${RS_pft} \
          ${prepared_esacci_file} ${temporary_dir}/tmp_in_file_renamed_pfts.nc

# tropical vs extra-tropical following the 1991-2020 Koppen-Geiger map
kg_class=${temporary_dir}/tmp_kg_class.nc
# extend the KG map such that land sea-mask differences can be accounted for
${cdo} -setmisstonn -setctomiss,0 -selvar,kg_class ${jsb_KG_map} ${kg_class}

# derive tropical mask as done in Georgievski & Hagemann (2019):
#   "Tropical WHERE (kgclim==1 .OR. kgclim==2 .OR. kgclim==3 .OR.kgclim==4 .OR. kgclim==6) else extra-tropical"
${cdo} -lec,6 ${kg_class} ${temporary_dir}/tmp_trop_mask_1.nc
${cdo} -eqc,5 ${kg_class} ${temporary_dir}/tmp_trop_mask_2.nc
${cdo} -ifthenelse  ${temporary_dir}/tmp_trop_mask_2.nc \
          ${temporary_dir}/tmp_zeros.nc ${temporary_dir}/tmp_trop_mask_1.nc ${temporary_dir}/tmp_trop_mask.nc
# and determine tropical share of evergreen and deciduous and raingreen share of shrubs
${cdo} -ifthenelse ${temporary_dir}/tmp_trop_mask.nc \
  -selvar,${TE_pft} ${temporary_dir}/tmp_in_file_renamed_pfts.nc ${temporary_dir}/tmp_zeros.nc ${temporary_dir}/tmp_TE.nc
${cdo} -ifthenelse ${temporary_dir}/tmp_trop_mask.nc \
  -selvar,${TD_pft} ${temporary_dir}/tmp_in_file_renamed_pfts.nc ${temporary_dir}/tmp_zeros.nc ${temporary_dir}/tmp_TD.nc
${cdo} -ifthenelse ${temporary_dir}/tmp_trop_mask.nc \
  -selvar,${RS_pft} ${temporary_dir}/tmp_in_file_renamed_pfts.nc ${temporary_dir}/tmp_zeros.nc ${temporary_dir}/tmp_RS.nc

# derive extra tropical
${cdo} -addc,1 -mulc,-1 -setmisstoc,0 ${temporary_dir}/tmp_trop_mask.nc ${temporary_dir}/tmp_extratrop_mask.nc
# and determine share of evergreen and deciduous
${cdo} -setname,${ETE_pft} -ifthenelse ${temporary_dir}/tmp_extratrop_mask.nc \
  -selvar,${TE_pft} ${temporary_dir}/tmp_in_file_renamed_pfts.nc ${temporary_dir}/tmp_zeros.nc ${temporary_dir}/tmp_ETE.nc
${cdo} -setname,${ETD_pft} -ifthenelse ${temporary_dir}/tmp_extratrop_mask.nc \
  -selvar,${TD_pft} ${temporary_dir}/tmp_in_file_renamed_pfts.nc ${temporary_dir}/tmp_zeros.nc ${temporary_dir}/tmp_ETD.nc
${cdo} -setname,${DS_pft} -ifthenelse ${temporary_dir}/tmp_extratrop_mask.nc \
  -selvar,${RS_pft} ${temporary_dir}/tmp_in_file_renamed_pfts.nc ${temporary_dir}/tmp_zeros.nc ${temporary_dir}/tmp_DS.nc

# replace tropical pfts
${cdo} -merge ${temporary_dir}/tmp_TE.nc ${temporary_dir}/tmp_TD.nc ${temporary_dir}/tmp_RS.nc ${temporary_dir}/tmp_new_tropical_pft_values.nc
${cdo} -replace ${temporary_dir}/tmp_in_file_renamed_pfts.nc \
  ${temporary_dir}/tmp_new_tropical_pft_values.nc ${temporary_dir}/tmp_in_file_updated_pft_values.nc

# add extra tropical pfts
mv ${prepared_esacci_file} ${temporary_dir}/tmp_prep_esacci_file_pre_split
${cdo} -merge ${temporary_dir}/tmp_in_file_updated_pft_values.nc \
  ${temporary_dir}/tmp_ETE.nc ${temporary_dir}/tmp_ETD.nc ${temporary_dir}/tmp_DS.nc ${prepared_esacci_file}

# assert new = old
# ${cdo} -sub -selvar,${E_source} ${in_file} -add ${temporary_dir}/tmp_TE.nc ${temporary_dir}/tmp_ETE.nc ${temporary_dir}/tmp_assert_zero_E.nc
# ${cdo} -sub -selvar,${D_source} ${in_file} -add ${temporary_dir}/tmp_TD.nc ${temporary_dir}/tmp_ETD.nc ${temporary_dir}/tmp_assert_zero_D.nc
echo ">> Tropical and extra-tropical pfts splitted"

