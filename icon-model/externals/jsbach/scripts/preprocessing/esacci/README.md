Derive a netcdf file (alternatively for JSBACH or for QUINCY) providing natural PFT distributions 
derived from ESACCI data for use in the icon-land init file scripts
(https://gitlab.dkrz.de/jsbach/jsbach/-/blob/land4icon-mpim/scripts/preprocessing/initial_files/create_icon-land_ini_files.sh)

## STEP1: map ESACCI to a set of specified PFTs
- The ESACCI tool can be downloaded after registering here: http://maps.elie.ucl.ac.be/CCI/viewer/download.php#usertool
- Note: the tool needs a certain java version (which was not available on levante nodes when testing, but e.g. on volhynia (BGC)).
- ESACCI maps can be downloaded (e.g. via wget) after registering on https://cds.climate.copernicus.eu 
  - see https://cds.climate.copernicus.eu/cdsapp#!/dataset/satellite-land-cover?tab=form
  - currently used: year 2022 and 300m
  - Note: during the meantime one requires a ECMWF account.
- The ESACCI tool requires a csv crosswalking (cw) table, a target directory, and a target grid.
- We use cw tables based on Harper et al. (2023), Table 4
  - but modified according to expert knowledge (Soenke + MPIM-group for jsbachs veg ratio max)
  - Note: having different target pfts, the cw tables differ for JSBACH and QUINCY
- Example call of ESACCI tool for QUINCY (you can find the QUINCY cross walking table in this repository next to this README):
```bash
../cf_ESA-CCI/lc-user-tools-4.3/bin/aggregate-map.sh -PgridName=GEOGRAPHIC_LAT_LON -PnumRows=720 -PuserPFTConversionTable="./cw_IQ_11-for-13_based-on-Tab4_20250829.csv" -PtargetDir="./IQ_11-for-13_lonlat025/" -PnumMajorityClasses=12 -PoutputLCCSClasses=false ../cf_ESA-CCI_032025/C3S-LC-L4-LCCS-Map-300m-P1Y-2022-v2.1.1.nc > out_esacci_2022.out 2>&1 &
```
- Example call of ESACCI tool for JSBACH (you can find the JSBACH cross walking table in this repository next to this README):
```bash
../cf_ESA-CCI/lc-user-tools-4.3/bin/aggregate-map.sh -PgridName=GEOGRAPHIC_LAT_LON -PnumRows=720 -PuserPFTConversionTable="./cw_jsb_based-on-Tab4_with_TabC1_bare_20250907.csv" -PtargetDir="./jsb_lonlat025/" -PnumMajorityClasses=11 -PoutputLCCSClasses=false ../cf_ESA-CCI_032025/C3S-LC-L4-LCCS-Map-300m-P1Y-2022-v2.1.1.nc > out_esacci_for_jsb_2022.out 2>&1 &
```
Note: the ESACCI tool throws a couple of warnings and then takes very very very long without any notifications!

## STEP2: derive the esacci file processed for jsbach/quincy for use in the ICON-Land init file creation.

#### how to run
- Login to levante.
- Check options in the config.txt file (model!, paths, ESACCI year, ...) and set your output directory (output_dir_name=./OUTPUTPATH).
- SET the path to your scripts in main.sh script (scriptDir=YOURSCRIPTSPATH/preprocessing/esacci/).
- If necessary change the account (#SBATCH --account=mj0143) in main.sh script.
- Have required climate files available.
- Run: sbatch ${scriptDir}/main.sh - slurm will create an o file for you.

#### what happens in the main routine
* There are several bash scripts called from the main routine on the way from the ESACCI tool output to the preprocessed file:
  * calculate_climatologies: derives climatological data required to split ESACCI PFTs into phenologically defined subtypes
  * split_into_c3_and_c4: splits c3 and c4 grasses according to avg temperature
  * split_broadleaved_pfts (only required for QUINCY which distinguishes phenologically defined forest types)
    * evergreen broadleaf is divided into rain- and xeric forest following thresholds for mean temperature and precipitation
    * deciduous broadleaf is divided into rain green and summer green following thresholds values for mean temperature
  * split_jsbach_woody_pfts (only required for JSBACH: splits tropical and extra-tropical forests and raingreen and deciduous shrubs)
  * further scripts / operations handling and formatting the data
    * remapping to target resolution (required to ensure that all files share the same orientation of latitudes)
    * setting file and variable attributes
    * ...

Note: only natural pfts (including bare land fraction) are used from ESACCI, crop and pasture will be replaced by those from LUH.
Note: for JSBACH the veg. ratio max (vrm) is derived from the bare fraction, since - as opposed to QUINCY - JSBACH currently still
      does not explicitly account for bare land but only implicitly using the vrm (which might change in the future).