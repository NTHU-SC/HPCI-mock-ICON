This folder contains scripts to pre-process land-use harmonisation (**LUH**, https://luh.umd.edu/data.shtml) data
- for use in the icon land init file preprocessing
- for use with quincy in icon land
Currently the scripts are called (submitted to the cluster) independently of each other.
 
#### aggregate_LUH_states.sh
Script to aggregate LUH type states for later use in the icon land init file preprocessing together with
esacci land cover maps.
**The script should not be called interactively but should be submitted with sbatch.**
The script can be setup in the header section (target/original grid, years, in- and output paths and files).
Please particularly check/set the in- and output paths (i.e. replace YOUROUTPUTFOLDER).
In this script the 13 LUH states (primf,primn,secdf,...,range,...,c4ann,c4per)
are aggregated to 5 different broader land cover types (nat,pastr,c3crops,c4crops,urban).
Rangelands are thereby treated depending on a variable called *fstnf* from the static data provided by LUH
(If fstnf=1, which means its a grid-cell with forests in LUH, there is a conversion to pasture;
if fstnf=0, i.e. a grid-cell without forests in LUH, no conversion happens and the fraction stays natural).
The script produces one file for each year in the demanded year range.
Note: lats in LUH data have another orientation as compared to the one used in the cdo grid "global_0.25".
Therefore, the data is remapped even if staying on 0.25deg, to ensure the same orientation in all files!

#### aggregate_LUH_land_use.sh
Script to aggregate different types of LUH land use data (here harvest and fertilisation data)
to be read in ICON-Land with QUINCY for the wood harvest and agriculture processes of QUINCY.
**The script should not be called interactively but should be submitted with sbatch.**
Please particularly check/set the in- and output paths (i.e. replace YOUROUTPUTFOLDER).
- In this script the LUH harvest types (primf_harv,primn_harv,secmf_harv,secyf_harv,secnf_harv,pltns_harv)
from the LUH transitions file are aggregated to one common harvest map.
- The fertilisation data is derived separately for c3crops and c4crops using the management and the state file.
The script produces one file for each year in the demanded year range.

#### calculate_mean_LUH_harvest.sh
Script to calculate mean LUH type harvest, to be read in ICON-Land for the QUINCY wood harvest process.
(Required for spin-up and S0 type simulations)
**The script should not be called interactively but should be submitted with sbatch.**
The script can be setup in the header section (grid, years, in- and output paths and files).
Please particularly check/set the in- and output paths (i.e. replace YOUROUTPUTFOLDER).
In this script the LUH harvest types (primf_harv,primn_harv,secmf_harv,secyf_harv,secnf_harv,pltns_harv)
are aggregated to one common harvest map.
The script produces one file with the mean of the demanded year range.

#### generate_static_LUH.sh
Script to prepare static land-use input files for the ICON-Land model with QUINCY.
It performs interpolation of crop types and forest harvest slash fractions to the selected grid 
and produces a merged NetCDF file for use in ICON-Land with QUINCY simulations.
**The script should not be called interactively but should be submitted with sbatch.**
The script can be setup in the header section (grid, in- and output paths and files).
Please particularly check/set the in- and output paths (i.e. replace YOUROUTPUTFOLDER).
This script interpolates LUH land-use input data (crop types and forest harvest slash fractions) to the ICON grid (R2B4 or R2B5).
The resulting variables are merged into a single NetCDF file used as static land-use input for ICON-Land with QUINCY simulations.
