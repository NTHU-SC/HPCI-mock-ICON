#!/bin/bash
#SBATCH --partition=shared     # Use compute node for R2B9 and higher
# #SBATCH --partition=compute  # Use in case of too little memory (OOM error)
# #SBATCH --mem=450G           #  ""
#SBATCH --account=mh0287
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8    # Use 128 for R2B9 and higher
#SBATCH --output=create_icon-land_ini_files.o%j
#SBATCH --time=08:00:00
#-----------------------------------------------------------------------------
#
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

#-----------------------------------------------------------------------------
# Master script to generate ICON-Land initial files
#
# It makes use of the existing scripts:
#  1) generate_fractional_mask.sh
#     - Generation of a fractional land sea mask for coupled setups:
#       Remapping of the ocean mask to the atmosphere grid
#
#  2) extpar4jsbach_mpim_icon.sh
#     - run extpar to generate soil texture data as well as albedo, roughness
#       length, forest fraction and LAI and vegetation fraction climatologies.
#
#  3) bc_files_from_extpar.sh
#     Generate ic/bc files containing the new expar data:
#     - Lake mask replaced (-> bc_land_frac)
#     - Albedo, roughness length, forest fraction, lai_clim and veg_fract
#       replaced.
#     - Rooting depth (and maxmoist) replaced, additional variables:
#       FR_SAND, FR_SILT, FR_CLAY, SUB_FR_SAND, SUB_FR_SILT and
#       SUB_FR_CLAY
#
#  4) adapt_nwp_extpar_file.sh
#     The extpar data file read by the NWP atmosphere contains several
#     variables, that are also included in the bc_land files. For consistency,
#     these variables are replaced by the respective bc_land file variables.
#
#  5) calc_hd_receive_mask.sh
#     Coupled setups with external HD model need a mask file indicating ocean
#     inflow cells. This mask depends on the fractional land sea mask. For
#     that reason it is part of this script suite.
#
#  6) generate_hdpara_file.sh
#     If the HD model is run within jsbach ("internal HD") it needs several
#     parameters on the ICON-Land grid. The script clones, configures and
#     runs the software in https://github.com/ThomasRiddick/DynamicHD.git
#
# This approach is meant to be preliminary. It documents the current process
# of initial data generation. The aim is however, to generate all initial data
# from extpar in the not so far future.
#-----------------------------------------------------------------------------
#
# Notes on the extpar file variables
#
# In this script we use two different variables defining extpar files although
# you generally need just one extpar file:
# - extpar_file: The extpar file we generate here using extpar4jsbach_mpim_icon.sh
#      (unless it is already available from /pool/data/JSBACH/icon/extpar4jsbach/).
#      This extpar file contains e.g. soil textures and other soil parameters.
# - nwp_extpar_file: The extpar file read by the NWP atmosphere at runtime.
#      It will be adapted to the newly generated ic/bc files in
#      adapt_nwp_extpar_file.sh. Unless a better extpar file is provided by DWD
#      (e.g. with additional variables needed for the run) you should here also
#      use the "extpar_file" from above.
#
#-----------------------------------------------------------------------------
set -e

# Which sections of this script should be run?
#  It is possible to run one task after the other, which might be useful in case of problems

run_extpar4jsbach_mpim_icon=false  # Generally available in /pool/data/JSBACH/icon/extpar4jsbach/
run_generate_fractional_mask=false # Only needed for coupled setups with new grid combination
run_bc_files_from_extpar=true
run_adapt_nwp_extpar_file=true     # Only needed for configurations with NWP atmosphere
run_calc_hd_receive_mask=false     # Only needed for coupled setups with external HD model
run_generate_hdpara_file=true      # Needed with internal HD model

#-----------------------------------------------------------------------------
#
# Settings - also needed in sub-scripts
# -------------------------------------
# Variables which typically need to be adapted to the target grid(s) and directory
#                   and file paths (from here up to "### End of settings ###" line:
# - work_dir_base and scratch_dir_base
# - output_root_dir and revision
# - mpi_grid
# - coupled
# - atmGridID and refinement
# - oceGridID (if coupled=true)
# - start_year and end_year or year_list
# - extpar_dir and extpar_file
# - extpar_source_dir (if extpar_file not yet existing)
# - nwp_extpar_file (for simulations with NWP atmosphere)
# - minimum/maximum fractions for grid cells that are not completely ocean/land (if coupled=false)
#
# --------------------------------------
# Base directories for constructing work and output directories (see below)
work_dir_base=/work/mh0287/${USER}    # levante
scratch_dir_base=/scratch/m/$USER     # levante
#scratch_dir_base=/work/mh0287/${USER} # breeze4

# --------------------------------------
# ICON grids used

mpim_grid=true     # Grid directory structure of MPIM grids
#mpim_grid=false   # Grid directory structure of DWD grids; use with icon-xpp (except 0012-0035)

if [[ $mpim_grid == true ]]; then
  icon_grid_rootdir=/pool/data/ICON/grids/public/mpim
  oce_grid_rootdir=${icon_grid_rootdir}
else
  icon_grid_rootdir=/pool/data/icon-xpp
  #icon_grid_rootdir=/pool/data/ICON/grids/public/edzw
  oce_grid_rootdir=/pool/data/ICON/grids/public/mpim
fi

export atmGridID=0049
#export atmGridID=0012 # icon-xpp
export atmRes=R02B04
export icon_grid=${icon_grid_rootdir}/${atmGridID}/icon_grid_${atmGridID}_${atmRes}_G.nc

export coupled=false       # land sea mask for coupled experiment?
if [[ ${coupled} == true ]]; then
  export oceGridID=0035      # Required for coupled configurations
  export oceRes=R02B06
  # Ocean grid/mask: We need the mask on the global (*_G.nc) grid.
  #    Using the reduced ocean grid (*_O.nc) technically works but leads to wrong results!
  export icon_grid_oce=${oce_grid_rootdir}/${oceGridID}/icon_mask_${oceGridID}_${oceRes}_G.nc
  export grid_label=$atmGridID-$oceGridID
else
  export grid_label=$atmGridID
fi

# --------------------------------------
# Output and working directories
export revision=r00xx      # ic/bc revision directory that will be generated
export work_dir=${scratch_dir_base}/ini_files/work/${grid_label}/${revision}  # temporary working directory
export output_root_dir=${work_dir_base}/ini_files     # root directory for the new ic/bc files
export bc_file_dir=${output_root_dir}/${grid_label}/land/${revision}
hd_rev=""                  # revision for external HD data; "" for no revision subdirectory
export hd_file_dir=${output_root_dir}/${grid_label}/hd/${hd_rev} # directory for external HD data

# --------------------------------------
# List of years:
# a) First and last year for a series of bc_land_frac files, e.g. for historical simulations
start_year=1850
end_year=1850
year_list=""
for ((yr=${start_year}; yr<=${end_year}; yr++)); do
  export year_list="${year_list} $yr"
done
#---
# b) Define list of years
#export year_list="1850 1979 1990 1992 2005 2015"

# For years after 2015: chose scenario (not available for all resolutions)
export scenario="ssp119"   # ssp119 / ssp126 / ssp245 / ssp370 / ssp434 / ssp460 / ssp534os / ssp585

# --------------------------------------
# extpar4jsbach file ("extpar_file"; compare above "Notes on the extpar files")
#
today=$(date +%Y%m%d) # time stamp for generated extpar file; set to a fixed YYYYMMDD to re-use file from different date
if [[ ${run_extpar4jsbach_mpim_icon} == false ]]; then
  # Use an already existing file and don't run extpar
  # -----
  if [[ ${mpim_grid} == true ]]; then
    export extpar_dir=/pool/data/JSBACH/icon/extpar4jsbach/mpim
  else
    export extpar_dir=/pool/data/JSBACH/icon/extpar4jsbach/dwd
  fi
  export extpar_file=icon_extpar4jsbach_${atmGridID}_20251007_tiles.nc   # Check for grid ID and corresponding date tag
  #export extpar_file=icon_extpar4jsbach_${atmGridID}_20251007_tiles.nc  # icon-xpp
  # -----
  #export extpar_dir=${output_root_dir}                                  # Use extpar file from previous run
  #export extpar_file=icon_extpar4jsbach_${atmGridID}_${today}_tiles.nc  #   with run_extpar4jsbach_mpim_icon=true
else
  # Run extpar and generate a new extpar4jsbach file
  # -----
  export extpar_version=v5.15.1p1 # extpar version if compiling extpar from scratch
  # Directory for extpar source code: either existing, to use a precompiled installation, or
  # writable by the user to clone and compile extpar anew - parent directory needs to exist!
  export extpar_source_dir=/work/mh0287/m212005/extpar/extpar4jsbach # on levante or breeze4
  # Directory and file name for generated extpar file
  export extpar_dir=${output_root_dir}
  export extpar_file=icon_extpar4jsbach_${atmGridID}_${today}_tiles.nc

  # Directory for extpar input data
  # If not running on levante or breeze4, add a case for your machine here
  case $(hostname) in
    levante*|*.lvt.dkrz.de)
      export extpar_data_dir=/work/pd1167/extpar-input-data
      ;;
    breeze4)
      export extpar_data_dir=/work/mh0287/icon-preprocessing/extpar-input-data
      ;;
    *)
      echo ERRO: Unknown host: $(hostname)
      exit 1
      ;;
  esac
fi

# --------------------------------------
# NWP Extpar file:  Only used in setups with NWP atmosphere  - compare above "Notes on the extpar files"
export nwp_extpar_file=${extpar_dir}/${extpar_file}       # Use extpar file defind above (to be generated or existing)

# --------------------------------------
# Fractional mask file
if [[ ${coupled} == true ]]; then
  if [[ ${run_generate_fractional_mask} == true ]]; then
    export fractional_mask=${output_root_dir}/${grid_label}/fractional_mask/fractional_lsm_${atmGridID}_${oceGridID}.nc
  else
    # Needs to exist
    export fractional_mask=${icon_grid_rootdir}/${grid_label}/fractional_mask/fractional_lsm_${atmGridID}_${oceGridID}.nc
    # ATTENTION: some XPP grids use "-" instead of "_"!
    # export fractional_mask=${icon_grid_rootdir}/${grid_label}/fractional_mask/fractional_lsm_${atmGridID}-${oceGridID}.nc
    if [[ ! -f ${fractional_mask} ]]; then
      echo "$0: ERROR: fractional mask file ${fractional_mask} does not exist!"
      exit 1
    fi
  fi
else
  export fractional_mask=${extpar_file}          # Fractional mask in extpar file (based on observations)
  # export fractional_mask=${icon_grid}          # Leads to integer LSM!
fi

# --------------------------------------
# Shared parameters
#
# minimum/maximum fractions for grid cells that are not completely ocean/land
#   0.000001 / 0.999999: Currently used; seams rather suitable for coupled setups, where water conservation is an issue.
#   0.001 / 0.999: Used in coupled and uncoupled setups until August 2023
#   0.25 / 0.75: Corresponding to coupled setup with matching grids, e.g. R2B5/R2B6
#   0.5 /0.5: No fractional grid cells
if [[ ${coupled} == true ]]; then
  export min_fract=0.000001 # minimum land grid cell fraction if not complete ocean
  export max_fract=0.999999 # maximum land grid cell fraction if not complete land
else
  #--- AMIP setup corresponding to coupled setup with matching grids, e.g. R2B5/R2B6
  export min_fract=0.25     # minimum land grid cell fraction if not complete ocean
  export max_fract=0.75     # maximum land grid cell fraction if not complete land
  #--- AMIP setup with integer land sea mask
  #export min_fract=0.5      # minimum land grid cell fraction if not complete ocean
  #export max_fract=0.5      # maximum land grid cell fraction if not complete land
fi

# --------------------------------------
# Definitions needed with HD input file generation
#
# --- Definitions needed with run_calc_hd_receive_mask (external HD)
# (existing) HD parameter file on 0.5 grid
export hdpara05_file=/pool/data/ICON/edzw-shadow/indepedent/hd/input/05deg/hdpara_vs1_12.nc
# HD receive mask file (to be generated)
export hd_receive=${hd_file_dir}/hd_receive_${grid_label}.nc

# --- Definitions needed with run_generate_hdpara_file (internal HD)
# Depending on min- and max_fract, the fractional mask in bc_land_frac ('notsea') might
# not be identical to the mask in ${fractional_mask}. For a consistent setup, we thus use
# the bc_land_frac fractions also for hd parameter generation. The defaults set here refer
# to a setup where land bc/ic files are generated anew (run_bc_files_from_extpar=true).
export keep_hd_output_dir=false   # keep additional output from HD parameter generation
                                  # for debugging etc.
export hd_fractional_mask=${bc_file_dir}/bc_land_frac.nc
export hd_orography_file=${bc_file_dir}/bc_land_sso.nc    # needed with R02B07 and higher

### End of settings ###

#-----------------------------------------------------------------------------
# Some preparations

function finish {
  if [[ -d ${work_dir} ]]; then
    echo "Cleaning up work dir ${work_dir}"
    rm -rf ${work_dir} >& /dev/null
  fi
}
trap finish EXIT

# Name and directory of this script
if [[ $SLURM_JOB_NAME == "" || $SLURM_JOB_NAME == "interactive" ]]; then
  # Interactive job
  scripts_dir=$(dirname $0)
  script_name=./${0##*/}
else
  # Batch job
  scripts_dir=$SLURM_SUBMIT_DIR
  script_name=$SLURM_JOB_NAME
fi
cd ${scripts_dir}
export scripts_dir=$PWD

# Check if output directory exists and you are allowed to write
if [[ ! -d ${output_root_dir} ]]; then
  echo "Directory 'output_root_dir' for output needs to exist (now: ${output_root_dir}). "
  exit 1
fi
if [[ ! -w ${output_root_dir} ]]; then
  echo "No write permission for output directory ${output_root_dir}."
  exit 1
fi

# Create directory to save all scripts used for this initial file revision
[[ -d ${bc_file_dir}/scripts ]] || mkdir -p ${bc_file_dir}/scripts
if [[ ${run_generate_fractional_mask} == true ]]; then
  [[ -d ${fractional_mask%/*}/scripts ]] || mkdir -p ${fractional_mask%/*}/scripts
fi
if [[ ${run_calc_hd_receive_mask} == true  || ${run_generate_hdpara_file} == true ]]; then
  [[ -d ${hd_file_dir}/scripts ]] || mkdir -p ${hd_file_dir}/scripts
fi

# Load the necessary modules
# Note: The modules might get overwritten below in case expar4jsbach needs to be recompiled!
set +u
. ${MODULESHOME}/init/bash
module purge
module load git
module load julia
case $(hostname) in
  levante* | *.lvt.dkrz.de)
    module load python3/2023.01-gcc-11.2.0
    module load cdo/2.5.0-gcc-11.2.0
    module load nco/5.0.6-gcc-11.2.0
    module load netcdf-c/4.8.1-gcc-11.2.0
    module load gcc/11.2.0-gcc-11.2.0
    ;;
  breeze4)
    module load python/3.10.4
    module load cdo/2.5.0
    module load nco/5.0.1
    module load netcdf-c/4.9.0
    module load gcc/11.2.0
    ;;
  *)
    echo ERROR: Unknown host: $(hostname)
    exit 1
    ;;
esac
module list
set -u

export cdo="cdo -s -P 8"

#-----------------------------------------------------------------------------
# 1. Generate fractional land sea mask
if [[ ${run_generate_fractional_mask} == true ]]; then
  if [[ ${coupled} != true ]]; then
    echo " Skipping the generation of fractional masks as this is only needed in coupled setups."
  else
    ${scripts_dir}/generate_fractional_mask.sh
    cp ${scripts_dir}/generate_fractional_mask.sh ${fractional_mask%/*}/scripts
    cp ${scripts_dir}/${script_name}              ${fractional_mask%/*}/scripts
  fi
fi

#-----------------------------------------------------------------------------
# 2. Run extpar4jsbach
if [[ ${run_extpar4jsbach_mpim_icon} == true ]]; then

  [[ ! -d ${extpar_dir} ]] && [[ ${extpar_dir} != ${output_root_dir} ]] && mkdir -p ${extpar_dir}

  if [[ ! -d ${extpar_data_dir} ]]; then
    echo "Directory with extpar input data doesn't exist:"
    echo "    ${extpar_data_dir}"
    echo "Preferably, since this is a very large data volume, this should be a clone of"
    echo "    https://gitlab.dkrz.de/extpar-data/extpar-input-data"
    echo "already existing in your file system. If you need to download this yourself"
    echo "contact jonas.jucker@c2sm.ethz.ch for access to the git lsf repository."
    exit 1
  fi

  # Clone extpar source from repository if necessary
  if [[ ! -d ${extpar_source_dir} ]]; then
    [[ ! -d ${extpar_source_dir%/*} ]] && mkdir -p ${extpar_source_dir%/*}
    cd ${extpar_source_dir%/*}
    if ! git clone --recursive git@github.com:C2SM-RCM/extpar.git ${extpar_source_dir##*/} >&/dev/null ;then
      if ! git clone --recursive git@gitlab.dkrz.de:m212005/extpar4jsbach.git ${extpar_source_dir##*/} >&/dev/null ;then
        echo "No valid git repository to clone extpar"
        echo "Make sure that you have password-less SSH read access to the extpar repository at"
        echo "   https://github.com/C2SM-RCM/extpar"
        echo "To get access, contact jonas.jucker@c2sm.ethz.ch"
        echo "Or check if there is already a compiled version of extpar on your computer by"
        echo "someone else and set variable 'extpar_source_dir' to point to it."
        exit 1
      else
        echo "Cloned extpar from git@gitlab.dkrz.de:m212005/extpar4jsbach.git"
      fi
    else
      echo "Cloned extpar from git@github.com:C2SM-RCM/extpar.git"
    fi
    cd ${extpar_source_dir}
    git checkout ${extpar_version}
    git submodule update
  fi

  # Configure and compile/make extpar if necessary assuming that the existence of the
  # modules.env file indicates a complete installation of extpar with the binaries
  # and python scripts already in the bin directory
  if [[ ! -f ${extpar_source_dir}/modules.env ]]; then
    if [[ ! -w ${extpar_source_dir} ]]; then
      echo "No write permission for extpar source directory."
      exit 1
    fi
    cd ${extpar_source_dir}
    # If not running on levante or breeze4, add a case for your machine here
    case $(hostname) in
      levante*|*.lvt.dkrz.de)
        ./configure.levante.gcc
        ;;
      breeze4)
        ./configure.breeze4.gcc
        ;;
      *)
        echo ERROR: Unknown host: $(hostname)
        exit 1
        ;;
    esac

    source modules.env   # This might changed the modules loaded above!
    make -j 4
  fi

  # Generate extpar data for Jsbach
  cd ${scripts_dir}
  ./extpar4jsbach_mpim_icon.sh
  cp ./extpar4jsbach_mpim_icon.sh ${bc_file_dir}/scripts
fi

#-----------------------------------------------------------------------------
# 3. Run bc_files_from_extpar
if [[ ${run_bc_files_from_extpar} == true ]]; then

  if [[ ! -f ${extpar_dir}/${extpar_file} ]]; then
    echo "${extpar_dir}/${extpar_file} does not exist."
    exit 1
  fi

  ${scripts_dir}/bc_files_from_extpar.sh
  cp ${scripts_dir}/bc_files_from_extpar.sh ${bc_file_dir}/scripts
fi

#-----------------------------------------------------------------------------
# 4. Run adapt_nwp_extpar_file
if [[ ${run_adapt_nwp_extpar_file} == true ]]; then
  ${scripts_dir}/adapt_nwp_extpar_file.sh
  cp ${scripts_dir}/adapt_nwp_extpar_file.sh ${bc_file_dir}/scripts
fi

#-----------------------------------------------------------------------------
# 5. Run calc_hd_receive_mask: generate receive mask for external HD model
if [[ ${run_calc_hd_receive_mask} == true ]]; then
  ${scripts_dir}/calc_hd_receive_mask.sh
  cp ${scripts_dir}/calc_hd_receive_mask.sh ${hd_file_dir}/scripts
  cp ${scripts_dir}/${script_name}          ${hd_file_dir}/scripts
fi

#-----------------------------------------------------------------------------
# 6. Generate the hdpara file: generate HD parameter file for internal HD
if [[ ${run_generate_hdpara_file} == true ]]; then
  ${scripts_dir}/generate_hdpara_file.sh
  cp ${scripts_dir}/generate_hdpara_file.sh ${hd_file_dir}/scripts
fi

cp ${scripts_dir}/${script_name}            ${bc_file_dir}/scripts

echo ""
echo "====  Initial and boundary file generation completed. ===="
echo ""

#-----------------------------------------------------------------------------
exit 0
