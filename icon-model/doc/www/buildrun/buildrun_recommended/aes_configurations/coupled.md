```{eval-rst}
:orphan:
```

(ref_buildrun_aes_coupled)=
# Coupled AES Configuration

Input dependencies:
:   `MODEL_DIR` is the checkout directory, `INPUT_ROOT` is `/pool/data/ICON` on levante; output for reduced resolution version
    ```
    $INPUT_ROOT/grids/public/mpim/0036/icon_grid_0036_R02B04_O.nc
    $INPUT_ROOT/grids/public/mpim/0043-0036/land/r0002/bc_land_frac.nc
    $INPUT_ROOT/grids/public/mpim/0043-0036/land/r0002/bc_land_phys.nc
    $INPUT_ROOT/grids/public/mpim/0043-0036/land/r0002/bc_land_soil.nc
    $INPUT_ROOT/grids/public/mpim/0043-0036/land/r0002/bc_land_sso.nc
    $INPUT_ROOT/grids/public/mpim/0043-0036/land/r0002/hdpara_icon_r2b4_013_using_grid_lsmask.nc
    $INPUT_ROOT/grids/public/mpim/0043-0036/land/r0002/hdrestart_R02B04_013_G_210120_1334_with_grid_file_mask.nc
    $INPUT_ROOT/grids/public/mpim/0043-0036/land/r0002/ic_land_soil.nc
    $INPUT_ROOT/grids/public/mpim/0043/aerosol_kinne/r0001/bc_aeropt_kinne_lw_b16_coa.nc
    $INPUT_ROOT/grids/public/mpim/0043/aerosol_kinne/r0001/bc_aeropt_kinne_sw_b14_coa.nc
    $INPUT_ROOT/grids/public/mpim/0043/aerosol_kinne/r0001/bc_aeropt_kinne_sw_b14_fin_2014.nc
    $INPUT_ROOT/grids/public/mpim/0043/icon_grid_0043_R02B04_G.nc
    $INPUT_ROOT/grids/public/mpim/0043/initial_condition/r0001/ifs2icon_1979010100_R02B04_G.nc
    $INPUT_ROOT/grids/public/mpim/0043/ozone/r0001/bc_ozone_historical_2014.nc
    $INPUT_ROOT/grids/public/mpim/0043/sst_and_seaice/r0001/bc_sic_1979_2016.nc
    $INPUT_ROOT/grids/public/mpim/0043/sst_and_seaice/r0001/bc_sst_1979_2016.nc
    $INPUT_ROOT/grids/public/mpim/independent/greenhouse_gases/greenhouse_ssp245.nc
    $INPUT_ROOT/grids/public/mpim/independent/solar_radiation/3.2/swflux_14band_cmip6_1850-2299-v3.2.nc
    $MODEL_DIR/data/ECHAM6_CldOptProps_rrtmgp_lw.nc
    $MODEL_DIR/data/ECHAM6_CldOptProps_rrtmgp_sw.nc
    $MODEL_DIR/data/rrtmgp-gas-lw-g128.nc
    $MODEL_DIR/data/rrtmgp-gas-sw-g112.nc
    $MODEL_DIR/externals/jsbach/data/lctlib_nlct21.def
    $MODEL_DIR/run/dict.iconam.mpim
    $INPUT_ROOT/grids/public/mpim/0036/ocean/restart/r0001/ler1166_restart_oce_21000101T000000Z.nc
    ```

Compatible machines and compilers:
:   &nbsp;
    - `levante_intel_hybrid`

Recommended resources:
: 4 CPU nodes

Estimated runtime (for resources indicated above):
: _The simplified test configuration (see [Description](#description)) is numerically stable only for a couple of days. It runs 6 model hours in about 40 s_

Sources:
:   &nbsp;
    - {{ '[Config]({}/run/mkexp/types/AES/coupled-R02B04L90.config)'.format(base_url) }} for `mkexp`
    - {{ '[test_aes_coupled.config]({}/run/checksuite.icon-dev/test_aes_coupled.config)'.format(base_url) }} - restart, nproma, openmp, mpi tests
    - {{ '[test_nextGEMS.config]({}/run/test_nextGEMS.config)'.format(base_url) }} - experiment workflow test

Analysis/postprocessing:
: For all models participating in the [nextGEMS](https://nextgems-h2020.eu/) project, examples, software and recipes for data processing have been collected on the [easyGEMS](https://easy.gems.dkrz.de/) site. See the ICON specific entries there for more information.


## Description

Unified coupled atmosphere/ocean configuration based on the AES atmosphere package as used for the [nextGEMS](https://nextgems-h2020.eu/) production experiments.
The AES package features a reduced set of sub-grid scale parametrizations due to the targeted km-scale model resolutions.
At the time of this writing, only a simple test model is implemented, using 160 km (R2B4) horizontal resolution with a comparably short timestep to compensate instabilities due to the reduced parametrization set. For scientific production, definitions for horizontal resolutions of at least 10 km atmospheric and 5 km oceanic (R2B8_R2B9) are currently under development.
