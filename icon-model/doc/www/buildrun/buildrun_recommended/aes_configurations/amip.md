```{eval-rst}
:orphan:
```

(ref_buildrun_amip)=
# AMIP Configuration

AMIP is a standard experimental protocol for global Atmospheric General Circulation Models (AGCMs) and prescribes realistic Sea Surface Temperature and Sea Ice from 1979 till near present. Virtually the entire international climate modeling community has participated in this project since its inception in 1990. This framework enables a diverse community of scientists to analyze AGCMs in a systematic fashion, a process which serves to facilitate model improvement. It is managed by PCMDI (Earth System Model Evaluation Project). See [AMIP](https://pcmdi.llnl.gov/mips/amip/home/overview.html).

Input dependencies:
: Files to set up the initial and boundary conditions are in `/pool/data/ICON/grids/public/mpim/`. Depending of the desired horizontal grid spacing, different folders are used. For instance, the file for the initial and boundary conditions for a grid spacing of ~40 km (R02B06)  are located in `/pool/data/ICON/grids/public/mpim//0052/`, and for a grid spacing of ~10 km (R02B08), in `/pool/data/ICON/grids/public/mpim/0054/`. Inside the corresponding folder (`/pool/data/ICON/grids/public/mpim/XXXX/`), there will be the grid file: `icon_grid_XXXX_R0NB0N_G.nc` and the folder containing initial conditions (`/pool/data/ICON/grids/public/mpim/XXXX/initial_conditions`), ozone (`/pool/data/ICON/grids/public/mpim/XXXX/ozone`), sea surface temperature and sea ice area (`/pool/data/ICON/grids/public/mpim/XXXX/sst_and_ice`), aerosols (`/pool/data/ICON/grids/public/mpim/XXXX/aerolsol_kine`), and land properties (`/pool/data/ICON/grids/public/mpim/XXXX/land`). The initial conditions, ozone, sea surface temperature, sea ice, aerosols and land properties are already interpolated to the used grid.

Compatible machines and compilers:
: Levante, CPU (Intel,GNU,NAG) and GPU (NVHPC)

Recommended resources:
: Depends of the horizontal resolution and if CPUs or GPUs are used. For CPUs, R02B06 run with 4 nodes, R02B08 with 32 nodes, R02B09 with ~500 nodes and R02B10 with ~600 nodes. For GPUs, R02B04 and R02B06 use one node with nblocks_c = 1 and nblocks_sub=4, R02B08 can run with 8, 16 and more than 32 nodes. However, the configuration of the nblocks_c and nblocks_sub differs. With 8 nodes, nblocks_c = 2 and nblocks_sub = 5. With 16 nodes, nblocks_c = 1 and nblocks_sub = 4. With 32 nodes and more, nblocks_c = 1 nblocks_sub = 1. R02B09 run with 24 nodes with nblocks_c = 1 and nproma_sub = 800.

Estimated runtime  (for resources indicated above):
: It takes 0.8 node hour for one hour simulation in R02B08. A similar performance is expected for GPUs.

Sources:
: {{ '[Config]({}/run/examples/amip-aes.config)'.format(base_url) }} for `mkexp` and CPU. For GPU, {{ '[Config]({}/run/examples/amip-aes_gpu.config)'.format(base_url) }} for `mkexp`. For CPUs and GPUs, the parameters can be modified in {{ '[default parameter settings]({}/run/exp.aes_amip)'.format(base_url) }}.

## Description
ICON can run in an AMIP configuration. This means solving the fluid dynamics equations in the atmopshere with three main paramaterization (microphysics, turbulence and radiation) on the entire globe. In this configuration, the atmosphere is coupled to a dynamical 1D land module. This means that land does not transport energy nor water horizontally. The ocean and the sea-ice are not dynamically active (uncoupled). In other words, sea surface temperature and sea-ice extension are prescribed.

In the default configuration, the atmosphere is divided in 90 vertical levels and the land uses 5 soil layers. The configuration of the grid is global and it can be used with horizontal grid spacing of 160 km and finer. The type of grid supported by `mkexp` are R2B4, R2B6, R2B8, R02B9, and R2B10. See [here](ref_buildrun_gridextpar) more information about grids.

The default configuration has the initial conditions of January 1st 2020. Other initial conditions can be found `/pool/data/ICON/grids/public/mpim/XXXX/initial_conditions/r0100`. Two dates are available to be used as initial conditions across the different grid configurations: 1979-01-01 and 2020-01-01. But certain grid configurations have more available dates, e.g., R02B08 has four: 1979-01-01, 1990-01-01, 2020-01-01, and 2020-08-01 (`/pool/data/ICON/grids/public/mpim/0054/initial_conditions/r0100`). If another date want to be used, use the script `/pool/data/ICON/grids/public/mpim/XXX/initial_conditions/r0100/11u-make-initial-data-ifs2icon-from-era5.sh` (only one script) to generate a new initial conditions from a pool of files in `/pool/data/ICON/grids/private/mpim/icon_preprocessing/source/ecmwf_initial_data/initial_conditions)`.

Aside from the initial conditions, the AMIP configuration needs time-dependent input data from aerosols, ozone, land properties, sea surface temperature, and sea ice extent. These inputs are automatically updated with the simulation time.

## Example
### Run AMIP with a R02B08 grid in Levante

- Follows steps for obtaining the code and configuring and building in [Quick Start](https://docs.icon-model.org/buildrun/buildrun_quickstart.html).

- To use GPU, use the following command in the folder containing the repository before using `make`

```
./config/dkrz/levante.gpu.nvhpc
```

- Then copy the example of the desired configuration in the run folder. `amip-aes.config` for CPUs and `amip-aes_gpu.config` for GPUs.
```
cd run/
cp ./example/amip-aes.config ./
```

- Open the amip-aes.config and choose the type of grid to use. Other grids are specified at the end of the file
```
EXP_TYPE = amip-aes-R2B8
```

- Change the account. The account starts with the letters mh or bb, following by numbers
```
ACCOUNT = mhXXXX
```

- Change the initial and last day of simulation as well as the interval for the restart file
```
INITIAL_DATE = 2020-01-01T00:00:00
FINAL_DATE = 2020-01-01T03:00:00
INTERVAL = PT3H
```

- Change the simulation ID to convenience
```
EXP_ID = amip-aes-R2B8
```


- Generate the run scripts using mkexp
```
../utils/mkexp/mkexp amip-aes.config
```

- This will generate a run script `../experiment/EXP_ID/scripts/EXP_ID.run_start`. Navigate to the folder and run the script
```
cd ../experiment/
sbatch EXP_ID.run_start
```
