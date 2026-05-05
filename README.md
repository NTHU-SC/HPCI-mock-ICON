# ICON

# Task description

[GitHub](https://github.com/NTHU-SC/HPCI-mock-ICON) (Git LFS is required)

This task requires the optimization of ICON application.

## Introduction

ICON, which obtains its name from the usage of spherical grids derived from the icosahedron (ICO) and the non-hydrostatic (N) dynamics, is a flexible, scalable, high-performance modeling framework for weather, climate and environmental prediction. ICON includes component models for the atmosphere, the ocean and the land, as well as chemical and biogeochemical cycles, all implemented on the basis of common data structures and sharing the same efficient technical infrastructure. ICON can be used in the most diverse resolutions and configurations to enable a whole range of applications — from global and regional weather forecasts and climate projections to very high-resolution digital twins of the Earth system. By providing actionable information for society and advances our understanding of the Earth's climate system. Thanks to these modeling capabilities, ICON is evolving into an internationally recognized and widely used modeling framework that advances our understanding of the Earth's climate system and provides actionable information for solving problems of high societal relevance.

## Testcase

* The test case provided is an aquaplanet configuration on the R02B04 grid.
* The test case simulation is based on a restart configuration that defines an atmosphere in equilibrium.
* This directory includes:
  * The test case configuration `ape_from_spinup_short.config`
  * The restart data directory `ape_from_spinup_short_restart_atm_19790701T000000Z.nc`
  * The restart date file `ape_from_spinup_short.date`
  * The restart run file `ape_from_spinup_short.run` (for reference)
  * Necessary input files (`pool/*`) (do not change the input file internal directory structure)

## Source code of ICON

You may use the icon-model source code in the provided package. There are lfs files in it, make sure you run `git lfs install` before you clone the repo.

## Dependencies

* Required
  * C, CXX and Fortran compilers
  * Software Libraries: MPI, ZLIB, HDF5, NetCDF (C and Fortran, with NetCDF-4 support), FYAML, XML2, LAPACK and BLAS
* Optional (feel free to use other tools you want)
  * CURL, GNU Make v3.81+, CMake v3.18+, Perl v5.10+
  * Python v3.9+, a Python environment with the 'six' and 'jinja2' packages installed

\*You may add any other dependency

## Build ICON

* Configure ICON options and compilation flags with a configuration wrapper file:
  * Navigate to the `icon-model` directory and create a new directory at the path `config/asc/<config-file-name>`.
  * Create your configuration wrapper file using as a reference the templates in the `config/` directory (e.g., `config/generic/gcc`). This configuration file might need to be adapted or adjusted to match your target hardware platform
  * Prepare the build by running your wrapper file from the `icon-model` directory
* Alternatively, you can also run your own configure command as long it succeeds.
* Build the source code with `make` (use flags like `-j32` to speedup)
* This step will generate an executable binary file `icon` at `bin/`

## Run Script Configuration

* Go to directory `icon-model`
* Copy `ape_from_spinup.config` config file to `run/`
  * adjust the `INPUT_ROOT` parameter to point to the provided input file directory `pool/data` if necessary.
  * (DO NOT MODIFY) The case settings are:
    * EXP_ID = ape_from_spinup_short
    * INITIAL_DATE = 1979-07-01
    * FINAL_DATE = 1979-07-03
    * ATMO_TIME_STEP = PT1M
    * INTERVAL = P2M
* Run `./utils/mkexp/mkexp` run/ape_from_spinup.config to create the experiment directory tree including scripts/, data/, work/, log/

## Run the test case

* Go to `icon-model/experiments/ape_from_spinup_short/scripts`, adjust run parameters in the `ape_from_spinup_short.run` file as needed based on your platform
* Run the run script `ape_from_spinup_short.run`
  * WARNING! The original run script uses Slurm, but you do not have Slurm installed on your cluster. You have to modify the script so that it can run successfully.
* Redirect the experiment log to your own `output.log` or find the experiment log in `experiments/ape_from_spinup_short/scripts/ape_from_spinup_short.run.log` (log archive can be found in `experiments/ape_from_spinup_short/log`) and change the file name to `output.log` when submitting. If the simulation completes successfully, you will see this in log file: 

```javascript
===========================
Script run successfully: OK
===========================
```

* The output files will be generated at `work/run_19790701-19790703`, we will ask you to provide a few files to validate correctness.

## Scoring

* (30%) Install dependencies
  * Show all the installation of all required dependencies
* (15%) Modify run script and start simulation
  * Show `ape_from_spinup_short.run.log`'s results, even if the result is incorrect.
* (15%) Correct simulation
  * Show `ape_from_spinup_short.run.log`'s results, the result has to be correct.
* (20%) Performance
  * Optimize the model and decrease your execution time.
  * Show your optimizations in the report.
* (10%) Profile
  * Use profiler to find the bottleneck.
  * Show your analysis in the report.

## Sample Output

```javascript
2026-05-05T22:22:43.805:  Timer report, ranks 0-7
2026-05-05T22:22:43.807:
2026-05-05T22:22:43.808:  -----------------------------------   -------   ------------   --------   ------------   ------------   --------   -------------   --------------   -------------   --------------   -------------   -----
2026-05-05T22:22:43.808:  name                                  # calls   t_min          min rank   t_avg          t_max          max rank   total min (s)   total min rank   total max (s)   total max rank   total avg (s)   # PEs
2026-05-05T22:22:43.808:  -----------------------------------   -------   ------------   --------   ------------   ------------   --------   -------------   --------------   -------------   --------------   -------------   -----
2026-05-05T22:22:43.808:
2026-05-05T22:22:43.808:  total                                 8               27m37s   [1]              27m37s         27m37s   [5]            1657.464    [1]                  1657.476    [5]                  1657.470    8
2026-05-05T22:22:43.808:   L wrt_output                         16            0.03164s   [4]            0.17526s       0.42813s   [0]               0.064    [5]                     0.799    [0]                     0.351    8
2026-05-05T22:22:43.808:   L integrate_nh                       7200           1.7595s   [0]             1.7998s        2.1973s   [4]            1615.363    [3]                  1624.735    [7]                  1619.793    8
2026-05-05T22:22:43.808:      L nh_solve                        36000         0.14650s   [0]            0.16283s       0.20595s   [6]             720.933    [0]                   738.711    [7]                   732.729    8
2026-05-05T22:22:43.808:         L nh_solve.veltend             43200         0.01110s   [4]            0.02351s       0.04523s   [0]             120.182    [6]                   134.786    [2]                   126.979    8
2026-05-05T22:22:43.808:         L nh_solve.cellcomp            72000         0.00234s   [3]            0.01006s       0.02564s   [1]              79.712    [7]                   104.119    [1]                    90.505    8
2026-05-05T22:22:43.808:         L nh_solve.edgecomp            72000         0.00206s   [0]            0.01501s       0.04550s   [0]             127.829    [7]                   143.666    [2]                   135.061    8
2026-05-05T22:22:43.808:         L nh_solve.vnupd               72000         0.00263s   [2]            0.01733s       0.04112s   [0]             151.254    [4]                   162.327    [2]                   155.951    8
2026-05-05T22:22:43.808:         L nh_solve.vimpl               72000         0.00550s   [6]            0.01487s       0.02590s   [1]             129.972    [2]                   138.779    [6]                   133.850    8
2026-05-05T22:22:43.808:         L nh_solve.exch                144000        0.00045s   [1]            0.00475s       0.05344s   [4]              55.454    [2]                   113.822    [5]                    85.582    8
2026-05-05T22:22:43.808:      L nh_hdiff                        7200          0.02948s   [1]            0.03724s       0.04898s   [0]              32.900    [7]                    33.999    [4]                    33.512    8
2026-05-05T22:22:43.808:      L transport                       7200          0.38681s   [2]            0.40509s       0.56098s   [7]             363.166    [2]                   366.263    [7]                   364.579    8
2026-05-05T22:22:43.808:         L adv_horiz                    7200          0.33250s   [2]            0.35029s       0.50324s   [6]             313.133    [0]                   317.649    [7]                   315.263    8
2026-05-05T22:22:43.809:            L adv_hflx                  7200          0.32234s   [0]            0.33791s       0.49237s   [6]             302.500    [0]                   305.830    [7]                   304.115    8
2026-05-05T22:22:43.809:               L back_traj              7200          0.00257s   [6]            0.00584s       0.00959s   [2]               3.678    [6]                     6.646    [2]                     5.260    8
2026-05-05T22:22:43.809:         L adv_vert                     7200          0.01430s   [3]            0.02683s       0.03702s   [7]              23.241    [4]                    24.923    [1]                    24.144    8
2026-05-05T22:22:43.809:            L adv_vflx                  7200          0.01206s   [4]            0.01682s       0.02725s   [7]              14.171    [4]                    15.819    [7]                    15.138    8
2026-05-05T22:22:43.809:      L iconam_aes                      7200          0.42235s   [3]            0.44905s       0.54976s   [1]             402.527    [1]                   406.346    [4]                   404.147    8
2026-05-05T22:22:43.809:         L dyn2phy                      7200          0.01288s   [2]            0.02892s       0.04088s   [0]              24.757    [3]                    26.864    [6]                    26.030    8
2026-05-05T22:22:43.809:            L d2p_sync                  14400         0.00060s   [6]            0.00400s       0.02165s   [0]               6.220    [2]                     7.954    [0]                     7.199    8
2026-05-05T22:22:43.809:         L aes_bcs                      7200          0.00001s   [4]            0.00006s       0.03299s   [7]               0.048    [5]                     0.054    [7]                     0.051    8
2026-05-05T22:22:43.809:         L aes_phy                      7200          0.34097s   [7]            0.36280s       0.46017s   [3]             325.111    [5]                   328.829    [2]                   326.522    8
2026-05-05T22:22:43.809:            L interface_aes_rad         576000        0.00000s   [0]            0.00000s       0.00085s   [0]               0.003    [1]                     0.004    [0]                     0.003    8
2026-05-05T22:22:43.809:            L interface_aes_rht         576000        0.00002s   [0]            0.00006s       0.00909s   [0]               3.160    [1]                     4.481    [3]                     4.013    8
2026-05-05T22:22:43.809:            L interface_aes_vdf         7200          0.25674s   [7]            0.28508s       0.38335s   [3]             251.881    [4]                   263.756    [1]                   256.575    8
2026-05-05T22:22:43.809:               L vdiff_down             7200          0.21109s   [4]            0.23011s       0.35301s   [3]             203.466    [4]                   212.773    [1]                   207.099    8
2026-05-05T22:22:43.809:               L update_surface         576000        0.00001s   [0]            0.00002s       0.00112s   [0]               1.406    [4]                     1.764    [2]                     1.579    8
2026-05-05T22:22:43.809:               L vdiff_up               576000        0.00016s   [1]            0.00037s       0.01761s   [3]              25.482    [7]                    28.111    [2]                    26.477    8
2026-05-05T22:22:43.809:            L interface_cloud_mig       576000        0.00020s   [2]            0.00069s       0.00693s   [3]              45.904    [1]                    53.957    [4]                    49.856    8
2026-05-05T22:22:43.809:               L cloud_mig              576000        0.00016s   [2]            0.00064s       0.00688s   [3]              42.510    [1]                    50.443    [4]                    46.403    8
2026-05-05T22:22:43.809:                  L satad               1152000       0.00003s   [5]            0.00007s       0.00616s   [3]              10.137    [6]                    10.741    [7]                    10.462    8
2026-05-05T22:22:43.809:                  L graupel             576000        0.00006s   [2]            0.00048s       0.00462s   [3]              30.780    [1]                    38.506    [4]                    34.641    8
2026-05-05T22:22:43.809:            L interface_aes_wmo         576000        0.00001s   [0]            0.00002s       0.00121s   [2]               1.244    [1]                     1.382    [2]                     1.302    8
2026-05-05T22:22:43.809:            L diagnose_cov              1152000       0.00000s   [0]            0.00001s       0.00238s   [6]               1.540    [4]                     1.787    [6]                     1.664    8
2026-05-05T22:22:43.809:         L phy2dyn                      7200          0.03170s   [7]            0.04277s       0.12362s   [4]              36.059    [1]                    41.335    [4]                    38.492    8
2026-05-05T22:22:43.809:            L p2d_sync                  14400         0.00061s   [1]            0.01084s       0.08713s   [5]              16.901    [2]                    22.353    [4]                    19.518    8
2026-05-05T22:22:43.809:         L diagnose_tcw                 576000        0.00000s   [0]            0.00000s       0.00009s   [2]               0.003    [3]                     0.004    [5]                     0.003    8
2026-05-05T22:22:43.809:         L diagnose_qvi                 576000        0.00002s   [7]            0.00008s       0.00701s   [4]               4.931    [1]                     6.209    [2]                     5.749    8
2026-05-05T22:22:43.809:         L diagnose_uvi                 1152000       0.00001s   [0]            0.00005s       0.00472s   [3]               6.637    [4]                     7.448    [3]                     7.008    8
2026-05-05T22:22:43.809:   L action                             7200          0.00001s   [3]            0.00008s       0.00150s   [1]               0.056    [3]                     0.085    [7]                     0.075    8
2026-05-05T22:22:43.809:  exch_data                             566736        0.00000s   [3]            0.00336s       0.39569s   [5]             150.596    [2]                   314.256    [6]                   238.383    8
2026-05-05T22:22:43.809:   L exch_data.wait                     566736        0.00000s   [0]            0.00270s       0.39523s   [5]             104.161    [2]                   264.883    [6]                   191.070    8
2026-05-05T22:22:43.809:  nh_diagnostics                        7216          0.01101s   [0]            0.02088s       0.04263s   [3]              14.829    [7]                    23.917    [3]                    18.837    8
2026-05-05T22:22:43.810:  diagnose_pres_temp                    28816         0.00504s   [4]            0.01120s       0.02860s   [3]              37.460    [4]                    42.837    [1]                    40.338    8
2026-05-05T22:22:43.810:  model_init                            24            0.01850s   [5]            0.74984s        1.4169s   [3]               2.242    [6]                     2.264    [0]                     2.250    8
2026-05-05T22:22:43.810:   L global_sum                         8             0.00011s   [3]            0.00100s       0.00189s   [7]               0.000    [3]                     0.002    [7]                     0.001    8
2026-05-05T22:22:43.810:   L prep_aes_phy                       8             0.03536s   [1]            0.03667s       0.04003s   [0]               0.035    [1]                     0.040    [0]                     0.037    8
2026-05-05T22:22:43.811:   L compute_domain_decomp              8             0.31877s   [5]            0.31891s       0.31969s   [0]               0.319    [5]                     0.320    [0]                     0.319    8
2026-05-05T22:22:43.812:      L ordglb_sum                      168           0.00006s   [3]            0.00015s       0.00018s   [1]               0.003    [3]                     0.003    [7]                     0.003    8
2026-05-05T22:22:43.812:   L compute_intp_coeffs                8              1.0575s   [4]             1.0577s        1.0583s   [3]               1.057    [4]                     1.058    [3]                     1.058    8
2026-05-05T22:22:43.812:   L init_ext_data                      8             0.00008s   [1]            0.00071s       0.00138s   [4]               0.000    [1]                     0.001    [4]                     0.001    8
2026-05-05T22:22:43.812:   L read_restart_files                 8             0.41706s   [2]            0.41956s       0.42396s   [0]               0.417    [2]                     0.424    [0]                     0.420    8
2026-05-05T22:22:43.812:      L load_restart                    8             0.41706s   [2]            0.41951s       0.42359s   [0]               0.417    [2]                     0.424    [0]                     0.420    8
2026-05-05T22:22:43.812:         L load_restart_io              19980         0.00000s   [0]            0.00001s       0.00838s   [0]               0.000    [3]                     0.102    [0]                     0.013    8
2026-05-05T22:22:43.812:         L load_restart_comm_setup      32            0.00000s   [0]            0.00433s       0.01115s   [5]               0.014    [0]                     0.018    [5]                     0.017    8
2026-05-05T22:22:43.812:         L load_restart_communication   19968         0.00000s   [0]            0.00015s       0.00551s   [7]               0.303    [0]                     0.397    [7]                     0.384    8
2026-05-05T22:22:43.812:         L load_restart_get_var_id      16            0.00000s   [0]            0.00020s       0.00323s   [0]               0.000    [1]                     0.003    [0]                     0.000    8
2026-05-05T22:22:43.812:  upper_atmosphere                      32            0.00000s   [7]            0.00006s       0.00036s   [5]               0.000    [3]                     0.000    [4]                     0.000    8
2026-05-05T22:22:43.812:   L upatmo_construction                16            0.00001s   [3]            0.00010s       0.00036s   [4]               0.000    [3]                     0.000    [4]                     0.000    8
2026-05-05T22:22:43.812:   L upatmo_destruction                 16            0.00000s   [0]            0.00001s       0.00017s   [0]               0.000    [7]                     0.000    [0]                     0.000    8
2026-05-05T22:22:43.812:  write_restart                         24            0.00324s   [0]            0.10351s       0.28942s   [5]               0.295    [0]                     0.316    [7]                     0.311    8
2026-05-05T22:22:43.812:   L write_restart_io                   25            0.00148s   [0]            0.03371s       0.05188s   [2]               0.099    [5]                     0.114    [2]                     0.105    8
2026-05-05T22:22:43.812:   L write_restart_communication        8             0.00000s   [1]            0.00000s       0.00000s   [2]               0.000    [1]                     0.000    [2]                     0.000    8
2026-05-05T22:22:43.812:   L write_restart_wait                 24            0.00000s   [0]            0.03195s       0.17490s   [3]               0.009    [4]                     0.175    [3]                     0.096    8
2026-05-05T22:22:43.812:   L write_restart_setup                16            0.00017s   [0]            0.00949s       0.02094s   [7]               0.003    [0]                     0.025    [5]                     0.019    8
2026-05-05T22:22:43.812:      L write_restart_collectors        760           0.00000s   [0]            0.00000s       0.00018s   [5]               0.000    [0]                     0.000    [5]                     0.000    8
2026-05-05T22:22:43.812:         L write_restart_indices        24            0.00001s   [0]            0.00004s       0.00009s   [5]               0.000    [0]                     0.000    [5]                     0.000    8
2026-05-05T22:22:43.812:  ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
2026-05-05T22:22:43.812:
2026-05-05T22:22:43.812:  mo_ext_data_state:destruct_ext_data: Destruction of data structure for external data started
2026-05-05T22:22:43.812:  mo_ext_data_state:destruct_ext_data: Destruction of data structure for external data finished
2026-05-05T22:22:43.812:  mo_interpolation:destruct_int_state: start to destruct int state
2026-05-05T22:22:43.812:  mo_interpolation:destruct_int_state: destruction of interpolation state finished
2026-05-05T22:22:43.812:  mo_atmo_model:destruct_atmo_model: destruct_2d_interpol_state is done
2026-05-05T22:22:43.813:  mo_alloc_patches:destruct_comm_patterns: start
2026-05-05T22:22:43.813:  mo_alloc_patches:destruct_comm_patterns: destruct_comm_patterns finished
2026-05-05T22:22:43.813:  mo_alloc_patches:destruct_patches: start
2026-05-05T22:22:43.813:  mo_alloc_patches:destruct_patches: destruct_patches finished
2026-05-05T22:22:43.813:  mo_alloc_patches:destruct_patches: start
2026-05-05T22:22:43.813:  mo_alloc_patches:destruct_patches: destruct_patches finished
2026-05-05T22:22:43.813:  mo_atmo_model:destruct_atmo_model: destruct_patches is done
2026-05-05T22:22:43.813:  mo_atmo_model:destruct_atmo_model: clean-up finished
2026-05-05T22:22:43.909: Tue May  5 10:22:43 PM UTC 2026
2026-05-05T22:22:43.910: quiz case finished; skip restart hand-over
2026-05-05T22:22:43.912: OK
2026-05-05T22:22:43.912: ============================
2026-05-05T22:22:43.912: Script run successfully: OK
2026-05-05T22:22:43.912: ============================
```

## Submission

Pack in `ICON_{YOUR_GROUP_NUMBER}.zip`

* `output.log`
* `nml.atmo.log`
* `ape_from_spinup_quiz_atm_2d_P1D_19790701T000000Z.nc`
* `ape_from_spinup_quiz_atm_3d_P1D_19790701T000000Z.nc`