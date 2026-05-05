```{eval-rst}
:orphan:
```

(ref_buildrun_omip)=
# OMIP Configuration

Input Files:
- R2B4 L40 grid:<br>`OceanOnly_Icos_0158km_etopo40.nc` -> `/pool/data/ICON/oes/input/r0002/OceanOnly_Icos_0158km_etopo40/OceanOnly_Icos_0158km_etopo40.nc`
- OMIP forcing:<br>`ocean-flux.nc` -> `/pool/data/ICON/oes/input/r0002/OceanOnly_Icos_0158km_etopo40/omipForcing-mpiomDaily-OceanOnly_Icos_0158km_etopo40.nc`
- PHC initial conditions:<br>`initial_state.nc` -> `/pool/data/ICON/oes/input/r0002/OceanOnly_Icos_0158km_etopo40/omipInitialState-AnnualAverage-OceanOnly_Icos_0158km_etopo40-40levels.nc`
- PHC salinity restoring:<br>`ocean-relax.nc` -> `/pool/data/ICON/oes/input/r0002/OceanOnly_Icos_0158km_etopo40/omipRelaxSurface-OceanOnly_Icos_0158km_etopo40.nc`

Compatible machines and compilers:
: Levante (CPU: nag, gfortran, ifort, and GPU: nvhpc)

Recommended resources:
: Single node

Estimated runtime (for resources indicated above):
: few minutes

Sources:
: {{ '[Config]({}/run/checksuite.ocean_internal/test_oes_omip.config)'.format(base_url) }} for `mkexp`.

Analysis/postprocessing:
: (under development).


## Description

This is a coarse configuration (r2b4 l40 which corresponds to 160km horizontal resolution) of the ICON ocean model component which is not coupled to the atmosphere component but it includes sea ice. The configuration is driven with OMIP (Ocean Model Intercomparison Project) forcing data and it uses a realistic bathymetry. It applies the so-called `zstar` vertical coordinate.
