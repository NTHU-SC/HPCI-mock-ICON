```{eval-rst}
:orphan:
```

(ref_buildrun_nest)=
# Nest Configuration

Input dependencies:
: All files required to run this case (initial/boundary conditions and grids) are stored in "/pool/data/ICON/grids/public/mpim/0049/lamnest/r0100/". These are: bc_aeropt_kinne_lw_b16_coa_DOM*, bc_aeropt_kinne_sw_b14_coa_DOM*, bc_aeropt_kinne_sw_b14_fin_DOM*, bc_land_frac_DOM*,  bc_land_phys_DOM*, bc_land_soil_DOM*, bc_land_sso_DOM*, bc_ozone_DOM*, bc_sic_DOM*, bc_sst_DOM*, ic_atmo_DOM*, ic_land_soil_DOM* with * going from domain 01 to 04. Finally, the meshes are global_160km_DOM01.nc (parent domain), Tropical_Atlantic_80km_DOM02.nc (first-level nested domain), Tropical_Atlantic_40km_DOM03.nc and Tropical_Atlantic_40km_DOM04.nc (second-level nested domains).

Compatible machines and compilers:
: Levante, CPU (Intel,GNU,NAG) and GPU (NVHPC)

Recommended resources:
: Single Node.

Estimated runtime (for resources indicated above):
: It takes 25 minutes on one Levante compute node (2xAMD 7763, 256 Gb main memory) to simulate 7 days, for a setup combining MPI for inter-node parallelisation and OpenMP for multi-threading. The same case was tested on one Levante GPU node (4xA100, 160Gb + 2x AMD 7763 CPU; 128 cores in total, 512 GB/1024 GB main memory) and completed 7 days in 16 minutes.

Sources:
: {{ '[Config]({}/run/checksuite.atm/test_aes_nest.config)'.format(base_url) }} for `mkexp`. It was build from the older script {{ '[Config]({}/run/exp.atm_cph_nest_R2B4)'.format(base_url) }} which used `make_runscripts`.

Analysis/postprocessing:
: (under development).


## Description

This test case consists of four domains nested within one another, starting from a parent global R2B4 mesh. The second domain includes parts of the central and southern Atlantic Ocean, spanning from 70°W to 40°E and 50°S to 38°N. It contains a second level of nested domains: one over the ITCZ (48°W to 8°W and 0°N to 20°N) and the other over the stratocumulus region off the Namibian coast (9°W to 17°E and 30°S to 4°S). The initial and boundary conditions files were created by interpolating the R2B4 case files for January 1st, 1988, to all nested domains. R2B4 is a global 160km-resolution mesh, it corresponds to grid ID 0049 as documented in http://icon-downloads.mpimet.mpg.de/mpim_grids.xml.

This configuration not only tests how the model handles multiple domains but also many key features used by limited-area models, such as nudging and interpolation to the boundary and halo of non-global domains. Each domain generates a series of its own _domXX output and restart files, which are all tested within the current development workflow. This test does not, however, cover the routines in charge of reading lateral boundary conditions.

Physically, this test offers a practical evaluation of the ITCZ and Namibian clouds' sensitivity to grid resolution. Indeed, as the reversed child-to-parent domain feedback is disabled, the coarser meshes are uninformed of the finer meshes' solutions. Therefore, projecting all domains' solutions onto the last-level children domains provides a mesh resolution sensitivity study by itself.
