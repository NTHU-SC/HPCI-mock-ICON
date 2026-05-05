# Welcome to the ICON Model documentation

:::topic
_It's the job that's never started as takes longest to finish - J.R.R. Tolkien_
:::

This website contains a collection of ICON documentation including references to documentation available at other places. We do not claim that this documentation is complete, but we hope you will still find it helpful.

```{toctree}
:hidden:
buildrun/buildrun_quickstart.md
tools/tools.md
atmosphere/atmosphere.md
ocean/ocean.md
waves/waves.md
land/land.md
infrastructure/infrastructure.md
literature/literature.md
```

:::{admonition} Release Information
:class: admonition-icontheme
ICON {{ '{}'.format(release) }} has been published and is available for download.
Information on the changes are available in the {{ '[**Release Notes**]({}/RELEASE_NOTES.md)'.format(base_url) }}.
:::

::::{grid} 1 2 2 3
:gutter: 1 1 1 2

:::{grid-item-card}
[**Getting Started**](ref_buildrun_quickstart)
^^^
[Building](ref_buildrun_introduction) & [Running](ref_buildrun_running)  
[](ref_buildrun_environments)  
[](ref_buildrun_gridextpar)  
[](ref_buildrun_icbc)  
[](ref_buildrun_recommconf)  
:::

:::{grid-item-card}
[**Interfaces & Tools**](ref_tools)
^^^
[Zonda](ref_tools_gridextpargui)  
[](ref_tools_yac)  
[](ref_tools_cdo)  
[](ref_tools_comin)
[Forward Operators](ref_tools_fwo)
:::

:::{grid-item-card}
[**Atmosphere**](ref_atmosphere)
^^^
[](ref_atmosphere_dycore)  
[](ref_atmosphere_physics)  
[](ref_atmosphere_art)
:::

:::{grid-item-card}
[**Ocean**](ref_ocean_overview)
^^^
[Sea-ice Model](ref_ocean_seaice)  
[Ocean Biogeochemistry](ref_ocean_biogeochem)  
:::

:::{grid-item-card}
[**Waves**](ref_waves_overview)
^^^
[Configuration](ref_waves_config)  
[Output Parameters](ref_waves_output)
:::

:::{grid-item-card}
[**Land**](ref_land)
^^^
[](ref_land_schemes)  
[](ref_land_cover_change)  
[](ref_land_biogeochem)  
:::

:::{grid-item-card}
[**Infrastructure**](ref_infrastructure)
^^^
[](ref_infrastructure_parallelization)  
[](ref_infrastructure_io)  
[](ref_infrastructure_testing)  
:::

:::{grid-item-card}
[**Literature**](ref_literature)
^^^
[](ref_literature_tutorials)  
[](ref_literature_technical)  
[](ref_literature_science)  
:::

::::
