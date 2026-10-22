# Umzimvubu_CRA

Climate Risk Assessment (CRA) for the Umzimvubu River Catchment and Alfred Nzo District, South Africa.

**Interactive map:** <https://vervis.github.io/Umzimvubu_CRA/>

The tool maps the risk that agricultural livelihoods will be lost to water insecurity as climate change shifts rainfall patterns and increases meteorological drought (sub-risk SR1). It is meant to help decide where nature-based solutions should be implemented, particularly invasive alien plant management. It was developed for the Umzimvubu Hub of the ALBATROSS project.

The tool builds on the climate risk assessment of the GCF-funded project [*Nature-based solutions for transforming smallholder farming systems*](https://www.greenclimate.fund/document/nature-based-solutions-transforming-smallholder-farming-systems-are-vulnerable-impacts), which covered eastern South Africa at coarse resolution. This tool improves on that assessment in two ways:

- It produces risk layers at 1 km resolution.
- It adapts the indicators to the local context, using input from Umzimvubu Hub stakeholders and finer, more recent local datasets where they are available.

> **Status:** draft research output. The layers are indicative and are not intended for site-level decisions.

## Approach

The tool follows the IPCC AR6 risk concept: **risk = hazard × vulnerability × exposure**. Each component is built from spatial indicators. Each indicator is scaled from 0 to 1, where 1 means higher risk. Scaling spans all time periods, so projected change stays visible. The scores will be combined into component scores and an overall risk score; this is stages 40–50 of the pipeline, which are in development.

| Component | Indicators (active) | Main sources |
|---|---|---|
| Hazard | Mean annual rainfall; inter-annual rainfall variability; consecutive dry days; aridity index (P/PET, Hargreaves) | NASA NEX-GDDP-CMIP6 v2: ensemble median, 1981–2014 baseline, SSP1-2.6 and SSP2-4.5 for 2021–2040, 2041–2070 and 2071–2100 |
| Vulnerability | Smallholder presence; reliance on agriculture; grazing capacity; reliance on rainfall-dependent water sources; livestock density; old fields; invasive alien trees; lower-income households | SANLC 2022 (DFFE); Census 2011 and 2022 (Stats SA via DataFirst); DAFF 2016 grazing capacity; FAO GLW4; MapWAPS (Stellenbosch University) |
| Exposure | Rural settlement households; Strategic Water Source Areas | DWS WSKS settlements; SANBI/CSIR SWSA |

The full indicator definitions, sources, licences and citations are in the "About & data sources" section of the interactive map. The processing settings are in `R/30_indicators/indicator_register.csv`.

## Repository structure

```
_targets.R                  pipeline definition ({targets})
R/
  00_setup/                 project folders
  10_download/              raw data downloads (101-112)
  20_prepare/               source-specific preparation (201-212)
  30_indicators/            indicator_register.csv: one row per indicator
  40_components/            component scores (in development)
  50_risk/                  risk scores (in development)
  60_validation/            validation (in development)
  70_visualisation/         701: interactive Leaflet map
  80_export/                exports
  functions/                shared functions loaded by tar_source()
docs/                       published map (GitHub Pages)
```

**Pipeline stages:**

- **10 download:** fetches the raw data.
- **20 prepare:** reads, crops, reprojects and derives one physical variable per source, on the source's own grid.
- **30 indicators:** a generic function resamples every layer to the 1 km template grid (SA Albers Equal Area), masks it to the catchment ∪ Alfred Nzo, and normalises it to 0–1.
- **70 visualisation:** builds the interactive map.

## Reproducing the analysis

The pipeline uses [{targets}](https://docs.ropensci.org/targets/) and [{renv}](https://rstudio.github.io/renv/).

```r
renv::restore()
targets::tar_visnetwork()
targets::tar_make()
```

The data (about 2 GB of rasters, NetCDF and shapefiles) are **not** in this repository. Only small provenance files are tracked: manifests, lookups and `.url` files. Run the scripts in `R/10_download/` to fetch the open datasets. Some inputs are not openly downloadable:

- Census tables via DataFirst, which need a login.
- DWS settlement data, which are available on request from DWS.

## Data use and citation

Each data source keeps its own licence. See the map's "About & data sources" section, and cite the original sources when you use derived layers. Basemaps are © CARTO / OpenStreetMap contributors and © Esri and its data providers.

## Licence

Code: MIT (see `LICENSE`). Third-party data remain under their original licences.
