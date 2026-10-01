# Pollen spatial distribution

R scripts to reproduce results from the study of the spatial distribution of airborne pollen across 25 sampling stations in Montreal.

**Publication:** [will be added upon publication]


## Scripts

| File | Description |
| --- | --- |
| `01_preparer_pollen.R` | Import the yearly parquet predictions, derive period/year/location, filter on classification confidence, write the cleaned CSVs |
| `02_carte_stations.R` | Map of the 25 sampling stations, by canopy cover and population density |
| `03_climat_montreal.R` | Select weather stations, download and merge daily climate data, plot climate trends for 2022 and 2023 |
| `04_article_analyses.R` | Main analysis: abundance models, PERMANOVA, PCoA, Shannon diversity, PAM clustering, Mantel tests, all figures |
| `05_article_tableaux.R` | Format the fitted models into publication-ready tables (markdown + Word) |

## Shared modules

Loaded all at once by `source("R/init.R")`, which every script except `03_climat_montreal.R` calls first.

| File | Description |
| --- | --- |
| `R/config.R` | Thresholds, CRS, paths, palettes, random seed and permutation count |
| `R/geo.R` | Station coordinates, MTM projection, distance matrices, canopy classes |
| `R/matrices.R` | True-zero grids and station/sample × taxon matrices |
| `R/diversity.R` | Shannon index |
| `R/pollen.R` | Pollen loading and the station × period × year × taxon count grid |
| `R/models.R` | glmmTMB validation bundle (type II Anova, R², DHARMa) |
| `R/plots.R` | Taxon palettes, envfit arrows, ordination plots, heatmaps, figure saving |
| `R/report.R` | Table builders and markdown/Word writers |

## Workflow

`01_preparer_pollen.R` → `04_article_analyses.R` → `05_article_tableaux.R`

`05_article_tableaux.R` reads the objects left in the environment by `04_article_analyses.R`, so run them in the same R session; otherwise it sources `04` itself. `02_carte_stations.R` and `03_climat_montreal.R` are independent and can be run at any time.

Outputs are written to `Outputs/article/` (figures and tables) and `Outputs/climatedata/`.

## Data

The `Data/` folder is not included in this repository: several of the pollen files exceed GitHub's 100 MB per-file limit. The scripts expect the following layout:

```
Data/
├── locations_traps.csv                                   # station ids, addresses, lon/lat (";" separated)
├── Gradients_plots.csv                                   # per-station canopy cover, population density, NDVI
├── LIM_ADMIN/                                            # Montreal administrative and land-extent shapefiles
└── Data_pollen/
    └── calibrated_predictions_{2021,2022,2023}.parquet   # raw classifier output, read by script 01
```

`01_preparer_pollen.R` writes `Data/Data_pollen/data{2021,2022,2023}.csv` from the parquet files; the other scripts read those CSVs.

## Requirements

R (≥ 4.x) with arrow, car, cluster, DHARMa, dplyr, emmeans, flextable, ggplot2, ggrepel, glmmTMB, multcomp, multcompView, officer, performance, reshape2, scales, scatterpie, sf, tmap, vegan, weathercan.

Run the scripts from the repository root, so the relative paths in `R/config.R` resolve.

## License

GNU General Public License v3.0 — see [LICENSE](LICENSE).
