# Single entry point for every analysis script: source("R/init.R") loads all
# modules, in dependency order, so a script never has to source them one by
# one or worry about ordering.

source("R/config.R")     # constants — no dependencies
source("R/geo.R")        # uses CRS_WGS84 / CRS_MTM from config.R
source("R/matrices.R")   # true-zero grids and station x taxon matrices
source("R/diversity.R")  # uses matrices.R output
source("R/pollen.R")     # uses matrices.R + config.R (CONF_MIN, DIR_DATA_POLLEN)
source("R/trees.R")      # uses matrices.R + config.R (RAYONS, DIR_DATA_TREES)
source("R/models.R")     # glmmTMB validation + coefficient extraction
source("R/plots.R")      # uses config.R (palettes, FLECHE_SCALE, ...)
source("R/report.R")     # markdown/Word writers, uses models.R's dependencies
