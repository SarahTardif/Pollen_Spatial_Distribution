# Single entry point for the analysis scripts: source("R/init.R") loads every
# module in dependency order.

source("R/config.R")     # constants
source("R/geo.R")        # coordinates, projections, distances
source("R/matrices.R")   # true zeros and station x taxon matrices
source("R/diversity.R")  # Shannon index
source("R/pollen.R")     # pollen loading and count grid
source("R/models.R")     # glmmTMB model validation
source("R/plots.R")      # shared palettes and figures
source("R/report.R")     # tables and markdown/Word writers
