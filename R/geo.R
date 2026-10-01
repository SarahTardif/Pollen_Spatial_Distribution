# Station coordinates, projection in metres and distance matrices.
# Every real-world distance goes through the same projection (EPSG 32188).

library(sf)

# Station metadata (trap, address, borough, lon/lat).
charger_stations <- function(fichier = "./Data/locations_traps.csv") {
  read.csv(fichier, sep = ";", header = TRUE)
}

# Projects a lon/lat data.frame to MTM (metres), as an sf object.
en_mtm <- function(df, coords = c("Longitude", "Latitude"), crs_origine = CRS_WGS84) {
  pts <- st_as_sf(df, coords = coords, crs = crs_origine)
  st_transform(pts, crs = CRS_MTM)
}

# Distance matrix (metres) between the given stations, as a dist object.
distances_stations <- function(locs_sf, stations, id_col = "trap") {
  idx <- match(stations, locs_sf[[id_col]])
  as.dist(matrix(as.numeric(st_distance(locs_sf[idx, ])),
                 nrow = length(idx), ncol = length(idx),
                 dimnames = list(stations, stations)))
}

# Socio-environmental gradients per station (NDVI, canopy cover, population density)
charger_gradients <- function(fichier = "./Data/Gradients_plots.csv") {
  gradients <- read.csv(fichier, header = TRUE)
  names(gradients)[1] <- "Plot"
  gradients$Plot <- as.character(gradients$Plot)
  gradients
}

# Canopy cover (%) -> 4-class factor (thresholds in R/config.R).
classer_canopy <- function(canopy_pct) {
  cut(canopy_pct, breaks = CANOPY_BREAKS, labels = CANOPY_LABELS, right = FALSE)
}
