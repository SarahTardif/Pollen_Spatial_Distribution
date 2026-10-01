# Geographic helpers: station coordinates, MTM (metres) projection, distance
# matrices. Every real-world distance in the codebase goes through the same
# projection (EPSG 32188) so results stay identical across scripts.

library(sf)

# Station metadata (trap id, address, borough, lon/lat). ";"-separated.
charger_stations <- function(fichier = "./Data/locations_traps.csv") {
  read.csv(fichier, sep = ";", header = TRUE)
}

# Projects a data.frame with lon/lat columns to MTM (metres), as an sf object.
en_mtm <- function(df, coords = c("Longitude", "Latitude"), crs_origine = CRS_WGS84) {
  pts <- st_as_sf(df, coords = coords, crs = crs_origine)
  st_transform(pts, crs = CRS_MTM)
}

# Full pairwise distance matrix (metres) between the given stations, as a
# dist object with station ids as dimnames. as.numeric() unwraps the units
# object returned by st_distance, which otherwise breaks downstream comparisons.
distances_stations <- function(locs_sf, stations, id_col = "trap") {
  idx <- match(stations, locs_sf[[id_col]])
  as.dist(matrix(as.numeric(st_distance(locs_sf[idx, ])),
                 nrow = length(idx), ncol = length(idx),
                 dimnames = list(stations, stations)))
}

# Socio-environmental gradients per station (NDVI, canopy cover, population
# density, etc.), one row per station ("Plot"). The header carries a UTF-8
# BOM, which read.csv turns into "X.U.FEFF.Plot" on the first column name; it
# is renamed back rather than relying on the encoding being handled upstream.
charger_gradients <- function(fichier = "./Data/Gradients_plots.csv") {
  gradients <- read.csv(fichier, header = TRUE)
  names(gradients)[1] <- "Plot"
  gradients$Plot <- as.character(gradients$Plot)
  gradients
}

# Canopy cover (%) -> 4-class factor, using the shared breaks/labels
# (CANOPY_BREAKS/CANOPY_LABELS, R/config.R) so the station map
# (map_25plots.r) and any other figure colouring stations by this class
# (e.g. fig2's PCoA) can never silently diverge.
classer_canopy <- function(canopy_pct) {
  cut(canopy_pct, breaks = CANOPY_BREAKS, labels = CANOPY_LABELS, right = FALSE)
}

# Point-to-point euclidean distance (metres) between two sets of projected
# coordinates, matched by idx — e.g. each tree's distance to its plot's trap.
# coords_obj / coords_ref are matrices as returned by sf::st_coordinates().
distance_au_piege <- function(coords_obj, coords_ref, idx) {
  sqrt((coords_obj[, 1] - coords_ref[idx, 1])^2 +
       (coords_obj[, 2] - coords_ref[idx, 2])^2)
}
