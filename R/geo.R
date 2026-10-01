# Coordonnees des stations, projection en metres et matrices de distances.
# Toutes les distances reelles passent par la meme projection (EPSG 32188).

library(sf)

# Metadonnees des stations (piege, adresse, arrondissement, lon/lat).
charger_stations <- function(fichier = "./Data/locations_traps.csv") {
  read.csv(fichier, sep = ";", header = TRUE)
}

# Projette un data.frame lon/lat en MTM (metres), en objet sf.
en_mtm <- function(df, coords = c("Longitude", "Latitude"), crs_origine = CRS_WGS84) {
  pts <- st_as_sf(df, coords = coords, crs = crs_origine)
  st_transform(pts, crs = CRS_MTM)
}

# Matrice des distances (metres) entre les stations donnees, en objet dist.
distances_stations <- function(locs_sf, stations, id_col = "trap") {
  idx <- match(stations, locs_sf[[id_col]])
  as.dist(matrix(as.numeric(st_distance(locs_sf[idx, ])),
                 nrow = length(idx), ncol = length(idx),
                 dimnames = list(stations, stations)))
}

# Gradients socio-environnementaux par station (NDVI, couvert, densite de
# population). L'entete porte un BOM UTF-8, d'ou le renommage de la 1re colonne.
charger_gradients <- function(fichier = "./Data/Gradients_plots.csv") {
  gradients <- read.csv(fichier, header = TRUE)
  names(gradients)[1] <- "Plot"
  gradients$Plot <- as.character(gradients$Plot)
  gradients
}

# Couvert forestier (%) -> facteur a 4 classes (seuils dans R/config.R).
classer_canopy <- function(canopy_pct) {
  cut(canopy_pct, breaks = CANOPY_BREAKS, labels = CANOPY_LABELS, right = FALSE)
}
