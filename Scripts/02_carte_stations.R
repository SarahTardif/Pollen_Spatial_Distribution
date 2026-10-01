# Map of the 25 sampling stations (fig_carte_stations): canopy cover as colour, population density as symbol size.

source("R/init.R")
library(sf)
library(tmap)

lim_map   <- st_read("./Data/LIM_ADMIN/limites_maps.shp")
lim_admin <- st_read("./Data/LIM_ADMIN/LIM_ADMIN.shp")

locations <- charger_stations()
gradients <- charger_gradients()
locations <- merge(locations, gradients, by.x = "trap", by.y = "Plot")

# Canopy cover and population density, 4 classes each. The canopy thresholds
# come from classer_canopy() (R/geo.R), shared with the PCoA of script 04 so
# the two figures can never diverge.
locations$canopy_classe <- classer_canopy(locations$Canopy.cover....)
locations$densite_classe <- cut(locations$Population.density..people.km.2.,
                                breaks = c(-Inf, 3000, 5000, 11000, Inf),
                                labels = c("less than 3000", " from 3000 to 4999",
                                           "from 5000 to 10999", "11000 and more"),
                                right = FALSE)

locations_points <- st_as_sf(locations, coords = c("Longitude", "Latitude"), crs = CRS_WGS84)
locations_points <- st_transform(locations_points, crs = CRS_MTM)   # same projection as lim_map

lim_admin <- st_make_valid(lim_admin)
lim_map   <- st_make_valid(lim_map)
lim_admin_terre <- st_intersection(lim_admin, st_union(lim_map))

# Zoom on the centre of the island (stations bbox + margin)
marge <- 4000 # metres
bbox_zoom <- st_bbox(locations_points)
bbox_zoom["xmin"] <- bbox_zoom["xmin"] - marge
bbox_zoom["xmax"] <- bbox_zoom["xmax"] + marge
bbox_zoom["ymin"] <- bbox_zoom["ymin"] - marge
bbox_zoom["ymax"] <- bbox_zoom["ymax"] + marge

# Label adjustments
locations_points$label_xmod <- 0.7
locations_points$label_ymod <- 1
locations_points$label_xmod[locations_points$trap == "12C"] <- -0.7

# North arrow and scale bar
map_samplers <- tm_shape(lim_map, bbox = bbox_zoom) + tm_fill() + tm_borders() +
  tm_shape(lim_admin_terre) + tm_borders(col = "grey60") +
  tm_shape(locations_points) +
  tm_symbols(fill = "canopy_classe",
             fill.scale = tm_scale_categorical(values = COULEURS_CANOPY),
             fill.legend = tm_legend(title = "Canopy cover (%)", text.size = 0.8, title.size = 0.8,
                                     item.height = 0.6, item.width = 0.6,
                                     position = tm_pos_in("left", "top")),
             size = "densite_classe",
             size.scale = tm_scale_categorical(values = c(0.6, 1.2, 1.8, 2.4)),
             size.legend = tm_legend(title = "Population density(hab/km²)", text.size = 0.8, title.size = 0.8,
                                     item.height = 0.6, item.width = 0.6,
                                     position = tm_pos_in("left", "top"))) +
  tm_text("trap", size = 0.8, xmod = "label_xmod", ymod = "label_ymod",
          bgcol = "white", bgcol_alpha = 0.7) +
  tm_compass(position = c("right", "bottom")) +
  tm_scalebar(position = c("right", "bottom"))
print(map_samplers)

# Same convention as the other figures (300 dpi, Outputs/article/figures), but through tmap_save()
dir.create(DIR_FIG_ARTICLE, recursive = TRUE, showWarnings = FALSE)

tmap_save(map_samplers, file.path(DIR_FIG_ARTICLE, "fig_carte_stations.png"),
          width = 8, height = 8, units = "in", dpi = 300)
