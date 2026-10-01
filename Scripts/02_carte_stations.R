# Carte des 25 stations d'echantillonnage (fig0 de l'article): couvert
# forestier en couleur, densite de population en taille de symbole.

source("R/init.R")
library(sf)
library(tmap)

lim_map   <- st_read("./Data/LIM_ADMIN/limites_maps.shp")
lim_admin <- st_read("./Data/LIM_ADMIN/LIM_ADMIN.shp")

locations <- charger_stations()
gradients <- charger_gradients()
locations <- merge(locations, gradients, by.x = "trap", by.y = "Plot")

# Couvert forestier et densite de population en 4 classes chacun. Les seuils de
# couvert viennent de classer_canopy() (R/geo.R), partages avec la PCoA du
# script 04 pour que les deux figures ne puissent pas diverger.
locations$canopy_classe <- classer_canopy(locations$Canopy.cover....)
locations$densite_classe <- cut(locations$Population.density..people.km.2.,
                                breaks = c(-Inf, 3000, 5000, 11000, Inf),
                                labels = c("less than 3000", " from 3000 to 4999",
                                           "from 5000 to 10999", "11000 and more"),
                                right = FALSE)

locations_points <- st_as_sf(locations, coords = c("Longitude", "Latitude"), crs = CRS_WGS84)
locations_points <- st_transform(locations_points, crs = CRS_MTM)   # meme projection que lim_map

# Les polygones d'arrondissement debordent sur le fleuve la ou la limite
# administrative suit le milieu du cours d'eau; on decoupe sur l'emprise
# terrestre pour que les traits pointilles ne traversent plus l'eau.
lim_admin <- st_make_valid(lim_admin)
lim_map   <- st_make_valid(lim_map)
lim_admin_terre <- st_intersection(lim_admin, st_union(lim_map))

# Emprise zoomee sur le centre de l'ile (bbox des stations + marge), plutot que
# l'emprise complete de lim_map qui deborde sur la Rive-Sud et Laval.
marge <- 4000 # metres
bbox_zoom <- st_bbox(locations_points)
bbox_zoom["xmin"] <- bbox_zoom["xmin"] - marge
bbox_zoom["xmax"] <- bbox_zoom["xmax"] + marge
bbox_zoom["ymin"] <- bbox_zoom["ymin"] - marge
bbox_zoom["ymax"] <- bbox_zoom["ymax"] + marge

# Etiquette en diagonale haut-droite du cercle; seule 12C, collee a 12B, passe
# a gauche pour ne pas se chevaucher.
locations_points$label_xmod <- 0.7
locations_points$label_ymod <- 1
locations_points$label_xmod[locations_points$trap == "12C"] <- -0.7

# Fleche nord et echelle: fonctions tmap natives, les equivalents ggspatial ne
# s'appliquent qu'aux objets ggplot2.
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

# Meme convention que les autres figures (300 dpi, Outputs/article/figures),
# mais via tmap_save(): sauver_figure() appelle ggsave(), qui ne fonctionne pas
# sur un objet tmap.
dir.create(DIR_FIG_ARTICLE, recursive = TRUE, showWarnings = FALSE)

tmap_save(map_samplers, file.path(DIR_FIG_ARTICLE, "fig0_carte_stations.png"),
          width = 8, height = 8, units = "in", dpi = 300)
