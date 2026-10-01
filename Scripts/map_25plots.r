#préparation des cartes
source("R/init.R") # pour sauver_figure() et DIR_FIG_ARTICLE (R/plots.R, R/config.R)
library(sf)
library(mapview)
library(tmap)
lim_map<-st_read("./Data/Data_utiles/LIM_ADMIN/limites_maps.shp")
lim_admin<-st_read("./Data/Data_utiles/LIM_ADMIN/LIM_ADMIN.shp")
locations<-read.csv("./Data/locations_traps.csv", sep=";", header=T)
gradients<-read.csv("./Data/Gradients_plots.csv", header=T)
locations<-merge(locations, gradients, by.x="trap", by.y="Plot")

# Canopy cover et densite de population regroupes en 4 classes chacun; seuils
# de canopee partages avec fig2 (classer_canopy(), R/geo.R, seuils dans
# R/config.R) pour que la carte et la PCoA ne puissent pas diverger.
locations$canopy_classe<-classer_canopy(locations$Canopy.cover....)
locations$densite_classe<-cut(locations$Population.density..people.km.2.,
                               breaks=c(-Inf, 3000, 5000, 11000, Inf),
                               labels=c("less than 3000", " from 3000 to 4999", "from 5000 to 10999", "11000 and more"),
                               right=FALSE)

locations_points<-st_as_sf(locations, coords= c("Longitude","Latitude"),crs=4326) # EPSG wgs84 google maps
locations_points<-st_transform(locations_points, crs=32188) #projection MTM comme lim_map

# les polygones d'arrondissement (lim_admin) debordent sur le fleuve/les plans
# d'eau la ou la limite administrative reelle suit le milieu du cours d'eau;
# on decoupe sur l'emprise terrestre (lim_map) pour que les traits pointilles
# ne traversent plus l'eau sur la carte
lim_admin<-st_make_valid(lim_admin)
lim_map<-st_make_valid(lim_map)
lim_admin_terre<-st_intersection(lim_admin, st_union(lim_map))

# emprise zoomee sur le centre de l'ile (bbox des stations + marge), plutot que
# l'emprise complete de lim_map qui deborde largement sur la Rive-Sud et Laval
marge<-4000 # metres
bbox_stations<-st_bbox(locations_points)
bbox_zoom<-bbox_stations
bbox_zoom["xmin"]<-bbox_zoom["xmin"]-marge
bbox_zoom["xmax"]<-bbox_zoom["xmax"]+marge
bbox_zoom["ymin"]<-bbox_zoom["ymin"]-marge
bbox_zoom["ymax"]<-bbox_zoom["ymax"]+marge

# etiquette placee en diagonale haut-droite du cercle par defaut (proche du
# cercle, contrairement au calcul par distance essaye plus tot qui envoyait
# les etiquettes trop loin); seules 12B et 12C, colles l'un contre l'autre,
# passent en haut-gauche/haut-droite opposes pour ne pas se chevaucher
locations_points$label_xmod<-0.7
locations_points$label_ymod<-1
locations_points$label_xmod[locations_points$trap=="12C"]<- -0.7

# fleche nord et echelle : fonctions tmap natives (tm_compass/tm_scalebar), pas
# les fonctions ggspatial (annotation_scale/annotation_north_arrow) qui ne
# s'appliquent qu'aux objets ggplot2, pas aux objets tmap
map_samplers<-tm_shape(lim_map, bbox=bbox_zoom)+ tm_fill()+ tm_borders()+
  tm_shape(lim_admin_terre)+ tm_borders(col="grey60")+
  tm_shape(locations_points)+
  tm_symbols(fill="canopy_classe",
             fill.scale=tm_scale_categorical(values=COULEURS_CANOPY),
             fill.legend=tm_legend(title="Canopy cover (%)", text.size=0.8, title.size=0.8,
                                    item.height=0.6, item.width=0.6,
                                    position=tm_pos_in("left", "top")),
             size="densite_classe",
             size.scale=tm_scale_categorical(values=c(0.6, 1.2, 1.8, 2.4)),
             size.legend=tm_legend(title="Population density(hab/km²)", text.size=0.8, title.size=0.8,
                                    item.height=0.6, item.width=0.6,
                                    position=tm_pos_in("left", "top")))+
  tm_text("trap", size=0.8, xmod="label_xmod", ymod="label_ymod", bgcol="white",bgcol_alpha=0.7)+
  tm_compass(position=c("right", "bottom"))+
  tm_scalebar(position=c("right", "bottom"))
print(map_samplers)

# meme convention que les autres figures de l'article (dpi = 300, dossier
# Outputs/article/figures via sauver_figure(), R/plots.R); ggsave() ne
# fonctionne pas sur un objet tmap donc on utilise tmap_save() directement
# plutot que sauver_figure() (qui appelle ggsave en interne)
dir.create(DIR_FIG_ARTICLE, recursive = TRUE, showWarnings = FALSE)

tmap_save(map_samplers, file.path(DIR_FIG_ARTICLE, "fig0_carte_stations.png"),
          width = 8, height = 8, units = "in", dpi = 300)
