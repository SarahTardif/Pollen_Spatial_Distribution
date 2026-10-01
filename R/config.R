# Constantes partagees par les scripts: seuils, CRS, chemins, palettes et
# constantes graphiques.

## Nettoyage des donnees
CONF_MIN <- 0.8   # seuil de confiance, applique dans 01_preparer_pollen.R
TAXATORM <- c("Gramineae", "Ambrosia")   # taxons non arborescents, retires

## CRS
CRS_WGS84 <- 4326   # lat/lon, format de locations_traps.csv
CRS_MTM   <- 32188  # MTM zone 8 (Quebec), en metres, pour les distances reelles

## Chemins
DIR_DATA_POLLEN <- "./Data/Data_pollen"
DIR_OUT_ARTICLE <- "./Outputs/article"
DIR_FIG_ARTICLE <- "./Outputs/article/figures"

## Sous-ensemble de taxons (palettes, figures)
N_TOP_TAXONS <- 10   # taille du sous-ensemble "top", voir top_taxons() dans R/matrices.R

## Palettes
# Taxons polliniques: 10 couleurs distinctes + un gris pour "Others".
COULEURS_TAXONS <- c("#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f",
                     "#edc948", "#b07aa1", "#ff9da7", "#9c755f", "#bab0ac")

COULEUR_AUTRES <- "grey80"   # couleur du groupe "Others"

# Classes de couvert forestier, partagees par la carte des stations (02) et la
# PCoA (04) pour que les deux figures ne puissent pas diverger.
CANOPY_BREAKS <- c(-Inf, 10, 20, 30, Inf)
CANOPY_LABELS <- c("less than 10", "from 10 to 19.9", "from 20 to 29.9", "30 and more")
COULEURS_CANOPY <- setNames(c("white", "#a1d99b", "#41ab5d", "#00441b"), CANOPY_LABELS)

COULEURS_ANNEES <- c("2022" = "#4e79a7", "2023" = "#e15759")

## Constantes graphiques des ordinations
FLECHE_SCALE <- 0.65   # longueur des fleches envfit, en fraction du demi-domaine
LABEL_OFFSET <- 1.12   # etiquette de taxon placee juste au-dela de la pointe

## Reproductibilite
GRAINE <- 42    # set.seed() pour toute permutation ou ordination
N_PERM <- 999   # permutations pour adonis2 / envfit / mantel
