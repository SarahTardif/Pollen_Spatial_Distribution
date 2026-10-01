# Constants shared across scripts: thresholds, CRS, paths, palettes and
# graphical constants.

## Data cleaning
CONF_MIN <- 0.8   # confidence threshold, applied in 01_preparer_pollen.R
TAXATORM <- c("Gramineae", "Ambrosia")   # non-tree taxa, removed

## CRS
CRS_WGS84 <- 4326   # lat/lon, as stored in locations_traps.csv
CRS_MTM   <- 32188  # MTM zone 8, in metres, for real-world distances

## Paths
DIR_DATA_POLLEN <- "./Data/Data_pollen"
DIR_OUT_ARTICLE <- "./Outputs/article"
DIR_FIG_ARTICLE <- "./Outputs/article/figures"

## Taxon subset (palettes, figures)
N_TOP_TAXONS <- 10 # see top_taxons() in R/matrices.R

## Palettes
# Pollen taxa: 10 distinct colours + a grey for "Others".
COULEURS_TAXONS <- c("#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f",
                     "#edc948", "#b07aa1", "#ff9da7", "#9c755f", "#bab0ac")

COULEUR_AUTRES <- "grey80"   # colour of the "Others" group

# Canopy cover classes, shared by the station map (02) and the PCoA (04) 
CANOPY_BREAKS <- c(-Inf, 10, 20, 30, Inf)
CANOPY_LABELS <- c("less than 10", "from 10 to 19.9", "from 20 to 29.9", "30 and more")
COULEURS_CANOPY <- setNames(c("white", "#a1d99b", "#41ab5d", "#00441b"), CANOPY_LABELS)

COULEURS_ANNEES <- c("2022" = "#4e79a7", "2023" = "#e15759")

## Graphical constants for ordination plots
FLECHE_SCALE <- 0.65   # envfit arrow length, as a fraction of the half-range
LABEL_OFFSET <- 1.12   # taxon label placed just beyond the arrow tip

## Reproducibility
GRAINE <- 42    # set.seed() for every permutation or ordination
N_PERM <- 999   # permutations for adonis2 / envfit / mantel
