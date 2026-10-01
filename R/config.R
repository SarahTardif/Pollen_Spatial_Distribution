# Constants shared across scripts: thresholds, plot radii, excluded plots, CRS,
# paths, palettes, and the graphical constants that used to be copy-pasted
# (and sometimes diverged) across files. Nothing here should change any model
# result — only figure appearance depends on some of these values.

## Data cleaning
CONF_MIN <- 0.8   # Confidence threshold applied once, in 01_preparer_pollen.R
TAXATORM<- c("Gramineae", "Ambrosia")

## Plot radii (tree analyses)
RAYONS    <- c(50, 100, 200)   # the 3 radii every tree analysis is repeated for
RAYON_REF <- 200               # radius used for figures when only one is shown

## Excluded / special plots — named separately because the two scripts that use
## them do NOT use the same set (Analyses_PollenVSTrees* excludes only 20C;
## exploratory_results.r excludes the 3 miniforest plots). Keeping two names
## instead of one avoids silently changing which plots a script drops.
PLOTS_EXCLUS_ARBRES <- c()
PLOTS_MINIFORET     <- c("19B", "20C", "23A")

## CRS
CRS_WGS84 <- 4326   # lat/lon, as stored in locations_traps.csv
CRS_MTM   <- 32188  # MTM zone 8 (Quebec), metres — used for real-world distances

## Paths
DIR_DATA_POLLEN <- "./Data/Data_pollen"
DIR_DATA_TREES  <- "./Data/Data_utiles/Data_trees"
DIR_OUT_ARBRES  <- "./Outputs/pollen_vs_arbres"
DIR_FIG_ARBRES  <- "./Outputs/pollen_vs_arbres/figures"
DIR_OUT_ARTICLE <- "./Outputs/article"
DIR_FIG_ARTICLE <- "./Outputs/article/figures"
DIR_OUT_ARTICLE_PERIODES <- "./Outputs/article_periodes"
DIR_FIG_ARTICLE_PERIODES <- "./Outputs/article_periodes/figures"
DIR_OUT_ARTICLE_EFFETS_FIXES <- "./Outputs/article_effets_fixes"
DIR_OUT_POLLINISATION <- "./Outputs/pollinisation_mode"
DIR_FIG_POLLINISATION <- "./Outputs/pollinisation_mode/figures"

## Taxon restriction (interaction models, envfit arrows, palettes)
N_TOP_TAXONS <- 10   # size of the "top" taxon subset; recomputed from whatever
                     # matrix/matrices are in scope, see top_taxons() in R/matrices.R

## Colour palettes
# Pollen taxa: 10 distinct colors + grey for "Others". Used by both the article
# pipeline and the exploratory script, which used to keep two identical copies.
COULEURS_TAXONS <- c("#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f",
                     "#edc948", "#b07aa1", "#ff9da7", "#9c755f", "#bab0ac")

# Tree genera: a different, colour-blind-friendly set so tree and pollen
# figures are never visually confused with each other.
COULEURS_ARBRES <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
                     "#e6ab02", "#a6761d", "#666666", "#a6cee3", "#1f78b4")

COULEUR_AUTRES <- "grey80"   # colour for the "Others"/"Other" bucket, both palettes

# Pollination mode / flower sex (24_pollinisation_mode.R): "unknown" reuses
# COULEUR_AUTRES so an unclassified genus reads the same way "Others" does
# elsewhere, instead of introducing a third grey.
COULEURS_POLLINISATION <- c(wind = "#f28e2b", animals = "#59a14f", unknown = COULEUR_AUTRES)
COULEURS_SEXE_FLEUR    <- c(monoecious = "#b07aa1", dioecious = "#edc948")

# Canopy cover classes (map_25plots.r's station map, and fig2's PCoA point
# colour): 4 classes, thresholds rounded at 10/20/30 % (25 stations, min
# 8.6 %, max 73.1 %). Shared so the two figures can never silently diverge.
CANOPY_BREAKS <- c(-Inf, 10, 20, 30, Inf)
CANOPY_LABELS <- c("less than 10", "from 10 to 19.9", "from 20 to 29.9", "30 and more")
COULEURS_CANOPY <- setNames(c("white", "#a1d99b", "#41ab5d", "#00441b"), CANOPY_LABELS)

# Years: Analyses_PollenVSTrees.r used firebrick/steelblue, Pipeline_article.r
# used #4e79a7/#e15759 — same information, two different hex sources. Unified
# to the article's palette (cosmetic only, does not change any model result).
COULEURS_ANNEES <- c("2022" = "#4e79a7", "2023" = "#e15759")

## Graphical constants for ordination plots (envfit arrows, pie charts)
FLECHE_SCALE  <- 0.65   # arrow length, as a fraction of the plot half-range
LABEL_OFFSET  <- 1.12   # taxon label placed slightly beyond the arrow tip
PIE_R         <- 0.04   # scatterpie radius, as a fraction of the plot range
NMDS_TRYMAX   <- 100    # metaMDS random restarts

## Reproducibility
GRAINE <- 42    # set.seed() value used everywhere a permutation/ordination is run
N_PERM <- 999   # permutations for adonis2 / envfit / mantel
