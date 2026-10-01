# Shannon diversity, computed one single way everywhere. The old
# Analyses_PollenVSTrees.r computed it by hand on raw taxon labels
# (p <- table(x) / length(x); -sum(p * log(p))); that is the exact same
# quantity as vegan::diversity() on a pre-tabulated count vector — table(x)/
# length(x) IS the proportion vector diversity() computes internally — just
# expressed from a different input shape. Confirmed numerically identical on
# the tree diversity data before migrating Analyses_PollenVSTrees.r over to
# this function (see the R/pollen_vs_arbres migration notes).

library(vegan)

# One Shannon value per row of a sample/station x taxon count matrix.
shannon <- function(mat, index = "shannon") {
  diversity(mat, index = index)
}
