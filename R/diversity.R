# Indice de Shannon, calcule d'une seule facon dans tout le projet.

library(vegan)

# Une valeur de Shannon par ligne d'une matrice de comptages taxon x echantillon.
shannon <- function(mat, index = "shannon") {
  diversity(mat, index = index)
}
