# Shannon index, computed one single way throughout the project.

library(vegan)

# One Shannon value per row of a sample x taxon count matrix.
shannon <- function(mat, index = "shannon") {
  diversity(mat, index = index)
}
