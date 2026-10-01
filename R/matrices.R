# True zeros and station/sample x taxon matrices.
#
# Two zero-fills coexist on purpose:
#   1. completer_zeros() fills the long-format grid that feeds the models,
#      keeping "taxon absent from a collected sample" (a true zero) distinct
#      from "sample never collected" (a missing row, no zero inserted).
#   2. acast(..., fill = 0) in both matrix functions is a purely structural
#      fill, applied once the true zeros are already in place.

library(reshape2)

# Joins `valeurs` onto the full grid `grille`, missing rows getting `defaut`
# (0) in `colonne_valeur`.
completer_zeros <- function(grille, valeurs, by, colonne_valeur, defaut = 0) {
  df <- merge(grille, valeurs, by = by, all.x = TRUE)
  df[[colonne_valeur]][is.na(df[[colonne_valeur]])] <- defaut
  df
}

# Station x taxon matrix, to be called on data already aggregated to that grain
# (one year at a time, for instance).
mat_station_taxon <- function(df, station_col, taxon_col, valeur_col = "count") {
  formule <- stats::as.formula(paste(station_col, "~", taxon_col))
  acast(df, formule, value.var = valeur_col, fill = 0)
}

# Sample-level matrix: one row per combination of `id_cols` (joined by "|"),
# one column per taxon. Rows summing to zero are dropped unless
# drop_vides = FALSE.
mat_echantillon_taxon <- function(df, id_cols, taxon_col, valeur_col = "count",
                                  sep = "|", drop_vides = TRUE) {
  row_id <- do.call(paste, c(df[id_cols], sep = sep))
  long   <- data.frame(row_id = row_id, taxon = df[[taxon_col]], valeur = df[[valeur_col]])
  mat    <- acast(long, row_id ~ taxon, value.var = "valeur", fill = 0)
  if (drop_vides) mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  mat
}

# Relative-abundance version of a count matrix (rows sum to 1).
en_relatif <- function(mat) mat / rowSums(mat)

# Names of the n most abundant taxa, pooled across whichever matrices are
# passed in. The caller passes exactly the matrix/matrices feeding the model or
# figure at hand, so the "top" always matches what is actually shown.
top_taxons <- function(..., n = N_TOP_TAXONS) {
  totaux <- colSums(do.call(rbind, list(...)))
  names(sort(totaux, decreasing = TRUE))[1:n]
}
