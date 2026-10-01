# True-zero grids and station/sample x taxon matrices.
#
# There are TWO different zero-fills in this codebase and both live here on
# purpose:
#   1. completer_zeros() fills the "long format" grid used to feed models —
#      merging the full set of (sample x category) combinations against what
#      was actually observed, so "this taxon/genus is absent from a sample
#      that was collected" (a true zero) stays distinct from "this sample was
#      never collected" (a missing row, no zero inserted). This is the pattern
#      from Pipeline_article.r's ab_taxon grid and from the old
#      Analyses_PollenVSTrees_par_genre.r's G_du_genre()/comptages_du_genre().
#   2. mat_station_taxon() / mat_echantillon_taxon() use acast(..., fill = 0),
#      a purely structural fill: once the long-format data already has one row
#      per (sample, category) — true zeros included via (1) — acast just needs
#      *some* value for combinations that still don't appear (there aren't
#      normally any left, it's a safety net, not where the "true zero" decision
#      is made).
# Do not "simplify" these into one function: (1) encodes a modelling decision
# (which absences are real), (2) is bookkeeping for a wide matrix.

library(reshape2)

# Left-joins `valeurs` onto the full `grille` on `by`, filling missing rows'
# `colonne_valeur` with `defaut` (0). `grille` is the set of rows that SHOULD
# exist (e.g. every collected sample x every taxon, or every plot for one
# genus); `valeurs` is what was actually observed.
completer_zeros <- function(grille, valeurs, by, colonne_valeur, defaut = 0) {
  df <- merge(grille, valeurs, by = by, all.x = TRUE)
  df[[colonne_valeur]][is.na(df[[colonne_valeur]])] <- defaut
  df
}

# Station (or Plot) x taxon count matrix — one row per level of `station_col`,
# meant to be called on data already aggregated to that grain (e.g. one year
# at a time). fill = 0 here is the structural kind, see file header.
mat_station_taxon <- function(df, station_col, taxon_col, valeur_col = "count") {
  formule <- stats::as.formula(paste(station_col, "~", taxon_col))
  acast(df, formule, value.var = valeur_col, fill = 0)
}

# Sample-level matrix: rows are one (station|period|year) sample (or any
# combination of `id_cols`, "|"-joined), columns are taxa. Used for Shannon
# diversity and for PERMANOVA/ordination on individual samples rather than
# annual totals. Drops rows that sum to zero (e.g. a sample containing only
# OTHER, once OTHER has been excluded upstream) unless drop_vides = FALSE.
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

# Names of the n taxa with the highest total abundance, pooled across
# whichever matrices are passed in. Callers pass exactly the matrix/matrices
# that feed the model or figure at hand (one year, one period, or several
# pooled together) so the "top" subset always matches what's actually being
# modeled/plotted, instead of one global list reused regardless of context.
top_taxons <- function(..., n = N_TOP_TAXONS) {
  totaux <- colSums(do.call(rbind, list(...)))
  names(sort(totaux, decreasing = TRUE))[1:n]
}
