# Vrais zeros et matrices station/echantillon x taxon.
#
# Deux remplissages de zeros coexistent volontairement:
#   1. completer_zeros() remplit la grille en format long qui alimente les
#      modeles, en distinguant "taxon absent d'un echantillon collecte" (vrai
#      zero) de "echantillon jamais collecte" (ligne absente, aucun zero).
#   2. acast(..., fill = 0) dans les deux fonctions de matrice est un simple
#      remplissage structurel, une fois les vrais zeros deja en place.

library(reshape2)

# Jointure de `valeurs` sur la grille complete `grille`, les lignes manquantes
# recevant `defaut` (0) dans `colonne_valeur`.
completer_zeros <- function(grille, valeurs, by, colonne_valeur, defaut = 0) {
  df <- merge(grille, valeurs, by = by, all.x = TRUE)
  df[[colonne_valeur]][is.na(df[[colonne_valeur]])] <- defaut
  df
}

# Matrice station x taxon, a appeler sur des donnees deja agregees a ce grain
# (par exemple une annee a la fois).
mat_station_taxon <- function(df, station_col, taxon_col, valeur_col = "count") {
  formule <- stats::as.formula(paste(station_col, "~", taxon_col))
  acast(df, formule, value.var = valeur_col, fill = 0)
}

# Matrice au grain de l'echantillon: une ligne par combinaison de `id_cols`
# (jointes par "|"), une colonne par taxon. Les lignes de somme nulle sont
# retirees sauf si drop_vides = FALSE.
mat_echantillon_taxon <- function(df, id_cols, taxon_col, valeur_col = "count",
                                  sep = "|", drop_vides = TRUE) {
  row_id <- do.call(paste, c(df[id_cols], sep = sep))
  long   <- data.frame(row_id = row_id, taxon = df[[taxon_col]], valeur = df[[valeur_col]])
  mat    <- acast(long, row_id ~ taxon, value.var = "valeur", fill = 0)
  if (drop_vides) mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  mat
}

# Version en abondances relatives d'une matrice de comptages (lignes = 1).
en_relatif <- function(mat) mat / rowSums(mat)

# Noms des n taxons les plus abondants, toutes matrices passees confondues.
# L'appelant passe exactement la ou les matrices qui alimentent le modele ou la
# figure en cours, pour que le "top" corresponde a ce qui est reellement montre.
top_taxons <- function(..., n = N_TOP_TAXONS) {
  totaux <- colSums(do.call(rbind, list(...)))
  names(sort(totaux, decreasing = TRUE))[1:n]
}
