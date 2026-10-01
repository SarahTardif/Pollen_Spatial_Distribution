# Chargement des donnees polliniques et grille station x periode x annee x
# taxon avec les vrais zeros. OTHER est conserve ici; il est exclu plus loin,
# dans les scripts, partout ou il s'agit de composition et non de charge totale.

# Lit et empile les CSV annuels nettoyes (ecrits par 01_preparer_pollen.R). Les
# colonnes utilisees comme facteurs sont forcees en caractere pour que les
# droplevels()/factor() en aval se comportent toujours pareil.
charger_pollen <- function(annees, dir = DIR_DATA_POLLEN) {
  data_pollen <- NULL
  for (an in annees) {
    morceau     <- read.csv(file.path(dir, paste0("data", an, ".csv")), header = TRUE)
    data_pollen <- rbind(data_pollen, morceau)
  }
  for (col in c("location", "period", "year", "Pred_Genus")) {
    data_pollen[[col]] <- as.character(data_pollen[[col]])
  }
  data_pollen
}

# Les echantillons qui existent reellement, soit la grille de reference de
# completer_zeros(). Le plan d'echantillonnage est incomplet (25 stations x 7
# periodes x 2 annees = 350 combinaisons possibles, moins collectees): sans
# cette table, "taxon absent" et "echantillon non collecte" seraient confondus.
echantillons <- function(df, cols = c("location", "period", "year")) {
  ech <- unique(df[, cols])
  rownames(ech) <- NULL
  ech
}

# Comptages par station x periode x annee x taxon, vrais zeros inclus.
comptages_pollen <- function(data_pollen, ech) {
  ab_counts <- as.data.frame(table(location   = data_pollen$location,
                                   period     = data_pollen$period,
                                   year       = data_pollen$year,
                                   Pred_Genus = data_pollen$Pred_Genus),
                             stringsAsFactors = FALSE)
  names(ab_counts)[names(ab_counts) == "Freq"] <- "count"
  ab_counts <- ab_counts[ab_counts$count > 0, ]

  taxons <- sort(unique(data_pollen$Pred_Genus))
  grille <- merge(ech, data.frame(Pred_Genus = taxons, stringsAsFactors = FALSE))

  completer_zeros(grille, ab_counts,
                  by = c("location", "period", "year", "Pred_Genus"),
                  colonne_valeur = "count")
}

# Charge pollinique totale par echantillon (tous taxons, OTHER inclus).
abondance_totale <- function(ab_taxon) {
  aggregate(count ~ location + period + year, data = ab_taxon, FUN = sum)
}
