# Pollen data loading and the station x period x year x taxon grid with true
# zeros. OTHER is kept here; it is excluded further down, in the scripts,

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


echantillons <- function(df, cols = c("location", "period", "year")) {
  ech <- unique(df[, cols])
  rownames(ech) <- NULL
  ech
}

# Counts per station x period x year x taxon, true zeros included.
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

# Total pollen load per sample (all taxa, OTHER included).
abondance_totale <- function(ab_taxon) {
  aggregate(count ~ location + period + year, data = ab_taxon, FUN = sum)
}
