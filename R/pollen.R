# Pollen data loading and the article's true-zero station x period x year x
# taxon grid (originally Pipeline_article.r:41-68). OTHER is kept here — it is
# excluded later, wherever composition (not total load) is involved, by
# whichever script needs that (Pred_Genus != "OTHER" is a one-line filter, not
# worth a function of its own).

# Reads and stacks the cleaned yearly CSVs (as written by 01_preparer_pollen.R
# / the old read_clean_data.r), coercing the columns used as model factors to
# character so downstream droplevels()/factor() calls behave the same
# regardless of how read.csv guessed their type.
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

# The samples that actually exist, i.e. the reference grid for completer_zeros().
# The sampling design is incomplete (25 stations x 7 periods x 2 years = 350
# possible combinations, fewer actually collected): without this, "taxon
# absent from a sample" and "sample never collected" would be confounded.
echantillons <- function(df, cols = c("location", "period", "year")) {
  ech <- unique(df[, cols])
  rownames(ech) <- NULL
  ech
}

# Counts per station x period x year x taxon, true zeros included. table()
# counts the grains actually observed; completer_zeros() puts back the zeros
# for taxa absent from a collected sample.
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

# Total pollen load per sample (all taxa summed, OTHER included).
abondance_totale <- function(ab_taxon) {
  aggregate(count ~ location + period + year, data = ab_taxon, FUN = sum)
}
