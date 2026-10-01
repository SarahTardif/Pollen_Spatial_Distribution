# Formats the models fitted by 04_article_analyses.R into publication-ready
# tables. Reads the objects left in the environment by the pipeline and writes
# the same set of tables to markdown and to Word. Run 04 first, or let the
# guard below do it.
# The formatting and writing functions live in R/report.R; what goes in which
# table (titles, legends) stays here, as editorial content.

source("R/init.R")
if (!exists("mod_ab_2022_noOTHER")) source("Scripts/04_article_analyses.R")

file_md   <- file.path(DIR_OUT_ARTICLE, "model_results.md")
file_docx <- file.path(DIR_OUT_ARTICLE, "model_results.docx")


#### BUILDING THE TABLES ####

# Every table is added to this list once, and both writers consume it, so the
# markdown and the Word document can never diverge.
tables <- list()

ajouter <- function(titre, legende, df) {
  tables[[length(tables) + 1]] <<- list(titre = titre, legende = legende, df = df)
}

## Table 1 - total abundance, location + period fixed, OTHER excluded
ajouter(
  "Table 1. Effect of location and period on total pollen abundance, OTHER excluded",
  paste("Negative binomial models (glmmTMB, nbinom2), one per year, with",
        "location and period both fitted as fixed effects: count ~ location +",
        "period, on the total load with OTHER excluded (identified grains",
        "only). Type II Wald chi-square tests."),
  rbind(tab_anova(mod_ab_2022_noOTHER, "2022"),
        tab_anova(mod_ab_2023_noOTHER, "2023")))

## Table 2 - estimated means per station
ajouter(
  "Table 2. Estimated mean abundance per station",
  paste("Estimated marginal means on the response scale, with Tukey-adjusted",
        "comparisons. Stations sharing a letter do not differ significantly."),
  rbind(tab_cld(cld_loc_2022_noOTHER, "2022"),
        tab_cld(cld_loc_2023_noOTHER, "2023")))

## Table 3 - abundance per taxon, one model per taxon (location + period
## fixed). Taxa that could not be fit for a given year (too few levels,
## all-zero counts) or whose Anova failed to converge are listed with a note
## instead of being silently dropped.

anova_abrel_taxon <- list()
for (an in annees) {
  mods <- get(paste0("mod_abrel_taxon_", an))
  skip <- get(paste0("skip_abrel_taxon_", an))

  for (taxon in names(mods)) {
    a <- try(tab_anova(mods[[taxon]], paste(an, taxon), garder_p_num = TRUE), silent = TRUE)
    if (!inherits(a, "try-error")) {
      a$Note <- ""
      a$Year <- an
      anova_abrel_taxon[[paste(an, taxon)]] <- a
    } else {
      skip[taxon] <- paste("Anova failed:", as.character(a))
    }
  }

  for (taxon in names(skip)) {
    nr <- note_row(paste(an, taxon), skip[[taxon]])
    nr$p_num <- NA_real_
    nr$Year  <- an
    anova_abrel_taxon[[paste(an, taxon, "skip")]] <- nr
  }
}
if (length(anova_abrel_taxon) > 0) {
  tab_abrel_taxon <- do.call(rbind, anova_abrel_taxon)

  ## Benjamini-Hochberg correction for multiple testing across taxa: one
  ## family per year, location and period p-values pooled together, 2022 and
  ## 2023 never mixed.
  tab_abrel_taxon$p_adj <- NA_real_
  for (an in annees) {
    idx <- tab_abrel_taxon$Year == an & !is.na(tab_abrel_taxon$p_num)
    tab_abrel_taxon$p_adj[idx] <- p.adjust(tab_abrel_taxon$p_num[idx], method = "BH")
  }
  tab_abrel_taxon$`p (BH-adjusted)` <- fmt_p(tab_abrel_taxon$p_adj)
  tab_abrel_taxon$p_num <- NULL
  tab_abrel_taxon$p_adj <- NULL
  tab_abrel_taxon$Year  <- NULL
  tab_abrel_taxon <- tab_abrel_taxon[, c("Model", "Term", "Chi2", "df", "p", "p (BH-adjusted)",
                                         "N obs", "R2 marginal", "R2 conditional", "R2 type", "Note")]

  ajouter(
    "Table 3. Effect of location and period on the abundance of each taxon, one model per taxon",
    paste("Negative binomial models (glmmTMB, nbinom2), one per taxon per year:",
          "count ~ location + period, fit only on that taxon's rows (all taxa,",
          "OTHER excluded, not restricted to top10). Type II Wald chi-square",
          "tests. p (BH-adjusted): Benjamini-Hochberg correction for multiple",
          "testing, applied per year across all taxa and both terms (location",
          "and period) together. Taxa that could not be fit for a given year",
          "(too few levels or all-zero counts) or whose Anova failed to",
          "converge are listed with a Note instead of a result."),
    tab_abrel_taxon)
}

## Table 4 - PERMANOVA
ajouter(
  "Table 4. PERMANOVA on pollen composition",
  paste("Bray-Curtis dissimilarities on the relative abundances per sample",
        "(station x period), OTHER excluded, 999 permutations. Period is",
        "fitted first so that the station effect is tested after the seasonal",
        "succession has been removed."),
  do.call(rbind, lapply(annees, function(an) tab_permanova(perm_comp[[an]], an))))

## Table 5 - PCoA axes
pct_tab <- do.call(rbind, pct_list)
ajouter(
  "Table 5. PCoA of the annual composition per station",
  paste("Principal coordinates analysis on the same Bray-Curtis matrix.",
        "Bray-Curtis is semi-metric, so negative eigenvalues appear; the",
        "Cailliez correction is applied when they exceed 10 % of the positive",
        "variance."),
  data.frame(
    Year       = pct_tab$year,
    Axis1      = paste0(fmt_num(pct_tab$axis1_pct), " %"),
    Axis2      = paste0(fmt_num(pct_tab$axis2_pct), " %"),
    NegEigen   = paste0(fmt_num(pct_tab$neg_eigen_pct), " %"),
    Correction = ifelse(pct_tab$cailliez, "Cailliez", "none"),
    stringsAsFactors = FALSE,
    row.names = NULL))
names(tables[[length(tables)]]$df) <- c("Year", "Axis 1 (% variance)",
                                        "Axis 2 (% variance)",
                                        "Negative eigenvalues", "Correction")

## Table 6 - envfit
envfit_tab <- NULL
for (an in annees) {
  ef  <- envfit_list[[an]]
  sig <- ef[ef$p < 0.05, ]
  sig <- sig[order(sig$p), ]
  if (nrow(sig) > 0) {
    envfit_tab <- rbind(envfit_tab,
                        data.frame(Year  = an,
                                   Taxon = sig$taxon,
                                   Axis1 = fmt_num(sig$Axe1),
                                   Axis2 = fmt_num(sig$Axe2),
                                   p     = fmt_p(sig$p),
                                   stringsAsFactors = FALSE))
  }
}
if (!is.null(envfit_tab)) {
  names(envfit_tab) <- c("Year", "Taxon", "Axis 1", "Axis 2", "p")
  ajouter(
    "Table 6. Taxa correlated with the PCoA axes",
    "envfit on the PCoA coordinates, 999 permutations, significant taxa only (p < 0.05).",
    envfit_tab)
}

## Table 7 - Shannon diversity
ajouter(
  "Table 7. Effect of the station on pollen diversity",
  paste("Gaussian mixed models (glmmTMB), one per year:",
        "Shannon ~ station + (1 | period). Shannon computed per sample on the",
        "composition matrix, OTHER excluded."),
  rbind(tab_anova(mod_compdiv_2022, "2022"),
        tab_anova(mod_compdiv_2023, "2023")))

## Table 8 - PAM partitioning
pam_list <- list("2022" = pam_2022, "2023" = pam_2023)
pam_tab  <- NULL
for (an in annees) {
  p <- pam_list[[an]]
  pam_tab <- rbind(pam_tab,
                   data.frame(Year       = an,
                              k          = fmt_int(p$courbe$k),
                              Silhouette = fmt_num(p$courbe$sil),
                              Optimal    = ifelse(p$courbe$k == p$k, "yes", ""),
                              stringsAsFactors = FALSE))
}
names(pam_tab) <- c("Year", "Number of groups (k)", "Average silhouette width",
                    "Retained")
ajouter(
  "Table 8. k-medoids partitioning of the stations",
  paste("PAM on the Bray-Curtis dissimilarities between annual station",
        "compositions, k swept from 2 to 10. An average silhouette width below",
        "0.25 indicates that no real group structure exists."),
  pam_tab)

## Table 9 - conclusion of the partitioning.
# k = 3 is added below the retained k for every year, even when its silhouette
# stays under 0.25, since it is the number of groups sometimes suggested by eye
# on the PCoA (fig2). A third block answers a different question again: reduced
# to a single station (k = 1, not part of the PAM sweep), which one is closest
# to the average composition ? Silhouette and structure do not apply to a lone
# station, hence the "-".
tab9_retenu <- do.call(rbind, lapply(annees, function(an) {
  p <- pam_list[[an]]
  data.frame(
    Year       = an,
    k          = fmt_int(p$k),
    Silhouette = fmt_num(p$sil),
    Structure  = ifelse(p$sil < 0.25, "no real structure", "structure detected"),
    Medoids    = ifelse(p$sil < 0.25, "-", paste(p$medoids, collapse = ", ")),
    Selection  = "retained (highest silhouette)",
    Distance   = "-",
    stringsAsFactors = FALSE)
}))
tab9_k3 <- do.call(rbind, lapply(annees, function(an) {
  p <- pam_list[[an]]
  data.frame(
    Year       = an,
    k          = "3",
    Silhouette = fmt_num(p$k3$sil),
    Structure  = ifelse(p$k3$sil < 0.25, "no real structure", "structure detected"),
    Medoids    = paste(p$k3$medoids, collapse = ", "),
    Selection  = "k = 3 (comparison)",
    Distance   = "-",
    stringsAsFactors = FALSE)
}))
tab9_moyenne <- do.call(rbind, lapply(annees, function(an) {
  s <- get(paste0("station_moy_", an))
  data.frame(
    Year       = an,
    k          = "1",
    Silhouette = "-",
    Structure  = "-",
    Medoids    = s$station,
    Selection  = "single station (closest to average composition)",
    Distance   = fmt_num(s$distance),
    stringsAsFactors = FALSE)
}))
ajouter(
  "Table 9. Optimal partition and medoid stations",
  paste("Medoid stations are only meaningful when the silhouette width reaches",
        "0.25. k = 3 is shown for every year for comparison, even when below",
        "threshold. The single-station row (k = 1) is not part of the PAM sweep:",
        "it is the station whose annual composition is closest (Bray-Curtis) to",
        "the average composition across all stations that year."),
  rbind(tab9_retenu, tab9_k3, tab9_moyenne))
names(tables[[length(tables)]]$df) <- c("Year", "k", "Silhouette width",
                                        "Interpretation", "Medoid stations",
                                        "Selection", "Distance to average")

## Table 10 - Mantel tests
ajouter(
  "Table 10. Mantel tests of distance decay",
  paste("Correlation between the geographic distance matrix (EPSG 32188,",
        "metres) and the pollen distance matrices, 999 permutations."),
  data.frame(Year       = mantel_res$year,
             Comparison = mantel_res$comparison,
             r          = fmt_num(mantel_res$mantel_r),
             p          = fmt_p(mantel_res$p_value),
             stringsAsFactors = FALSE))

## Table 11 - estimated means per period.
# Period's Anova is already in Table 1 (same two-factor model); this table only
# adds the period Tukey letters.
ajouter(
  "Table 11. Estimated mean abundance per period",
  paste("Estimated marginal means on the response scale, with Tukey-adjusted",
        "comparisons. Periods sharing a letter do not differ significantly."),
  rbind(tab_cld(cld_per_2022_noOTHER, "2022", groupe_col = "period", nom_groupe = "Period"),
        tab_cld(cld_per_2023_noOTHER, "2023", groupe_col = "period", nom_groupe = "Period")))

## Table S1 - model diagnostics.
# No row for the per-taxon models: none of them goes through valider(), so
# there is no single DHARMa object to report for them.
ajouter(
  "Table S1. Model diagnostics",
  paste("DHARMa tests on the scaled residuals (1000 simulations). A p-value",
        "above 0.05 means the assumption is not rejected."),
  rbind(tab_diag(res_ab_2022_noOTHER,  "Abundance 2022 (OTHER excluded)"),
        tab_diag(res_ab_2023_noOTHER,  "Abundance 2023 (OTHER excluded)"),
        tab_diag(res_compdiv_2022, "Diversity 2022"),
        tab_diag(res_compdiv_2023, "Diversity 2023")))
names(tables[[length(tables)]]$df) <- c("Model", "Uniformity (KS)",
                                        "Dispersion", "Outliers")

## Table S2 - estimated mean diversity (Shannon) per station
ajouter(
  "Table S2. Estimated mean pollen diversity (Shannon) per station",
  paste("Estimated marginal means on the response scale, with Tukey-adjusted",
        "comparisons. Stations sharing a letter do not differ significantly."),
  rbind(tab_cld(cld_compdiv_2022, "2022"),
        tab_cld(cld_compdiv_2023, "2023")))


#### FIGURES ####

# The figures themselves are written by 04_article_analyses.R at 300 dpi. They
# are only listed here so the document states which file goes with which
# figure.
figures <- data.frame(
  fichier = c("fig1_abondance_totale_par_station.png",
              "fig1b_abondance_totale_par_period.png",
              "fig1c_heatmap_abondance_station_periode.png",
              "fig1d_abondance_par_taxon_station_2022.png",
              "fig1d_abondance_par_taxon_station_2023.png",
              "fig2_pcoa_composition_2022.png",
              "fig2_pcoa_composition_2023.png",
              "fig3_shannon_par_station.png",
              "fig4_composition_par_station.png",
              "fig5_silhouette_pam.png",
              "fig6_carte_groupes_stations.png"),
  legende = c("Total pollen abundance per sample and per station, OTHER excluded, log10 scale, by year. See Table 2 for Tukey-adjusted comparisons between stations.",
              "Total pollen abundance per sample and per period, OTHER excluded, log10 scale, by year. See Table 11 for Tukey-adjusted comparisons between periods.",
              "Total pollen abundance per station x period, OTHER excluded, log10 scale, by year, coloured by deviation from that period's own median (red = above, blue = below). Descriptive complement to Tables 2/11: no model is fit separately per period, since a single location x year value per period leaves no replication.",
              "Pollen abundance per station and per taxon, 2022: one panel per taxon, showing the 12 most abundant taxa of that year in decreasing order of abundance. Each horizontal boxplot is one station, with its individual samples jittered on top; the dashed line is that taxon's mean across all stations. The x-axis is free between panels, since taxon counts span several orders of magnitude. Descriptive counterpart to Table 3, which tests the station effect taxon by taxon.",
              "Pollen abundance per station and per taxon, 2023: same as the previous figure, for 2023 (that year's own 12 most abundant taxa).",
              "PCoA of the annual pollen composition per station, 2022 (Bray-Curtis), with the taxa correlated with the axes.",
              "PCoA of the annual pollen composition per station, 2023 (Bray-Curtis), with the taxa correlated with the axes.",
              "Shannon diversity per sample and per station, by year, with Tukey letters (Table S2).",
              "Annual relative abundance of the 10 most abundant taxa per station, 2022 and 2023 (each year's own top 10 taxa; a taxon common to both years' top 10 keeps the same colour in both panels, a taxon that is top 10 in only one year gets its own additional colour).",
              "Average silhouette width as a function of the number of groups k.",
              "Map of the station groups and of the medoid stations, by year."),
  stringsAsFactors = FALSE)


#### WRITE ####

ecrire_md(tables, figures, file_md)
ecrire_docx(tables, figures, file_docx)

cat("\nWritten:\n", file_md, "\n", file_docx, "\n")
