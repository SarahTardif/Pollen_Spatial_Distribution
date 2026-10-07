#write tables from script04 results to markdown and to Word. 
# formatting and writing functions are in R/report.R

source("R/init.R")
if (!exists("mod_ab_noOTHER")) source("Scripts/04_article_analyses.R")

file_md   <- file.path(DIR_OUT_ARTICLE, "model_results.md")
file_docx <- file.path(DIR_OUT_ARTICLE, "model_results.docx")


#### BUILDING THE TABLES ####

# Every table is added to this list once
tables <- list()

# the legend is optional: ajouter(titre, df) or ajouter(titre, legende, df)
ajouter <- function(titre, legende, df) {
  if (missing(df)) { df <- legende; legende <- "" }
  tables[[length(tables) + 1]] <<- list(titre = titre, legende = legende, df = df)
}

## Table 1 - total abundance, location + period fixed, OTHER excluded
ajouter(
  "Table 1. Effect of location and period on total pollen abundance, OTHER excluded",
  tab_anova(mod_ab_noOTHER, "Total abundance"))

## Table 2 - estimated means per station
ajouter(
  "Table 2. Estimated mean abundance per station",
  tab_cld(cld_loc_noOTHER))

## Table 3 - abundance per taxon, one model per taxon (location + period
## fixed). 

anova_abrel_taxon <- list()
skip <- skip_abrel_taxon

for (taxon in names(mod_abrel_taxon)) {
  a <- try(tab_anova(mod_abrel_taxon[[taxon]], taxon, garder_p_num = TRUE), silent = TRUE)
  if (!inherits(a, "try-error")) {
    a$Note <- ""
    anova_abrel_taxon[[taxon]] <- a
  } else {
    skip[taxon] <- paste("Anova failed:", as.character(a))
  }
}

for (taxon in names(skip)) {
  nr <- note_row(taxon, skip[[taxon]])
  nr$p_num <- NA_real_
  anova_abrel_taxon[[paste(taxon, "skip")]] <- nr
}

if (length(anova_abrel_taxon) > 0) {
  tab_abrel_taxon <- do.call(rbind, anova_abrel_taxon)

  ## Benjamini-Hochberg correction for multiple testing across taxa
  tab_abrel_taxon$p_adj <- NA_real_
  idx <- !is.na(tab_abrel_taxon$p_num)
  tab_abrel_taxon$p_adj[idx] <- p.adjust(tab_abrel_taxon$p_num[idx], method = "BH")
  tab_abrel_taxon$`p (BH-adjusted)` <- fmt_p(tab_abrel_taxon$p_adj)
  tab_abrel_taxon$p_num <- NULL
  tab_abrel_taxon$p_adj <- NULL
  tab_abrel_taxon <- tab_abrel_taxon[, c("Model", "Term", "Chi2", "df", "p", "p (BH-adjusted)",
                                         "N obs", "R2 marginal", "R2 conditional", "R2 type", "Note")]

  ajouter(
    "Table 3. Effect of location and period on the abundance of each taxon, one model per taxon",
    tab_abrel_taxon)
}

## Table 4 - PERMANOVA
ajouter(
  "Table 4. PERMANOVA on pollen composition",
  tab_permanova(perm_comp))

## Table 5 - PCoA axes
ajouter(
  "Table 5. PCoA of the seasonal composition per station",
  data.frame(
    Axis1      = paste0(fmt_num(pct_res$axis1_pct), " %"),
    Axis2      = paste0(fmt_num(pct_res$axis2_pct), " %"),
    NegEigen   = paste0(fmt_num(pct_res$neg_eigen_pct), " %"),
    Correction = ifelse(pct_res$cailliez, "Cailliez", "none"),
    stringsAsFactors = FALSE,
    row.names = NULL))
names(tables[[length(tables)]]$df) <- c("Axis 1 (% variance)",
                                        "Axis 2 (% variance)",
                                        "Negative eigenvalues", "Correction")

## Table 6 - envfit
sig <- envfit_res[envfit_res$p < 0.05, ]
sig <- sig[order(sig$p), ]
if (nrow(sig) > 0) {
  envfit_tab <- data.frame(Taxon = sig$taxon,
                           Axis1 = fmt_num(sig$Axe1),
                           Axis2 = fmt_num(sig$Axe2),
                           p     = fmt_p(sig$p),
                           stringsAsFactors = FALSE)
  names(envfit_tab) <- c("Taxon", "Axis 1", "Axis 2", "p")
  ajouter(
    "Table 6. Taxa correlated with the PCoA axes",
    "envfit on the PCoA coordinates, 999 permutations, significant taxa only (p < 0.05).",
    envfit_tab)
}

## Table 7 - Shannon diversity
ajouter(
  "Table 7. Effect of the station on pollen diversity",
  tab_anova(mod_compdiv, "Shannon diversity"))

## Table 8 - Mantel tests
ajouter(
  "Table 8. Mantel tests of distance decay",
  data.frame(Comparison = mantel_res$comparison,
             r          = fmt_num(mantel_res$mantel_r),
             p          = fmt_p(mantel_res$p_value),
             stringsAsFactors = FALSE))

## Table 9 - estimated means per period.

ajouter(
  "Table 9. Estimated mean abundance per period",
  tab_cld(cld_per_noOTHER, groupe_col = "period", nom_groupe = "Period"))

## Table S1 - model diagnostics.

ajouter(
  "Table S1. Model diagnostics",
  rbind(tab_diag(res_ab_noOTHER, "Abundance (OTHER excluded)"),
        tab_diag(res_compdiv,    "Diversity")))
names(tables[[length(tables)]]$df) <- c("Model", "Uniformity (KS)",
                                        "Dispersion", "Outliers")

## Table S2 - estimated mean diversity (Shannon) per station
ajouter(
  "Table S2. Estimated mean pollen diversity (Shannon) per station",
  tab_cld(cld_compdiv))


#### FIGURES ####

# The figures themselves are written by 04_article_analyses.R listed here to know which file goes with which figure.
figures <- data.frame(
  fichier = c("fig_abondance_totale_par_station.png",
              "fig_abondance_totale_par_period.png",
              "fig_heatmap_abondance_station_periode.png",
              "fig_abondance_par_taxon_station.png",
              "fig_pcoa_composition.png",
              "fig_shannon_par_station.png",
              "fig_composition_par_station.png"),
  legende = c("Total pollen abundance per sample and per station, OTHER excluded, log10 scale. See Table 2 for Tukey-adjusted comparisons between stations.",
              "Total pollen abundance per sample and per period, OTHER excluded, log10 scale. See Table 9 for Tukey-adjusted comparisons between periods.",
              "Total pollen abundance per station x period, OTHER excluded, log10 scale, coloured by deviation from that period's own median (red = above, blue = below). Descriptive complement to Tables 2/9: no model is fit separately per period, since a single value per station and period leaves no replication.",
              "Pollen abundance per station and per taxon: one panel per taxon, showing the 12 most abundant taxa in decreasing order of abundance. Each horizontal boxplot is one station, with its individual samples jittered on top; the dashed line is that taxon's mean across all stations. The x-axis is free between panels, since taxon counts span several orders of magnitude. Descriptive counterpart to Table 3, which tests the station effect taxon by taxon.",
              "PCoA of the seasonal pollen composition per station (Bray-Curtis), with the taxa correlated with the axes.",
              "Shannon diversity per sample and per station (Table S2 for Tukey-adjusted comparisons).",
              "Seasonal relative abundance of the 10 most abundant taxa per station, all periods combined."),
  stringsAsFactors = FALSE)


#### WRITE ####

ecrire_md(tables, figures, file_md)
ecrire_docx(tables, figures, file_docx)

cat("\nWritten:\n", file_md, "\n", file_docx, "\n")
