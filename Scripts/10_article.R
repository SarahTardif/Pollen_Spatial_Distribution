# Analysis pipeline for the article: spatial distribution of pollen in Montreal
# Three questions:
#   1. does pollen abundance vary between stations ?
#   2. does pollen composition vary between stations ?
#   3. how many stations, and which ones, are needed for a representative network ?
# 2022 and 2023 are always modelled separately.
# ex-Pipeline_article.r

source("R/init.R")
library(glmmTMB)
library(vegan)
library(cluster)
library(emmeans)
library(multcompView)
library(multcomp)   # provides cld(); emmeans only registers the method for it

dir.create(DIR_FIG_ARTICLE, recursive = TRUE, showWarnings = FALSE)

#### DATA PREPARATION ####

## pollen data, already filtered at Confidence >= 0.8 by 01_preparer_pollen.R
data_pollen <- charger_pollen(c(2022, 2023))
locations   <- charger_stations()

# canopy cover class per station, same 4 classes/colours as the station map
# (map_25plots.r) -- used to colour the PCoA points below (fig2).
gradients <- charger_gradients()
locations$canopy_classe <- classer_canopy(
  gradients$Canopy.cover....[match(locations$trap, gradients$Plot)])

annees <- sort(unique(data_pollen$year))

## samples that actually exist
# The sampling design is not complete: 25 stations x 7 periods x 2 years = 350
# possible combinations, but fewer were actually collected. This table is the
# reference used below to fill in the true zeros: without it, "taxon absent from
# a sample" and "sample never collected" would be confounded.
data_samples <- echantillons(data_pollen)
cat("Samples collected:", nrow(data_samples), "\n")

## counts per station x period x year x taxon, zeros included
# OTHER is kept here because it is part of the total pollen load; it is
# removed further down wherever composition is involved.
ab_taxon <- comptages_pollen(data_pollen, data_samples)

## total pollen load (Ab_tot), OTHER included
ab_tot      <- abondance_totale(ab_taxon)
ab_tot_year <- aggregate(count ~ location + year, data = ab_tot, FUN = sum)

## composition tables, OTHER excluded
ab_taxon_comp <- ab_taxon[ab_taxon$Pred_Genus != "OTHER", ]

## total pollen load without OTHER, for the sensitivity check in section 1
ab_tot_noOTHER <- abondance_totale(ab_taxon_comp)

# annual matrices, one per year, rows = stations
ab_year_taxon <- aggregate(count ~ location + year + Pred_Genus,
                           data = ab_taxon_comp, FUN = sum)

mat_2022 <- mat_station_taxon(ab_year_taxon[ab_year_taxon$year == "2022", ], "location", "Pred_Genus")
mat_2023 <- mat_station_taxon(ab_year_taxon[ab_year_taxon$year == "2023", ], "location", "Pred_Genus")

mat_rel_2022 <- en_relatif(mat_2022)
mat_rel_2023 <- en_relatif(mat_2023)

# sample-level matrix (station|period|year x taxon), used for Shannon and PERMANOVA.
# Built from ab_taxon_comp, i.e. every taxon (only OTHER excluded) — PERMANOVA
# below is NOT restricted to top10, unlike the taxon:station/taxon:period
# interaction GLMMs further down.
mat_samples     <- mat_echantillon_taxon(ab_taxon_comp, c("location", "period", "year"), "Pred_Genus")
mat_rel_samples <- en_relatif(mat_samples)

env_samples <- data.frame(
  row_id   = rownames(mat_samples),
  location = sapply(strsplit(rownames(mat_samples), "\\|"), `[`, 1),
  period   = sapply(strsplit(rownames(mat_samples), "\\|"), `[`, 2),
  year     = sapply(strsplit(rownames(mat_samples), "\\|"), `[`, 3),
  stringsAsFactors = FALSE
)

cat("Annual rows:", nrow(mat_2022), "+", nrow(mat_2023),
    "| sample rows:", nrow(mat_samples),
    "| taxa (without OTHER):", ncol(mat_samples), "\n")

## per-year top10 taxa — computed separately per year (not pooled), so each
## year's models/figures reflect that year's own dominant taxa, matching
## top_taxons()'s own contract (R/matrices.R): pass exactly the matrix that
## feeds the model/figure at hand, not a list pooled across years.
top10_2022 <- top_taxons(mat_2022)
top10_2023 <- top_taxons(mat_2023)
top10_list <- list("2022" = top10_2022, "2023" = top10_2023)



#### 1. SPATIAL VARIATION IN POLLEN ABUNDANCE ####
## OTHER excluded only: this section used to carry a with/without-OTHER pair
## (mod_abglob_*/mod_abglob_*_noOTHER, etc.); now only the OTHER-excluded
## analysis remains. location and period are both fixed effects
## (count ~ location + period, no random term) instead of the previous two
## complementary mixed models (location fixed + (1|period), and period fixed
## + (1|location)): one model now supplies both the station and the period
## Tukey letters (reported in Table 2 / Table 11, not on the figures — the
## boxplots below are log10, unlettered). Same formula as
## 14_article_effets_fixes.R's mod_ab_*_noOTHER, refit locally here (not
## sourced) so the two pipelines stay independent.

ab_tot_noOTHER$location <- factor(ab_tot_noOTHER$location)
ab_tot_noOTHER$period   <- factor(ab_tot_noOTHER$period)

ab_tot_noOTHER_2022 <- droplevels(ab_tot_noOTHER[ab_tot_noOTHER$year == "2022", ])
ab_tot_noOTHER_2023 <- droplevels(ab_tot_noOTHER[ab_tot_noOTHER$year == "2023", ])

## total abundance model, one per year
mod_ab_2022_noOTHER <- glmmTMB(count ~ location + period,
                               data = ab_tot_noOTHER_2022, family = nbinom2)
res_ab_2022_noOTHER <- valider(mod_ab_2022_noOTHER, "Abundance 2022 (OTHER excluded, location + period fixed)", tracer = FALSE, verbose = FALSE)
emm_loc_2022_noOTHER <- emmeans(mod_ab_2022_noOTHER, ~ location, type = "response")
cld_loc_2022_noOTHER <- cld(emm_loc_2022_noOTHER, adjust = "tukey", Letters = letters)
emm_per_2022_noOTHER <- emmeans(mod_ab_2022_noOTHER, ~ period, type = "response")
cld_per_2022_noOTHER <- cld(emm_per_2022_noOTHER, adjust = "tukey", Letters = letters)

mod_ab_2023_noOTHER <- glmmTMB(count ~ location + period,
                               data = ab_tot_noOTHER_2023, family = nbinom2)
res_ab_2023_noOTHER <- valider(mod_ab_2023_noOTHER, "Abundance 2023 (OTHER excluded, location + period fixed)", tracer = FALSE, verbose = FALSE)
emm_loc_2023_noOTHER <- emmeans(mod_ab_2023_noOTHER, ~ location, type = "response")
cld_loc_2023_noOTHER <- cld(emm_loc_2023_noOTHER, adjust = "tukey", Letters = letters)
emm_per_2023_noOTHER <- emmeans(mod_ab_2023_noOTHER, ~ period, type = "response")
cld_per_2023_noOTHER <- cld(emm_per_2023_noOTHER, adjust = "tukey", Letters = letters)


## descriptive heatmap of total abundance per station x period (log scale),
## one panel per year. A model fit separately per period would collapse to a
## single location x year value per period (no replication left to test
## anything, same reasoning as 12_article_periodes.R) so this heatmap, not a
## model, is how the location x period pattern is inspected directly.
## Coloured by deviation from that period x year's OWN median (log10 scale,
## via ave()) rather than a single median pooled across all periods/years:
## this isolates station-to-station variation within each period instead of
## being dominated by the seasonal (period-to-period) signal, and the
## deviation is signed so the legend runs negative (below that period's
## median) to positive (above), instead of always-positive raw log10 values.
ab_tot_noOTHER$log_count         <- log10(ab_tot_noOTHER$count)
ab_tot_noOTHER$mediane_log_count <- ave(ab_tot_noOTHER$log_count,
                                        ab_tot_noOTHER$period, ab_tot_noOTHER$year,
                                        FUN = median)
ab_tot_noOTHER$ecart_mediane     <- ab_tot_noOTHER$log_count - ab_tot_noOTHER$mediane_log_count

heatmap_ab_station_period <- heatmap_station_periode(ab_tot_noOTHER,
                                                      x = "period", y = "location", fill = "ecart_mediane",
                                                      facet = "year", log10 = FALSE, divergent = TRUE, midpoint = 0,
                                                      xlab = "Period", ylab = "Station",
                                                      legend_lab = "Deviation from period median\n(log10)")

sauver_figure("fig1c_heatmap_abondance_station_periode.png", heatmap_ab_station_period,
              largeur = 8, hauteur = 10)


## per-taxon model, one count ~ location + period per taxon (all taxa, OTHER
## excluded), mirroring 14_article_effets_fixes.R's fit_taxon_models() —
## duplicated locally (not sourced) for the same independence reason as
## above. Replaces the old top10-restricted Pred_Genus:location /
## Pred_Genus:period interaction models. Their omnibus tests are reported in
## 11_article_tableaux.R's Table 3 (mirroring
## 15_article_effets_fixes_tableaux.R's Table 16); those tests say whether a
## taxon's abundance varies between stations but not how, so fig1d below
## shows the per-station distribution behind them, taxon by taxon — fig1 and
## fig1c only ever show it pooled over taxa.
fit_taxon_models <- function(dat) {

  taxa <- unique(as.character(dat$Pred_Genus))
  models  <- vector("list", length(taxa)); names(models) <- taxa
  skipped <- character(0)

  for (taxon in taxa) {

    dat_t <- droplevels(dat[dat$Pred_Genus == taxon, ])
    dat_t$location <- factor(dat_t$location)
    dat_t$period   <- factor(dat_t$period)

    if (nlevels(dat_t$location) < 2 || nlevels(dat_t$period) < 2 ||
        sum(dat_t$count) == 0) {
      skipped[taxon] <- "insufficient levels or all-zero counts - not fit"
      next
    }

    mod <- tryCatch(
      glmmTMB(count ~ location + period, data = dat_t, family = nbinom2),
      error = function(e) e
    )
    if (inherits(mod, "error")) {
      skipped[taxon] <- paste("fit error:", conditionMessage(mod))
      next
    }

    models[[taxon]] <- mod
  }

  list(models  = models[!vapply(models, is.null, logical(1))],
       skipped = skipped)
}

ab_taxon_2022 <- droplevels(ab_taxon_comp[ab_taxon_comp$year == "2022", ])
ab_taxon_2022$location   <- factor(ab_taxon_2022$location)
ab_taxon_2022$period     <- factor(ab_taxon_2022$period)
ab_taxon_2022$Pred_Genus <- factor(ab_taxon_2022$Pred_Genus)

ab_taxon_2023 <- droplevels(ab_taxon_comp[ab_taxon_comp$year == "2023", ])
ab_taxon_2023$location   <- factor(ab_taxon_2023$location)
ab_taxon_2023$period     <- factor(ab_taxon_2023$period)
ab_taxon_2023$Pred_Genus <- factor(ab_taxon_2023$Pred_Genus)

out_abrel_taxon_2022 <- fit_taxon_models(ab_taxon_2022)
mod_abrel_taxon_2022  <- out_abrel_taxon_2022$models    # named list, 1 glmmTMB per taxon
skip_abrel_taxon_2022 <- out_abrel_taxon_2022$skipped   # named char vector, taxon -> reason not fit

out_abrel_taxon_2023 <- fit_taxon_models(ab_taxon_2023)
mod_abrel_taxon_2023  <- out_abrel_taxon_2023$models
skip_abrel_taxon_2023 <- out_abrel_taxon_2023$skipped


## taxa shown on fig1d = the 12 most abundant of each year, each year having
## its own list (same contract as top10_list above: the figure of a year
## shows that year's own dominant taxa). All ~33 taxa x 25 stations on one
## figure would be unreadable, and the tail taxa are near-zero everywhere.
top12_list <- list("2022" = top_taxons(mat_2022, n = 12),
                   "2023" = top_taxons(mat_2023, n = 12))


## fig1d: one panel per taxon, horizontal station boxplots (stations on the
## y-axis) with the individual per-sample points jittered on top — same
## reading as fig1, but per taxon instead of pooled. Raw counts, i.e. the
## distribution each per-taxon model above was fit on; the models' own tests
## stay in Table 3. Points are plain grey, not coloured by period: with 7
## periods x 25 stations already dense per panel, a period legend added
## clutter without making any station-level pattern easier to read.
## scales = "free_x" is required, not cosmetic: taxon counts span several
## orders of magnitude (Betula in the thousands, minor taxa in single
## digits), so a shared x-axis would flatten every panel but the largest.
## A straight log10 x-axis is not usable here — a taxon absent from a
## collected sample is a true zero (completer_zeros()), and log10(0) would
## drop exactly the points that show a station has none of that taxon.
## Using a natural-log transform instead, via log1p (ln(1+x), scales::
## log1p_trans, bundled with ggplot2): it stays defined and equal to 0 at
## x = 0, so the true zeros stay visible while still compressing the
## multi-order-of-magnitude spread within each taxon panel. Breaks are fixed
## explicitly (rather than left to the default per-panel algorithm), which
## with free_x facets produced overlapping, unreadable tick labels.
## Station order is the canonical one (locations$trap), NOT reordered by
## median: with a different order per panel, comparing where a station sits
## across taxa — the whole point of the figure — would be impossible.
boxplot_taxon_station <- function(df, taxons, ordre_stations, titre) {

  df <- df[df$Pred_Genus %in% taxons, ]
  df$Pred_Genus <- factor(df$Pred_Genus, levels = taxons)   # most abundant taxon first
  df$location   <- factor(df$location, levels = rev(ordre_stations))  # rev: first station on top

  # per-panel mean line, pre-aggregated and passed through geom_vline()'s own
  # `data` for the same reason as the mean line in boxplot_station_periodes()
  # below: a stat_summary() layer would summarise per y-category (one line
  # per station) instead of per panel.
  moy <- aggregate(count ~ Pred_Genus, data = df, FUN = mean)

  ggplot(df, aes(x = count, y = location)) +
    geom_vline(data = moy, aes(xintercept = count),
               linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    # outlier.shape = NA: every point is already drawn individually below
    geom_boxplot(fill = "grey90", colour = "grey40", outlier.shape = NA, linewidth = 0.3,
                orientation = "y") +
    geom_jitter(height = 0.15, width = 0, size = 0.7, alpha = 0.5, colour = "grey30") +
    facet_wrap(~ Pred_Genus, ncol = 3, scales = "free_x") +
    scale_x_continuous(trans  = scales::log1p_trans(),
                        breaks = c(0, 10, 100, 1000, 10000),
                        labels = scales::label_number(big.mark = " ")) +
    theme_bw() +
    theme(axis.text.y      = element_text(size = 6),
          axis.text.x      = element_text(size = 6, angle = 45, hjust = 1),
          strip.text       = element_text(face = "italic"),
          panel.grid.minor = element_blank()) +
    labs(x = "Pollen count per sample (ln scale)", y = "Station", title = titre)
}

for (an in annees) {

  dat_an     <- get(paste0("ab_taxon_", an))
  taxons_fig <- top12_list[[an]]

  fig_taxon_station <- boxplot_taxon_station(dat_an, taxons_fig, locations$trap,
                                             titre = paste("Pollen abundance per station and taxon,", an))

  sauver_figure(paste0("fig1d_abondance_par_taxon_station_", an, ".png"), fig_taxon_station,
                largeur = 11,
                hauteur = min(2.4 * ceiling(length(taxons_fig) / 3) + 1, 16))
}


## Station boxplot with the 7 per-period points jittered on top and a global
## (all-stations) mean line per year panel, used for fig1 and fig3. Not added
## to boxplot_station() (R/plots.R) itself — same reasoning as
## boxplot_station_annee() in 12_article_periodes.R: that function is also
## used by box_ab_period below (x = "period", so per-period points don't
## apply) and by 14_article_effets_fixes.R, which keep the plain-grey look;
## this styling stays local to the figures that need it.
# outlier.shape = NA avoids double-plotting outliers, since every point is
# now shown individually via geom_jitter. The mean line is pre-aggregated
# (one row per facet level, e.g. per year) rather than left to stat_summary's
# implicit grouping — stat_summary(aes(group = 1), ...) alongside a discrete
# x still summarised per x-category here (one hline per station instead of
# one per panel), so the mean is computed by hand and passed to geom_hline()
# via its own `data` instead.
boxplot_station_periodes <- function(df, x, y, xlab = "Station", ylab, titre,
                                     facet = NULL, log10 = FALSE, angle_x = 45,
                                     lettres = NULL, nudge = 1.1, ordre = "median") {
  df[[x]] <- if (identical(ordre, "alpha")) factor(df[[x]]) else reorder(df[[x]], df[[y]], FUN = median)

  moyenne_globale <- if (!is.null(facet)) {
    moy <- aggregate(df[[y]], by = list(df[[facet]]), FUN = mean)
    names(moy) <- c(facet, "y")
    moy
  } else {
    data.frame(y = mean(df[[y]]))
  }

  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]])) +
    geom_boxplot(fill = "grey80", colour = "grey30", outlier.shape = NA) +
    geom_jitter(width = 0.15, size = 1, alpha = 0.5, colour = "grey20") +
    geom_hline(data = moyenne_globale, aes(yintercept = y),
               colour = "grey30", linetype = "dashed", linewidth = 0.6)

  if (!is.null(facet)) p <- p + facet_wrap(stats::as.formula(paste0("~", facet)), ncol = 1)
  if (log10) p <- p + scale_y_log10()
  if (!is.null(lettres)) p <- p + lettres_layer(df, x, y, lettres, facet, nudge)

  p + theme_bw() +
    theme(axis.text.x = element_text(angle = angle_x, hjust = 1),
          panel.grid = element_blank()) +
    labs(x = xlab, y = ylab, title = titre)
}


## boxplot of the total abundance per station
# log10 scale (counts span several orders of magnitude), x axis in natural
# station order rather than sorted by median, no Tukey letters on the figure
# itself — station means/letters are reported in Table 2 instead.
box_ab_station <- boxplot_station_periodes(ab_tot_noOTHER, x = "location", y = "count",
                                  ylab = "Total abundance (number of pollen grains, log10)",
                                  titre = "",
                                  facet = "year", log10 = TRUE, ordre = "alpha")

sauver_figure("fig1_abondance_totale_par_station.png", box_ab_station, largeur = 9, hauteur = 8)


## boxplot of the total abundance per period
box_ab_period <- boxplot_station_periodes(ab_tot_noOTHER, x = "period", y = "count", xlab = "Period",
                                 ylab = "Total abundance (number of pollen grains, log10)",
                                 titre = "",
                                 facet = "year", log10 = TRUE, ordre = "alpha")

sauver_figure("fig1b_abondance_totale_par_period.png", box_ab_period, largeur = 9, hauteur = 8)


#### 2. SPATIAL VARIATION IN POLLEN COMPOSITION ####

## Bray-Curtis dissimilarity + PERMANOVA, one per year
# The test is run on the sample-level matrix and not on the annual one: with 25
# rows and 25 station levels the annual model would be saturated.
# period comes first in the model because the seasonal succession of taxa is by
# far the strongest signal and would otherwise soak up the station effect. With
# by = "terms" adonis2 permutes the residuals of the reduced model, so location
# is tested AFTER the period effect has been removed, which is exactly what we
# want here.
# Do not add blocks = period: under permutation within blocks the period design
# matrix is invariant, its permuted F always equals the observed F, and vegan
# reuses that same permutation set for every term. Both terms then come out with
# the same meaningless p-value (checked: 0.234 / 0.234 in 2022).

perm_comp <- list()
for (an in annees) {

  lignes <- env_samples$row_id[env_samples$year == an]
  mat_rel_samples_an <- mat_rel_samples[lignes, , drop = FALSE]
  env_samples_an     <- env_samples[env_samples$year == an, ]

  set.seed(GRAINE)
  perm_comp[[an]] <- adonis2(mat_rel_samples_an ~ period + location,
                             data = env_samples_an, method = "bray",
                             permutations = N_PERM, by = "terms")

  cat("\nPERMANOVA", an, "\n")
}


## PCoA on the annual composition, with taxa arrows
# Not the same matrix as the PERMANOVA: same taxa/exclusions and same
# Bray-Curtis method, but aggregated at station level (mat_rel_2022/2023,
# one row per station) instead of sample level (mat_rel_samples, one row per
# station|period|year) — a station-level PCoA plot is what's readable/
# interpretable here, whereas PERMANOVA needs the sample-level replication to
# test location at all (see comment above). Bray-Curtis is semi-metric, so
# the PCoA produces negative eigenvalues: they are quantified, and the
# Cailliez correction is applied if they exceed 10 % of the positive
# variance.

mat_rel_list <- list("2022" = mat_rel_2022, "2023" = mat_rel_2023)
pcoa_list    <- list()
envfit_list  <- list()   # ef_df of each year, kept for 11_article_tableaux.R
pct_list     <- list()   # axis percentages of each year, same reason
dominants    <- NULL     # dominant taxon per station, accumulated across years below
comp_plot_all <- NULL   # per-year composition, accumulated for the merged fig4 below

for (an in annees) {

  # drop the taxa absent that year (null columns make envfit fail)
  mat_an <- mat_rel_list[[an]]
  mat_an <- mat_an[, colSums(mat_an) > 0, drop = FALSE]

  d_an    <- vegdist(mat_an, method = "bray")
  pcoa_an <- cmdscale(d_an, k = 2, eig = TRUE)

  eig     <- pcoa_an$eig
  var_pos <- sum(eig[eig > 0])
  var_neg <- abs(sum(eig[eig < 0]))
  # share of negative eigenvalues, computed before any correction
  pct_neg <- var_neg / var_pos * 100
  cat("\nPCoA", an, ": negative eigenvalues =",
      round(pct_neg, 1), "% of the positive variance\n")

  cailliez <- pct_neg > 10
  if (cailliez) {
    cat("Cailliez correction applied\n")
    pcoa_an <- wcmdscale(d_an, k = 2, eig = TRUE, add = "cailliez")
    eig     <- pcoa_an$eig
    var_pos <- sum(eig[eig > 0])
  }

  coord <- pcoa_an$points
  colnames(coord) <- c("Axe1", "Axe2")
  pct <- eig[1:2] / var_pos * 100
  cat("Axis 1 =", round(pct[1], 1), "% | Axis 2 =", round(pct[2], 1), "%\n")

  pcoa_list[[an]] <- pcoa_an

  sc_p          <- as.data.frame(coord)
  sc_p$location <- rownames(coord)
  sc_p$row_id   <- paste0(sc_p$location, "|", an)
  sc_p$canopy_classe <- locations$canopy_classe[match(sc_p$location, locations$trap)]

  # envfit: taxa correlated with the axes. Of the significant ones, only the
  # 10 with the highest envfit r2 (best-fitting to the axes, not necessarily
  # the most abundant) are drawn, otherwise almost all 33 taxa come out
  # significant and the plot becomes unreadable.
  fl     <- fleches_envfit(coord, mat_an, sc_p, axes = c("Axe1", "Axe2"), top_n = 15)
  ef_df  <- fl$complet
  ef_sig <- fl$sig

  # Stations as labelled points only: the pies made the plot unreadable once the
  # taxa arrows were added. Points are coloured by the station's canopy cover
  # class (same 4 classes/colours as the map_25plots.r station map) to check
  # by eye whether composition groups line up with canopy cover.
  pcoa_arrows <- plot_ordination(
    sc_p, axes = c("Axe1", "Axe2"), id_col = "location", pies = FALSE,
    fleches = ef_sig,
    point_colour_col = "canopy_classe", point_palette = COULEURS_CANOPY,
    point_legend = "Canopy cover (%)",
    titre = paste("PCoA of the pollen composition per station -", an, "(Bray-Curtis)"),
    sous_titre = "Stations and taxa correlated with the axes (envfit, p < 0.05)",
    xlab = sprintf("Axis 1 (%.1f %%)", pct[1]),
    ylab = sprintf("Axis 2 (%.1f %%)", pct[2]))

  sauver_figure(paste0("fig2_pcoa_composition_", an, ".png"), pcoa_arrows, largeur = 7, hauteur = 7)

  cat("\nTaxa correlated with the PCoA axes", an, "(envfit, p < 0.05):\n")

  # kept for 11_article_tableaux.R, which would otherwise only see the last year
  envfit_list[[an]] <- ef_df
  pct_list[[an]]    <- data.frame(year = an, axis1_pct = pct[1], axis2_pct = pct[2],
                                  neg_eigen_pct = pct_neg, cailliez = cailliez)

  ## composition of the year's pollen counts, kept per year (this year's own
  # top10 and its own "Others" bucket) but accumulated below so both years
  # can be drawn on one merged figure with a single, taxon-consistent palette
  # (see after the loop) instead of two separate files with independent,
  # rank-based colours that could put the same taxon under a different colour
  # in each year, or a different taxon under the same colour.
  top10_an <- top10_list[[an]]

  comp_long_an <- data.frame(
    location = rep(rownames(mat_an), times = ncol(mat_an)),
    taxon    = rep(colnames(mat_an), each = nrow(mat_an)),
    percent  = as.vector(mat_an) * 100,
    stringsAsFactors = FALSE
  )

  # taxa outside this year's top10 are lumped into "Others"
  comp_long_an$taxon_grp <- ifelse(comp_long_an$taxon %in% top10_an, comp_long_an$taxon, "Others")
  comp_plot_an <- aggregate(percent ~ location + taxon_grp, data = comp_long_an, FUN = sum)
  comp_plot_an$taxon_grp <- factor(comp_plot_an$taxon_grp, levels = c(top10_an, "Others"))

  # dominant taxon of each station, this year
  dominants_an <- comp_plot_an[comp_plot_an$taxon_grp != "Others", ]
  dominants_an <- dominants_an[order(dominants_an$location, -dominants_an$percent), ]
  dominants_an <- dominants_an[!duplicated(dominants_an$location), ]
  dominants_an$year <- an
  dominants <- rbind(dominants, dominants_an)

  comp_plot_an$year <- an
  comp_plot_all <- rbind(comp_plot_all, comp_plot_an)
}


## merged stacked bar of the annual composition per station, both years on
# one figure (facetted by year) instead of two separate files. Each year
# keeps its own top10/"Others" split computed above; what's unified here is
# only the colour scheme, via the union of the two years' top10 taxa, ordered
# by their combined abundance so the legend reads dominant-to-rare. A taxon
# common to both years' top10 gets the same colour in both panels; a taxon
# that is top10 in only one year still gets its own distinct colour rather
# than colliding with an unrelated taxon's colour or falling back to "Others"
# grey — palette_taxons() (R/plots.R) appends extra colours beyond the base
# 10 automatically when the union exceeds 10 taxa.
taxa_union   <- Reduce(union, top10_list)
totaux_union <- colSums(mat_2022[, taxa_union, drop = FALSE]) +
               colSums(mat_2023[, taxa_union, drop = FALSE])
taxa_union   <- names(sort(totaux_union, decreasing = TRUE))

pal_taxa_union <- palette_taxons(taxa_union, COULEURS_TAXONS)

comp_plot_all$taxon_grp <- factor(comp_plot_all$taxon_grp, levels = c(taxa_union, "Others"))

bar_comp_station <- ggplot(comp_plot_all, aes(x = location, y = percent, fill = taxon_grp)) +
  geom_bar(stat = "identity", colour = "grey30", linewidth = 0.15) +
  scale_fill_manual(values = pal_taxa_union, name = NULL) +
  facet_wrap(~ year, ncol = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(face = "italic")) +
  labs(x = "Station", y = "Annual relative abundance (%)",
       title = "Annual pollen composition per station")

sauver_figure("fig4_composition_par_station.png", bar_comp_station, largeur = 10, hauteur = 14)


## diversity model, one per year
# Shannon computed per sample on the composition matrix (OTHER excluded), so
# each station has up to 7 periods of replication within a year.

shannon_vals <- shannon(mat_samples)

div_samples <- data.frame(
  location = factor(env_samples$location),
  period   = factor(env_samples$period),
  year     = env_samples$year,
  shannon  = as.numeric(shannon_vals[env_samples$row_id]),
  stringsAsFactors = FALSE
)

div_2022 <- droplevels(div_samples[div_samples$year == "2022", ])
div_2023 <- droplevels(div_samples[div_samples$year == "2023", ])

# mod div 2022
mod_compdiv_2022 <- glmmTMB(shannon ~ location + (1 | period),
                            data = div_2022, family = gaussian)
res_compdiv_2022 <- valider(mod_compdiv_2022, "Diversity 2022", tracer = FALSE, verbose = FALSE)
emm_compdiv_2022 <- emmeans(mod_compdiv_2022, ~ location, type = "response")
cld_compdiv_2022 <- cld(emm_compdiv_2022, adjust = "tukey", Letters = letters)

# mod div 2023
mod_compdiv_2023 <- glmmTMB(shannon ~ location + (1 | period),
                            data = div_2023, family = gaussian)
res_compdiv_2023 <- valider(mod_compdiv_2023, "Diversity 2023", tracer = FALSE, verbose = FALSE)
emm_compdiv_2023 <- emmeans(mod_compdiv_2023, ~ location, type = "response")
cld_compdiv_2023 <- cld(emm_compdiv_2023, adjust = "tukey", Letters = letters)



## boxplot of the Shannon index per station
box_div_station <- boxplot_station_periodes(div_samples, x = "location", y = "shannon",
                                   ylab = "Shannon index",
                                   titre = "",
                                   facet = "year", log10 = FALSE)

sauver_figure("fig3_shannon_par_station.png", box_div_station, largeur = 9, hauteur = 8)



#### 3. HOW MANY STATIONS, AND WHICH ONES ####

## k-medoids (PAM) on the annual composition
# Sweep k from 2 to 10 and keep the k with the highest average silhouette
# width. Below 0.25 the silhouette indicates there is no real group structure.
# Runs on the relative-abundance matrix (mat_rel_2022/2023), the same one the
# PCoA above is built from: on raw counts, Bray-Curtis is not scale-invariant,
# so stations end up separated by total pollen load (which differs a lot
# between stations, see question 1) rather than by composition, and the
# resulting groups no longer match what the PCoA plot is showing.
run_pam <- function(mat, etiquette) {
    set.seed(GRAINE)
  d <- vegdist(mat, method = "bray")
  k_range <- 2:10
  sils <- numeric(length(k_range))
  for (i in seq_along(k_range)) {
    sils[i] <- pam(d, k = k_range[i], diss = TRUE)$silinfo$avg.width
  }
  k_opt   <- k_range[which.max(sils)]
  sil_max <- max(sils)
  pam_fin <- pam(d, k = k_opt, diss = TRUE)
  medoids <- rownames(mat)[pam_fin$id.med]

  # k = 3 kept alongside the optimal k, even when below the 0.25 threshold:
  # it is the number of groups sometimes suggested by eye on the PCoA (fig2),
  # so it is reported for comparison in Table 9 regardless of the retained k.
  pam_k3     <- pam(d, k = 3, diss = TRUE)
  medoids_k3 <- rownames(mat)[pam_k3$id.med]
  sil_k3     <- sils[k_range == 3]

  # k = 3 kept alongside the optimal k, even when below the 0.25 threshold:
  # it is the number of groups sometimes suggested by eye on the PCoA (fig2),
  # so it is reported for comparison in Table 9 regardless of the retained k.
  pam_k3     <- pam(d, k = 3, diss = TRUE)
  medoids_k3 <- rownames(mat)[pam_k3$id.med]
  sil_k3     <- sils[k_range == 3]

  cat("\nPAM -", etiquette, "\n")
  cat("optimal k =", k_opt, "(silhouette =", round(sil_max, 3), ")\n")
  if (sil_max < 0.25) {
    cat("silhouette < 0.25: no real group structure, the stations do not",
        "separate into clear groups\n")
  } else {
    cat("recommended medoid stations:", paste(medoids, collapse = ", "), "\n")
  }

  list(k = k_opt, sil = sil_max, medoids = medoids,
       clustering = pam_fin$clustering,
       k3 = list(sil = sil_k3, medoids = medoids_k3, clustering = pam_k3$clustering),
       k3 = list(sil = sil_k3, medoids = medoids_k3, clustering = pam_k3$clustering),
       courbe = data.frame(k = k_range, sil = sils))
}

pam_2022 <- run_pam(mat_rel_2022, "2022")
pam_2023 <- run_pam(mat_rel_2023, "2023")
pam_2022 <- run_pam(mat_rel_2022, "2022")
pam_2023 <- run_pam(mat_rel_2023, "2023")

## single "average" station
# Different question from PAM above: PAM groups stations and reports one
# medoid per group; this asks what a single-station network would look like —
# which one station's annual composition is closest to the average
# composition across all stations that year ? Same Bray-Curtis metric as
# PCoA/PAM/PERMANOVA above, obtained by adding the average as an extra
# "station" row to the matrix rather than switching to Euclidean distance
# here. Always run separately per year, like everything else in this script.
station_moyenne <- function(mat, etiquette) {

  moyenne <- colMeans(mat)   # mean relative composition, one station = one weight

  mat_avec_moyenne <- rbind(mat, moyenne = moyenne)
  d         <- as.matrix(vegdist(mat_avec_moyenne, method = "bray"))
  d_moyenne <- sort(d["moyenne", rownames(mat)])   # distance to the average, closest first

  station_rep <- names(d_moyenne)[1]

  cat("\nStation closest to the average composition -", etiquette, ":",
      station_rep, "(Bray-Curtis distance to the average =",
      round(d_moyenne[1], 3), ")\n")

  list(station = station_rep, distance = d_moyenne[1], distances = d_moyenne)
}

station_moy_2022 <- station_moyenne(mat_rel_2022, "2022")
station_moy_2023 <- station_moyenne(mat_rel_2023, "2023")

## silhouette curves
courbes <- rbind(data.frame(pam_2022$courbe, year = "2022"),
                 data.frame(pam_2023$courbe, year = "2023"))

plot_silhouette <- ggplot(courbes, aes(x = k, y = sil, colour = year)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.25, linetype = "dashed", colour = "grey40") +
  scale_x_continuous(breaks = 2:10) +
  scale_colour_manual(values = COULEURS_ANNEES, name = "Year") +
  theme_bw() +
  labs(x = "Number of groups (k)", y = "Average silhouette width",
       title = "How many distinct groups of stations ?",
       subtitle = "Dashed line = 0.25 threshold below which there is no real structure")

sauver_figure("fig5_silhouette_pam.png", plot_silhouette, largeur = 7, hauteur = 5)

## map of the groups and medoids
carte_2022 <- data.frame(location = names(pam_2022$clustering),
                         cluster  = factor(pam_2022$clustering),
                         year     = "2022",
                         is_medoid = names(pam_2022$clustering) %in% pam_2022$medoids,
                         stringsAsFactors = FALSE)
carte_2023 <- data.frame(location = names(pam_2023$clustering),
                         cluster  = factor(pam_2023$clustering),
                         year     = "2023",
                         is_medoid = names(pam_2023$clustering) %in% pam_2023$medoids,
                         stringsAsFactors = FALSE)

carte_clusters <- rbind(carte_2022, carte_2023)
carte_clusters <- merge(carte_clusters, locations[, c("trap", "Longitude", "Latitude")],
                        by.x = "location", by.y = "trap")

map_clusters <- ggplot(carte_clusters, aes(x = Longitude, y = Latitude, colour = cluster)) +
  geom_point(aes(size = is_medoid, shape = is_medoid)) +
  geom_text(aes(label = location), vjust = -1, size = 3, colour = "black") +
  facet_wrap(~year) +
  scale_size_manual(values = c(`FALSE` = 2.5, `TRUE` = 5), guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude", colour = "Group",
       title = "Groups of stations based on the annual pollen composition",
       subtitle = "Triangles = recommended medoid stations")

sauver_figure("fig6_carte_groupes_stations.png", map_clusters, largeur = 11, hauteur = 6)


## Mantel test: geographic distance vs pollen distance
# Coordinates are projected to EPSG 32188 (MTM zone 8, metres) before computing
# the distances.
locs_sf <- en_mtm(locations[, c("trap", "Longitude", "Latitude")])

mantel_res <- NULL   # accumulated across years, kept for 11_article_tableaux.R

for (an in annees) {

  mat_an   <- mat_rel_list[[an]]
  stations <- rownames(mat_an)

  dist_geo <- distances_stations(locs_sf, stations)

  # composition: Bray-Curtis between stations
  dist_comp <- vegdist(mat_an, method = "bray")
  set.seed(GRAINE)
  mantel_comp <- mantel(dist_comp, dist_geo, permutations = N_PERM)

  # abundance: difference in total annual load, on a log10 scale
  ab_an <- setNames(ab_tot_year$count[ab_tot_year$year == an],
                    ab_tot_year$location[ab_tot_year$year == an])
  dist_ab <- dist(log10(ab_an[stations]))
  set.seed(GRAINE)
  mantel_ab <- mantel(dist_ab, dist_geo, permutations = N_PERM)

  cat("\nMantel test", an, "- composition vs geography\n")
  cat("\nMantel test", an, "- abundance vs geography\n")

  mantel_res <- rbind(mantel_res,
                      data.frame(year       = an,
                                 comparison = c("Composition (Bray-Curtis) vs geography",
                                                "Total abundance (log10) vs geography"),
                                 mantel_r   = c(mantel_comp$statistic, mantel_ab$statistic),
                                 p_value    = c(mantel_comp$signif, mantel_ab$signif),
                                 stringsAsFactors = FALSE))
}
