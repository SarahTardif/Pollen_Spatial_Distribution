# Analysis pipeline for the article: spatial distribution of pollen in Montreal.
# Two questions: does pollen abundance vary between stations, and does
# composition vary between stations ? A single sampling season (Data_Maya).

source("R/init.R")
library(glmmTMB)
library(vegan)
library(emmeans)
library(multcompView)
library(multcomp)

dir.create(DIR_FIG_ARTICLE, recursive = TRUE, showWarnings = FALSE)

#### DATA PREPARATION ####

data_pollen <- charger_pollen()
locations   <- charger_stations()

# canopy cover class per station, same classes and colours as the map (script 02)
gradients <- charger_gradients()
locations$canopy_classe <- classer_canopy(
  gradients$Canopy.cover....[match(locations$trap, gradients$Plot)])

## samples actually collected, the reference grid for the true zeros
data_samples <- echantillons(data_pollen)
cat("Samples collected:", nrow(data_samples), "\n")

## counts per station x period x taxon, zeros included
ab_taxon <- comptages_pollen(data_pollen, data_samples)

## OTHER (unidentified grains) is excluded from every analysis below
ab_taxon_comp <- ab_taxon[ab_taxon$Pred_Genus != "OTHER", ]

ab_tot_noOTHER <- abondance_totale(ab_taxon_comp)
ab_tot_station <- aggregate(count ~ location, data = ab_tot_noOTHER, FUN = sum)

## seasonal matrix, rows = stations, all periods combined
ab_station_taxon <- aggregate(count ~ location + Pred_Genus,
                              data = ab_taxon_comp, FUN = sum)

mat_station     <- mat_station_taxon(ab_station_taxon, "location", "Pred_Genus")
mat_rel_station <- en_relatif(mat_station)

## sample-level matrix (station|period x taxon), for Shannon and the PERMANOVA
mat_samples     <- mat_echantillon_taxon(ab_taxon_comp, c("location", "period"), "Pred_Genus")
mat_rel_samples <- en_relatif(mat_samples)

env_samples <- data.frame(
  row_id   = rownames(mat_samples),
  location = sapply(strsplit(rownames(mat_samples), "\\|"), `[`, 1),
  period   = sapply(strsplit(rownames(mat_samples), "\\|"), `[`, 2),
  stringsAsFactors = FALSE
)

cat("Station rows:", nrow(mat_station),
    "| sample rows:", nrow(mat_samples),
    "| taxa (without OTHER):", ncol(mat_samples), "\n")

top10 <- top_taxons(mat_station)



#### 1. SPATIAL VARIATION IN ABUNDANCE ####

## location and period are both fixed effects

ab_tot_noOTHER$location <- factor(ab_tot_noOTHER$location)
ab_tot_noOTHER$period   <- factor(ab_tot_noOTHER$period, levels = levels(data_pollen$period))

## total abundance model
mod_ab_noOTHER <- glmmTMB(count ~ location + period,
                          data = ab_tot_noOTHER, family = nbinom2)
res_ab_noOTHER <- valider(mod_ab_noOTHER, "Abundance (OTHER excluded, location + period fixed)", tracer = TRUE, verbose = FALSE)
emm_loc_noOTHER <- emmeans(mod_ab_noOTHER, ~ location, type = "response")
cld_loc_noOTHER <- cld(emm_loc_noOTHER, adjust = "tukey", Letters = letters)
emm_per_noOTHER <- emmeans(mod_ab_noOTHER, ~ period, type = "response")
cld_per_noOTHER <- cld(emm_per_noOTHER, adjust = "tukey", Letters = letters)


## descriptive heatmap of the total abundance per station x period.
# Coloured by the deviation from its own period median

ab_tot_noOTHER$log_count         <- log10(ab_tot_noOTHER$count)
ab_tot_noOTHER$mediane_log_count <- ave(ab_tot_noOTHER$log_count,
                                        ab_tot_noOTHER$period,
                                        FUN = median)
ab_tot_noOTHER$ecart_mediane     <- ab_tot_noOTHER$log_count - ab_tot_noOTHER$mediane_log_count

heatmap_ab_station_period <- heatmap_station_periode(ab_tot_noOTHER,
                                                      x = "period", y = "location", fill = "ecart_mediane",
                                                      log10 = FALSE, divergent = TRUE, midpoint = 0,
                                                      xlab = "Period", ylab = "Station",
                                                      legend_lab = "Deviation from period median\n(log10)")

sauver_figure("fig_heatmap_abondance_station_periode.png", heatmap_ab_station_period,
              largeur = 8, hauteur = 6)


## one count ~ location + period model per taxon, reported in Table 3 by script 05
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

ab_taxon_mod <- droplevels(ab_taxon_comp)
ab_taxon_mod$location   <- factor(ab_taxon_mod$location)
ab_taxon_mod$period     <- factor(ab_taxon_mod$period, levels = levels(data_pollen$period))
ab_taxon_mod$Pred_Genus <- factor(ab_taxon_mod$Pred_Genus)

out_abrel_taxon  <- fit_taxon_models(ab_taxon_mod)
mod_abrel_taxon  <- out_abrel_taxon$models
skip_abrel_taxon <- out_abrel_taxon$skipped


## taxa shown on fig_abondance_par_taxon_station: the 12 most abundant
top12 <- top_taxons(mat_station, n = 12)


## one panel per taxon, horizontal station boxplots with the samples as points.

boxplot_taxon_station <- function(df, taxons, ordre_stations, titre) {

  df <- df[df$Pred_Genus %in% taxons, ]
  df$Pred_Genus <- factor(df$Pred_Genus, levels = taxons)
  df$location   <- factor(df$location, levels = rev(ordre_stations))  # rev: first station on top

  # per-panel mean, pre-aggregated: stat_summary() would summarise per station
  moy <- aggregate(count ~ Pred_Genus, data = df, FUN = mean)

  ggplot(df, aes(x = count, y = location)) +
    geom_vline(data = moy, aes(xintercept = count),
               linetype = "dashed", colour = "grey40", linewidth = 0.3) +
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

fig_taxon_station <- boxplot_taxon_station(ab_taxon_mod, top12, locations$trap,
                                           titre = "Pollen abundance per station and taxon")

sauver_figure("fig_abondance_par_taxon_station.png", fig_taxon_station,
              largeur = 11,
              hauteur = min(2.4 * ceiling(length(top12) / 3) + 1, 16))


## station boxplot with the period points on top and a global mean line per panel
# The mean is pre-aggregated: stat_summary() would draw one line per station.
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
# the means and letters are in Table 2.
box_ab_station <- boxplot_station_periodes(ab_tot_noOTHER, x = "location", y = "count",
                                  ylab = "Total abundance (number of pollen grains, log10)",
                                  titre = "",
                                  log10 = TRUE, ordre = "alpha")

sauver_figure("fig_abondance_totale_par_station.png", box_ab_station, largeur = 9, hauteur = 5)


## boxplot of the total abundance per period (period is already a factor in numeric order)
box_ab_period <- boxplot_station_periodes(ab_tot_noOTHER, x = "period", y = "count", xlab = "Period",
                                 ylab = "Total abundance (number of pollen grains, log10)",
                                 titre = "",
                                 log10 = TRUE, ordre = "alpha")

sauver_figure("fig_abondance_totale_par_period.png", box_ab_period, largeur = 9, hauteur = 5)


#### 2. SPATIAL VARIATION IN COMPOSITION ####

## Bray-Curtis dissimilarity + PERMANOVA

set.seed(GRAINE)
perm_comp <- adonis2(mat_rel_samples ~ period + location,
                     data = env_samples, method = "bray",
                     permutations = N_PERM, by = "terms")


## PCoA of the seasonal composition, with taxa arrows.
# Same taxa and metric as the PERMANOVA but aggregated at station level all periods combined

# drop the absent taxa (null columns make envfit fail)
mat_pcoa <- mat_rel_station[, colSums(mat_rel_station) > 0, drop = FALSE]

d_pcoa <- vegdist(mat_pcoa, method = "bray")
pcoa   <- cmdscale(d_pcoa, k = 2, eig = TRUE)

eig     <- pcoa$eig
var_pos <- sum(eig[eig > 0])
var_neg <- abs(sum(eig[eig < 0]))
pct_neg <- var_neg / var_pos * 100   # computed before any correction
cat("\nPCoA: negative eigenvalues =",
    round(pct_neg, 1), "% of the positive variance\n")

cailliez <- pct_neg > 10
if (cailliez) {
  cat("Cailliez correction applied\n")
  pcoa    <- wcmdscale(d_pcoa, k = 2, eig = TRUE, add = "cailliez")
  eig     <- pcoa$eig
  var_pos <- sum(eig[eig > 0])
}

coord <- pcoa$points
colnames(coord) <- c("Axe1", "Axe2")
pct <- eig[1:2] / var_pos * 100
cat("Axis 1 =", round(pct[1], 1), "% | Axis 2 =", round(pct[2], 1), "%\n")

sc_p          <- as.data.frame(coord)
sc_p$location <- rownames(coord)
sc_p$canopy_classe <- locations$canopy_classe[match(sc_p$location, locations$trap)]

# envfit: taxa correlated with the axes. Only the best-fitting significant
# ones are drawn, otherwise almost every taxon comes out significant.
fl     <- fleches_envfit(coord, mat_pcoa, sc_p, axes = c("Axe1", "Axe2"), top_n = 15)
ef_sig <- fl$sig

# stations as labelled points, coloured by canopy cover class
pcoa_arrows <- plot_ordination(
  sc_p, axes = c("Axe1", "Axe2"), id_col = "location", pies = FALSE,
  fleches = ef_sig,
  point_colour_col = "canopy_classe", point_palette = COULEURS_CANOPY,
  point_legend = "Canopy cover (%)",
  titre = "PCoA of the pollen composition per station (Bray-Curtis)",
  sous_titre = "Stations and taxa correlated with the axes (envfit, p < 0.05)",
  xlab = sprintf("Axis 1 (%.1f %%)", pct[1]),
  ylab = sprintf("Axis 2 (%.1f %%)", pct[2]))

sauver_figure("fig_pcoa_composition.png", pcoa_arrows, largeur = 7, hauteur = 7)

envfit_res <- fl$complet   # kept for script 05
pct_res    <- data.frame(axis1_pct = pct[1], axis2_pct = pct[2],   # kept for script 05
                         neg_eigen_pct = pct_neg, cailliez = cailliez)

## seasonal composition, top10 + "Others"
comp_long <- data.frame(
  location = rep(rownames(mat_pcoa), times = ncol(mat_pcoa)),
  taxon    = rep(colnames(mat_pcoa), each = nrow(mat_pcoa)),
  percent  = as.vector(mat_pcoa) * 100,
  stringsAsFactors = FALSE
)

comp_long$taxon_grp <- ifelse(comp_long$taxon %in% top10, comp_long$taxon, "Others")
comp_plot <- aggregate(percent ~ location + taxon_grp, data = comp_long, FUN = sum)
comp_plot$taxon_grp <- factor(comp_plot$taxon_grp, levels = c(top10, "Others"))

# dominant taxon of each station
dominants <- comp_plot[comp_plot$taxon_grp != "Others", ]
dominants <- dominants[order(dominants$location, -dominants$percent), ]
dominants <- dominants[!duplicated(dominants$location), ]


## stacked bars of the seasonal composition per station
pal_taxa <- palette_taxons(top10, COULEURS_TAXONS)

bar_comp_station <- ggplot(comp_plot, aes(x = location, y = percent, fill = taxon_grp)) +
  geom_bar(stat = "identity", colour = "grey30", linewidth = 0.15) +
  scale_fill_manual(values = pal_taxa, name = NULL) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(face = "italic")) +
  labs(x = "Station", y = "Seasonal relative abundance (%)",
       title = "Seasonal pollen composition per station")

sauver_figure("fig_composition_par_station.png", bar_comp_station, largeur = 10, hauteur = 7)


## diversity model
# Shannon computed per sample

shannon_vals <- shannon(mat_samples)

div_samples <- data.frame(
  location = factor(env_samples$location),
  period   = factor(env_samples$period),
  shannon  = as.numeric(shannon_vals[env_samples$row_id]),
  stringsAsFactors = FALSE
)

mod_compdiv <- glmmTMB(shannon ~ location + (1 | period),
                       data = div_samples, family = gaussian)
res_compdiv <- valider(mod_compdiv, "Diversity", tracer = FALSE, verbose = FALSE)
emm_compdiv <- emmeans(mod_compdiv, ~ location, type = "response")
cld_compdiv <- cld(emm_compdiv, adjust = "tukey", Letters = letters)


## boxplot of the Shannon index per station
box_div_station <- boxplot_station_periodes(div_samples, x = "location", y = "shannon",
                                   ylab = "Shannon index",
                                   titre = "",
                                   log10 = FALSE)

sauver_figure("fig_shannon_par_station.png", box_div_station, largeur = 9, hauteur = 5)



## Mantel test: geographic distance vs pollen distance.
# Coordinates projected to EPSG 32188 (MTM zone 8, metres)
locs_sf <- en_mtm(locations[, c("trap", "Longitude", "Latitude")])

stations <- rownames(mat_rel_station)
dist_geo <- distances_stations(locs_sf, stations)

# composition: Bray-Curtis between stations
dist_comp <- vegdist(mat_rel_station, method = "bray")
set.seed(GRAINE)
mantel_comp <- mantel(dist_comp, dist_geo, permutations = N_PERM)

# abundance: difference in total seasonal load, on a log10 scale
ab_saison <- setNames(ab_tot_station$count, ab_tot_station$location)
dist_ab   <- dist(log10(ab_saison[stations]))
set.seed(GRAINE)
mantel_ab <- mantel(dist_ab, dist_geo, permutations = N_PERM)

mantel_res <- data.frame(comparison = c("Composition (Bray-Curtis) vs geography",   # kept for script 05
                                        "Total abundance (log10) vs geography"),
                         mantel_r   = c(mantel_comp$statistic, mantel_ab$statistic),
                         p_value    = c(mantel_comp$signif, mantel_ab$signif),
                         stringsAsFactors = FALSE)
