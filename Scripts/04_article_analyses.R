# Analysis pipeline for the article: spatial distribution of pollen in Montreal.
# Two questions: does pollen abundance vary between stations, and does
# composition vary between stations ? 2022 and 2023 are always modelled separately.

source("R/init.R")
library(glmmTMB)
library(vegan)
library(emmeans)
library(multcompView)
library(multcomp)

dir.create(DIR_FIG_ARTICLE, recursive = TRUE, showWarnings = FALSE)

#### DATA PREPARATION ####

data_pollen <- charger_pollen(c(2022, 2023))
locations   <- charger_stations()

# canopy cover class per station, same classes and colours as the map (script 02)
gradients <- charger_gradients()
locations$canopy_classe <- classer_canopy(
  gradients$Canopy.cover....[match(locations$trap, gradients$Plot)])

annees <- sort(unique(data_pollen$year))

## samples actually collected, the reference grid for the true zeros
data_samples <- echantillons(data_pollen)
cat("Samples collected:", nrow(data_samples), "\n")

## counts per station x period x year x taxon, zeros included
ab_taxon <- comptages_pollen(data_pollen, data_samples)

## OTHER (unidentified grains) is excluded from every analysis below
ab_taxon_comp <- ab_taxon[ab_taxon$Pred_Genus != "OTHER", ]

ab_tot_noOTHER <- abondance_totale(ab_taxon_comp)
ab_tot_year    <- aggregate(count ~ location + year, data = ab_tot_noOTHER, FUN = sum)

## annual matrices, one per year, rows = stations
ab_year_taxon <- aggregate(count ~ location + year + Pred_Genus,
                           data = ab_taxon_comp, FUN = sum)

mat_2022 <- mat_station_taxon(ab_year_taxon[ab_year_taxon$year == "2022", ], "location", "Pred_Genus")
mat_2023 <- mat_station_taxon(ab_year_taxon[ab_year_taxon$year == "2023", ], "location", "Pred_Genus")

mat_rel_2022 <- en_relatif(mat_2022)
mat_rel_2023 <- en_relatif(mat_2023)

## sample-level matrix (station|period|year x taxon), for Shannon and the PERMANOVA
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

## top10 computed year by year, not pooled
top10_list <- list("2022" = top_taxons(mat_2022),
                   "2023" = top_taxons(mat_2023))



#### 1. SPATIAL VARIATION IN ABUNDANCE ####

## location and period are both fixed effects

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


## descriptive heatmap of the total abundance per station x period.
# Coloured by the deviation from its own period x year median 

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

sauver_figure("fig_heatmap_abondance_station_periode.png", heatmap_ab_station_period,
              largeur = 8, hauteur = 10)


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

ab_taxon_2022 <- droplevels(ab_taxon_comp[ab_taxon_comp$year == "2022", ])
ab_taxon_2022$location   <- factor(ab_taxon_2022$location)
ab_taxon_2022$period     <- factor(ab_taxon_2022$period)
ab_taxon_2022$Pred_Genus <- factor(ab_taxon_2022$Pred_Genus)

ab_taxon_2023 <- droplevels(ab_taxon_comp[ab_taxon_comp$year == "2023", ])
ab_taxon_2023$location   <- factor(ab_taxon_2023$location)
ab_taxon_2023$period     <- factor(ab_taxon_2023$period)
ab_taxon_2023$Pred_Genus <- factor(ab_taxon_2023$Pred_Genus)

out_abrel_taxon_2022 <- fit_taxon_models(ab_taxon_2022)
mod_abrel_taxon_2022  <- out_abrel_taxon_2022$models
skip_abrel_taxon_2022 <- out_abrel_taxon_2022$skipped

out_abrel_taxon_2023 <- fit_taxon_models(ab_taxon_2023)
mod_abrel_taxon_2023  <- out_abrel_taxon_2023$models
skip_abrel_taxon_2023 <- out_abrel_taxon_2023$skipped


## taxa shown on fig_abondance_par_taxon_station: the 12 most abundant of each year
top12_list <- list("2022" = top_taxons(mat_2022, n = 12),
                   "2023" = top_taxons(mat_2023, n = 12))


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

for (an in annees) {

  dat_an     <- get(paste0("ab_taxon_", an))
  taxons_fig <- top12_list[[an]]

  fig_taxon_station <- boxplot_taxon_station(dat_an, taxons_fig, locations$trap,
                                             titre = paste("Pollen abundance per station and taxon,", an))

  sauver_figure(paste0("fig_abondance_par_taxon_station_", an, ".png"), fig_taxon_station,
                largeur = 11,
                hauteur = min(2.4 * ceiling(length(taxons_fig) / 3) + 1, 16))
}


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
                                  facet = "year", log10 = TRUE, ordre = "alpha")

sauver_figure("fig_abondance_totale_par_station.png", box_ab_station, largeur = 9, hauteur = 8)


## boxplot of the total abundance per period
box_ab_period <- boxplot_station_periodes(ab_tot_noOTHER, x = "period", y = "count", xlab = "Period",
                                 ylab = "Total abundance (number of pollen grains, log10)",
                                 titre = "",
                                 facet = "year", log10 = TRUE, ordre = "alpha")

sauver_figure("fig_abondance_totale_par_period.png", box_ab_period, largeur = 9, hauteur = 8)


#### 2. SPATIAL VARIATION IN COMPOSITION ####

## Bray-Curtis dissimilarity + PERMANOVA, one per year.

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


## PCoA of the annual composition, with taxa arrows.
# Same taxa and metric as the PERMANOVA but aggregated at station level all periods combined

mat_rel_list <- list("2022" = mat_rel_2022, "2023" = mat_rel_2023)
envfit_list  <- list()   # kept for script 05
pct_list     <- list()   # kept for script 05
dominants    <- NULL
comp_plot_all <- NULL    # per-year composition, merged after the loop

for (an in annees) {

  # drop the taxa absent that year (null columns make envfit fail)
  mat_an <- mat_rel_list[[an]]
  mat_an <- mat_an[, colSums(mat_an) > 0, drop = FALSE]

  d_an    <- vegdist(mat_an, method = "bray")
  pcoa_an <- cmdscale(d_an, k = 2, eig = TRUE)

  eig     <- pcoa_an$eig
  var_pos <- sum(eig[eig > 0])
  var_neg <- abs(sum(eig[eig < 0]))
  pct_neg <- var_neg / var_pos * 100   # computed before any correction
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

  sc_p          <- as.data.frame(coord)
  sc_p$location <- rownames(coord)
  sc_p$row_id   <- paste0(sc_p$location, "|", an)
  sc_p$canopy_classe <- locations$canopy_classe[match(sc_p$location, locations$trap)]

  # envfit: taxa correlated with the axes. Only the best-fitting significant
  # ones are drawn, otherwise almost every taxon comes out significant.
  fl     <- fleches_envfit(coord, mat_an, sc_p, axes = c("Axe1", "Axe2"), top_n = 15)
  ef_df  <- fl$complet
  ef_sig <- fl$sig

  # stations as labelled points, coloured by canopy cover class
  pcoa_arrows <- plot_ordination(
    sc_p, axes = c("Axe1", "Axe2"), id_col = "location", pies = FALSE,
    fleches = ef_sig,
    point_colour_col = "canopy_classe", point_palette = COULEURS_CANOPY,
    point_legend = "Canopy cover (%)",
    titre = paste("PCoA of the pollen composition per station -", an, "(Bray-Curtis)"),
    sous_titre = "Stations and taxa correlated with the axes (envfit, p < 0.05)",
    xlab = sprintf("Axis 1 (%.1f %%)", pct[1]),
    ylab = sprintf("Axis 2 (%.1f %%)", pct[2]))

  sauver_figure(paste0("fig_pcoa_composition_", an, ".png"), pcoa_arrows, largeur = 7, hauteur = 7)

  envfit_list[[an]] <- ef_df
  pct_list[[an]]    <- data.frame(year = an, axis1_pct = pct[1], axis2_pct = pct[2],
                                  neg_eigen_pct = pct_neg, cailliez = cailliez)

  ## the year's composition, with its own top10 and its own "Others" group
  top10_an <- top10_list[[an]]

  comp_long_an <- data.frame(
    location = rep(rownames(mat_an), times = ncol(mat_an)),
    taxon    = rep(colnames(mat_an), each = nrow(mat_an)),
    percent  = as.vector(mat_an) * 100,
    stringsAsFactors = FALSE
  )

  comp_long_an$taxon_grp <- ifelse(comp_long_an$taxon %in% top10_an, comp_long_an$taxon, "Others")
  comp_plot_an <- aggregate(percent ~ location + taxon_grp, data = comp_long_an, FUN = sum)
  comp_plot_an$taxon_grp <- factor(comp_plot_an$taxon_grp, levels = c(top10_an, "Others"))

  # dominant taxon of each station, that year
  dominants_an <- comp_plot_an[comp_plot_an$taxon_grp != "Others", ]
  dominants_an <- dominants_an[order(dominants_an$location, -dominants_an$percent), ]
  dominants_an <- dominants_an[!duplicated(dominants_an$location), ]
  dominants_an$year <- an
  dominants <- rbind(dominants, dominants_an)

  comp_plot_an$year <- an
  comp_plot_all <- rbind(comp_plot_all, comp_plot_an)
}


## stacked bars of the annual composition per station, both years on one figure.
# Each year keeps its own top10/"Others" 
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

sauver_figure("fig_composition_par_station.png", bar_comp_station, largeur = 10, hauteur = 14)


## diversity model, one per year.
# Shannon computed per sample

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

mod_compdiv_2022 <- glmmTMB(shannon ~ location + (1 | period),
                            data = div_2022, family = gaussian)
res_compdiv_2022 <- valider(mod_compdiv_2022, "Diversity 2022", tracer = FALSE, verbose = FALSE)
emm_compdiv_2022 <- emmeans(mod_compdiv_2022, ~ location, type = "response")
cld_compdiv_2022 <- cld(emm_compdiv_2022, adjust = "tukey", Letters = letters)

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

sauver_figure("fig_shannon_par_station.png", box_div_station, largeur = 9, hauteur = 8)



## Mantel test: geographic distance vs pollen distance.
# Coordinates projected to EPSG 32188 (MTM zone 8, metres) 
locs_sf <- en_mtm(locations[, c("trap", "Longitude", "Latitude")])

mantel_res <- NULL   # accumulated across years, kept for script 05

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
