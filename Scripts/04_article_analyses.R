# Pipeline d'analyse de l'article: distribution spatiale du pollen a Montreal.
# Trois questions:
#   1. l'abondance pollinique varie-t-elle entre les stations ?
#   2. la composition pollinique varie-t-elle entre les stations ?
#   3. combien de stations, et lesquelles, pour un reseau representatif ?
# 2022 et 2023 sont toujours modelisees separement.

source("R/init.R")
library(glmmTMB)
library(vegan)
library(cluster)
library(emmeans)
library(multcompView)
library(multcomp)   # fournit cld(); emmeans n'enregistre que la methode

dir.create(DIR_FIG_ARTICLE, recursive = TRUE, showWarnings = FALSE)

#### PREPARATION DES DONNEES ####

## pollen, deja filtre a Confidence >= CONF_MIN par 01_preparer_pollen.R
data_pollen <- charger_pollen(c(2022, 2023))
locations   <- charger_stations()

# classe de couvert forestier par station, memes 4 classes et couleurs que la
# carte des stations (script 02) -- sert a colorer les points de la PCoA (fig2)
gradients <- charger_gradients()
locations$canopy_classe <- classer_canopy(
  gradients$Canopy.cover....[match(locations$trap, gradients$Plot)])

annees <- sort(unique(data_pollen$year))

## echantillons reellement collectes, grille de reference pour les vrais zeros
data_samples <- echantillons(data_pollen)
cat("Samples collected:", nrow(data_samples), "\n")

## comptages par station x periode x annee x taxon, zeros inclus.
# OTHER est conserve ici car il fait partie de la charge pollinique totale; il
# est retire plus bas partout ou il s'agit de composition.
ab_taxon <- comptages_pollen(data_pollen, data_samples)

## charge pollinique totale, OTHER inclus
ab_tot      <- abondance_totale(ab_taxon)
ab_tot_year <- aggregate(count ~ location + year, data = ab_tot, FUN = sum)

## tables de composition, OTHER exclu
ab_taxon_comp <- ab_taxon[ab_taxon$Pred_Genus != "OTHER", ]

## charge totale sans OTHER, utilisee dans la section 1
ab_tot_noOTHER <- abondance_totale(ab_taxon_comp)

# matrices annuelles, une par annee, lignes = stations
ab_year_taxon <- aggregate(count ~ location + year + Pred_Genus,
                           data = ab_taxon_comp, FUN = sum)

mat_2022 <- mat_station_taxon(ab_year_taxon[ab_year_taxon$year == "2022", ], "location", "Pred_Genus")
mat_2023 <- mat_station_taxon(ab_year_taxon[ab_year_taxon$year == "2023", ], "location", "Pred_Genus")

mat_rel_2022 <- en_relatif(mat_2022)
mat_rel_2023 <- en_relatif(mat_2023)

# matrice au grain de l'echantillon (station|periode|annee x taxon), pour
# Shannon et la PERMANOVA. Tous les taxons sauf OTHER: la PERMANOVA n'est pas
# restreinte au top10, contrairement aux modeles par taxon plus bas.
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

## top10 par annee, calcule annee par annee (non groupe) pour que les modeles
## et figures de chaque annee refletent ses propres taxons dominants
top10_2022 <- top_taxons(mat_2022)
top10_2023 <- top_taxons(mat_2023)
top10_list <- list("2022" = top10_2022, "2023" = top10_2023)



#### 1. VARIATION SPATIALE DE L'ABONDANCE ####
## OTHER exclu. location et period sont tous deux des effets fixes
## (count ~ location + period, sans terme aleatoire): un seul modele fournit
## ainsi les lettres de Tukey des stations et celles des periodes.

ab_tot_noOTHER$location <- factor(ab_tot_noOTHER$location)
ab_tot_noOTHER$period   <- factor(ab_tot_noOTHER$period)

ab_tot_noOTHER_2022 <- droplevels(ab_tot_noOTHER[ab_tot_noOTHER$year == "2022", ])
ab_tot_noOTHER_2023 <- droplevels(ab_tot_noOTHER[ab_tot_noOTHER$year == "2023", ])

## modele d'abondance totale, un par annee
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


## heatmap descriptive de l'abondance totale par station x periode.
## Un modele ajuste periode par periode se reduirait a une seule valeur par
## station x annee (plus aucune replication a tester), d'ou cette heatmap
## plutot qu'un modele. Coloree par l'ecart a la mediane de SA periode x annee
## (via ave()) plutot qu'a une mediane globale: cela isole la variation entre
## stations au lieu d'etre dominee par le signal saisonnier.
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


## un modele count ~ location + period par taxon (tous les taxons, OTHER
## exclu). Leurs tests omnibus sont rapportes au tableau 3 par le script 05;
## ces tests disent si l'abondance d'un taxon varie entre stations mais pas
## comment, d'ou la fig1d plus bas qui montre la distribution par station,
## taxon par taxon.
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
mod_abrel_taxon_2022  <- out_abrel_taxon_2022$models    # liste nommee, 1 glmmTMB par taxon
skip_abrel_taxon_2022 <- out_abrel_taxon_2022$skipped   # taxon -> raison de non-ajustement

out_abrel_taxon_2023 <- fit_taxon_models(ab_taxon_2023)
mod_abrel_taxon_2023  <- out_abrel_taxon_2023$models
skip_abrel_taxon_2023 <- out_abrel_taxon_2023$skipped


## taxons montres sur la fig1d: les 12 plus abondants de chaque annee. Les ~33
## taxons x 25 stations sur une figure seraient illisibles.
top12_list <- list("2022" = top_taxons(mat_2022, n = 12),
                   "2023" = top_taxons(mat_2023, n = 12))


## fig1d: un panneau par taxon, boxplots horizontaux par station avec les
## echantillons individuels en points. Comptages bruts, soit la distribution
## sur laquelle chaque modele par taxon a ete ajuste.
## scales = "free_x" est necessaire: les comptages par taxon couvrent plusieurs
## ordres de grandeur, un axe commun aplatirait tous les panneaux sauf un.
## Transformation log1p (ln(1+x)) et non log10: un taxon absent d'un
## echantillon collecte est un vrai zero, et log10(0) supprimerait justement
## les points montrant qu'une station n'a aucun grain de ce taxon.
## L'ordre des stations est l'ordre canonique (locations$trap) et non un tri
## par mediane: avec un ordre different par panneau, comparer la position d'une
## station d'un taxon a l'autre serait impossible.
boxplot_taxon_station <- function(df, taxons, ordre_stations, titre) {

  df <- df[df$Pred_Genus %in% taxons, ]
  df$Pred_Genus <- factor(df$Pred_Genus, levels = taxons)   # taxon le plus abondant en premier
  df$location   <- factor(df$location, levels = rev(ordre_stations))  # rev: 1re station en haut

  # moyenne par panneau, pre-agregee et passee via le `data` de geom_vline():
  # une couche stat_summary() resumerait par categorie y (une ligne par
  # station) au lieu d'une par panneau.
  moy <- aggregate(count ~ Pred_Genus, data = df, FUN = mean)

  ggplot(df, aes(x = count, y = location)) +
    geom_vline(data = moy, aes(xintercept = count),
               linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    # outlier.shape = NA: chaque point est deja trace individuellement plus bas
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


## Boxplot par station avec les 7 points de periode en surimpression et une
## ligne de moyenne globale par panneau d'annee (fig1 et fig3).
# La moyenne est pre-agregee (une ligne par niveau de facette) plutot que
# laissee au regroupement implicite de stat_summary(), qui avec un x discret
# resumait par categorie x (une ligne par station au lieu d'une par panneau).
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


## boxplot de l'abondance totale par station
# echelle log10 (les comptages couvrent plusieurs ordres de grandeur), axe x en
# ordre naturel, pas de lettres de Tukey sur la figure: les moyennes et lettres
# par station sont au tableau 2.
box_ab_station <- boxplot_station_periodes(ab_tot_noOTHER, x = "location", y = "count",
                                  ylab = "Total abundance (number of pollen grains, log10)",
                                  titre = "",
                                  facet = "year", log10 = TRUE, ordre = "alpha")

sauver_figure("fig1_abondance_totale_par_station.png", box_ab_station, largeur = 9, hauteur = 8)


## boxplot de l'abondance totale par periode
box_ab_period <- boxplot_station_periodes(ab_tot_noOTHER, x = "period", y = "count", xlab = "Period",
                                 ylab = "Total abundance (number of pollen grains, log10)",
                                 titre = "",
                                 facet = "year", log10 = TRUE, ordre = "alpha")

sauver_figure("fig1b_abondance_totale_par_period.png", box_ab_period, largeur = 9, hauteur = 8)


#### 2. VARIATION SPATIALE DE LA COMPOSITION ####

## dissimilarite de Bray-Curtis + PERMANOVA, une par annee.
# Le test tourne sur la matrice au grain de l'echantillon et non sur
# l'annuelle: avec 25 lignes et 25 niveaux de station, le modele annuel serait
# sature. period vient en premier dans le modele car la succession saisonniere
# est de loin le signal le plus fort et absorberait sinon l'effet station; avec
# by = "terms", location est teste APRES retrait de l'effet periode.
# Ne pas ajouter blocks = period: sous permutation intra-bloc, la matrice de
# design de period est invariante et les deux termes sortent avec la meme
# p-value, denuee de sens.

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


## PCoA de la composition annuelle, avec fleches de taxons.
# Memes taxons et meme methode Bray-Curtis que la PERMANOVA, mais agregee au
# grain de la station (une ligne par station): c'est ce qui est lisible sur un
# graphique, alors que la PERMANOVA a besoin de la replication au grain de
# l'echantillon pour tester location. Bray-Curtis est semi-metrique, donc la
# PCoA produit des valeurs propres negatives: elles sont quantifiees, et la
# correction de Cailliez est appliquee si elles depassent 10 % de la variance
# positive.

mat_rel_list <- list("2022" = mat_rel_2022, "2023" = mat_rel_2023)
pcoa_list    <- list()
envfit_list  <- list()   # ef_df de chaque annee, conserve pour le script 05
pct_list     <- list()   # pourcentages d'axes de chaque annee, meme raison
dominants    <- NULL     # taxon dominant par station, accumule sur les annees
comp_plot_all <- NULL    # composition par annee, accumulee pour la fig4 fusionnee

for (an in annees) {

  # retirer les taxons absents cette annee-la (les colonnes nulles font echouer envfit)
  mat_an <- mat_rel_list[[an]]
  mat_an <- mat_an[, colSums(mat_an) > 0, drop = FALSE]

  d_an    <- vegdist(mat_an, method = "bray")
  pcoa_an <- cmdscale(d_an, k = 2, eig = TRUE)

  eig     <- pcoa_an$eig
  var_pos <- sum(eig[eig > 0])
  var_neg <- abs(sum(eig[eig < 0]))
  # part des valeurs propres negatives, calculee avant toute correction
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

  # envfit: taxons correles aux axes. Parmi les significatifs, seuls les mieux
  # ajustes aux axes (r2 envfit le plus eleve) sont traces, sinon presque tous
  # les taxons sortent significatifs et le graphique devient illisible.
  fl     <- fleches_envfit(coord, mat_an, sc_p, axes = c("Axe1", "Axe2"), top_n = 15)
  ef_df  <- fl$complet
  ef_sig <- fl$sig

  # Stations en points etiquetes: les camemberts rendaient le graphique
  # illisible une fois les fleches de taxons ajoutees. Les points sont colores
  # par classe de couvert forestier (memes classes que la carte du script 02)
  # pour verifier a l'oeil si la composition s'aligne sur le couvert.
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

  # conserve pour le script 05, qui ne verrait sinon que la derniere annee
  envfit_list[[an]] <- ef_df
  pct_list[[an]]    <- data.frame(year = an, axis1_pct = pct[1], axis2_pct = pct[2],
                                  neg_eigen_pct = pct_neg, cailliez = cailliez)

  ## composition de l'annee, gardee par annee (son propre top10 et son propre
  # groupe "Others") mais accumulee pour que les deux annees soient tracees sur
  # une figure fusionnee avec une palette coherente par taxon (voir apres la
  # boucle), au lieu de deux fichiers aux couleurs independantes.
  top10_an <- top10_list[[an]]

  comp_long_an <- data.frame(
    location = rep(rownames(mat_an), times = ncol(mat_an)),
    taxon    = rep(colnames(mat_an), each = nrow(mat_an)),
    percent  = as.vector(mat_an) * 100,
    stringsAsFactors = FALSE
  )

  # les taxons hors du top10 de l'annee sont regroupes dans "Others"
  comp_long_an$taxon_grp <- ifelse(comp_long_an$taxon %in% top10_an, comp_long_an$taxon, "Others")
  comp_plot_an <- aggregate(percent ~ location + taxon_grp, data = comp_long_an, FUN = sum)
  comp_plot_an$taxon_grp <- factor(comp_plot_an$taxon_grp, levels = c(top10_an, "Others"))

  # taxon dominant de chaque station, cette annee-la
  dominants_an <- comp_plot_an[comp_plot_an$taxon_grp != "Others", ]
  dominants_an <- dominants_an[order(dominants_an$location, -dominants_an$percent), ]
  dominants_an <- dominants_an[!duplicated(dominants_an$location), ]
  dominants_an$year <- an
  dominants <- rbind(dominants, dominants_an)

  comp_plot_an$year <- an
  comp_plot_all <- rbind(comp_plot_all, comp_plot_an)
}


## barres empilees de la composition annuelle par station, les deux annees sur
# une figure (facettes par annee). Chaque annee garde son propre decoupage
# top10/"Others" calcule plus haut; ce qui est unifie ici est seulement la
# palette, via l'union des top10 des deux annees, ordonnee par abondance
# combinee pour que la legende se lise du dominant au rare.
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


## modele de diversite, un par annee.
# Shannon calcule par echantillon sur la matrice de composition (OTHER exclu),
# donc chaque station a jusqu'a 7 periodes de replication dans une annee.

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


## boxplot de l'indice de Shannon par station
box_div_station <- boxplot_station_periodes(div_samples, x = "location", y = "shannon",
                                   ylab = "Shannon index",
                                   titre = "",
                                   facet = "year", log10 = FALSE)

sauver_figure("fig3_shannon_par_station.png", box_div_station, largeur = 9, hauteur = 8)



#### 3. COMBIEN DE STATIONS, ET LESQUELLES ####

## k-medoides (PAM) sur la composition annuelle.
# k balaye de 2 a 10, on garde le k de plus forte silhouette moyenne. Sous
# 0.25, la silhouette indique qu'il n'y a pas de vraie structure de groupes.
# Tourne sur la matrice d'abondances relatives, la meme que la PCoA: sur des
# comptages bruts, Bray-Curtis n'est pas invariant a l'echelle et les stations
# se separeraient par charge pollinique totale plutot que par composition.
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

  # k = 3 conserve a cote du k optimal, meme sous le seuil de 0.25: c'est le
  # nombre de groupes parfois suggere a l'oeil sur la PCoA (fig2), donc
  # rapporte au tableau 9 pour comparaison quel que soit le k retenu.
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
       courbe = data.frame(k = k_range, sil = sils))
}

pam_2022 <- run_pam(mat_rel_2022, "2022")
pam_2023 <- run_pam(mat_rel_2023, "2023")

## station "moyenne" unique.
# Question differente de PAM: PAM groupe les stations et donne un medoide par
# groupe; ici on demande a quoi ressemblerait un reseau d'une seule station,
# soit laquelle a la composition annuelle la plus proche de la composition
# moyenne de toutes les stations cette annee-la. Meme metrique Bray-Curtis,
# obtenue en ajoutant la moyenne comme ligne supplementaire a la matrice.
station_moyenne <- function(mat, etiquette) {

  moyenne <- colMeans(mat)   # composition relative moyenne, une station = un poids

  mat_avec_moyenne <- rbind(mat, moyenne = moyenne)
  d         <- as.matrix(vegdist(mat_avec_moyenne, method = "bray"))
  d_moyenne <- sort(d["moyenne", rownames(mat)])   # distance a la moyenne, plus proche d'abord

  station_rep <- names(d_moyenne)[1]

  cat("\nStation closest to the average composition -", etiquette, ":",
      station_rep, "(Bray-Curtis distance to the average =",
      round(d_moyenne[1], 3), ")\n")

  list(station = station_rep, distance = d_moyenne[1], distances = d_moyenne)
}

station_moy_2022 <- station_moyenne(mat_rel_2022, "2022")
station_moy_2023 <- station_moyenne(mat_rel_2023, "2023")

## courbes de silhouette
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

## carte des groupes et des medoides
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


## test de Mantel: distance geographique vs distance pollinique.
# Les coordonnees sont projetees en EPSG 32188 (MTM zone 8, metres) avant le
# calcul des distances.
locs_sf <- en_mtm(locations[, c("trap", "Longitude", "Latitude")])

mantel_res <- NULL   # accumule sur les annees, conserve pour le script 05

for (an in annees) {

  mat_an   <- mat_rel_list[[an]]
  stations <- rownames(mat_an)

  dist_geo <- distances_stations(locs_sf, stations)

  # composition: Bray-Curtis entre stations
  dist_comp <- vegdist(mat_an, method = "bray")
  set.seed(GRAINE)
  mantel_comp <- mantel(dist_comp, dist_geo, permutations = N_PERM)

  # abondance: ecart de charge annuelle totale, en echelle log10
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
