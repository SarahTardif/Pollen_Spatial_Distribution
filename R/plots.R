# Shared plotting functions: taxon palettes, envfit arrows, ordinations,
# heatmaps, figure saving.

library(ggplot2)
library(scatterpie)
library(vegan)
library(ggrepel)

# Named colour vector for the "top" taxa + a catch-all group. 
palette_taxons <- function(top, couleurs = COULEURS_TAXONS,
                           autres_label = "Others", autres_couleur = COULEUR_AUTRES) {
  if (length(top) > length(couleurs)) {
    couleurs <- c(couleurs, grDevices::hcl.colors(length(top) - length(couleurs), palette = "Dark 3"))
  }
  setNames(c(couleurs[seq_along(top)], autres_couleur), c(top, autres_label))
}

# envfit of the taxa onto the ordination axes

fleches_envfit <- function(ord, mat, sc, axes, seuil_p = 0.05, taxons_filtre = NULL,
                           top_n = NULL, permutations = N_PERM, graine = GRAINE) {
  set.seed(graine)
  ef    <- envfit(ord, mat, permutations = permutations)
  ef_df <- as.data.frame(scores(ef, display = "vectors"))
  ef_df$p     <- ef$vectors$pvals
  ef_df$r2    <- ef$vectors$r  # vegan names it "r" but it is already the squared r
  ef_df$taxon <- rownames(ef_df)

  garder <- ef_df$p < seuil_p
  if (!is.null(taxons_filtre)) garder <- garder & ef_df$taxon %in% taxons_filtre
  ef_sig <- ef_df[garder, ]

  # top_n keeps the taxa that best fit the axes (highest envfit r2)
  if (!is.null(top_n) && nrow(ef_sig) > top_n) {
    ef_sig <- ef_sig[order(-ef_sig$r2), ][seq_len(top_n), ]
  }

  if (nrow(ef_sig) > 0) {
    ax1 <- axes[1]; ax2 <- axes[2]
    rng   <- max(diff(range(sc[[ax1]])), diff(range(sc[[ax2]]))) / 2
    a_scl <- rng * FLECHE_SCALE / max(sqrt(ef_sig[[ax1]]^2 + ef_sig[[ax2]]^2))
    ef_sig$x1 <- ef_sig[[ax1]] * a_scl
    ef_sig$y1 <- ef_sig[[ax2]] * a_scl
  }

  list(complet = ef_df, sig = ef_sig)
}

# Ordination plot
plot_ordination <- function(sc, axes, id_col, pies = FALSE,
                            pie_cols = NULL, pie_r_col = "r", palette = NULL,
                            fleches = NULL, label_offset = LABEL_OFFSET,
                            label_y_offset = 0.007,
                            titre = NULL, sous_titre = NULL, stress = NULL,
                            xlab = axes[1], ylab = axes[2],
                            point_colour_col = NULL, point_palette = NULL,
                            point_legend = point_colour_col) {

  ax1 <- axes[1]; ax2 <- axes[2]

  p <- ggplot() +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey75", linewidth = 0.3)

  # Station labels 
  station_labels <- NULL

  if (pies) {
    p <- p +
      geom_scatterpie(data = sc, aes(x = .data[[ax1]], y = .data[[ax2]], r = .data[[pie_r_col]]),
                       cols = pie_cols, colour = "grey25", linewidth = 0.25) +
      geom_label_repel(data = sc, aes(x = .data[[ax1]], y = .data[[ax2]] + .data[[pie_r_col]] + label_y_offset,
                                label = .data[[id_col]]),
                 size = 3, label.padding = unit(0.12, "lines"), label.size = 0.2, fill = "white",
                 seed = GRAINE, segment.color = NA) +
      scale_fill_manual(values = palette, name = NULL)
  } else if (!is.null(point_colour_col)) {
    p <- p +
      geom_point(data = sc, aes(x = .data[[ax1]], y = .data[[ax2]], fill = .data[[point_colour_col]]),
                shape = 21, size = 2.6, colour = "grey30", stroke = 0.4) +
      scale_fill_manual(values = point_palette, name = point_legend)
    station_labels <- data.frame(x = sc[[ax1]], y = sc[[ax2]], label = sc[[id_col]],
                                 fontface = "plain", colour = "black", size = 2.8,
                                 stringsAsFactors = FALSE)
  } else {
    p <- p +
      geom_point(data = sc, aes(x = .data[[ax1]], y = .data[[ax2]]), size = 1.4, colour = "grey40")
    station_labels <- data.frame(x = sc[[ax1]], y = sc[[ax2]], label = sc[[id_col]],
                                 fontface = "plain", colour = "black", size = 2.8,
                                 stringsAsFactors = FALSE)
  }

  if (!is.null(fleches) && nrow(fleches) > 0) {
    arrow_len <- if (pies) 0.22 else 0.2
    arrow_lwd <- if (pies) 0.7  else 0.6
    txt_size  <- if (pies) 4    else 3.5
    p <- p +
      geom_segment(data = fleches, aes(x = 0, y = 0, xend = x1, yend = y1),
                   arrow = arrow(length = unit(arrow_len, "cm"), type = "open"),
                   colour = "grey35", linewidth = arrow_lwd)

    if (is.null(station_labels)) {
      p <- p +
        geom_text_repel(data = fleches, aes(x = x1 * label_offset, y = y1 * label_offset, label = taxon),
                  size = txt_size, fontface = "italic", colour = "grey10",
                  seed = GRAINE, segment.color = NA)
    } else {
      taxon_labels <- data.frame(x = fleches$x1 * label_offset, y = fleches$y1 * label_offset,
                                 label = fleches$taxon, fontface = "italic", colour = "grey10",
                                 size = txt_size, stringsAsFactors = FALSE)

      
      t_interior <- c(0.25, 0.5, 0.75, 0.9)
      arrow_obstacles <- data.frame(
        x = as.vector(outer(t_interior, fleches$x1)),
        y = as.vector(outer(t_interior, fleches$y1)))

      label_data <- rbind(station_labels, taxon_labels)

      p <- p +
        geom_point(data = arrow_obstacles, aes(x = x, y = y), alpha = 0, size = 0.1) +
        geom_text_repel(data = label_data,
                  aes(x = x, y = y, label = label, fontface = fontface, colour = colour, size = size),
                  seed = GRAINE, segment.color = NA, max.overlaps = Inf) +
        scale_colour_identity() + scale_size_identity()
    }
  } else if (!is.null(station_labels)) {
    p <- p +
      geom_text_repel(data = station_labels,
                aes(x = x, y = y, label = label, fontface = fontface, colour = colour, size = size),
                seed = GRAINE, segment.color = NA, max.overlaps = Inf) +
      scale_colour_identity() + scale_size_identity()
  }

  if (!is.null(stress)) {
    p <- p + annotate("text", x = Inf, y = Inf,
                      label = paste0("Stress = ", round(stress, 2)),
                      hjust = 1.05, vjust = 1.5, size = 4.5)
  }

  p <- p + coord_fixed()

  if (pies) {
    p <- p + theme_bw(base_size = 13) +
      theme(legend.text     = element_text(face = "italic", size = 11),
            legend.key.size = unit(0.55, "cm"),
            panel.grid      = element_blank(),
            axis.title      = element_text(size = 13))
  } else {
    p <- p + theme_bw() + theme(panel.grid = element_blank())
  }

  p + labs(x = xlab, y = ylab, title = titre, subtitle = sous_titre)
}


lettres_layer <- function(df, x, y, lettres, facet = NULL, nudge = 1.1) {
  by_cols <- if (!is.null(facet)) c(x, facet) else x
  y_max <- aggregate(df[y], by = as.list(df[by_cols]), FUN = max)
  names(y_max) <- c(by_cols, "y_pos")
  y_max$y_pos <- y_max$y_pos * nudge
  pos <- merge(lettres, y_max, by = by_cols)
  geom_text(data = pos, aes(x = .data[[x]], y = y_pos, label = label),
            inherit.aes = FALSE, vjust = 0)
}

# Heatmap 
heatmap_station_periode <- function(df, x = "period", y = "location", fill = "count",
                                     facet = NULL, facet_nrow = 1, log10 = FALSE, divergent = FALSE,
                                     midpoint = NULL, angle_x = 45, xlab = "Period", ylab = "Station",
                                     titre = NULL, legend_lab = NULL) {
  df$valeur_remplissage <- if (log10) log10(df[[fill]]) else df[[fill]]

  if (is.null(legend_lab)) legend_lab <- if (log10) paste0(fill, " (log10)") else fill
  if (is.null(midpoint))   midpoint   <- median(df$valeur_remplissage, na.rm = TRUE)

  p <- ggplot(df, aes(x = .data[[x]], y = .data[[y]], fill = valeur_remplissage)) +
    geom_tile(colour = "white")

  p <- if (divergent) {
    p + scale_fill_gradient2(name = legend_lab,
                      low = "#2166AC", mid = "grey95", high = "#B2182B",
                      midpoint = midpoint)
  } else {
    p + scale_fill_viridis_c(name = legend_lab)
  }

  if (!is.null(facet)) p <- p + facet_wrap(stats::as.formula(paste0("~", facet)), nrow = facet_nrow)

  p + coord_fixed(ratio = 1) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = angle_x, hjust = 1)) +
    labs(x = xlab, y = ylab, title = titre)
}

# ggsave wrapper: creates the output directory if needed, fixes dpi = 300.
sauver_figure <- function(nom, plot, dir = DIR_FIG_ARTICLE, largeur = 8, hauteur = 8, dpi = 300) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(dir, nom), plot, width = largeur, height = hauteur, dpi = dpi)
}
