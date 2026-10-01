# Mise en forme, construction des tableaux et ecriture des sorties
# publiables (markdown + Word). Le contenu editorial (quel tableau, avec quel
# titre et quelle legende) reste dans le script qui appelle ces fonctions.

library(officer)
library(flextable)
library(car)          # Anova()
library(performance)  # r2()
library(DHARMa)       # testUniformity/testDispersion/testOutliers


#### MISE EN FORME ####

# p-values: sous 0.001 la valeur exacte n'apporte rien, on rapporte le seuil.
fmt_p <- function(p) {
  out <- character(length(p))
  for (i in seq_along(p)) {
    if (is.na(p[i])) {
      out[i] <- ""
    } else if (p[i] < 0.001) {
      out[i] <- "< 0.001"
    } else {
      out[i] <- formatC(p[i], format = "f", digits = 3)
    }
  }
  out
}

# Nombres: 3 chiffres significatifs, jamais de notation scientifique. Les
# comptages au-dela du millier sont arrondis et separes par des espaces, sinon
# les abondances polliniques sortent en 3.45e+04.
fmt_num <- function(x, digits = 3) {
  out <- character(length(x))
  for (i in seq_along(x)) {
    if (is.na(x[i])) {
      out[i] <- ""
    } else if (x[i] != 0 && abs(x[i]) >= 1000) {
      out[i] <- format(round(x[i]), big.mark = " ", scientific = FALSE, trim = TRUE)
    } else {
      out[i] <- format(signif(x[i], digits), scientific = FALSE, trim = TRUE)
    }
  }
  out
}

# Entiers (degres de liberte, k, nombre de stations).
fmt_int <- function(x) {
  out <- character(length(x))
  for (i in seq_along(x)) {
    out[i] <- if (is.na(x[i])) "" else format(round(x[i]), trim = TRUE)
  }
  out
}

# R2 marginal et conditionnel d'un modele glmmTMB, en chaines formatees, plus
# le type de R2 reellement utilise. Meme repli sur McFadden que valider()
# (R/models.R) quand r2() echoue; McFadden n'a pas de volet conditionnel.
r2_vals <- function(mod) {
  r <- try(r2(mod), silent = TRUE)
  if (!inherits(r, "try-error") && !is.null(r)) {
    return(c(fmt_num(as.numeric(r$R2_marginal)), fmt_num(as.numeric(r$R2_conditional)), "R2"))
  }
  rm <- try(r2_mcfadden_local(mod), silent = TRUE)
  if (inherits(rm, "try-error") || is.null(rm) || !is.finite(rm)) return(c("", "", ""))
  c(fmt_num(as.numeric(rm[1])), "", "McFadden")
}


#### CONSTRUCTION DES TABLEAUX ####

# Anova type II d'un modele glmmTMB + son R2, sur un bloc de lignes.
# L'etiquette du modele et le R2 ne sont ecrits que sur la 1re ligne, pour
# qu'un tableau empilant plusieurs modeles reste lisible.
tab_anova <- function(mod, etiquette, garder_p_num = FALSE) {
  a <- as.data.frame(Anova(mod, type = "II"))
  r <- r2_vals(mod)
  n <- nrow(a)
  # N obs: pour les modeles de comptage, les vrais zeros (completer_zeros(),
  # R/matrices.R) sont des lignes valides, donc nobs(mod) surestimerait le
  # nombre d'echantillons ayant reellement du pollen de ce taxon; on rapporte
  # les comptages non nuls. Pour les modeles gaussiens (Shannon), 0 est une
  # vraie valeur de diversite et nobs(mod) est rapporte tel quel.
  # La reponse passe par model.response(model.frame(mod)) et non
  # insight::get_response(), qui reevalue l'appel stocke et renvoie donc les
  # donnees de la derniere iteration pour des modeles ajustes en boucle.
  n_obs <- if (family(mod)$family == "gaussian") {
    nobs(mod)
  } else {
    sum(model.response(model.frame(mod)) != 0, na.rm = TRUE)
  }
  out <- data.frame(
    Model  = c(etiquette, rep("", n - 1)),
    Term   = rownames(a),
    Chi2   = fmt_num(a[["Chisq"]]),
    df     = fmt_int(a[["Df"]]),
    p      = fmt_p(a[["Pr(>Chisq)"]]),
    Nobs   = c(fmt_int(n_obs), rep("", n - 1)),
    R2m    = c(r[1], rep("", n - 1)),
    R2c    = c(r[2], rep("", n - 1)),
    R2type = c(r[3], rep("", n - 1)),
    stringsAsFactors = FALSE)
  names(out) <- c("Model", "Term", "Chi2", "df", "p", "N obs",
                  "R2 marginal", "R2 conditional", "R2 type")
  # p-value brute a cote de la colonne "p" deja formatee, pour les appelants
  # qui corrigent pour tests multiples sur plusieurs tab_anova() empiles.
  if (garder_p_num) out$p_num <- a[["Pr(>Chisq)"]]
  out
}

# emmeans + lettres de Tukey d'un modele. Les groupes partageant une lettre ne
# different pas significativement. `groupe_col` doit etre une colonne de
# as.data.frame(cld_obj), soit le terme sur lequel emmeans() a tourne.
tab_cld <- function(cld_obj, an, groupe_col = "location", nom_groupe = "Station") {
  d <- as.data.frame(cld_obj)
  # les noms de colonnes dependent de la famille du modele et de type = "response"
  est <- if ("response"  %in% names(d)) d$response  else d$emmean
  lcl <- if ("asymp.LCL" %in% names(d)) d$asymp.LCL else d$lower.CL
  ucl <- if ("asymp.UCL" %in% names(d)) d$asymp.UCL else d$upper.CL
  out <- data.frame(
    Year     = an,
    Groupe   = as.character(d[[groupe_col]]),
    Estimate = fmt_num(est),
    SE       = fmt_num(d$SE),
    Low      = fmt_num(lcl),
    High     = fmt_num(ucl),
    Group    = trimws(as.character(d$.group)),
    stringsAsFactors = FALSE)
  names(out) <- c("Year", nom_groupe, "Estimated mean (grains)", "SE",
                  "95% CI lower", "95% CI upper", "Tukey group")
  out
}

# Sortie adonis2 d'une annee.
tab_permanova <- function(res, an) {
  d <- as.data.frame(res)
  out <- data.frame(
    Year = c(an, rep("", nrow(d) - 1)),
    Term = rownames(d),
    Df   = fmt_int(d[["Df"]]),
    SS   = fmt_num(d[["SumOfSqs"]]),
    R2   = fmt_num(d[["R2"]]),
    F    = fmt_num(d[["F"]]),
    p    = fmt_p(d[["Pr(>F)"]]),
    stringsAsFactors = FALSE)
  names(out) <- c("Year", "Term", "df", "Sum of squares", "R2", "F", "p")
  out
}

# Diagnostics DHARMa d'un modele, sur une seule ligne.
tab_diag <- function(res, etiquette) {
  data.frame(
    Model      = etiquette,
    Uniformity = fmt_p(testUniformity(res, plot = FALSE)$p.value),
    Dispersion = fmt_p(testDispersion(res, plot = FALSE)$p.value),
    Outliers   = fmt_p(testOutliers(res, plot = FALSE)$p.value),
    stringsAsFactors = FALSE)
}

# Ligne de remplacement, au format de tab_anova(), pour un modele qui n'a pas pu
# etre ajuste: la ligne reste visible avec une note au lieu d'etre supprimee.
note_row <- function(etiquette, note) {
  data.frame(Model = etiquette, Term = "", Chi2 = "", df = "", p = "",
             `N obs` = "", `R2 marginal` = "", `R2 conditional` = "", `R2 type` = "",
             Note = note, check.names = FALSE, stringsAsFactors = FALSE)
}


#### ECRITURE ####

# markdown: une section par tableau, syntaxe pipe
ecrire_md <- function(tables, figures, fichier) {

  cat("# Model results - spatial distribution of pollen in Montreal\n\n",
      "_Generated on ", format(Sys.time(), "%Y-%m-%d %H:%M"), "_\n\n",
      sep = "", file = fichier)

  for (t in tables) {
    cat("\n## ", t$titre, "\n\n", sep = "", file = fichier, append = TRUE)
    cat(t$legende, "\n\n", sep = "", file = fichier, append = TRUE)

    d <- t$df
    cat("| ", paste(names(d), collapse = " | "), " |\n",
        "| ", paste(rep("---", ncol(d)), collapse = " | "), " |\n",
        sep = "", file = fichier, append = TRUE)
    for (i in seq_len(nrow(d))) {
      valeurs <- as.character(unlist(d[i, ]))
      valeurs[is.na(valeurs)] <- ""
      cat("| ", paste(valeurs, collapse = " | "), " |\n",
          sep = "", file = fichier, append = TRUE)
    }
  }

  if (!is.null(figures)) {
    cat("\n## Figures\n\n", sep = "", file = fichier, append = TRUE)
    for (i in seq_len(nrow(figures))) {
      cat("**Figure ", i, ".** ", figures$legende[i], "\n\n",
          "![Figure ", i, "](figures/", figures$fichier[i], ")\n\n",
          sep = "", file = fichier, append = TRUE)
    }
  }
}

# Word: officer assemble le document, flextable met en forme chaque tableau
ecrire_docx <- function(tables, figures, fichier) {

  doc <- read_docx()
  doc <- body_add_par(doc, "Model results - spatial distribution of pollen in Montreal",
                      style = "heading 1")
  doc <- body_add_par(doc, paste("Generated on", format(Sys.time(), "%Y-%m-%d %H:%M")),
                      style = "Normal")

  for (t in tables) {
    doc <- body_add_par(doc, t$titre, style = "heading 2")
    doc <- body_add_par(doc, t$legende, style = "Normal")

    ft <- flextable(t$df)
    ft <- theme_booktabs(ft)
    ft <- fontsize(ft, size = 9, part = "all")
    ft <- bold(ft, part = "header")
    ft <- autofit(ft)
    doc <- body_add_flextable(doc, ft)
    doc <- body_add_par(doc, "", style = "Normal")
  }

  if (!is.null(figures)) {
    # les figures sont livrees en fichiers separes a 300 dpi; le document les liste
    doc <- body_add_par(doc, "Figures", style = "heading 2")
    for (i in seq_len(nrow(figures))) {
      doc <- body_add_par(doc, paste0("Figure ", i, ". ", figures$legende[i],
                                      " (", figures$fichier[i], ")"),
                          style = "Normal")
    }
  }

  print(doc, target = fichier)
}
