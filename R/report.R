# Formatting, table building and writing of the publication-ready output
# (markdown + Word). The editorial content (which table, with which title and
# legend) stays in the script that calls these functions.

library(officer)
library(flextable)
library(car)          # Anova()
library(performance)  # r2()
library(DHARMa)       # testUniformity/testDispersion/testOutliers


#### FORMATTING ####

# p-values: below 0.001 the exact value carries no information, so the
# threshold is reported instead.
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

# Numbers: 3 significant digits, never scientific notation. Counts above a
# thousand are rounded and thousand-separated instead, otherwise the pollen
# abundances come out as 3.45e+04.
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

# Integers (degrees of freedom, k, number of stations).
fmt_int <- function(x) {
  out <- character(length(x))
  for (i in seq_along(x)) {
    out[i] <- if (is.na(x[i])) "" else format(round(x[i]), trim = TRUE)
  }
  out
}

# Marginal and conditional R2 of a glmmTMB model, as formatted strings, plus
# the R2 type actually used. Same McFadden fallback as valider() (R/models.R)
# when r2() fails; McFadden has no conditional counterpart.
r2_vals <- function(mod) {
  r <- try(r2(mod), silent = TRUE)
  if (!inherits(r, "try-error") && !is.null(r)) {
    return(c(fmt_num(as.numeric(r$R2_marginal)), fmt_num(as.numeric(r$R2_conditional)), "R2"))
  }
  rm <- try(r2_mcfadden_local(mod), silent = TRUE)
  if (inherits(rm, "try-error") || is.null(rm) || !is.finite(rm)) return(c("", "", ""))
  c(fmt_num(as.numeric(rm[1])), "", "McFadden")
}


#### TABLE BUILDERS ####

# Type II Anova of a glmmTMB model + its R2, as one block of rows. The model
# label and the R2 are only written on the first row so that a table stacking
# several models stays readable.
tab_anova <- function(mod, etiquette, garder_p_num = FALSE) {
  a <- as.data.frame(Anova(mod, type = "II"))
  r <- r2_vals(mod)
  n <- nrow(a)
  # N obs: for count models the true zeros (completer_zeros(), R/matrices.R)
  # are valid rows, so nobs(mod) would overstate how many samples actually had
  # pollen of that taxon; the non-zero counts are reported instead. For
  # gaussian models (Shannon), 0 is a real diversity value and nobs(mod) is
  # reported as-is.
  # The response goes through model.response(model.frame(mod)) and not
  # insight::get_response(), which re-evaluates the stored call and therefore
  # returns the last iteration's data for models fit in a loop.
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
  # raw p-value alongside the already-formatted "p" column, for callers that
  # correct for multiple testing across several stacked tab_anova() blocks.
  if (garder_p_num) out$p_num <- a[["Pr(>Chisq)"]]
  out
}

# emmeans + Tukey letters of one model. Groups sharing a letter do not differ
# significantly. `groupe_col` must be a column of as.data.frame(cld_obj), i.e.
# the term the emmeans() call was run over.
tab_cld <- function(cld_obj, an, groupe_col = "location", nom_groupe = "Station") {
  d <- as.data.frame(cld_obj)
  # the column names depend on the model family and on type = "response"
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

# adonis2 output of one year.
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

# DHARMa diagnostics of one model, as a single row.
tab_diag <- function(res, etiquette) {
  data.frame(
    Model      = etiquette,
    Uniformity = fmt_p(testUniformity(res, plot = FALSE)$p.value),
    Dispersion = fmt_p(testDispersion(res, plot = FALSE)$p.value),
    Outliers   = fmt_p(testOutliers(res, plot = FALSE)$p.value),
    stringsAsFactors = FALSE)
}

# Placeholder row, in tab_anova()'s shape, for a model that could not be fit:
# the row stays visible with a note instead of being dropped.
note_row <- function(etiquette, note) {
  data.frame(Model = etiquette, Term = "", Chi2 = "", df = "", p = "",
             `N obs` = "", `R2 marginal` = "", `R2 conditional` = "", `R2 type` = "",
             Note = note, check.names = FALSE, stringsAsFactors = FALSE)
}


#### WRITERS ####

# markdown: one section per table, pipe syntax
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

# Word: officer assembles the document, flextable formats each table
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
    # figures are delivered as separate 300 dpi files; the document lists them
    doc <- body_add_par(doc, "Figures", style = "heading 2")
    for (i in seq_len(nrow(figures))) {
      doc <- body_add_par(doc, paste0("Figure ", i, ". ", figures$legende[i],
                                      " (", figures$fichier[i], ")"),
                          style = "Normal")
    }
  }

  print(doc, target = fichier)
}
