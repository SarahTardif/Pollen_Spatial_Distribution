# Formatting helpers, table builders and writers for publication-ready output
# (markdown + Word). Moved here almost unchanged from the old model_results.r
# (its lines 21-140 and 339-403) because this is formatting, not model
# fitting (doesn't belong in models.R) and not orchestration (doesn't belong
# in a script). What table goes in which document, with which title and
# legend, is editorial content and stays in the analysis script that calls
# these functions (11_article_tableaux.R).

library(officer)
library(flextable)
library(car)          # Anova()
library(performance)  # r2()
library(DHARMa)       # testUniformity/testDispersion/testOutliers


#### FORMATTING HELPERS ####

# p-values: below 0.001 the exact value carries no information, so it is
# reported as a threshold, as is standard in the literature.
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

# numbers: 3 significant digits, never scientific notation. Counts in the
# thousands are rounded and thousand-separated instead, otherwise the pollen
# abundances come out as 3.45e+04 and become unreadable.
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

# integers (degrees of freedom, k, number of stations)
fmt_int <- function(x) {
  out <- character(length(x))
  for (i in seq_along(x)) {
    out[i] <- if (is.na(x[i])) "" else format(round(x[i]), trim = TRUE)
  }
  out
}

# R2 marginal / conditional of a glmmTMB model, as two formatted strings, plus
# the R2 type actually used. r2() errors outright on glmmTMB models with no
# random effects at all (e.g. the location + period fixed-effects models in
# 14_article_effets_fixes.R) - same McFadden pseudo-R2 fallback as valider()
# (R/models.R), so the two stay consistent. McFadden is a single value with
# no marginal/conditional split, so it is only put in the marginal slot.
r2_vals <- function(mod) {
  r <- try(r2(mod), silent = TRUE)
  if (!inherits(r, "try-error") && !is.null(r)) {
    return(c(fmt_num(as.numeric(r$R2_marginal)), fmt_num(as.numeric(r$R2_conditional)), "R2"))
  }
  # r2_mcfadden_local() (R/models.R), never performance::r2_mcfadden(): the
  # latter refits the null model by re-evaluating the model's stored call, so
  # for the per-taxon models fit in a loop it used the last taxon's data for
  # every null and produced meaningless (sometimes negative) values.
  rm <- try(r2_mcfadden_local(mod), silent = TRUE)
  if (inherits(rm, "try-error") || is.null(rm) || !is.finite(rm)) return(c("", "", ""))
  c(fmt_num(as.numeric(rm[1])), "", "McFadden")
}


#### TABLE BUILDERS ####

# Anova type II of a glmmTMB model + the R2 of that model, on one block of rows.
# The model label and the R2 are only written on the first row so that a table
# stacking several models stays readable.
tab_anova <- function(mod, etiquette, garder_p_num = FALSE) {
  a <- as.data.frame(Anova(mod, type = "II"))
  r <- r2_vals(mod)
  n <- nrow(a)
  # N obs: for count models (abundance), the true-zero pattern
  # (completer_zeros(), R/matrices.R) keeps "taxon absent from a collected
  # sample" as a valid row with count = 0, so nobs(mod) alone would overstate
  # how many samples actually had a non-zero count for this response - report
  # the non-zero count instead. For gaussian models (Shannon diversity), 0 is
  # a real diversity value (a sample with a single taxon), not a placeholder,
  # so nobs(mod) is reported as-is.
  # Response extracted via model.response(model.frame(mod)), NOT
  # insight::get_response(mod): the latter re-evaluates the model's stored
  # call against the environment it was fit in, so for models fit in a loop
  # that reuses one variable name across iterations (fit_taxon_models() in
  # 14_article_effets_fixes.R uses `dat_t` for every taxon), it silently
  # returns the LAST iteration's data for every model - confirmed with a
  # 3-taxon repro where get_response() returned identical values for all
  # three models. model.frame(mod) uses glmmTMB's own stored data snapshot,
  # not a re-evaluated call, so it stays correct per model.
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
  # raw, unrounded p-value alongside the formatted "p" column -- needed by
  # callers that adjust for multiple testing across many stacked tab_anova()
  # rows (e.g. one model per taxon) after fmt_p() has already turned "p" into
  # a display string. FALSE by default so every other caller keeps today's
  # exact column set.
  if (garder_p_num) out$p_num <- a[["Pr(>Chisq)"]]
  out
}

# emmeans + Tukey letters of one abundance model. Groups (stations by default,
# but any emmeans grouping factor, e.g. period) sharing a letter are not
# significantly different. groupe_col must be a column of as.data.frame(cld_obj),
# i.e. the term the emmeans() call was run over.
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

# adonis2 output of one year
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

# DHARMa diagnostics of one model, as a single row. plot = FALSE keeps the
# tests silent; only the p-values are reported.
tab_diag <- function(res, etiquette) {
  data.frame(
    Model      = etiquette,
    Uniformity = fmt_p(testUniformity(res, plot = FALSE)$p.value),
    Dispersion = fmt_p(testDispersion(res, plot = FALSE)$p.value),
    Outliers   = fmt_p(testOutliers(res, plot = FALSE)$p.value),
    stringsAsFactors = FALSE)
}

# Placeholder row for a tab_anova()-shaped table, in place of a model/taxon/
# period that could not be fit or whose Anova() failed - keeps the row visible
# with a Note instead of silently dropping it. Used by 11_article_tableaux.R
# (per-taxon models) and 13_article_periodes_tableaux.R (per-period models).
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

    # raw preformatted text (e.g. capture.output(summary(mod))) instead of a
    # table: some model summaries don't fit the tab_*() column shapes, so they
    # are kept verbatim in a fenced code block rather than forced into a table
    if (!is.null(t$texte)) {
      cat("```\n", paste(t$texte, collapse = "\n"), "\n```\n\n",
          sep = "", file = fichier, append = TRUE)
      next
    }

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

    # raw preformatted text: one monospace paragraph per line, so alignment
    # (e.g. glmmTMB's summary() coefficient columns) survives instead of being
    # collapsed by a table cell
    if (!is.null(t$texte)) {
      police_brute <- fp_text(font.family = "Courier New", font.size = 8)
      for (ligne in t$texte) {
        doc <- body_add_fpar(doc, fpar(ftext(if (ligne == "") " " else ligne, police_brute)))
      }
      doc <- body_add_par(doc, "", style = "Normal")
      next
    }

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
