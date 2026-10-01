# Fitting-side helpers: the standard validation bundle run under every model
# in this codebase (Anova type II + r2 + DHARMa), and the two ways models get
# turned into a results table (one row per term of one model: extraire_coefs;
# one term across many models: extraire_terme).
#
# The model formula itself is NOT hidden in here — glmmTMB() stays a visible
# line in the analysis script, valider() only runs the checks that would
# otherwise be copy-pasted under each model:
#   mod_abglob_2022 <- glmmTMB(count ~ location + (1 | period), data = ab_tot_2022, family = nbinom2)
#   res_abglob_2022 <- valider(mod_abglob_2022, "Abondance 2022")

library(car)          # Anova()
library(performance)  # r2()
library(DHARMa)       # simulateResiduals()
library(glmmTMB)      # r2_mcfadden_local() refits the null model itself

# McFadden pseudo-R2 = 1 - logLik(model) / logLik(intercept-only model), with
# the null model refit HERE rather than through performance::r2_mcfadden().
#
# r2_mcfadden() delegates the null refit to insight::null_model(), which
# re-evaluates the model's stored call in the environment it was fit in. For
# models fit in a loop that reuses one variable name for the data
# (fit_taxon_models() in 10_article.R / 14_article_effets_fixes.R uses `dat_t`
# for every taxon), that environment only holds the LAST iteration's data by
# the time the R2 is computed, so every taxon's null model was silently fit on
# the last taxon's counts. The R2 was then a ratio of two log-likelihoods
# computed on different datasets, and went negative whenever the last taxon's
# null fit happened to have a higher logLik than the current taxon's full fit
# (a rare taxon compared against an abundant one — that's where Table 3's
# negative R2 came from). Confirmed on a 2-taxon repro: taxon A got 0.571
# through r2_mcfadden() vs 0.002 against its own null, and only the last taxon
# of each loop was correct. Same failure mode as the one already documented
# under tab_anova() (R/report.R) for insight::get_response().
#
# model.frame(mod) is glmmTMB's own stored data snapshot, so the response
# below is always this model's, never a re-evaluated call's. Returns NA if the
# null model cannot be fit.
r2_mcfadden_local <- function(mod) {
  y  <- model.response(model.frame(mod))
  d0 <- data.frame(.y = as.numeric(y))
  n0 <- try(glmmTMB(.y ~ 1, data = d0, family = family(mod)), silent = TRUE)
  if (inherits(n0, "try-error")) return(NA_real_)
  ll0 <- as.numeric(logLik(n0))
  if (!is.finite(ll0) || ll0 == 0) return(NA_real_)
  1 - as.numeric(logLik(mod)) / ll0
}

# Anova (type II) + r2 + DHARMa residual simulation/plot for one glmmTMB
# model. Returns the DHARMa simulation object (as the diagnostics table
# builders in report.R need it). anova_try = TRUE for models whose Anova can
# fail to converge (e.g. the taxon:station interaction models) — the error is
# then printed instead of stopping the script, matching the original
# print(try(Anova(...))) calls.
valider <- function(mod, etiquette, n_sim = 1000, anova_try = FALSE, tracer = TRUE, verbose = TRUE) {
  cat("\n===", etiquette, "===\n")
  a <- if (anova_try) try(Anova(mod, type = "II")) else Anova(mod, type = "II")
  if (verbose) print(a)
  # r2() errors outright on glmmTMB models with no random effects at all (a
  # fully-fixed-effects specification, which is a valid model here) combined
  # with the nbinom2 family — no method falls back to a GLM-style pseudo R2.
  # r2_nagelkerke() runs without erroring in that case but is not trustworthy
  # for glmmTMB + nbinom2: it's built on deviance(), and deviance.glmmTMB()
  # (sum(residuals(type = "deviance")^2)) does not match -2*logLik for this
  # family — checked directly on mod_ab_2022: deviance() = 186 vs -2*logLik =
  # 3472, off by more than an order of magnitude, so the resulting Cox-Snell/
  # Nagelkerke R2 is built on a broken number. McFadden is logLik-ratio based
  # instead (not deviance()) and matches a manual logLik comparison against the
  # refit null model exactly, so it's the fallback here rather than
  # r2_nagelkerke(). It goes through r2_mcfadden_local() above, NOT
  # performance::r2_mcfadden() — see that function's comment for why.
  r2_res <- try(r2(mod), silent = TRUE)
  mcfadden_utilise <- inherits(r2_res, "try-error")
  if (mcfadden_utilise) r2_res <- r2_mcfadden_local(mod)
  if (verbose) {
    cat("\n---", if (mcfadden_utilise) "Pseudo-R2 (McFadden, fallback)" else "R2", "\n")
    print(r2_res)
  }
  res <- simulateResiduals(mod, n = n_sim)
  if (tracer) plot(res)
  res
}

# All fixed-effect coefficients of one model, as a data.frame row per term —
# used to build a summary table across several models/radii (e.g.
# resume_modeles in 20_pollen_vs_arbres.R).
extraire_coefs <- function(mod, nom_modele, rayon, reponse) {
  coefs <- summary(mod)$coefficients$cond
  data.frame(
    modele    = nom_modele,
    reponse   = reponse,
    rayon_m   = rayon,
    terme     = rownames(coefs),
    estime    = round(coefs[, "Estimate"], 4),
    err_type  = round(coefs[, "Std. Error"], 4),
    z         = round(coefs[, "z value"], 3),
    p         = round(coefs[, ncol(coefs)], 4),
    AIC       = round(AIC(mod), 1),
    n_obs     = nobs(mod),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

# One named term of one model, as a named vector — used when looping over
# many separate models (one per genus) and collecting a single effect each
# time. Returns NA everywhere if the model failed to fit or didn't converge
# (non-finite AIC), so the caller's loop doesn't need its own try/catch.
extraire_terme <- function(mod, terme) {
  if (inherits(mod, "try-error") || !is.finite(AIC(mod))) {
    return(c(estime = NA, err_type = NA, z = NA, p = NA, AIC = NA))
  }
  coefs <- summary(mod)$coefficients$cond
  c(estime   = round(coefs[terme, "Estimate"], 4),
    err_type = round(coefs[terme, "Std. Error"], 4),
    z        = round(coefs[terme, "z value"], 3),
    p        = round(coefs[terme, ncol(coefs)], 4),
    AIC      = round(AIC(mod), 1))
}

# Significance stars from a p-value (or vector of p-values). NA -> "".
etoiles <- function(p) {
  ifelse(is.na(p), "",
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "**",
  ifelse(p < 0.05,  "*",
  ifelse(p < 0.1,   ".", "")))))
}
