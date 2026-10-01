# Model validation: the block of checks (type II Anova + R2 + DHARMa) run under
# every model in the project. The model formula stays visible in the analysis
# script, valider() only runs the checks.

library(car)          # Anova()
library(performance)  # r2()
library(DHARMa)       # simulateResiduals()
library(glmmTMB)      # null model refit

# McFadden pseudo-R2 = 1 - logLik(model) / logLik(null model).
# Do not replace with performance::r2_mcfadden(): it refits the null model by
# re-evaluating the stored call in its original environment, so for models fit
# in a loop that reuses one variable name for the data, every null model ends
# up fit on the last iteration's data (wrong R2, sometimes negative).
# model.frame(mod) uses the data snapshot stored by glmmTMB, hence always this
# model's own data.
r2_mcfadden_local <- function(mod) {
  y  <- model.response(model.frame(mod))
  d0 <- data.frame(.y = as.numeric(y))
  n0 <- try(glmmTMB(.y ~ 1, data = d0, family = family(mod)), silent = TRUE)
  if (inherits(n0, "try-error")) return(NA_real_)
  ll0 <- as.numeric(logLik(n0))
  if (!is.finite(ll0) || ll0 == 0) return(NA_real_)
  1 - as.numeric(logLik(mod)) / ll0
}

# Type II Anova + R2 + DHARMa residual simulation for one glmmTMB model.
# Returns the DHARMa object, which the diagnostics tables need.
# anova_try = TRUE for models whose Anova may fail to converge: the error is
# then printed instead of stopping the script.
# Falls back to McFadden when r2() fails, which happens for glmmTMB models with
# no random effect at all combined with the nbinom2 family. r2_nagelkerke() is
# not used as the fallback because it builds on deviance(), which for this
# family does not match -2*logLik.
valider <- function(mod, etiquette, n_sim = 1000, anova_try = FALSE, tracer = TRUE, verbose = TRUE) {
  cat("\n===", etiquette, "===\n")
  a <- if (anova_try) try(Anova(mod, type = "II")) else Anova(mod, type = "II")
  if (verbose) print(a)
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
