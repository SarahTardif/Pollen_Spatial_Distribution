# Model validation

library(car)          # Anova()
library(performance)  # r2()
library(DHARMa)       # simulateResiduals()
library(glmmTMB)      # null model refit

# McFadden pseudo-R2 = 1 - logLik(model) / logLik(null model).

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
