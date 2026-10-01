# Validation des modeles: le bloc de verifications (Anova type II + R2 +
# DHARMa) lance sous chaque modele du projet. La formule du modele reste
# visible dans le script d'analyse, valider() ne fait que les verifications.

library(car)          # Anova()
library(performance)  # r2()
library(DHARMa)       # simulateResiduals()
library(glmmTMB)      # refit du modele nul

# Pseudo-R2 de McFadden = 1 - logLik(modele) / logLik(modele nul).
# Ne pas remplacer par performance::r2_mcfadden(): celle-ci refait le modele nul
# en reevaluant l'appel stocke dans son environnement d'origine, donc pour des
# modeles ajustes dans une boucle qui reutilise le meme nom de variable pour les
# donnees, tous les modeles nuls se retrouvent ajustes sur les donnees de la
# derniere iteration (R2 faux, parfois negatif). model.frame(mod) utilise la
# copie des donnees stockee par glmmTMB, donc toujours celles de ce modele-ci.
r2_mcfadden_local <- function(mod) {
  y  <- model.response(model.frame(mod))
  d0 <- data.frame(.y = as.numeric(y))
  n0 <- try(glmmTMB(.y ~ 1, data = d0, family = family(mod)), silent = TRUE)
  if (inherits(n0, "try-error")) return(NA_real_)
  ll0 <- as.numeric(logLik(n0))
  if (!is.finite(ll0) || ll0 == 0) return(NA_real_)
  1 - as.numeric(logLik(mod)) / ll0
}

# Anova (type II) + R2 + simulation des residus DHARMa pour un modele glmmTMB.
# Retourne l'objet DHARMa, dont les tableaux de diagnostics ont besoin.
# anova_try = TRUE pour les modeles dont l'Anova peut ne pas converger: l'erreur
# est alors affichee au lieu d'arreter le script.
# Repli sur McFadden quand r2() echoue: c'est le cas des modeles glmmTMB sans
# aucun effet aleatoire combines a la famille nbinom2. r2_nagelkerke() n'est pas
# utilise comme repli car il repose sur deviance(), qui pour cette famille ne
# correspond pas a -2*logLik.
valider <- function(mod, etiquette, n_sim = 1000, anova_try = FALSE, tracer = TRUE, verbose = TRUE) {
  cat("\n===", etiquette, "===\n")
  a <- if (anova_try) try(Anova(mod, type = "II")) else Anova(mod, type = "II")
  if (verbose) print(a)
  r2_res <- try(r2(mod), silent = TRUE)
  mcfadden_utilise <- inherits(r2_res, "try-error")
  if (mcfadden_utilise) r2_res <- r2_mcfadden_local(mod)
  if (verbose) {
    cat("\n---", if (mcfadden_utilise) "Pseudo-R2 (McFadden, repli)" else "R2", "\n")
    print(r2_res)
  }
  res <- simulateResiduals(mod, n = n_sim)
  if (tracer) plot(res)
  res
}
