# Point d'entree unique des scripts d'analyse: source("R/init.R") charge tous
# les modules dans l'ordre de leurs dependances.

source("R/config.R")     # constantes
source("R/geo.R")        # coordonnees, projections, distances
source("R/matrices.R")   # vrais zeros et matrices station x taxon
source("R/diversity.R")  # indice de Shannon
source("R/pollen.R")     # chargement du pollen et grille de comptages
source("R/models.R")     # validation des modeles glmmTMB
source("R/plots.R")      # palettes et figures partagees
source("R/report.R")     # tableaux et ecriture markdown/Word
