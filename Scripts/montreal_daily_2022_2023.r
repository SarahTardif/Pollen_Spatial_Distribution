# ============================================================================
# Script autonome: Données climatiques journalières de Montréal (2022-2023)
# ============================================================================
# Télécharge les données de Montréal et crée 2 graphiques (2022 et 2023)
# avec T, vent, orientation, pluie
#
# Ce script est totalement autonome et peut être copié dans un autre repo
# ============================================================================

# ── 1. Installation et chargement des packages ────────────────────────────

# Décommenter pour installer si nécessaire:
# install.packages("weathercan",
#                  repos = c("https://ropensci.r-universe.dev",
#                            "https://cloud.r-project.org"))
# install.packages("dplyr")

library(weathercan)
library(dplyr)

# ── 2. Créer dossier de sortie ──────────────────────────────────────────

if (!dir.exists("./Outputs/climatedata")) {
  dir.create("./Outputs/climatedata", showWarnings = FALSE)
}

# ── 3. Fonction haversine pour calculer distance entre deux points ─────────

haversine_km <- function(lat1, lon1, lat2, lon2) {
  R    <- 6371
  dlat <- (lat2 - lat1) * pi / 180
  dlon <- (lon2 - lon1) * pi / 180
  a    <- sin(dlat/2)^2 + cos(lat1*pi/180) * cos(lat2*pi/180) * sin(dlon/2)^2
  2 * R * asin(sqrt(a))
}

# ── 4. Télécharger toutes les stations et chercher la meilleure pour Montréal

cat("Téléchargement de la liste des stations...\n")
all_stations <- stations()

# Coordonnées de Montréal
montreal_lat <- 45.50
montreal_lon <- -73.57

# Filtrer stations journalières disponibles - PRIORISER LES AÉROPORTS
stations_candidates <- all_stations %>%
  filter(
    interval == "day",
    !is.na(lat), !is.na(lon),
    end >= 2022  # Doit être actif en 2022 ou après
  ) %>%
  mutate(
    distance_km = haversine_km(lat, lon, montreal_lat, montreal_lon),
    # Prioriser les aéroports (contiennent "INTL", "AIRPORT", ou noms spécifiques)
    is_airport = grepl("INTL|AIRPORT|TRUDEAU", station_name, ignore.case = TRUE),
    priority = if_else(is_airport, 1, 2)
  ) %>%
  filter(distance_km <= 100)

# Afficher stations aéroportuaires disponibles
airport_stations <- stations_candidates %>%
  filter(is_airport) %>%
  arrange(distance_km)

cat("=== Stations aéroportuaires disponibles ===\n")
print(airport_stations %>%
  select(station_name, station_id, start, end, lat, lon, distance_km))

# Sélectionner les meilleures stations aéroportuaires
selected_stations <- airport_stations %>%
  slice_head(n = min(2, nrow(airport_stations)))

if (nrow(selected_stations) == 0) {
  # Fallback aux stations les plus proches
  selected_stations <- stations_candidates %>%
    arrange(distance_km) %>%
    slice_head(n = 2)
  cat("\nPas de stations aéroportuaires. Utilisation des 2 stations les plus proches.\n")
}

cat("\n=== Stations sélectionnées pour fusion ===\n")
print(selected_stations %>%
  select(station_name, station_id, lat, lon, distance_km))

station_ids <- selected_stations$station_id

# ── 5. Télécharger les données journalières pour 2022-2023 ───────────────

cat("\nTéléchargement des données journalières de", length(station_ids), "station(s)...\n")

# Télécharger et combiner les données de toutes les stations
all_data <- list()
for (i in seq_along(station_ids)) {
  cat("Téléchargement station", i, "...\n")
  data_station <- weather_dl(
    station_ids = station_ids[i],
    start       = "2022-01-01",
    end         = "2023-12-31",
    interval    = "day",
    quiet       = TRUE
  )
  all_data[[i]] <- data_station
}

# Combiner les données
donnees <- bind_rows(all_data)

cat("Données téléchargées :", nrow(donnees), "lignes\n")
cat("Nombre de stations avec données:", n_distinct(donnees$station_id), "\n")

# ── 6. Préparer les données ─────────────────────────────────────────────

# Convertir la date en format Date et garder les colonnes essentielles
donnees <- donnees %>%
  mutate(
    date = as.Date(date),
    year = as.integer(year),
    month = as.integer(month),
    day_of_year = as.integer(day)
  ) %>%
  select(
    date, year, month, day_of_year, station_name, station_id, lat, lon,
    min_temp, max_temp, mean_temp,
    spd_max_gust, dir_max_gust,
    total_rain
  ) %>%
  arrange(date)

# Combiner intelligemment les données de plusieurs stations:
# Pour chaque date, prendre la première valeur non-NA pour chaque variable
donnees_combined <- donnees %>%
  group_by(date) %>%
  summarise(
    year = first(year),
    month = first(month),
    day_of_year = first(day_of_year),
    # Pour chaque variable, prendre la première valeur non-NA
    min_temp = first(na.omit(min_temp)),
    max_temp = first(na.omit(max_temp)),
    mean_temp = first(na.omit(mean_temp)),
    spd_max_gust = first(na.omit(spd_max_gust)),
    dir_max_gust = first(na.omit(dir_max_gust)),
    total_rain = first(na.omit(total_rain)),
    # Tracer quelle station a fourni les données principales
    stations_used = paste(unique(na.omit(station_name)), collapse = " + "),
    .groups = "drop"
  ) %>%
  arrange(date)

# Remplacer les valeurs 0 potentiellement fausses par NA pour total_rain
donnees_combined <- donnees_combined %>%
  mutate(
    total_rain = case_when(
      is.na(total_rain) ~ NA_real_,
      total_rain == 0 ~ NA_real_,
      TRUE ~ total_rain
    )
  )

# Utiliser les données combinées
donnees <- donnees_combined

# Vérifier les données
cat("\nRésumé des données:\n")
print(summary(donnees))

# Sauvegarder les données brutes
write.csv(donnees, "./Outputs/climatedata/montreal_daily_2022_2023.csv", row.names = FALSE)
cat("\nDonnées sauvegardées: ./Outputs/climatedata/montreal_daily_2022_2023.csv\n")

# ── 7. Définir les semaines de phénologie pour chaque année ─────────────

weeks_2022 <- data.frame(
  week = c("W1", "W2", "W3", "W4", "W5", "W6", "W7"),
  start = as.Date(c("2022-04-05", "2022-04-20", "2022-05-03", "2022-05-17",
                     "2022-05-31", "2022-06-14", "2022-06-28")),
  end = as.Date(c("2022-04-20", "2022-05-03", "2022-05-17", "2022-05-31",
                   "2022-06-14", "2022-06-28", "2022-07-12"))
)

weeks_2023 <- data.frame(
  week = c("W1", "W2", "W3", "W4", "W5", "W6", "W7"),
  start = as.Date(c("2023-03-28", "2023-04-11", "2023-04-25", "2023-05-09",
                     "2023-05-23", "2023-06-06", "2023-06-20")),
  end = as.Date(c("2023-04-11", "2023-04-25", "2023-05-09", "2023-05-23",
                   "2023-06-06", "2023-06-20", "2023-07-04"))
)

# Convertit un nombre de lignes de marge (comme line= dans mtext()) en
# coordonnée Y "user", pour pouvoir positionner segments()/text() à une
# distance constante de l'axe (indépendante de ylim), tout en gardant le
# contrôle sur la taille du texte et le dessin des accolades (mtext() seul
# ne permet pas de tracer les traits des accolades).
line_to_user_y <- function(line) {
  ligne_pouces <- par("cin")[2] * par("cex") * par("lheight")
  decalage <- diff(grconvertY(c(0, ligne_pouces), from = "inches", to = "user"))
  par("usr")[3] - line * decalage
}

# Fonction pour dessiner les accolades des semaines
draw_week_brackets <- function(weeks_data) {
  usr <- par()$usr
  y_max <- usr[4]

  # L'accolade (line 2.6) et son label (line 3.7) restent proches l'un de
  # l'autre, à ~1 ligne d'écart, tout en étant sous les labels de mois de R
  # (~ligne 0.7-1.7) et sous "Date" (ligne 2.5 sur le panel du bas).
  y_bracket <- line_to_user_y(2.6)
  y_text <- line_to_user_y(3.7)
  patte <- 0.02 * (y_max - usr[3])

  for (i in 1:nrow(weeks_data)) {
    x_start <- as.numeric(weeks_data$start[i])
    x_end <- as.numeric(weeks_data$end[i])
    x_mid <- (x_start + x_end) / 2

    # Dessiner une accolade simple (trait horizontal avec trait vertical)
    segments(x_start, y_bracket, x_end, y_bracket, col = "grey30", lwd = 1.5)
    segments(x_start, y_bracket, x_start, y_bracket + patte, col = "grey30", lwd = 1.5)
    segments(x_end, y_bracket, x_end, y_bracket + patte, col = "grey30", lwd = 1.5)

    # Label W1, W2, etc., plus petit que les labels d'axe (cex.axis = 0.85)
    # pour rester discret sous l'accolade
    text(x_mid, y_text, weeks_data$week[i], cex = 0.7, col = "black", font = 1)
  }
}

# ── 8. Créer une figure combinée: 2022 à gauche, 2023 à droite ──────────

# Préparer les données pour une année (début avril à mi-juillet)
preparer_annee <- function(data_year, year_value) {
  data_year %>%
    filter(year == year_value) %>%
    mutate(month_day = paste0(sprintf("%02d", month), "-", sprintf("%02d", as.integer(format(date, "%d"))))) %>%
    filter(month_day >= "03-28" & month_day <= "07-15") %>%  # Début avril à début août
    arrange(date)
}

d_2022 <- preparer_annee(donnees, 2022)
d_2023 <- preparer_annee(donnees, 2023)

# Limites Y calculées sur les deux années combinées (pas année par année
# comme dans une version précédente) pour que les deux colonnes de la
# figure soient directement comparables à l'œil. Mêmes vérifications
# Inf/-Inf qu'avant.
temp_vals <- c(d_2022$min_temp, d_2022$max_temp, d_2023$min_temp, d_2023$max_temp)
if (sum(!is.na(temp_vals)) > 0) {
  ylim_temp <- range(temp_vals, na.rm = TRUE) + c(-2, 2)
} else {
  ylim_temp <- c(-30, 30)
}

# ylim doit partir de 0, pas de min(valeurs): type="h" dessine chaque
# barre depuis y=0 jusqu'à la valeur, donc si l'axe démarre au-dessus de
# 0 (ex. la vitesse de vent minimale observée), la base des barres tombe
# sous le bas du panel — avec xpd=NA elles se dessinent quand même, mais
# débordent hors du panel au lieu d'y rester.
wind_vals <- c(d_2022$spd_max_gust, d_2023$spd_max_gust)
if (sum(!is.na(wind_vals)) > 0) {
  ylim_wind_speed <- c(0, max(wind_vals, na.rm = TRUE) + 5)
} else {
  ylim_wind_speed <- c(0, 50)
}
ylim_wind_dir <- c(0, 360)

rain_vals <- c(d_2022$total_rain, d_2023$total_rain)
if (sum(!is.na(rain_vals)) > 0) {
  ylim_rain <- c(0, max(rain_vals, na.rm = TRUE) + 5)
} else {
  ylim_rain <- c(0, 50)
}

filename_combine <- "./Outputs/climatedata/montreal_daily_2022_2023.png"

# Dessine les 4 panels d'une année dans la colonne courante du layout
# (appelée une fois par année; par(mfcol=...) avance de colonne après 4
# plots). Les ylim_* sont pris dans l'environnement englobant (calculés
# une fois pour les deux années ci-dessus, communs aux deux colonnes).
dessiner_colonne_annee <- function(d, year_value, weeks_current) {

  if (nrow(d) == 0) {
    cat("Pas de données pour l'année", year_value, "\n")
    return()
  }

  # ── Panel 1: Température (Tmin, Tmoy, Tmax) ────────────────────────────
  par(mar = c(5, 4, 2, 1))
  plot(d$date, d$mean_temp,
       type = "l", col = "grey50", lwd = 1.5,
       ylim = ylim_temp,
       xlab = "", ylab = "Température (°C)",
       main = paste("Montréal - daily climate variations -", year_value),
       cex.main = 1.2, cex.lab = 0.9, cex.axis = 0.85)

  # Ajouter Tmin et Tmax
  lines(d$date, d$min_temp, col = "steelblue", lwd = 1, lty = 2)
  lines(d$date, d$max_temp, col = "firebrick", lwd = 1, lty = 2)

  # Légende pour T
  legend("topleft",
         legend = c("Tmin", "Tmoy", "Tmax"),
         col = c("steelblue", "grey50", "firebrick"),
         lty = c(2, 1, 2), lwd = c(1, 1.5, 1),
         cex = 0.8, bty = "n")

  # Ajouter les accolades des semaines
  draw_week_brackets(weeks_current)

  # ── Panel 2: Vitesse du vent ────────────────────────────────────────────
  par(mar = c(5, 4, 1, 1))
  plot(d$date, d$spd_max_gust,
       type = "h", col = "navy", lwd = 1.5,
       ylim = ylim_wind_speed,
       xlab = "", ylab = "Vent (km/h)",
       cex.lab = 0.9, cex.axis = 0.85)

  # Ajouter les accolades des semaines
  draw_week_brackets(weeks_current)

  # ── Panel 3: Orientation du vent ────────────────────────────────────────
  # yaxt="n" + axis() manuel: les ticks par défaut de R pour ylim=c(0,360)
  # tombent sur 0/50/100/.../350, qui ne correspondent à aucune des 4
  # directions cardinales — les lignes de référence pointillées à
  # 0/90/180/270/360 se retrouvaient donc sans tick ni label en face
  # d'elles. Les ticks sont fixés directement sur les cardinaux et étiquetés
  # en points cardinaux (plus lisible que les degrés pour ce panel).
  par(mar = c(5, 4, 1, 1))
  plot(d$date, d$dir_max_gust,
       type = "p", col = "grey40", pch = 16, cex = 0.4,
       ylim = ylim_wind_dir, yaxt = "n",
       xlab = "", ylab = "Orientation",
       cex.lab = 0.9, cex.axis = 0.85)
  axis(2, at = c(0, 90, 180, 270, 360), labels = c("N", "E", "S", "O", "N"),
       cex.axis = 0.85)

  # Ajouter des lignes de référence pour les directions cardinales
  abline(h = c(0, 90, 180, 270, 360), col = "lightgrey", lty = 3, lwd = 0.5)

  # Ajouter les accolades des semaines
  draw_week_brackets(weeks_current)

  # ── Panel 4: Pluie totale ───────────────────────────────────────────────
  par(mar = c(5, 4, 1, 1))
  plot(d$date, d$total_rain,
       type = "h", col = "steelblue", lwd = 1.5,
       ylim = ylim_rain,
       xlab = "", ylab = "Pluie (mm)",
       cex.lab = 0.9, cex.axis = 0.85)

  # Ajouter les accolades des semaines
  draw_week_brackets(weeks_current)
}

# Créer le fichier PNG combiné. Canevas élargi pour les 2 colonnes et
# pointsize augmenté (10 -> 13) par rapport à l'ancienne figure par année:
# tous les cex.* ci-dessus sont relatifs à pointsize, donc ce seul
# changement fait grossir tout le texte (titres, axes, légende, labels de
# semaines) proportionnellement, sans retoucher chaque cex individuellement.
png(filename_combine, width = 18, height = 18, units = "cm", res = 600, pointsize = 13)

# par(mfcol=...) remplit par colonne (les 4 panels d'une année d'abord,
# puis la colonne suivante) plutôt que par ligne (par(mfrow=...) donnerait
# 2022/2023 alternés ligne par ligne) — c'est ce qui met 2022 entièrement
# à gauche et 2023 entièrement à droite. mar reste fixé panel par panel
# dans dessiner_colonne_annee(), comme dans l'ancienne version par année.
par(mfcol = c(4, 2),
    mgp = c(2.5, 0.7, 0),
    xpd = NA)

dessiner_colonne_annee(d_2022, 2022, weeks_2022)
dessiner_colonne_annee(d_2023, 2023, weeks_2023)

dev.off()

cat("Graphique combiné sauvegardé:", filename_combine, "\n")

# ── Résumé final ────────────────────────────────────────────────────────

cat("\n========== COMPLET ==========\n")
cat("Stations utilisées:\n")
for (i in seq_len(nrow(selected_stations))) {
  s <- selected_stations[i, ]
  cat("  -", s$station_name, "\n")
  cat("    ID:", s$station_id, "| Lat:", s$lat, "| Lon:", s$lon, "\n")
  cat("    Distance Montréal centre:", round(s$distance_km, 2), "km\n")
}
cat("\nPériode: 2022-01-01 à 2023-12-31\n")
cat("\nFichiers créés:\n")
cat("  - ./Outputs/climatedata/montreal_daily_2022_2023.csv (données brutes)\n")
cat("  - ./Outputs/climatedata/montreal_daily_2022_2023.png (graphique combiné 2022/2023)\n")
cat("=============================\n")
