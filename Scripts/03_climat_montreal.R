# Script autonome: donnees climatiques journalieres de Montreal (2022-2023).
# Telecharge les donnees de la station meteo la plus proche et produit une
# figure a 2 colonnes (2022 et 2023) avec temperature, vent, orientation du
# vent et pluie. Ne depend pas de R/init.R: copiable tel quel dans un autre
# projet.

library(weathercan)
library(dplyr)

if (!dir.exists("./Outputs/climatedata")) {
  dir.create("./Outputs/climatedata", recursive = TRUE, showWarnings = FALSE)
}

# Distance haversine entre deux points, en km.
haversine_km <- function(lat1, lon1, lat2, lon2) {
  R    <- 6371
  dlat <- (lat2 - lat1) * pi / 180
  dlon <- (lon2 - lon1) * pi / 180
  a    <- sin(dlat/2)^2 + cos(lat1*pi/180) * cos(lat2*pi/180) * sin(dlon/2)^2
  2 * R * asin(sqrt(a))
}

## Chercher la meilleure station pour Montreal, en priorisant les aeroports

cat("Telechargement de la liste des stations...\n")
all_stations <- stations()

montreal_lat <- 45.50
montreal_lon <- -73.57

stations_candidates <- all_stations %>%
  filter(
    interval == "day",
    !is.na(lat), !is.na(lon),
    end >= 2022
  ) %>%
  mutate(
    distance_km = haversine_km(lat, lon, montreal_lat, montreal_lon),
    is_airport = grepl("INTL|AIRPORT|TRUDEAU", station_name, ignore.case = TRUE),
    priority = if_else(is_airport, 1, 2)
  ) %>%
  filter(distance_km <= 100)

airport_stations <- stations_candidates %>%
  filter(is_airport) %>%
  arrange(distance_km)

cat("=== Stations aeroportuaires disponibles ===\n")
print(airport_stations %>%
  select(station_name, station_id, start, end, lat, lon, distance_km))

selected_stations <- airport_stations %>%
  slice_head(n = min(2, nrow(airport_stations)))

if (nrow(selected_stations) == 0) {
  # repli sur les 2 stations les plus proches
  selected_stations <- stations_candidates %>%
    arrange(distance_km) %>%
    slice_head(n = 2)
  cat("\nPas de stations aeroportuaires. Utilisation des 2 stations les plus proches.\n")
}

cat("\n=== Stations selectionnees pour fusion ===\n")
print(selected_stations %>%
  select(station_name, station_id, lat, lon, distance_km))

station_ids <- selected_stations$station_id

## Telecharger les donnees journalieres 2022-2023

cat("\nTelechargement des donnees journalieres de", length(station_ids), "station(s)...\n")

all_data <- list()
for (i in seq_along(station_ids)) {
  cat("Telechargement station", i, "...\n")
  data_station <- weather_dl(
    station_ids = station_ids[i],
    start       = "2022-01-01",
    end         = "2023-12-31",
    interval    = "day",
    quiet       = TRUE
  )
  all_data[[i]] <- data_station
}

donnees <- bind_rows(all_data)

cat("Donnees telechargees :", nrow(donnees), "lignes\n")
cat("Nombre de stations avec donnees:", n_distinct(donnees$station_id), "\n")

## Preparer les donnees

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

# Fusion des stations: pour chaque date, la 1re valeur non-NA de chaque variable
donnees_combined <- donnees %>%
  group_by(date) %>%
  summarise(
    year = first(year),
    month = first(month),
    day_of_year = first(day_of_year),
    min_temp = first(na.omit(min_temp)),
    max_temp = first(na.omit(max_temp)),
    mean_temp = first(na.omit(mean_temp)),
    spd_max_gust = first(na.omit(spd_max_gust)),
    dir_max_gust = first(na.omit(dir_max_gust)),
    total_rain = first(na.omit(total_rain)),
    stations_used = paste(unique(na.omit(station_name)), collapse = " + "),
    .groups = "drop"
  ) %>%
  arrange(date)

# les 0 de pluie sont potentiellement des valeurs manquantes
donnees_combined <- donnees_combined %>%
  mutate(
    total_rain = case_when(
      is.na(total_rain) ~ NA_real_,
      total_rain == 0 ~ NA_real_,
      TRUE ~ total_rain
    )
  )

donnees <- donnees_combined

cat("\nResume des donnees:\n")
print(summary(donnees))

write.csv(donnees, "./Outputs/climatedata/montreal_daily_2022_2023.csv", row.names = FALSE)
cat("\nDonnees sauvegardees: ./Outputs/climatedata/montreal_daily_2022_2023.csv\n")

## Semaines de phenologie de chaque annee

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
# coordonnee Y "user", pour positionner les accolades a distance constante de
# l'axe quel que soit ylim. mtext() seul ne permet pas de tracer les traits.
line_to_user_y <- function(line) {
  ligne_pouces <- par("cin")[2] * par("cex") * par("lheight")
  decalage <- diff(grconvertY(c(0, ligne_pouces), from = "inches", to = "user"))
  par("usr")[3] - line * decalage
}

# Accolades des semaines sous l'axe des x.
draw_week_brackets <- function(weeks_data) {
  usr <- par()$usr
  y_max <- usr[4]

  y_bracket <- line_to_user_y(2.6)
  y_text <- line_to_user_y(3.7)
  patte <- 0.02 * (y_max - usr[3])

  for (i in 1:nrow(weeks_data)) {
    x_start <- as.numeric(weeks_data$start[i])
    x_end <- as.numeric(weeks_data$end[i])
    x_mid <- (x_start + x_end) / 2

    segments(x_start, y_bracket, x_end, y_bracket, col = "grey30", lwd = 1.5)
    segments(x_start, y_bracket, x_start, y_bracket + patte, col = "grey30", lwd = 1.5)
    segments(x_end, y_bracket, x_end, y_bracket + patte, col = "grey30", lwd = 1.5)

    text(x_mid, y_text, weeks_data$week[i], cex = 0.7, col = "black", font = 1)
  }
}

## Figure combinee: 2022 a gauche, 2023 a droite

# Donnees d'une annee, de fin mars a mi-juillet.
preparer_annee <- function(data_year, year_value) {
  data_year %>%
    filter(year == year_value) %>%
    mutate(month_day = paste0(sprintf("%02d", month), "-", sprintf("%02d", as.integer(format(date, "%d"))))) %>%
    filter(month_day >= "03-28" & month_day <= "07-15") %>%
    arrange(date)
}

d_2022 <- preparer_annee(donnees, 2022)
d_2023 <- preparer_annee(donnees, 2023)

# Limites Y calculees sur les deux annees combinees pour que les deux colonnes
# de la figure soient directement comparables a l'oeil.
temp_vals <- c(d_2022$min_temp, d_2022$max_temp, d_2023$min_temp, d_2023$max_temp)
if (sum(!is.na(temp_vals)) > 0) {
  ylim_temp <- range(temp_vals, na.rm = TRUE) + c(-2, 2)
} else {
  ylim_temp <- c(-30, 30)
}

# ylim doit partir de 0: type="h" trace chaque barre depuis y=0, donc si l'axe
# demarre au-dessus de 0 la base des barres sort du panneau.
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

# Trace les 4 panneaux d'une annee dans la colonne courante du layout. Les
# ylim_* sont pris dans l'environnement englobant (calcules une fois pour les
# deux annees ci-dessus).
dessiner_colonne_annee <- function(d, year_value, weeks_current) {

  if (nrow(d) == 0) {
    cat("Pas de donnees pour l'annee", year_value, "\n")
    return()
  }

  # Panneau 1: temperature (Tmin, Tmoy, Tmax)
  par(mar = c(5, 4, 2, 1))
  plot(d$date, d$mean_temp,
       type = "l", col = "grey50", lwd = 1.5,
       ylim = ylim_temp,
       xlab = "", ylab = "Température (°C)",
       main = paste("Montréal - daily climate variations -", year_value),
       cex.main = 1.2, cex.lab = 0.9, cex.axis = 0.85)

  lines(d$date, d$min_temp, col = "steelblue", lwd = 1, lty = 2)
  lines(d$date, d$max_temp, col = "firebrick", lwd = 1, lty = 2)

  legend("topleft",
         legend = c("Tmin", "Tmoy", "Tmax"),
         col = c("steelblue", "grey50", "firebrick"),
         lty = c(2, 1, 2), lwd = c(1, 1.5, 1),
         cex = 0.8, bty = "n")

  draw_week_brackets(weeks_current)

  # Panneau 2: vitesse du vent
  par(mar = c(5, 4, 1, 1))
  plot(d$date, d$spd_max_gust,
       type = "h", col = "navy", lwd = 1.5,
       ylim = ylim_wind_speed,
       xlab = "", ylab = "Vent (km/h)",
       cex.lab = 0.9, cex.axis = 0.85)

  draw_week_brackets(weeks_current)

  # Panneau 3: orientation du vent. yaxt="n" + axis() manuel: les ticks par
  # defaut pour ylim=c(0,360) tombent sur 0/50/100/... et ne correspondent a
  # aucune direction cardinale.
  par(mar = c(5, 4, 1, 1))
  plot(d$date, d$dir_max_gust,
       type = "p", col = "grey40", pch = 16, cex = 0.4,
       ylim = ylim_wind_dir, yaxt = "n",
       xlab = "", ylab = "Orientation",
       cex.lab = 0.9, cex.axis = 0.85)
  axis(2, at = c(0, 90, 180, 270, 360), labels = c("N", "E", "S", "O", "N"),
       cex.axis = 0.85)

  abline(h = c(0, 90, 180, 270, 360), col = "lightgrey", lty = 3, lwd = 0.5)

  draw_week_brackets(weeks_current)

  # Panneau 4: pluie totale
  par(mar = c(5, 4, 1, 1))
  plot(d$date, d$total_rain,
       type = "h", col = "steelblue", lwd = 1.5,
       ylim = ylim_rain,
       xlab = "", ylab = "Pluie (mm)",
       cex.lab = 0.9, cex.axis = 0.85)

  draw_week_brackets(weeks_current)
}

# pointsize = 13: tous les cex.* ci-dessus sont relatifs a pointsize, donc ce
# seul reglage fait grossir tout le texte proportionnellement.
png(filename_combine, width = 18, height = 18, units = "cm", res = 600, pointsize = 13)

# par(mfcol=...) remplit par colonne (les 4 panneaux d'une annee d'abord), ce
# qui met 2022 entierement a gauche et 2023 entierement a droite.
par(mfcol = c(4, 2),
    mgp = c(2.5, 0.7, 0),
    xpd = NA)

dessiner_colonne_annee(d_2022, 2022, weeks_2022)
dessiner_colonne_annee(d_2023, 2023, weeks_2023)

dev.off()

cat("Graphique combine sauvegarde:", filename_combine, "\n")

cat("\n========== COMPLET ==========\n")
cat("Stations utilisees:\n")
for (i in seq_len(nrow(selected_stations))) {
  s <- selected_stations[i, ]
  cat("  -", s$station_name, "\n")
  cat("    ID:", s$station_id, "| Lat:", s$lat, "| Lon:", s$lon, "\n")
  cat("    Distance Montreal centre:", round(s$distance_km, 2), "km\n")
}
cat("\nPeriode: 2022-01-01 a 2023-12-31\n")
cat("\nFichiers crees:\n")
cat("  - ./Outputs/climatedata/montreal_daily_2022_2023.csv (donnees brutes)\n")
cat("  - ./Outputs/climatedata/montreal_daily_2022_2023.png (graphique combine)\n")
cat("=============================\n")
