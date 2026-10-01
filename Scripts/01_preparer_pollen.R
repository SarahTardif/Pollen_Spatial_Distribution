# Ingestion step, run ONCE (not sourced by the other scripts): reads the 3
# yearly parquet files, derives period/year/location/comments from h5_key,
# filters at Confidence >= CONF_MIN and writes the CSVs the rest of the project
# reads through charger_pollen(). Re-run after any change to the cleaning logic
# below, otherwise the downstream scripts keep reading the stale CSVs.

source("R/init.R")
library(arrow)

data2021_brut <- read_parquet(file.path(DIR_DATA_POLLEN, "calibrated_predictions_2021.parquet"))
data2022_brut <- read_parquet(file.path(DIR_DATA_POLLEN, "calibrated_predictions_2022.parquet"))
data2023_brut <- read_parquet(file.path(DIR_DATA_POLLEN, "calibrated_predictions_2023.parquet"))

data2023_brut$h5_key[data2023_brut$h5_key == "23_W7_22A1"] <- "23_W7_22A" # identification error on one sample

# year, period and location columns
for (df_name in c("data2021_brut", "data2022_brut", "data2023_brut")) {
  df <- get(df_name)
  df$period <- as.factor(df$window)
  df$year   <- as.factor(paste0("20", substr(df$h5_key, 1, 2)))
  raw_loc <- sub(".*_([0-9]{1,2}[A-Za-z]).*", "\\1", df$h5_key)
  raw_loc <- sub("^([0-9])([A-Za-z])", "0\\1\\2", raw_loc) # pad the location name with a leading 0 where missing
  df$location <- as.factor(raw_loc)
  suffix <- sub(".*[0-9]{1,2}[A-Za-z]", "", df$h5_key)
  df$comments <- ifelse(suffix == "",     NA,
                        ifelse(suffix == "toit", "rooftop",
                               ifelse(suffix == "-1", "new_cyto", suffix)))
  assign(df_name, df)
}

# keep only the identifications above the confidence threshold
data2021 <- data2021_brut[data2021_brut$Confidence >= CONF_MIN, ]
data2022 <- data2022_brut[data2022_brut$Confidence >= CONF_MIN, ]
data2023 <- data2023_brut[data2023_brut$Confidence >= CONF_MIN, ]

# remove gramineae and ambrosia: keep trees only
data2021 <- data2021[!data2021$Pred_Genus %in% TAXATORM, ]
data2022 <- data2022[!data2022$Pred_Genus %in% TAXATORM, ]
data2023 <- data2023[!data2023$Pred_Genus %in% TAXATORM, ]

write.csv(data2021, file.path(DIR_DATA_POLLEN, "data2021.csv"), row.names = FALSE)
write.csv(data2022, file.path(DIR_DATA_POLLEN, "data2022.csv"), row.names = FALSE)
write.csv(data2023, file.path(DIR_DATA_POLLEN, "data2023.csv"), row.names = FALSE)
