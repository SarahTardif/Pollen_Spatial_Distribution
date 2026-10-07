# preparation step, run ONCE (not sourced by the other scripts): reads the 3
# yearly parquet files, derives period/year/location/comments from h5_key,
# filters at Confidence >= CONF_MIN and writes the CSVs the rest of the project
# reads through charger_pollen(). Re-run after any change to the cleaning logic
# below, otherwise the downstream scripts keep reading the output CSVs.

source("R/init.R")
library(arrow)

data_brut <- read_parquet(file.path(DIR_DATA_POLLEN, "Data_pollen_Maya.parquet"))

# year, period and location columns
data_brut$period <- as.factor(data_brut$window)

raw_loc <- sub("^[0-9]{1,2}_([A-Za-z]).*", "\\1", data_brut$h5_key)
data_brut$location <- as.factor(raw_loc)

suffix <- sub("^[0-9]{1,2}_[A-Za-z]", "", data_brut$h5_key)
data_brut$comments <- dplyr::case_when(
  suffix == ""     ~ NA_character_,
  suffix == "_erreur_manip" ~ "erreur_manip",
  TRUE             ~ suffix
)

# keep only the identifications above the confidence threshold
dataMaya <- data_brut[data_brut$Confidence >= CONF_MIN, ]


# remove gramineae and ambrosia: keep trees only
dataMaya <- dataMaya[!dataMaya$Pred_Genus %in% TAXATORM, ]


write.csv(dataMaya, file.path(DIR_DATA_POLLEN, "dataMaya.csv"), row.names = FALSE)

