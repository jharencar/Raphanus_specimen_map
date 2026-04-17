#!/usr/bin/env Rscript
# Herbarium_tissue_collected.csv rows whose specimen is NOT in
# Raphanus_lab_database_sequencing_plan.csv after the same Qubit filter as the map script.
# Join Raphanus_specimens_all_2026.csv by catalogNumber (then by id) for locality & coordinates.
# Note: all_2026 has no verbatimLocality field; locality + verbatimCoordinates are exported.

suppressPackageStartupMessages(library(dplyr))

proj <- "/Users/admin-jgharenc/Documents/Github/Raphanus_specimen_map"
setwd(proj)
source("sequencing_plan_filters.R")

norm_cat <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "#N/A", "0", "0.0")] <- NA_character_
  x
}

plan <- read.csv(
  "Raphanus_lab_database_sequencing_plan.csv",
  header = TRUE,
  na.strings = c("", "NA", "#N/A")
)
qc <- qbit_col(plan)
plan_f <- plan %>%
  mutate(
    id_key = suppressWarnings(as.numeric(as.character(id))),
    cat_norm = dplyr::coalesce(
      norm_cat(as.character(catalogNumber)),
      norm_cat(as.character(cat_id))
    )
  ) %>%
  filter(qbit_keep(.data[[qc]]))

plan_id_set <- unique(stats::na.omit(plan_f$id_key))
plan_cat_set <- unique(plan_f$cat_norm[!is.na(plan_f$cat_norm)])

herb <- read.csv(
  "Herbarium_tissue_collected.csv",
  header = TRUE,
  check.names = FALSE,
  na.strings = c("", "NA", "#N/A")
)
names(herb) <- ifelse(
  is.na(names(herb)) | names(herb) == "",
  paste0("V", seq_along(names(herb))),
  names(herb)
)

herb2 <- herb %>%
  mutate(
    tissue_sheet_row = dplyr::row_number(),
    id_key = suppressWarnings(as.numeric(as.character(id))),
    cat_norm = norm_cat(as.character(catalogNumber))
  )

represented <- (herb2$id_key %in% plan_id_set) |
  (!is.na(herb2$cat_norm) & herb2$cat_norm %in% plan_cat_set)

missing <- herb2 %>% filter(!represented)

spec <- read.csv(
  "Raphanus_specimens_all_2026.csv",
  header = TRUE,
  na.strings = c("", "NA")
)

spec_cat <- spec %>%
  mutate(cat_norm = norm_cat(as.character(catalogNumber))) %>%
  filter(!is.na(cat_norm)) %>%
  arrange(id) %>%
  group_by(cat_norm) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(
    cat_norm,
    specimen_id_match = id,
    locality = locality,
    verbatimCoordinates = verbatimCoordinates,
    decimalLatitude = decimalLatitude,
    decimalLongitude = decimalLongitude,
    county = county,
    municipality = municipality,
    institutionCode_spec = institutionCode
  )

spec_id <- spec %>%
  transmute(
    id_key = id,
    specimen_id_match_id = id,
    locality_id = locality,
    verbatimCoordinates_id = verbatimCoordinates,
    decimalLatitude_id = decimalLatitude,
    decimalLongitude_id = decimalLongitude,
    county_id = county,
    municipality_id = municipality,
    institutionCode_spec_id = institutionCode
  )

out <- missing %>%
  left_join(spec_cat, by = "cat_norm") %>%
  left_join(spec_id, by = "id_key") %>%
  mutate(
    specimen_id = dplyr::coalesce(specimen_id_match, specimen_id_match_id),
    institutionCode_specimen = dplyr::coalesce(institutionCode_spec, institutionCode_spec_id),
    locality = dplyr::coalesce(locality, locality_id),
    verbatimCoordinates = dplyr::coalesce(verbatimCoordinates, verbatimCoordinates_id),
    verbatimLocality = NA_character_,
    decimalLatitude = dplyr::coalesce(decimalLatitude, decimalLatitude_id),
    decimalLongitude = dplyr::coalesce(decimalLongitude, decimalLongitude_id),
    county = dplyr::coalesce(county, county_id),
    municipality = dplyr::coalesce(municipality, municipality_id)
  ) %>%
  select(
    tissue_sheet_row,
    id_tissue = id,
    id_key,
    dplyr::any_of("institutionCode"),
    dplyr::any_of("catalogNumber"),
    cat_norm,
    dplyr::any_of("year"),
    specimen_id,
    institutionCode_specimen,
    locality,
    verbatimLocality,
    verbatimCoordinates,
    decimalLatitude,
    decimalLongitude,
    county,
    municipality,
    dplyr::any_of(c("DNA tissue", "leaf_tissue", "notes"))
  )

out_csv <- "tissue_not_in_sequencing_plan_after_qubit.csv"
write.csv(out, out_csv, row.names = FALSE, na = "")

message(
  "Tissue rows: ", nrow(herb2),
  " | In Qubit-filtered sequencing plan (id or catalog): ", sum(represented),
  " | NOT in plan: ", nrow(out),
  " -> ", normalizePath(out_csv)
)
message(
  "Note: verbatimLocality is not in Raphanus_specimens_all_2026.csv; column is blank; ",
  "use locality and verbatimCoordinates from the specimen record."
)
