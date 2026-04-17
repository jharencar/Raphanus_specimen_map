#!/usr/bin/env Rscript
# Rows from Raphanus_lab_database_sequencing_plan.csv that still lack coordinates
# after joining Raphanus_specimens_all_2026.csv (same logic as map_sequencing_plan_by_decade.R).
# Applies Qubit filter: drops "low" and numeric values < 0.01 (same as mapping script).

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

proj <- "/Users/admin-jgharenc/Documents/Github/Raphanus_specimen_map"
setwd(proj)
source("sequencing_plan_filters.R")

norm_cat <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "#N/A", "0", "0.0")] <- NA_character_
  x
}

lab <- read.csv(
  "Raphanus_lab_database_sequencing_plan.csv",
  header = TRUE,
  na.strings = c("", "NA", "#N/A")
)

spec <- read.csv(
  "Raphanus_specimens_all_2026.csv",
  header = TRUE,
  na.strings = c("", "NA")
)

spec_id <- spec %>%
  transmute(
    join_id = id,
    dec_lat_1 = decimalLatitude,
    dec_lon_1 = decimalLongitude,
    specificEpithet_1 = specificEpithet,
    year_1 = year,
    county_1 = county,
    institutionCode_1 = institutionCode,
    catalogNumber_1 = as.character(catalogNumber)
  )

spec_cat <- spec %>%
  transmute(
    cat_norm = norm_cat(catalogNumber),
    dec_lat_2 = decimalLatitude,
    dec_lon_2 = decimalLongitude,
    specificEpithet_2 = specificEpithet,
    year_2 = year,
    county_2 = county,
    institutionCode_2 = institutionCode,
    join_id_2 = id
  ) %>%
  filter(!is.na(cat_norm))

qc <- qbit_col(lab)
lab2 <- lab %>%
  mutate(
    lab_row = dplyr::row_number(),
    id_join = suppressWarnings(as.numeric(as.character(id))),
    cat_norm = dplyr::coalesce(
      norm_cat(as.character(catalogNumber)),
      norm_cat(as.character(cat_id))
    )
  ) %>%
  filter(qbit_keep(.data[[qc]]))

out <- lab2 %>%
  left_join(spec_id, by = c("id_join" = "join_id")) %>%
  left_join(spec_cat, by = "cat_norm", relationship = "many-to-one") %>%
  mutate(
    decimalLatitude = dplyr::coalesce(dec_lat_1, dec_lat_2),
    decimalLongitude = dplyr::coalesce(dec_lon_1, dec_lon_2),
    specimen_year = dplyr::coalesce(year_1, year_2),
    county_spec = dplyr::coalesce(county_1, county_2),
    institutionCode_spec = dplyr::coalesce(institutionCode_1, institutionCode_2),
    specimen_id = dplyr::case_when(
      !is.na(dec_lat_1) ~ id_join,
      !is.na(dec_lat_2) ~ join_id_2,
      TRUE ~ NA_real_
    )
  )

missing <- out %>%
  filter(is.na(decimalLatitude) | is.na(decimalLongitude)) %>%
  mutate(
    join_note = dplyr::case_when(
      !is.na(year_1) ~ "matched_by_occurrence_id_DB_missing_coordinates",
      !is.na(year_2) ~ "matched_by_catalog_DB_missing_coordinates",
      !is.na(id_join) ~ "occurrence_id_not_in_Raphanus_specimens_all_2026",
      !is.na(cat_norm) ~ "catalog_not_in_Raphanus_specimens_all_2026",
      TRUE ~ "no_usable_occurrence_id_or_catalog"
    )
  ) %>%
  dplyr::select(
    lab_row,
    join_note,
    region,
    decade,
    id,
    cat_id,
    dplyr::any_of("catalogNumber"),
    cat_norm,
    id_join,
    dplyr::any_of("institutionCode"),
    dplyr::any_of("year"),
    dplyr::any_of("county"),
    dplyr::any_of("locality"),
    dplyr::any_of("lab_id"),
    dplyr::any_of("sampled."),
    dplyr::any_of("extracted."),
    dplyr::any_of("notes"),
    specimen_id,
    institutionCode_spec,
    specimen_year,
    county_spec,
    dplyr::any_of(c("DNA_tissue", "leaf_tissue", "extr_tissue")),
    dplyr::any_of(qc),
    decimalLatitude,
    decimalLongitude
  )

out_csv <- "sequencing_plan_missing_coordinates.csv"
write.csv(missing, out_csv, row.names = FALSE, na = "")

n_after_qubit <- nrow(lab2)
message(
  "Missing-coords rows (after Qubit filter, ", n_after_qubit, " plan rows): ",
  nrow(missing), " -> ", file.path(getwd(), out_csv)
)

dec_all <- lab2 %>%
  dplyr::group_by(decade) %>%
  dplyr::summarise(
    yr = suppressWarnings(min(as.numeric(as.character(year)), na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(yr = dplyr::if_else(is.finite(yr), yr, NA_real_)) %>%
  dplyr::arrange(yr, decade) %>%
  dplyr::pull(decade)

reg_lv <- c("BAYS", "SDLA", "VALY", "no region")

missing_summary_long <- missing %>%
  dplyr::mutate(region_group = region_group(region)) %>%
  dplyr::count(region_group, decade, name = "n") %>%
  tidyr::complete(
    region_group = reg_lv,
    decade = dec_all,
    fill = list(n = 0)
  )

missing_summary_wide <- missing_summary_long %>%
  tidyr::pivot_wider(names_from = decade, values_from = n, values_fill = 0)

dec_cols_present <- dec_all[dec_all %in% names(missing_summary_wide)]
missing_summary_wide <- missing_summary_wide %>%
  dplyr::mutate(region_group = factor(region_group, levels = reg_lv)) %>%
  dplyr::arrange(region_group) %>%
  dplyr::mutate(region_group = as.character(region_group)) %>%
  dplyr::select(region_group, dplyr::all_of(dec_cols_present))

num_dec_cols <- dec_cols_present
missing_summary_wide$row_total <- rowSums(missing_summary_wide[num_dec_cols], na.rm = TRUE)
total_row <- missing_summary_wide %>%
  dplyr::summarise(
    region_group = "TOTAL",
    dplyr::across(dplyr::all_of(num_dec_cols), sum)
  ) %>%
  dplyr::mutate(row_total = rowSums(dplyr::pick(dplyr::all_of(num_dec_cols))))
missing_summary_wide <- dplyr::bind_rows(missing_summary_wide, total_row)

sum_csv <- "sequencing_plan_missing_coords_summary_decade_by_region.csv"
write.csv(missing_summary_wide, sum_csv, row.names = FALSE)
message("Wrote summary table: ", file.path(getwd(), sum_csv))

if (requireNamespace("writexl", quietly = TRUE)) {
  out_xlsx <- "sequencing_plan_missing_coordinates.xlsx"
  writexl::write_xlsx(
    list(
      missing_coordinates = missing,
      summary_decade_by_region = missing_summary_wide,
      summary_long = missing_summary_long
    ),
    out_xlsx
  )
  message("Wrote ", out_xlsx)
}
