# Shared helpers for sequencing plan CSVs (sourced by map + export scripts).

qbit_col <- function(df) {
  nm <- names(df)[grepl("^Qbit", names(df), ignore.case = TRUE)][1]
  nm
}

#' Vector: TRUE = keep row. Always keep "x" (case-insensitive). Otherwise drop "low"
#' (case-insensitive) and numeric values < 0.01; keep NA/blank and other non-numeric text.
qbit_keep <- function(x) {
  v <- trimws(as.character(x))
  is_x <- tolower(v) == "x"
  is_low <- tolower(v) == "low"
  num <- suppressWarnings(as.numeric(v))
  is_small <- !is.na(num) & num < 0.01
  drop <- !is_x & (is_low | is_small)
  !drop | is.na(drop)
}

#' Reason a value is excluded by qbit_keep(): literal "low" vs numeric < 0.01.
#' (Only meaningful for rows where !qbit_keep(x); otherwise may return "other".)
qubit_exclude_reason <- function(x) {
  v <- trimws(as.character(x))
  is_low <- tolower(v) == "low"
  num <- suppressWarnings(as.numeric(v))
  is_small <- !is.na(num) & num < 0.01
  dplyr::case_when(
    is_low ~ "text_low",
    is_small ~ "numeric_lt_0.01",
    TRUE ~ "other"
  )
}

#' SDLA, BAYS, VALY preserved; everything else (NA, blank, humboldt, …) -> "no region"
region_group <- function(region_vec) {
  r <- trimws(as.character(region_vec))
  r[r %in% c("", "NA")] <- NA_character_
  dplyr::case_when(
    r %in% c("SDLA", "BAYS", "VALY") ~ r,
    TRUE ~ "no region"
  )
}
