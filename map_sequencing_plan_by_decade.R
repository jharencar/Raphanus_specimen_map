#!/usr/bin/env Rscript
# Map specimens from Raphanus_lab_database_sequencing_plan.csv, colored by decade.
# Join Raphanus_specimens_all_2026.csv by CCH2 id (Herbarium_ID) and/or catalog number.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(leaflet)
  library(leaflet.extras)
  library(htmlwidgets)
  library(RColorBrewer)
})

proj <- "/Users/admin-jgharenc/Documents/Github/Raphanus_specimen_map"
setwd(proj)
source("sequencing_plan_filters.R")

#' Numeric ordering key for decade labels (e.g. pre.1890, 1890s) when year is missing.
decade_sort_key <- function(decade_chr) {
  d <- trimws(as.character(decade_chr))
  key <- rep(NA_real_, length(d))
  for (i in seq_along(d)) {
    di <- d[i]
    if (is.na(di) || !nzchar(di)) next
    if (grepl("^pre\\.", di, ignore.case = TRUE)) {
      y <- suppressWarnings(as.numeric(sub("^pre\\.(\\d{4}).*$", "\\1", di, perl = TRUE)))
      key[i] <- if (is.finite(y)) y - 0.5 else NA_real_
    } else {
      m <- regmatches(di, regexpr("[0-9]{4}", di))
      key[i] <- if (length(m)) suppressWarnings(as.numeric(m)) else NA_real_
    }
  }
  key
}

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
    specificEpithet = dplyr::coalesce(specificEpithet_1, specificEpithet_2),
    specimen_year = dplyr::coalesce(year_1, year_2),
    county = dplyr::coalesce(county_1, county_2),
    institutionCode_spec = dplyr::coalesce(institutionCode_1, institutionCode_2),
    specimen_id = dplyr::case_when(
      !is.na(dec_lat_1) ~ id_join,
      !is.na(dec_lat_2) ~ join_id_2,
      TRUE ~ NA_real_
    )
  )

map_df <- out %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
  distinct(lab_row, .keep_all = TRUE) %>%
  mutate(qubit_ng_ul = .data[[qc]])

n_lab_raw <- nrow(lab)
n_lab <- nrow(lab2)
n_map <- nrow(map_df)
message(
  "After Qubit filter: ", n_lab, " / ", n_lab_raw, " plan rows | ",
  "Plotted with coordinates: ", n_map, " | Missing coords: ", n_lab - n_map
)

# Decade order: chronological (min plan year per decade, else parsed label)
dec_rank <- lab2 %>%
  group_by(decade) %>%
  summarise(
    yr_min = suppressWarnings(min(as.numeric(as.character(year)), na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    yr_min = if_else(is.finite(yr_min), yr_min, NA_real_),
    sort_key = dplyr::coalesce(yr_min, decade_sort_key(decade))
  ) %>%
  arrange(sort_key, decade)

dec_levels <- dec_rank$decade
map_df <- map_df %>%
  mutate(decade_f = factor(as.character(decade), levels = dec_levels))

n_dec <- length(levels(map_df$decade_f))
# Discrete Spectral (same palette family as year maps); one shade per decade in time order
spectral_11 <- brewer.pal(11, "Spectral")
pal_cols <- if (n_dec <= 1L) {
  spectral_11[6L]
} else {
  colorRampPalette(spectral_11)(n_dec)
}

dec_pal <- colorFactor(palette = pal_cols, domain = levels(map_df$decade_f))

icon_size <- 12
make_triangle_uri <- function(hex_color, size = icon_size, stroke = "black", stroke_width = 1.5) {
  col_enc <- gsub("#", "%23", hex_color)
  sprintf(
    "data:image/svg+xml;utf8,<svg xmlns='http://www.w3.org/2000/svg' width='%d' height='%d' viewBox='0 0 20 20'><polygon points='10,3 18,17 2,17' fill='%s' stroke='%s' stroke-width='%.2f'/></svg>",
    size, size, col_enc, stroke, stroke_width
  )
}

icon_urls <- vapply(
  dec_pal(map_df$decade_f),
  make_triangle_uri,
  character(1),
  size = icon_size
)

tri_icons <- icons(
  iconUrl = icon_urls,
  iconWidth = icon_size,
  iconHeight = icon_size,
  iconAnchorX = round(icon_size / 2),
  iconAnchorY = round(icon_size / 2)
)

m <- leaflet(map_df) %>%
  addTiles() %>%
  addMarkers(
    lng = ~decimalLongitude,
    lat = ~decimalLatitude,
    icon = tri_icons,
    popup = ~paste0(
      "<b>Decade:</b> ", decade,
      "<br><b>Region:</b> ", region,
      "<br><b>Specimen ID:</b> ", ifelse(is.na(specimen_id), "", as.character(specimen_id)),
      "<br><b>Catalog:</b> ", ifelse(is.na(cat_norm), "", cat_norm),
      "<br><b>Herbarium (plan):</b> ", ifelse(is.na(institutionCode), "", as.character(institutionCode)),
      "<br><b>Herbarium (specimen):</b> ", ifelse(is.na(institutionCode_spec), "", as.character(institutionCode_spec)),
      "<br><b>Year (plan):</b> ", ifelse(is.na(year), "", as.character(year)),
      "<br><b>Year (database):</b> ", ifelse(is.na(specimen_year), "", as.character(specimen_year)),
      "<br><b>Taxon:</b> ", ifelse(is.na(specificEpithet), "", as.character(specificEpithet)),
      "<br><b>Lab ID:</b> ", ifelse(is.na(lab_id), "", as.character(lab_id)),
      "<br><b>Qubit ng/µL:</b> ", ifelse(is.na(qubit_ng_ul), "", as.character(qubit_ng_ul))
    ),
    label = ~paste(decade, lab_id),
    group = "seqPlan",
    options = markerOptions(riseOnHover = TRUE)
  ) %>%
  setView(lng = -119.42, lat = 36.78, zoom = 6) %>%
  addLegend(
    position = "bottomright",
    pal = dec_pal,
    values = map_df$decade_f,
    title = "Decade",
    opacity = 1
  ) %>%
  addSearchFeatures(
    targetGroups = "seqPlan",
    options = searchFeaturesOptions(
      zoom = 10,
      openPopup = TRUE,
      textPlaceholder = "Search map..."
    )
  )

out_html <- "sequencing_plan_map_by_decade.html"
has_pandoc <- nzchar(Sys.which("pandoc")) ||
  (requireNamespace("rmarkdown", quietly = TRUE) && rmarkdown::pandoc_available())
if (has_pandoc) {
  saveWidget(m, out_html, selfcontained = TRUE)
} else {
  saveWidget(m, out_html, selfcontained = FALSE)
  message("Pandoc not found; wrote ", out_html, " (not self-contained; may create *_files/).")
}

message("Wrote ", file.path(proj, out_html))

# Histogram: all Qubit-passing plan rows, stacked by primary regions (+ no region)
hist_df <- lab2 %>%
  mutate(
    decade_f = factor(as.character(decade), levels = dec_levels),
    region_f = factor(
      region_group(region),
      levels = c("BAYS", "SDLA", "VALY", "no region")
    )
  )

region_fill <- c(
  "BAYS" = "#1b9e77",
  "SDLA" = "#d95f02",
  "VALY" = "#7570b3",
  "no region" = "#cccccc"
)

p <- ggplot(hist_df, aes(x = decade_f, fill = region_f)) +
  geom_bar(
    position = position_dodge(width = 0.88, preserve = "single"),
    width = 0.8,
    color = "grey30",
    linewidth = 0.15
  ) +
  scale_x_discrete(limits = dec_levels, drop = FALSE) +
  scale_fill_manual(values = region_fill, name = "Region", drop = FALSE) +
  labs(
    title = "Sequencing plan: specimens per decade (Qubit filter applied; all plan rows)",
    subtitle = "Side-by-side bars by BAYS, SDLA, VALY; other/blank regions as 'no region'",
    x = NULL,
    y = "Count"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

ggsave("sequencing_plan_counts_by_decade.pdf", p, width = 11, height = 5.5, units = "in")
message("Wrote sequencing_plan_counts_by_decade.pdf")
