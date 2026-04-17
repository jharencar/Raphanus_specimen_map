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

dec_rank <- map_df %>%
  group_by(decade) %>%
  summarise(yr = suppressWarnings(min(as.numeric(as.character(year)), na.rm = TRUE)), .groups = "drop") %>%
  mutate(yr = ifelse(is.finite(yr), yr, NA_real_)) %>%
  arrange(yr, decade)

dec_levels <- dec_rank$decade
map_df <- map_df %>%
  mutate(decade_f = factor(as.character(decade), levels = unique(dec_levels)))

n_dec <- length(levels(map_df$decade_f))
pal_cols <- if (n_dec <= 12) {
  brewer.pal(max(3, n_dec), "Set3")[seq_len(n_dec)]
} else {
  colorRampPalette(brewer.pal(12, "Set3"))(n_dec)
}

dec_pal <- colorFactor(palette = pal_cols, domain = map_df$decade_f)

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

# Optional: quick static check
p <- ggplot(map_df, aes(x = decade_f, fill = decade_f)) +
  geom_bar() +
  scale_fill_manual(values = setNames(pal_cols, levels(map_df$decade_f)), guide = "none") +
  labs(
    title = "Sequencing plan: specimens per decade (Qubit >= 0.01, excluding 'low'; with coordinates)",
    x = NULL,
    y = "Count"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("sequencing_plan_counts_by_decade.pdf", p, width = 10, height = 5, units = "in")
message("Wrote sequencing_plan_counts_by_decade.pdf")
