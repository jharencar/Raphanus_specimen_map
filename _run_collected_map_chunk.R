# One-off: run "##Collected specimen maps" chunk from Raphanus_specimen_map.qmd
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(htmlwidgets)
  library(leaflet)
  library(leaflet.extras)
})

setwd("/Users/admin-jgharenc/Documents/Github/Raphanus_specimen_map")

# --- chunk: Collected specimen maps ---
CA_rad_spec <- read.csv("Raphanus_specimens_all_2026.csv", header = TRUE)

Nat_coords_path <- "/Users/juliaharencar/Documents/Github/Raphanus_Phenology/Output_Files/completed_specimen_data.csv"
if (file.exists(Nat_coords_path)) {
  Nat_coords <- read.csv(Nat_coords_path, header = TRUE)
  CA_rad_spec <- CA_rad_spec %>%
    rows_update(
      Nat_coords %>% dplyr::select(id, decimalLatitude, decimalLongitude),
      by = "id",
      unmatched = "ignore"
    )
}
if (file.exists("CA_rad_spec_1915to40_missing.coords.csv")) {
  more_coords <- read.csv("CA_rad_spec_1915to40_missing.coords.csv", header = TRUE)
  more_coords$decimalLatitude <- as.numeric(more_coords$decimalLatitude)
  CA_rad_spec <- CA_rad_spec %>%
    rows_update(
      more_coords %>% dplyr::select(id, decimalLatitude, decimalLongitude),
      by = "id",
      unmatched = "ignore"
    )
}
if (file.exists("CA_rad_spec_sampled.to.251114_missing.coords_w_some_filled.csv")) {
  yet_more <- read.csv("CA_rad_spec_sampled.to.251114_missing.coords_w_some_filled.csv", header = TRUE)
  yet_more$decimalLatitude <- as.numeric(yet_more$decimalLatitude)
  CA_rad_spec <- CA_rad_spec %>%
    rows_update(
      yet_more %>% dplyr::select(id, decimalLatitude, decimalLongitude),
      by = "id",
      unmatched = "ignore"
    )
}

herb <- read.csv("Herbarium_tissue_collected_260409.csv", header = TRUE, check.names = FALSE)
names(herb) <- ifelse(is.na(names(herb)) | names(herb) == "", paste0("V", seq_along(names(herb))), names(herb))

lookup <- CA_rad_spec %>%
  dplyr::select(id, institutionCode, catalogNumber, decimalLatitude, decimalLongitude,
                specificEpithet, year, county) %>%
  filter(!is.na(catalogNumber), catalogNumber != "")

herb_filled <- herb %>%
  mutate(
    id_char = as.character(id),
    id_filled = if_else(
      is.na(id) | id_char == "" | id_char == "NA",
      NA_character_,
      id_char
    )
  ) %>%
  left_join(
    lookup %>% dplyr::select(id, catalogNumber, institutionCode) %>% rename(id_from_db = id),
    by = c("catalogNumber" = "catalogNumber", "institutionCode" = "institutionCode")
  ) %>%
  mutate(
    id_filled = if_else(
      is.na(id_filled) | id_filled == "",
      as.character(id_from_db),
      id_filled
    )
  )

still_missing <- herb_filled %>% filter(is.na(id_filled) | id_filled == "")
if (nrow(still_missing) > 0) {
  by_cat <- still_missing %>%
    dplyr::select(catalogNumber) %>%
    distinct() %>%
    left_join(
      lookup %>% dplyr::select(id, catalogNumber) %>% distinct(catalogNumber, .keep_all = TRUE),
      by = "catalogNumber"
    )
  herb_filled <- herb_filled %>%
    left_join(by_cat %>% rename(id_by_cat = id), by = "catalogNumber") %>%
    mutate(
      id_filled = case_when(
        !is.na(id_filled) & id_filled != "" ~ id_filled,
        !is.na(id_by_cat) ~ as.character(id_by_cat),
        TRUE ~ id_filled
      )
    )
}

herb_filled <- herb_filled %>%
  mutate(id_resolved = suppressWarnings(as.numeric(id_filled))) %>%
  filter(!is.na(id_resolved))

selected <- unique(herb_filled$id_resolved)

samp_list_w.data <- CA_rad_spec %>%
  filter(id %in% selected) %>%
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
  distinct(id, .keep_all = TRUE)

message("Sampled specimens with coordinates: ", nrow(samp_list_w.data),
        " (rows in tissue list with resolved id: ", nrow(herb_filled), ")")

yr_max <- max(samp_list_w.data$year, na.rm = TRUE)
yr_min <- min(samp_list_w.data$year, na.rm = TRUE)

p_hist <- ggplot(samp_list_w.data, aes(x = year)) +
  geom_histogram(binwidth = 5, col = "black", fill = "darkolivegreen3") +
  scale_x_continuous(breaks = seq(floor(yr_min / 5) * 5, ceiling(yr_max / 5) * 5, by = 5)) +
  labs(title = "Herbarium tissue samples: collection counts by year", x = "Collection year", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave("Herbarium_sampled_histogram_by_year_260409.pdf", p_hist, width = 10, height = 5, units = "in")
message("Wrote Herbarium_sampled_histogram_by_year_260409.pdf")

make_triangle_uri <- function(hex_color, size = 10, stroke = "black", stroke_width = 2) {
  col_enc <- gsub("#", "%23", hex_color)
  sprintf(
    "data:image/svg+xml;utf8,
    <svg xmlns='http://www.w3.org/2000/svg' width='%d' height='%d' viewBox='0 0 20 20'>
      <polygon points='10,3 18,17 2,17' fill='%s' stroke='%s' stroke-width='%.2f'/>
    </svg>",
    size, size, col_enc, stroke, stroke_width
  )
}

selected_df <- samp_list_w.data %>%
  mutate(year = as.numeric(year))

colors <- colorNumeric("Spectral", domain = c(yr_min, yr_max))

icon_size <- 12
icon_urls <- vapply(
  colors(selected_df$year),
  FUN = make_triangle_uri,
  FUN.VALUE = character(1),
  size = icon_size
)

tri_icons <- icons(
  iconUrl = icon_urls,
  iconWidth = icon_size,
  iconHeight = icon_size,
  iconAnchorX = round(icon_size / 2),
  iconAnchorY = round(icon_size / 2)
)

full_specimen_map_w_sampling <- leaflet(selected_df) %>%
  addTiles() %>%
  addMarkers(
    lat = ~decimalLatitude,
    lng = ~decimalLongitude,
    icon = tri_icons,
    popup = ~paste(
      "Scientific Name:", specificEpithet,
      "<br>Year:", year,
      "<br>Herbarium:", institutionCode,
      "<br>ID:", id,
      "<br>Catalog:", catalogNumber
    ),
    label = ~paste(specificEpithet, institutionCode, year, id),
    group = "sampledPoints",
    options = markerOptions(riseOnHover = TRUE)
  ) %>%
  setView(lng = -119.4179, lat = 36.7783, zoom = 6) %>%
  addLegend(
    position = "bottomright",
    pal = colors,
    values = seq(yr_min, yr_max, length.out = 5),
    title = "Collection year",
    opacity = 1
  ) %>%
  addSearchFeatures(
    targetGroups = "sampledPoints",
    options = searchFeaturesOptions(
      zoom = 10,
      openPopup = TRUE,
      textPlaceholder = "Search map..."
    )
  )

saveWidget(full_specimen_map_w_sampling, "triangle_sampled_CA_specimen_color_by_year_260409.html", selfcontained = TRUE)
message("Wrote triangle_sampled_CA_specimen_color_by_year_260409.html")
