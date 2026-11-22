suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# Figure contract
# Core conclusion: MDR-profile richness declined after 2017 in both hosts.
# Archetype: single quantitative validation panel (accumulation curves).
# Backend: R only.
# Source: latest 5,299-isolate phenotype-linked class matrix in Source Data Fig. 2.
# Reviewer risk: MDR and non-MDR profiles must not be mixed; blank class values are
# not treated as a separate profile state; existing profile codes must remain stable.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
root_dir <- if (length(script_arg)) {
  dirname(sub("^--file=", "", script_arg[1]))
} else {
  "."
}
if (!nzchar(root_dir)) root_dir <- "."
input_dir <- file.path(root_dir, "data", "publication_inputs")
tmp_dir <- file.path(root_dir, "outputs", "fig2_profile_update")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
source_file <- file.path(input_dir, "Source_Data_Fig_02.xlsx")
old_table_file <- file.path(input_dir, "Supplementary_Table_04.xlsx")
out_table_file <- file.path(tmp_dir, "Supplementary_Table_04_Revised.xlsx")
out_source_file <- file.path(tmp_dir, "Fig_2c_Revised_source_data.xlsx")
out_fig2_source_file <- file.path(tmp_dir, "Source_Data_Fig_02_Revised.xlsx")
out_prefix <- file.path(tmp_dir, "Fig_2c_Revised")

years <- 2008:2023
class_cols_source <- c("Pen", "BLAI", "Cep", "Car", "Ami", "Tet", "Phe", "Sul", "Flu", "Pol")
# Preserve the class order used by the existing Supplementary Table 4 codes.
class_order_display <- c("Pen", "Lac", "Cep", "Ami", "Tet", "Phe", "Sul", "Flu", "Pol", "Car")

canonical_profile <- function(values_named) {
  values_named[is.na(values_named)] <- 0
  present <- names(values_named)[as.numeric(values_named) == 1]
  present[present == "BLAI"] <- "Lac"
  ordered <- class_order_display[class_order_display %in% present]
  if (length(ordered) == 0) "Pansusceptible" else paste(ordered, collapse = "-")
}

canonical_old_text <- function(x) {
  if (is.na(x)) return(NA_character_)
  if (x == "Pansusceptible") return(x)
  pieces <- strsplit(gsub("BLAI", "Lac", x, fixed = TRUE), "-", fixed = TRUE)[[1]]
  paste(class_order_display[class_order_display %in% pieces], collapse = "-")
}

class_data <- read_excel(source_file, sheet = "c")
mdr_data <- read_excel(source_file, sheet = "b") %>%
  transmute(`Isolate ID` = Identity_Code, MDR = as.integer(MDR))

stopifnot(nrow(class_data) == 5299L)
stopifnot(n_distinct(class_data$`Isolate ID`) == 5299L)

dat <- class_data %>%
  left_join(mdr_data, by = "Isolate ID")
stopifnot(!any(is.na(dat$MDR)))

dat$profile <- apply(dat[, class_cols_source], 1, function(z) {
  names(z) <- class_cols_source
  canonical_profile(z)
})

# Read and audit the historical code map.
old_s1_raw <- read_excel(old_table_file, sheet = "Sheet1", col_names = FALSE)
old_s2_raw <- read_excel(old_table_file, sheet = "Sheet2", col_names = FALSE)

old_map <- tibble(
  profile_id = as.character(old_s1_raw[[1]][3:nrow(old_s1_raw)]),
  profile_display_old = as.character(old_s1_raw[[2]][3:nrow(old_s1_raw)])
) %>%
  filter(!is.na(profile_display_old)) %>%
  mutate(profile = vapply(profile_display_old, canonical_old_text, character(1)))

# Confirm that the old Sheet2 was the cross-host shared-profile subset.
old_numeric <- old_s1_raw[3:nrow(old_s1_raw), ]
old_ids <- as.character(old_numeric[[1]])
old_poultry_matrix <- as.data.frame(lapply(old_numeric[, 3:18], function(x) suppressWarnings(as.numeric(x))))
old_swine_matrix <- as.data.frame(lapply(old_numeric[, 20:35], function(x) suppressWarnings(as.numeric(x))))
old_poultry_total <- rowSums(old_poultry_matrix, na.rm = TRUE)
old_swine_total <- rowSums(old_swine_matrix, na.rm = TRUE)
expected_old_shared <- old_ids[!is.na(old_ids) & old_poultry_total > 0 & old_swine_total > 0]
actual_old_shared <- na.omit(as.character(old_s2_raw[[1]][3:nrow(old_s2_raw)]))
stopifnot(setequal(expected_old_shared, actual_old_shared))

current_profiles <- unique(dat$profile)
known_profiles <- old_map$profile
new_profiles <- setdiff(current_profiles, known_profiles)
max_old_id <- max(as.integer(gsub("Profile ", "", old_map$profile_id[grepl("^Profile [0-9]+$", old_map$profile_id)])))

new_map <- tibble(
  profile_id = paste("Profile", seq.int(max_old_id + 1L, length.out = length(new_profiles))),
  profile_display_old = new_profiles,
  profile = new_profiles
)

profile_map <- bind_rows(old_map, new_map) %>%
  distinct(profile, .keep_all = TRUE)

dat <- dat %>% left_join(profile_map %>% select(profile, profile_id), by = "profile")
stopifnot(all(dat$profile == "Pansusceptible" | !is.na(dat$profile_id)))

# Keep the pansusceptible row unnumbered, matching the original workbook.
dat$profile_id[dat$profile == "Pansusceptible"] <- NA_character_

profile_counts <- dat %>%
  count(profile_id, profile, Host, Year, name = "n") %>%
  complete(nesting(profile_id, profile), Host = c("Poultry", "Swine"), Year = years,
           fill = list(n = 0L))

current_map <- dat %>%
  distinct(profile_id, profile) %>%
  mutate(order_id = ifelse(profile == "Pansusceptible", 0L,
                           as.integer(gsub("Profile ", "", profile_id)))) %>%
  arrange(order_id)

host_profile_presence <- dat %>% distinct(Host, profile)
shared_profiles <- intersect(
  host_profile_presence$profile[host_profile_presence$Host == "Poultry"],
  host_profile_presence$profile[host_profile_presence$Host == "Swine"]
)

make_table_matrix <- function(map_df, sheet2 = FALSE) {
  nrows <- nrow(map_df)
  out <- matrix(NA, nrow = nrows, ncol = 35)
  out[, 1] <- map_df$profile_id
  out[, 2] <- map_df$profile
  for (i in seq_len(nrows)) {
    for (j in seq_along(years)) {
      pv <- profile_counts %>% filter(profile == map_df$profile[i], Host == "Poultry", Year == years[j]) %>% pull(n)
      sv <- profile_counts %>% filter(profile == map_df$profile[i], Host == "Swine", Year == years[j]) %>% pull(n)
      out[i, 2 + j] <- ifelse(length(pv) == 0 || pv == 0, NA, pv)
      out[i, 19 + j] <- ifelse(length(sv) == 0 || sv == 0, NA, sv)
    }
  }
  as.data.frame(out, stringsAsFactors = FALSE)
}

sheet1_map <- current_map
sheet2_map <- current_map %>% filter(profile %in% shared_profiles, profile != "Pansusceptible")
sheet1_data <- make_table_matrix(sheet1_map)
sheet2_data <- make_table_matrix(sheet2_map, sheet2 = TRUE)

write_profile_sheet <- function(wb, sheet, data_df, sheet2 = FALSE) {
  addWorksheet(wb, sheet, gridLines = TRUE)
  writeData(wb, sheet, x = "ID", startCol = 1, startRow = 1, colNames = FALSE)
  writeData(wb, sheet, x = "Profile composition", startCol = 2, startRow = 1, colNames = FALSE)
  writeData(wb, sheet, x = if (sheet2) "No. poultry isolates" else "Number of poultry isolates",
            startCol = 3, startRow = 1, colNames = FALSE)
  writeData(wb, sheet, x = if (sheet2) "No. swine isolates" else "Number of swine isolates",
            startCol = 20, startRow = 1, colNames = FALSE)
  if (sheet2) {
    writeData(wb, sheet, x = "Antibiogram", startCol = 2, startRow = 2, colNames = FALSE)
  }
  writeData(wb, sheet, x = years, startCol = 3, startRow = 2, colNames = FALSE)
  writeData(wb, sheet, x = years, startCol = 20, startRow = 2, colNames = FALSE)
  writeData(wb, sheet, x = data_df, startCol = 1, startRow = 3, colNames = FALSE, keepNA = FALSE)

  mergeCells(wb, sheet, cols = 1, rows = 1:2)
  if (!sheet2) mergeCells(wb, sheet, cols = 2, rows = 1:2)
  mergeCells(wb, sheet, cols = 3:18, rows = 1)
  mergeCells(wb, sheet, cols = 20:35, rows = 1)

  header_style <- createStyle(fontName = "Arial", fontSize = 10, textDecoration = "bold",
                              halign = "center", valign = "center", border = "TopBottomLeftRight",
                              borderStyle = "thin", wrapText = TRUE)
  body_style <- createStyle(fontName = "Arial", fontSize = 9, halign = "center", valign = "center",
                            border = "TopBottomLeftRight", borderStyle = "thin")
  profile_style <- createStyle(fontName = "Arial", fontSize = 9, halign = "left", valign = "center",
                               border = "TopBottomLeftRight", borderStyle = "thin")
  spacer_style <- createStyle(border = "TopBottomLeftRight", borderStyle = "none")

  addStyle(wb, sheet, header_style, rows = 1:2, cols = c(1:18, 20:35), gridExpand = TRUE, stack = TRUE)
  if (nrow(data_df) > 0) {
    addStyle(wb, sheet, body_style, rows = 3:(nrow(data_df) + 2), cols = c(1, 3:18, 20:35), gridExpand = TRUE)
    addStyle(wb, sheet, profile_style, rows = 3:(nrow(data_df) + 2), cols = 2, gridExpand = TRUE)
    addStyle(wb, sheet, spacer_style, rows = 1:(nrow(data_df) + 2), cols = 19, gridExpand = TRUE)
  }
  setColWidths(wb, sheet, cols = 1, widths = if (sheet2) 13.33 else 12.33)
  setColWidths(wb, sheet, cols = 2, widths = if (sheet2) 38.11 else 36.11)
  setColWidths(wb, sheet, cols = 3:18, widths = if (sheet2) 4.89 else 4.78)
  setColWidths(wb, sheet, cols = 19, widths = if (sheet2) 1.44 else 1.78)
  setColWidths(wb, sheet, cols = 20:35, widths = if (sheet2) 5.78 else 4.78)
  setRowHeights(wb, sheet, rows = 1:2, heights = 22)
  freezePane(wb, sheet, firstActiveRow = 3, firstActiveCol = 3)
}

wb <- createWorkbook()
write_profile_sheet(wb, "Sheet1", sheet1_data, sheet2 = FALSE)
write_profile_sheet(wb, "Sheet2", sheet2_data, sheet2 = TRUE)
saveWorkbook(wb, out_table_file, overwrite = TRUE)

# Synchronize the seven named profiles in Source Data Fig. 2d while preserving
# the existing workbook layout and the historical profile codes.
key_ids_compact <- c("Profile152", "Profile211", "Profile238", "Profile239",
                     "Profile255", "Profile258", "Profile263")
key_ids_spaced <- sub("Profile", "Profile ", key_ids_compact, fixed = TRUE)
fig2_wb <- loadWorkbook(source_file)
fig2_d <- as.data.frame(read_excel(source_file, sheet = "d", col_names = FALSE))
host_year_n <- dat %>% count(Host, Year, name = "denominator")

for (i in seq_along(key_ids_compact)) {
  profile_i <- profile_map$profile[match(key_ids_spaced[i], profile_map$profile_id)]
  stopifnot(length(profile_i) == 1L, !is.na(profile_i))
  poultry_counts <- profile_counts %>%
    filter(profile == profile_i, Host == "Poultry") %>%
    arrange(match(Year, years)) %>% pull(n)
  swine_counts <- profile_counts %>%
    filter(profile == profile_i, Host == "Swine") %>%
    arrange(match(Year, years)) %>% pull(n)
  poultry_den <- host_year_n %>% filter(Host == "Poultry") %>%
    arrange(match(Year, years)) %>% pull(denominator)
  swine_den <- host_year_n %>% filter(Host == "Swine") %>%
    arrange(match(Year, years)) %>% pull(denominator)

  quantity_row <- 3L + i
  percentage_row <- 11L + i
  writeData(fig2_wb, "d", t(poultry_counts), startRow = quantity_row, startCol = 3,
            colNames = FALSE, rowNames = FALSE)
  writeData(fig2_wb, "d", sum(poultry_counts), startRow = quantity_row, startCol = 19,
            colNames = FALSE)
  writeData(fig2_wb, "d", t(swine_counts), startRow = quantity_row, startCol = 21,
            colNames = FALSE, rowNames = FALSE)
  writeData(fig2_wb, "d", sum(swine_counts), startRow = quantity_row, startCol = 37,
            colNames = FALSE)
  writeData(fig2_wb, "d", t(100 * poultry_counts / poultry_den),
            startRow = percentage_row, startCol = 3, colNames = FALSE, rowNames = FALSE)
  writeData(fig2_wb, "d", 100 * sum(poultry_counts) / sum(poultry_den),
            startRow = percentage_row, startCol = 19, colNames = FALSE)
  writeData(fig2_wb, "d", t(100 * swine_counts / swine_den),
            startRow = percentage_row, startCol = 21, colNames = FALSE, rowNames = FALSE)
  writeData(fig2_wb, "d", 100 * sum(swine_counts) / sum(swine_den),
            startRow = percentage_row, startCol = 37, colNames = FALSE)
}
saveWorkbook(fig2_wb, out_fig2_source_file, overwrite = TRUE)

# MDR-only accumulation curves for Fig. 2c.
mdr <- dat %>%
  filter(MDR == 1L) %>%
  mutate(Group = paste(Host, Stage, sep = "_"))

group_levels <- c("Poultry_Pre-2017", "Poultry_Post-2017", "Swine_Pre-2017", "Swine_Post-2017")
mdr$Group <- factor(mdr$Group, levels = group_levels)

accumulation_curve <- function(profile_vec, n_repeats = 100L, seed = 20261122L) {
  set.seed(seed)
  n <- length(profile_vec)
  results <- matrix(NA_real_, nrow = n, ncol = n_repeats)
  for (i in seq_len(n_repeats)) {
    shuffled <- sample(profile_vec, size = n, replace = FALSE)
    results[, i] <- cumsum(!duplicated(shuffled))
  }
  tibble(
    SampleSize = seq_len(n),
    Mean = rowMeans(results),
    Lower = apply(results, 1, quantile, probs = 0.025, names = FALSE),
    Upper = apply(results, 1, quantile, probs = 0.975, names = FALSE)
  )
}

curve_list <- lapply(seq_along(group_levels), function(i) {
  g <- group_levels[i]
  x <- mdr %>% filter(Group == g)
  accumulation_curve(x$profile, n_repeats = 100L, seed = 20261122L + i) %>% mutate(Group = g)
})
curve_df <- bind_rows(curve_list) %>%
  mutate(Group = factor(Group, levels = group_levels))

richness_summary <- mdr %>%
  group_by(Host, Stage, Group) %>%
  summarise(n_isolates = n(), profile_richness = n_distinct(profile), .groups = "drop") %>%
  arrange(Group)

expected <- c(162L, 87L, 118L, 88L)
observed <- richness_summary$profile_richness[match(group_levels, as.character(richness_summary$Group))]
stopifnot(identical(as.integer(observed), expected))

all_amr_richness <- dat %>%
  group_by(Host, Stage) %>%
  summarise(n_isolates = n(), profile_richness = n_distinct(profile), .groups = "drop") %>%
  arrange(factor(Host, c("Poultry", "Swine")), factor(Stage, c("Pre-2017", "Post-2017")))
stopifnot(identical(as.integer(all_amr_richness$profile_richness), c(178L, 106L, 130L, 108L)))

group_labels <- c(
  "Poultry_Pre-2017" = "Poultry pre-2017",
  "Poultry_Post-2017" = "Poultry post-2017",
  "Swine_Pre-2017" = "Swine pre-2017",
  "Swine_Post-2017" = "Swine post-2017"
)
cols <- c(
  "Poultry_Pre-2017" = "#9BD7F0",
  "Poultry_Post-2017" = "#1E88E5",
  "Swine_Pre-2017" = "#F6C77A",
  "Swine_Post-2017" = "#F28E2B"
)

label_df <- curve_df %>%
  group_by(Group) %>% slice_tail(n = 1) %>% ungroup() %>%
  mutate(
    label = group_labels[as.character(Group)],
    label_offset = c(
      "Poultry_Pre-2017" = 4,
      "Poultry_Post-2017" = 8,
      "Swine_Pre-2017" = 4,
      "Swine_Post-2017" = -7
    )[as.character(Group)],
    label_y = Mean + label_offset
  )

p <- ggplot(curve_df, aes(x = SampleSize, y = Mean, colour = Group, fill = Group)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.18, colour = NA) +
  geom_line(linewidth = 0.65) +
  geom_text(data = label_df, aes(y = label_y, label = label), hjust = 1.03, vjust = 0,
            size = 2.15, show.legend = FALSE) +
  scale_colour_manual(values = cols, guide = "none") +
  scale_fill_manual(values = cols, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.08))) +
  labs(x = "Sample size", y = "Unique MDR profiles") +
  theme_minimal(base_size = 7, base_family = "Arial") +
  theme(
    panel.grid.major = element_line(linewidth = 0.25, colour = "#E6E6E6"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 7),
    axis.text = element_text(size = 6.4, colour = "black"),
    plot.margin = margin(3, 4, 3, 3, unit = "mm")
  )

width_in <- 90 / 25.4
height_in <- 68 / 25.4
svglite::svglite(paste0(out_prefix, ".svg"), width = width_in, height = height_in)
print(p)
dev.off()
grDevices::cairo_pdf(paste0(out_prefix, ".pdf"), width = width_in, height = height_in, family = "Arial")
print(p)
dev.off()
ragg::agg_tiff(paste0(out_prefix, ".tiff"), width = width_in, height = height_in,
               units = "in", res = 600, background = "white")
print(p)
dev.off()
ragg::agg_png(paste0(out_prefix, ".png"), width = width_in, height = height_in,
              units = "in", res = 300, background = "white")
print(p)
dev.off()

source_wb <- createWorkbook()
addWorksheet(source_wb, "MDR_accumulation_curves")
writeData(source_wb, "MDR_accumulation_curves", curve_df)
addWorksheet(source_wb, "Richness_summary")
writeData(source_wb, "Richness_summary", richness_summary)
addWorksheet(source_wb, "All_AMR_richness_check")
writeData(source_wb, "All_AMR_richness_check", all_amr_richness)
addWorksheet(source_wb, "Profile_code_map")
writeData(source_wb, "Profile_code_map", current_map %>% select(profile_id, profile))
addWorksheet(source_wb, "New_profile_codes")
writeData(source_wb, "New_profile_codes", new_map %>% select(profile_id, profile))
saveWorkbook(source_wb, out_source_file, overwrite = TRUE)

key_ids <- c("Profile 152", "Profile 211", "Profile 238", "Profile 239", "Profile 255", "Profile 258", "Profile 263")
key_check <- current_map %>% filter(profile_id %in% key_ids) %>% arrange(match(profile_id, key_ids))
stopifnot(identical(key_check$profile_id, key_ids))

summary_lines <- c(
  "Fig. 2c update completed using MDR isolates only.",
  paste0("Poultry: ", observed[1], " pre-2017 profiles; ", observed[2], " post-2017 profiles."),
  paste0("Swine: ", observed[3], " pre-2017 profiles; ", observed[4], " post-2017 profiles."),
  "Independent all-profile check (including pansusceptible): poultry 178 to 106; swine 130 to 108.",
  paste0("All observed profiles in Supplementary Table 4 Sheet1: ", nrow(sheet1_map), " (including pansusceptible)."),
  paste0("Cross-host shared resistant profiles in Sheet2: ", nrow(sheet2_map), "."),
  if (nrow(new_map)) {
    paste0("New profile codes assigned: ", paste(new_map$profile_id, collapse = ", "), ".")
  } else {
    "No new profile codes were required; the canonical code map was retained."
  },
  "Key profile codes preserved:",
  paste(capture.output(print(key_check %>% select(profile_id, profile), n = Inf)), collapse = "\n")
)
writeLines(summary_lines, file.path(tmp_dir, "Fig2c_Table4_update_summary.txt"), useBytes = TRUE)

cat(paste(summary_lines, collapse = "\n"), "\n")
