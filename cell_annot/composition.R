## ============================================================
## Panel G Equivalent: Per-Cluster Retention Overlap
##
## Stacked barplot showing per ddqc cluster what fraction
## of cells was retained by:
##   - All three methods (ddqc + miQC + standard cutoff)
##   - ddqc only
##   - miQC only
##   - Standard cutoff only
##   - ddqc + miQC only
##   - ddqc + standard cutoff only
##   - miQC + standard cutoff only
##   - None (removed by all)
##
## Mirrors Fig 3G from Subramanian et al. Genome Biology 2022
## but extended to three methods including miQC
##
## Input:
##   - compare_annot/ddqc_classification_mad_results_res_1.3.csv
##   - compare_annot/miqc_classification_mad_results.csv
##   - compare_annot/std_classification_mad_results.csv
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})


# ---- PATHS --------------------------------------------------------------------

BASE_DIR <- "/projectnb/ds596/projects/Team 9/scQC_project"
comp_dir <- file.path(BASE_DIR, "compare_annot")
out_dir  <- file.path(BASE_DIR, "compare_annot/panel_g")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---- HELPERS -----------------------------------------------------------------

bool_fix <- function(df) {
  cols <- grep("_passed_qc$|^keep_", colnames(df), value = TRUE)
  for (col in cols) {
    if (is.character(df[[col]])) df[[col]] <- df[[col]] == "True"
  }
  df
}


# ---- 1. LOAD DATA ------------------------------------------------------------

message("Loading MAD results...")

ddqc_path <- file.path(comp_dir, "ddqc_classification_mad_results_res_1.3.csv")
miqc_path <- file.path(comp_dir, "miqc_classification_mad_results.csv")
std_path  <- file.path(comp_dir, "std_classification_mad_results.csv")

ddqc_df <- read.csv(ddqc_path, row.names = 1) %>% bool_fix()
miqc_df <- read.csv(miqc_path, row.names = 1) %>% bool_fix()
std_df  <- read.csv(std_path,  row.names = 1) %>% bool_fix()

# Rename cluster columns
if ("cluster" %in% colnames(ddqc_df) && !"ddqc_cluster" %in% colnames(ddqc_df))
  ddqc_df <- ddqc_df %>% rename(ddqc_cluster = cluster)
if ("cluster" %in% colnames(std_df) && !"std_cluster" %in% colnames(std_df))
  std_df <- std_df %>% rename(std_cluster = cluster)

message(sprintf("ddqc: %d cells | miQC: %d cells | Std: %d cells",
                nrow(ddqc_df), nrow(miqc_df), nrow(std_df)))


# ---- 2. BUILD OVERLAP TABLE --------------------------------------------------
# Join all three methods by cell barcode
# ddqc_cluster is used as the reference cluster assignment
# (all cells that passed basic filter have a ddqc_cluster)

# Base: ddqc cells (all cells post basic filter)
base_df <- ddqc_df %>%
  select(ddqc_cluster, keep_ddqc) %>%
  mutate(cell = rownames(ddqc_df))

# Join miQC
miqc_keep <- data.frame(
  cell       = rownames(miqc_df),
  keep_miqc  = miqc_df$keep_miqc,
  stringsAsFactors = FALSE
)

# Join standard cutoff
std_keep <- data.frame(
  cell      = rownames(std_df),
  keep_std  = std_df$keep_std,
  stringsAsFactors = FALSE
)

combined <- base_df %>%
  left_join(miqc_keep, by = "cell") %>%
  left_join(std_keep,  by = "cell") %>%
  mutate(
    keep_ddqc = ifelse(is.na(keep_ddqc), FALSE, keep_ddqc),
    keep_miqc = ifelse(is.na(keep_miqc), FALSE, keep_miqc),
    keep_std  = ifelse(is.na(keep_std),  FALSE, keep_std)
  )


# ---- 3. CLASSIFY EACH CELL INTO RETENTION CATEGORY --------------------------

combined <- combined %>%
  mutate(
    retention_category = case_when(
      keep_ddqc  & keep_miqc  & keep_std  ~ "All three",
      keep_ddqc  & keep_miqc  & !keep_std ~ "ddqc + miQC",
      keep_ddqc  & !keep_miqc & keep_std  ~ "ddqc + Standard",
      !keep_ddqc & keep_miqc  & keep_std  ~ "miQC + Standard",
      keep_ddqc  & !keep_miqc & !keep_std ~ "ddqc only",
      !keep_ddqc & keep_miqc  & !keep_std ~ "miQC only",
      !keep_ddqc & !keep_miqc & keep_std  ~ "Standard only",
      TRUE                                ~ "None retained"
    )
  )

# Save per-cell table
write.csv(combined,
          file.path(out_dir, "per_cell_retention_categories.csv"),
          row.names = FALSE)


# ---- 4. COMPUTE PER-CLUSTER PROPORTIONS --------------------------------------

cluster_props <- combined %>%
  group_by(ddqc_cluster, retention_category) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(ddqc_cluster) %>%
  mutate(
    total = sum(n),
    pct   = 100 * n / total
  ) %>%
  ungroup()

# Order clusters numerically
cluster_order <- as.character(
  sort(as.numeric(unique(cluster_props$ddqc_cluster)))
)
cluster_props$ddqc_cluster <- factor(cluster_props$ddqc_cluster,
                                      levels = cluster_order)

# Color scheme matching paper (Fig 3G) extended for miQC
category_colors <- c(
  "All three"       = "#2c7bb6",   # blue  — retained by all
  "ddqc + miQC"     = "#74add1",   # light blue
  "ddqc + Standard" = "#abd9e9",   # lighter blue
  "miQC + Standard" = "#74c476",   # green
  "ddqc only"       = "#d9ef8b",   # yellow-green — ddqc unique
  "miQC only"       = "#fee08b",   # yellow
  "Standard only"   = "#fc8d59",   # orange — standard cutoff unique
  "None retained"   = "#d7d7d7"    # grey  — removed by all
)

category_order <- c(
  "All three", "ddqc + miQC", "ddqc + Standard", "miQC + Standard",
  "ddqc only", "miQC only", "Standard only", "None retained"
)

cluster_props$retention_category <- factor(
  cluster_props$retention_category,
  levels = category_order
)


# ---- 5. PANEL G: STACKED BARPLOT ---------------------------------------------

message("Plotting Panel G stacked barplot...")

n_clusters <- length(unique(cluster_props$ddqc_cluster))

p_bar <- ggplot(cluster_props,
                aes(x = ddqc_cluster, y = pct,
                    fill = retention_category)) +
  geom_bar(stat = "identity", width = 0.85) +
  scale_fill_manual(
    values = category_colors,
    breaks = category_order,
    name   = "Retention Status"
  ) +
  scale_y_continuous(
    limits = c(0, 101),
    expand = expansion(mult = c(0, 0)),
    labels = function(x) paste0(x, "%")
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.caption    = element_text(size = 8, color = "grey50", hjust = 0),
    axis.text.x     = element_text(angle = 45, hjust = 1,
                                    size = max(6, min(10, 200 / n_clusters))),
    axis.text.y     = element_text(size = 10),
    axis.title.x    = element_blank(),
    axis.title.y    = element_text(size = 11),
    legend.title    = element_text(size = 10, face = "bold"),
    legend.text     = element_text(size = 9),
    legend.key.size = unit(0.4, "cm"),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title   = "Proportion of Cells Retained by QC Method per ddqc Cluster",
    y       = "% Cells",
    caption = paste0(
      "ddqc resolution = 1.3 | miQC probabilistic mixture model | ",
      "Standard cutoff: nFeature_RNA \u2265 200, percent_mito < 10%\n",
      "Cluster IDs from ddqc (res=1.3) graph-based clustering"
    )
  )

ggsave(
  file.path(out_dir, "panel_g_retention_barplot.png"),
  p_bar,
  width  = max(10, n_clusters * 0.35),
  height = 6,
  dpi    = 450,
  limitsize = FALSE
)
message("Saved: panel_g_retention_barplot.png")


# ---- 6. SUMMARY: OVERALL RETENTION COUNTS ------------------------------------

overall_summary <- combined %>%
  count(retention_category, name = "n_cells") %>%
  mutate(
    pct = round(100 * n_cells / nrow(combined), 2),
    retention_category = factor(retention_category, levels = category_order)
  ) %>%
  arrange(retention_category)

message("\nOverall retention summary:")
print(overall_summary)

write.csv(overall_summary,
          file.path(out_dir, "overall_retention_summary.csv"),
          row.names = FALSE)


# ---- 7. SIMPLIFIED VERSION (matching paper Fig 3G style) --------------------
# Collapse into: retained by all, ddqc only, standard only, miQC only, none

combined_simple <- combined %>%
  mutate(
    simple_category = case_when(
      keep_ddqc & keep_miqc & keep_std  ~ "Retained by all",
      keep_ddqc & !keep_miqc & !keep_std ~ "ddqc only",
      !keep_ddqc & !keep_miqc & keep_std ~ "Standard cutoff only",
      !keep_ddqc & keep_miqc & !keep_std ~ "miQC only",
      TRUE ~ "Other / None"
    )
  )

cluster_simple <- combined_simple %>%
  group_by(ddqc_cluster, simple_category) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(ddqc_cluster) %>%
  mutate(pct = 100 * n / sum(n)) %>%
  ungroup() %>%
  mutate(
    ddqc_cluster    = factor(ddqc_cluster, levels = cluster_order),
    simple_category = factor(simple_category, levels = c(
      "Retained by all", "ddqc only",
      "Standard cutoff only", "miQC only", "Other / None"
    ))
  )

simple_colors <- c(
  "Retained by all"      = "#2c7bb6",
  "ddqc only"            = "#d9ef8b",
  "Standard cutoff only" = "#fc8d59",
  "miQC only"            = "#74c476",
  "Other / None"         = "#d7d7d7"
)

p_simple <- ggplot(cluster_simple,
                   aes(x = ddqc_cluster, y = pct,
                       fill = simple_category)) +
  geom_bar(stat = "identity", width = 0.85) +
  scale_fill_manual(values = simple_colors, name = "Retention Status") +
  scale_y_continuous(
    limits = c(0, 101),
    expand = expansion(mult = c(0, 0)),
    labels = function(x) paste0(x, "%")
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title         = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.caption       = element_text(size = 8, color = "grey50", hjust = 0),
    axis.text.x        = element_text(angle = 45, hjust = 1,
                                       size = max(6, min(10, 200 / n_clusters))),
    axis.text.y        = element_text(size = 10),
    axis.title.x       = element_blank(),
    axis.title.y       = element_text(size = 11),
    legend.title       = element_text(size = 10, face = "bold"),
    legend.text        = element_text(size = 9),
    legend.key.size    = unit(0.4, "cm"),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title   = "Proportion of Cells Retained by QC Method per ddqc Cluster",
    y       = "% Cells",
    caption = paste0(
      "ddqc res=1.3 | miQC | Standard cutoff (nFeature_RNA \u2265 200, percent_mito < 10%)\n",
      "Mirrors Fig 3G — Subramanian et al. Genome Biology 23:267 (2022)"
    )
  )

ggsave(
  file.path(out_dir, "panel_g_simple_barplot.png"),
  p_simple,
  width  = max(10, n_clusters * 0.35),
  height = 6,
  dpi    = 450,
  limitsize = FALSE
)
message("Saved: panel_g_simple_barplot.png")


# ---- SESSION INFO -------------------------------------------------------------

sink(file.path(out_dir, "session_info.txt"))
sessionInfo()
sink()

message("\n=== Done. Outputs in: ", out_dir, " ===")
message("Outputs:")
message("  panel_g_retention_barplot.png  — full 8-category stacked barplot")
message("  panel_g_simple_barplot.png     — simplified 5-category (mirrors paper)")
message("  per_cell_retention_categories.csv")
message("  overall_retention_summary.csv")