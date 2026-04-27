## ============================================================
## Cross-tabulation Heatmap: Cluster ID vs CellTypist/SingleR
##
## For each condition:
##   - ddqc resolutions (0.4, 0.8, 1.3, 1.6, 2.0)
##   - Standard cutoff (res=1.4)
##
## Rows    = cluster IDs (from ddqc or standard cutoff)
## Columns = CellTypist / SingleR labels
## Values  = % of cluster cells with that label
##
## High diagonal dominance = clusters are biologically coherent
##
## Inputs:
##   - compare_annot/ddqc_classification_mad_results_res_{res}.csv
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
out_dir  <- file.path(BASE_DIR, "compare_annot/crosstab_heatmap")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

resolutions <- c(0.4, 0.8, 1.3, 1.6, 2)


# ---- HELPERS -----------------------------------------------------------------

bool_fix <- function(df) {
  cols <- grep("_passed_qc$|^keep_", colnames(df), value = TRUE)
  for (col in cols) {
    if (is.character(df[[col]])) df[[col]] <- df[[col]] == "True"
  }
  df
}

#' Build row-normalised cross-tabulation matrix
#' rows = cluster IDs, columns = annotation labels
#' values = % of cluster cells with that annotation label

build_crosstab <- function(df, cluster_col, annot_col, keep_col) {

  sub <- df %>%
    filter(!is.na(.data[[keep_col]]),
           .data[[keep_col]] == TRUE,
           !is.na(.data[[cluster_col]]),
           !is.na(.data[[annot_col]]),
           .data[[annot_col]] != "",
           .data[[annot_col]] != "Unknown")

  if (nrow(sub) == 0) return(NULL)

  # Raw counts
  tab <- table(
    cluster   = sub[[cluster_col]],
    cell_type = sub[[annot_col]]
  )

  # Row-normalise to % within each cluster
  tab_pct <- prop.table(tab, margin = 1) * 100

  as.data.frame(tab_pct) %>%
    rename(cluster = cluster, cell_type = cell_type, pct = Freq)
}

#' Order clusters numerically if possible, otherwise alphabetically
order_clusters <- function(x) {
  nums <- suppressWarnings(as.numeric(as.character(x)))
  if (!any(is.na(nums))) {
    factor(x, levels = as.character(sort(unique(nums))))
  } else {
    factor(x, levels = sort(unique(as.character(x))))
  }
}

#' Order cell types by their peak cluster (column with max value)
order_celltypes <- function(mat_df) {
  peak <- mat_df %>%
    group_by(cell_type) %>%
    slice_max(pct, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(as.numeric(as.character(cluster))) %>%
    pull(cell_type)
  factor(mat_df$cell_type, levels = unique(peak))
}


#' Plot one cross-tabulation heatmap
#'
#' @param mat_df    data.frame with columns: cluster, cell_type, pct
#' @param title     plot title
#' @param n_cells   total number of cells used (for subtitle)

plot_crosstab <- function(mat_df, title, n_cells = NULL) {

  mat_df$cluster   <- order_clusters(mat_df$cluster)
  mat_df$cell_type <- order_celltypes(mat_df)

  subtitle <- if (!is.null(n_cells))
    paste0("n = ", scales::comma(n_cells), " retained cells") else NULL

  n_clusters  <- length(unique(mat_df$cluster))
  n_celltypes <- length(unique(mat_df$cell_type))

  # Scale all text and geometry proportionally to plot dimensions
  # physical width ~ n_celltypes * 0.28, height ~ n_clusters * 0.22 (inches)
  phys_w    <- max(5, n_celltypes * 0.28)
  phys_h    <- max(3, n_clusters  * 0.22)

  # Base font size scales with smaller dimension
  base_size  <- max(7, min(11, min(phys_w, phys_h) * 1.1))

  # Axis tick font: shrink further when many labels
  xtick_size <- max(6, min(base_size, phys_w * 0.9 / n_celltypes * 10))
  ytick_size <- max(6, min(base_size, phys_h * 0.9 / n_clusters  * 10))

  # In-tile label size (ggplot size units ~ pt / 2.845)
  label_size <- max(1.8, min(3.0, min(phys_w / n_celltypes,
                                       phys_h / n_clusters) * 3.5))

  # Tile border width proportional to tile size
  tile_lwd   <- max(0.1, min(0.4, min(phys_w / n_celltypes,
                                       phys_h / n_clusters) * 0.5))

  # Only show % labels for cells > 5% to avoid clutter
  mat_df$label <- ifelse(mat_df$pct >= 5,
                         sprintf("%.0f", mat_df$pct), "")

  ggplot(mat_df, aes(x = cell_type, y = cluster, fill = pct)) +
    geom_tile(color = "white", linewidth = tile_lwd) +
    geom_text(aes(label = label),
              size  = label_size,
              color = ifelse(mat_df$pct > 60, "white", "grey20")) +
    scale_fill_gradient(
      low    = "#f7fbff",
      high   = "#2171b5",
      name   = "% of cluster",
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100)
    ) +
    scale_x_discrete(position = "bottom") +
    theme_bw(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5,
                                   size = base_size * 1.2),
      plot.subtitle = element_text(hjust = 0.5, size = base_size * 0.9,
                                   color = "grey40"),
      axis.text.x   = element_text(angle = 45, hjust = 1,
                                   size  = xtick_size),
      axis.text.y   = element_text(size  = ytick_size),
      axis.title    = element_blank(),
      legend.title  = element_text(size = base_size * 0.9, face = "bold"),
      legend.text   = element_text(size = base_size * 0.85),
      legend.key.height = unit(phys_h * 0.15, "inches"),
      panel.grid    = element_blank()
    ) +
    labs(title = title, subtitle = subtitle)
}


# ============================================================
# LOOP: ddqc resolutions
# ============================================================

message("=== ddqc cross-tabulation heatmaps ===")

for (res in resolutions) {

  message(paste0("\n--- Resolution: ", res, " ---"))

  path <- file.path(comp_dir,
                    paste0("ddqc_classification_mad_results_res_", res, ".csv"))
  if (!file.exists(path)) {
    message("  Not found — skipping")
    next
  }

  df <- read.csv(path, row.names = 1) %>% bool_fix()

  # Rename cluster column if needed
  if ("cluster" %in% colnames(df) && !"ddqc_cluster" %in% colnames(df))
    df <- df %>% rename(ddqc_cluster = cluster)

  n_retained <- sum(df$keep_ddqc == TRUE, na.rm = TRUE)
  res_out    <- file.path(out_dir, paste0("res_", res))
  dir.create(res_out, showWarnings = FALSE)

  for (annot in c("cell_typist", "single_r")) {

    annot_label <- if (annot == "cell_typist") "CellTypist" else "SingleR"

    if (!annot %in% colnames(df)) {
      message("  Missing column: ", annot, " — skipping")
      next
    }

    mat <- build_crosstab(df,
                          cluster_col = "ddqc_cluster",
                          annot_col   = annot,
                          keep_col    = "keep_ddqc")

    if (is.null(mat)) next

    p <- plot_crosstab(
      mat_df  = mat,
      title   = paste0("ddqc (res=", res, ") clusters vs ", annot_label),
      n_cells = n_retained
    )

    n_clusters  <- length(unique(mat$cluster))
    n_celltypes <- length(unique(mat$cell_type))
    h <- max(3, n_clusters  * 0.22)
    w <- max(5, n_celltypes * 0.28)

    ggsave(
      file.path(res_out,
                paste0("crosstab_ddqc_res", res, "_vs_", annot, ".png")),
      p, width = w, height = h, dpi = 450, limitsize = FALSE
    )

    # Also save raw CSV
    write.csv(
      mat,
      file.path(res_out,
                paste0("crosstab_ddqc_res", res, "_vs_", annot, ".csv")),
      row.names = FALSE
    )

    message("  Saved: ", annot_label, " heatmap")
  }
}


# ============================================================
# STANDARD CUTOFF
# ============================================================

message("\n=== Standard cutoff cross-tabulation heatmap ===")

std_path <- file.path(comp_dir, "std_classification_mad_results.csv")

if (file.exists(std_path)) {

  df <- read.csv(std_path, row.names = 1) %>% bool_fix()

  if ("cluster" %in% colnames(df) && !"std_cluster" %in% colnames(df))
    df <- df %>% rename(std_cluster = cluster)

  n_retained <- sum(!is.na(df$keep_std) & df$keep_std == TRUE)
  std_out    <- file.path(out_dir, "standard_cutoff")
  dir.create(std_out, showWarnings = FALSE)

  for (annot in c("cell_typist", "single_r")) {

    annot_label <- if (annot == "cell_typist") "CellTypist" else "SingleR"

    if (!annot %in% colnames(df)) {
      message("  Missing column: ", annot, " — skipping")
      next
    }

    mat <- build_crosstab(df,
                          cluster_col = "std_cluster",
                          annot_col   = annot,
                          keep_col    = "keep_std")

    if (is.null(mat)) next

    p <- plot_crosstab(
      mat_df  = mat,
      title   = paste0("Standard cutoff clusters vs ", annot_label),
      n_cells = n_retained
    )

    n_clusters  <- length(unique(mat$cluster))
    n_celltypes <- length(unique(mat$cell_type))
    h <- max(3, n_clusters  * 0.22)
    w <- max(5, n_celltypes * 0.28)

    ggsave(
      file.path(std_out,
                paste0("crosstab_std_vs_", annot, ".png")),
      p, width = w, height = h, dpi = 450, limitsize = FALSE
    )

    write.csv(
      mat,
      file.path(std_out,
                paste0("crosstab_std_vs_", annot, ".csv")),
      row.names = FALSE
    )

    message("  Saved: Standard cutoff vs ", annot_label)
  }

} else {
  message("Standard cutoff MAD results not found — skipping")
}


# ============================================================
# COMBINED SUMMARY: cluster purity across all conditions
# ============================================================

message("\n=== Computing cluster purity summary ===")

# For each condition, compute dominant cell type % per cluster
# (= the max % in each row of the heatmap)
# High purity = clusters are biologically coherent

purity_rows <- list()

compute_purity <- function(mat_df, condition_label) {
  mat_df %>%
    group_by(cluster) %>%
    summarise(
      dominant_celltype = cell_type[which.max(pct)],
      purity_pct        = max(pct),
      .groups           = "drop"
    ) %>%
    mutate(condition = condition_label)
}

# ddqc resolutions
for (res in resolutions) {
  for (annot in c("cell_typist", "single_r")) {
    csv_path <- file.path(out_dir, paste0("res_", res),
                          paste0("crosstab_ddqc_res", res, "_vs_", annot, ".csv"))
    if (!file.exists(csv_path)) next
    mat <- read.csv(csv_path)
    annot_label <- if (annot == "cell_typist") "CellTypist" else "SingleR"
    purity_rows[[paste(res, annot)]] <- compute_purity(
      mat, paste0("ddqc res=", res, " (", annot_label, ")")
    )
  }
}

# Standard cutoff
for (annot in c("cell_typist", "single_r")) {
  csv_path <- file.path(out_dir, "standard_cutoff",
                        paste0("crosstab_std_vs_", annot, ".csv"))
  if (!file.exists(csv_path)) next
  mat <- read.csv(csv_path)
  annot_label <- if (annot == "cell_typist") "CellTypist" else "SingleR"
  purity_rows[[paste("std", annot)]] <- compute_purity(
    mat, paste0("Standard cutoff (", annot_label, ")")
  )
}

if (length(purity_rows) > 0) {

  purity_df <- bind_rows(purity_rows)

  write.csv(purity_df,
            file.path(out_dir, "cluster_purity_summary.csv"),
            row.names = FALSE)

  # Summary: mean purity per condition
  purity_summary <- purity_df %>%
    group_by(condition) %>%
    summarise(
      n_clusters      = n(),
      mean_purity     = round(mean(purity_pct), 1),
      median_purity   = round(median(purity_pct), 1),
      pct_above_80    = round(100 * mean(purity_pct >= 80), 1),
      .groups         = "drop"
    )

  message("\nCluster purity summary:")
  print(purity_summary)

  write.csv(purity_summary,
            file.path(out_dir, "cluster_purity_by_condition.csv"),
            row.names = FALSE)

  # Purity distribution plot
  purity_df$condition <- factor(purity_df$condition,
                                 levels = unique(purity_df$condition))

  p_purity <- ggplot(purity_df,
                     aes(x = condition, y = purity_pct, fill = condition)) +
    geom_violin(alpha = 0.7, color = NA) +
    geom_boxplot(width = 0.15, fill = "white",
                 outlier.size = 1, color = "grey30") +
    geom_hline(yintercept = 80, linetype = "dashed",
               color = "red", linewidth = 0.7) +
    annotate("text", x = 0.6, y = 82,
             label = "80% purity threshold",
             color = "red", size = 3.5, hjust = 0) +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(limits = c(0, 105),
                       labels = function(x) paste0(x, "%")) +
    theme_bw(base_size = 12) +
    theme(
      plot.title   = element_text(face = "bold", hjust = 0.5, size = 13),
      axis.text.x  = element_text(angle = 35, hjust = 1, size = 10),
      axis.title.x = element_blank(),
      legend.position = "none"
    ) +
    labs(
      title = "Cluster Purity Distribution by Condition",
      y     = "% Dominant Cell Type per Cluster"
    )

  ggsave(file.path(out_dir, "cluster_purity_distribution.png"),
         p_purity,
         width  = max(6, length(unique(purity_df$condition)) * 0.9),
         height = 4, dpi = 450)

  message("Saved: cluster_purity_distribution.png")
}


# ---- SESSION INFO -------------------------------------------------------------

sink(file.path(out_dir, "session_info.txt"))
sessionInfo()
sink()

message("\n=== Done. Outputs in: ", out_dir, " ===")
message("Structure:")
message("  crosstab_heatmap/")
message("    res_0.4/ ... res_2.0/   — one heatmap per resolution per annotation")
message("    standard_cutoff/        — one heatmap per annotation")
message("    cluster_purity_summary.csv")
message("    cluster_purity_by_condition.csv")
message("    cluster_purity_distribution.png")