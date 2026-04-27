## ============================================================
## Cell Type Composition Comparison Across QC Methods
##
## Shows whether ddqc, miQC, and standard cutoff retain
## the same biological cell type proportions — proving
## that QC method choice does not introduce biological bias
##
## Inputs:
##   - cells_initial_res_{res}.csv  (ddqc, one per resolution)
##   - recreation/miQC/cells_initial.csv
##   - recreation/standard_cutoff/cells_initial.csv
##   - cell_classification/colon_immune_celltypist.csv
##   - cell_classification/SingleR_colon.csv
##
## Outputs:
##   - composition_barplot_celltypist.png
##   - composition_barplot_singler.png
##   - composition_heatmap_celltypist.png
##   - composition_heatmap_singler.png
##   - cell_type_composition_table.csv
##   - retention_rate_by_celltype.png  (per cell type, per method)
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})


# ---- PATHS --------------------------------------------------------------------

BASE_DIR    <- "/projectnb/ds596/projects/Team 9/scQC_project"
ddqc_dir    <- file.path(BASE_DIR, "ddqc_annot")
miqc_dir    <- file.path(BASE_DIR, "recreation/miQC")
std_dir     <- file.path(BASE_DIR, "recreation/standard_cutoff")
annot_dir   <- file.path(BASE_DIR, "cell_annot")
out_dir     <- file.path(BASE_DIR, "compare_annot/composition")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

celltypist_csv <- file.path(annot_dir, "colon_immune_celltypist.csv")
singler_csv    <- file.path(annot_dir, "SingleR_colon_immuneRef.csv")

# ddqc resolution to use for main comparison (use 1.3 as default — paper default)
# All resolutions will also be compared in a separate plot
focal_resolution <- 1.3
resolutions      <- c(0.4, 0.8, 1.3, 1.6, 2)


# ---- 1. LOAD ANNOTATION LABELS -----------------------------------------------

message("Loading annotation labels...")

ct_df <- read.csv(celltypist_csv, row.names = 1)
# Find the majority_voting or predicted_labels column
ct_col <- if ("majority_voting" %in% colnames(ct_df)) "majority_voting" else colnames(ct_df)[1]
ct_df  <- data.frame(cell = rownames(ct_df), cell_typist = ct_df[[ct_col]],
                     stringsAsFactors = FALSE)

sr_df <- read.csv(singler_csv, row.names = 1)
sr_df <- data.frame(cell = rownames(sr_df), single_r = sr_df[[1]],
                    stringsAsFactors = FALSE)

message(sprintf("CellTypist: %d cells, %d types",
                nrow(ct_df), length(unique(ct_df$cell_typist))))
message(sprintf("SingleR: %d cells, %d types",
                nrow(sr_df), length(unique(sr_df$single_r))))


# ---- 2. HELPER: LOAD CELLS + JOIN ANNOTATIONS --------------------------------

load_and_join <- function(cells_path, keep_col) {
  if (!file.exists(cells_path)) {
    message("  Not found: ", cells_path)
    return(NULL)
  }
  cells <- read.csv(cells_path, stringsAsFactors = FALSE)

  # Convert boolean if needed
  if (is.character(cells[[keep_col]])) {
    cells[[keep_col]] <- cells[[keep_col]] == "True"
  }

  # Join annotations
  cells <- cells %>%
    left_join(ct_df, by = "cell") %>%
    left_join(sr_df, by = "cell")

  cells
}


# ---- 3. COMPUTE CELL TYPE COMPOSITION PER METHOD -----------------------------

#' For a given cells dataframe, compute proportion of each cell type
#' among retained cells
#'
#' @param cells     per-cell dataframe with keep_col and annotation columns
#' @param keep_col  column indicating which cells are retained
#' @param annot_col annotation column (cell_typist or single_r)
#' @param method    label for this method (e.g. "ddqc_res1.3")

get_composition <- function(cells, keep_col, annot_col, method) {
  retained <- cells %>%
    filter(.data[[keep_col]] == TRUE, !is.na(.data[[annot_col]]))

  total_retained <- nrow(retained)
  if (total_retained == 0) return(NULL)

  retained %>%
    count(.data[[annot_col]], name = "n") %>%
    rename(cell_type = 1) %>%
    mutate(
      proportion = n / total_retained,
      pct        = proportion * 100,
      method     = method,
      n_total    = total_retained
    )
}

#' Compute retention RATE per cell type:
#' what fraction of each cell type survives this QC method

get_retention_rate <- function(cells, keep_col, annot_col, method) {
  cells %>%
    filter(!is.na(.data[[annot_col]])) %>%
    group_by(cell_type = .data[[annot_col]]) %>%
    summarise(
      n_total    = n(),
      n_retained = sum(.data[[keep_col]] == TRUE, na.rm = TRUE),
      retention_rate = n_retained / n_total,
      .groups = "drop"
    ) %>%
    mutate(method = method)
}


# ---- 4. BUILD COMPOSITION TABLE ACROSS METHODS --------------------------------

message("Computing cell type compositions...")

composition_list    <- list()
retention_rate_list <- list()

for (annot_col in c("cell_typist", "single_r")) {

  # -- ddqc: focal resolution --------------------------------------------------
  cells_path <- file.path(ddqc_dir,
                          paste0("cells_initial_res_", focal_resolution, ".csv"))
  cells <- load_and_join(cells_path, "keep_ddqc")
  if (!is.null(cells)) {
    method_label <- paste0("ddqc_res", focal_resolution)
    composition_list[[paste(method_label, annot_col)]] <-
      get_composition(cells, "keep_ddqc", annot_col, method_label)
    retention_rate_list[[paste(method_label, annot_col)]] <-
      get_retention_rate(cells, "keep_ddqc", annot_col, method_label)
  }

  # -- ddqc: all resolutions ---------------------------------------------------
  for (res in resolutions) {
    cells_path <- file.path(ddqc_dir,
                            paste0("cells_initial_res_", res, ".csv"))
    cells <- load_and_join(cells_path, "keep_ddqc")
    if (!is.null(cells)) {
      method_label <- paste0("ddqc_res", res)
      composition_list[[paste(method_label, annot_col)]] <-
        get_composition(cells, "keep_ddqc", annot_col, method_label)
    }
  }

  # -- miQC --------------------------------------------------------------------
  cells_path <- file.path(miqc_dir, "cells_initial_miQC.csv")
  cells <- load_and_join(cells_path, "keep_miqc")
  if (!is.null(cells)) {
    composition_list[[paste("miQC", annot_col)]] <-
      get_composition(cells, "keep_miqc", annot_col, "miQC")
    retention_rate_list[[paste("miQC", annot_col)]] <-
      get_retention_rate(cells, "keep_miqc", annot_col, "miQC")
  }

  # -- Standard cutoff ---------------------------------------------------------
  cells_path <- file.path(std_dir, "stdCutoff_cells_initial.csv")
  cells <- load_and_join(cells_path, "keep_std")
  if (!is.null(cells)) {
    composition_list[[paste("standard_cutoff", annot_col)]] <-
      get_composition(cells, "keep_std", annot_col, "standard_cutoff")
    retention_rate_list[[paste("standard_cutoff", annot_col)]] <-
      get_retention_rate(cells, "keep_std", annot_col, "standard_cutoff")
  }
}


# ---- 5. PLOTS ----------------------------------------------------------------

theme_comp <- theme_bw(base_size = 13) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y     = element_text(size = 11),
    legend.title    = element_text(size = 12, face = "bold"),
    legend.text     = element_text(size = 10),
    legend.key.size = unit(0.4, "cm"),
    strip.text      = element_text(face = "bold", size = 11)
  )


for (annot_col in c("cell_typist", "single_r")) {

  annot_label <- if (annot_col == "cell_typist") "CellTypist" else "SingleR"

  # -- 5a. Composition barplot: three main methods side by side ----------------

  methods_main <- c(paste0("ddqc_res", focal_resolution), "miQC", "standard_cutoff")
  keys_main    <- paste(methods_main, annot_col)

  comp_main <- bind_rows(composition_list[keys_main]) %>%
    filter(!is.na(cell_type))

  if (nrow(comp_main) > 0) {

    # Order cell types by mean proportion across methods
    ct_order <- comp_main %>%
      group_by(cell_type) %>%
      summarise(mean_pct = mean(pct), .groups = "drop") %>%
      arrange(desc(mean_pct)) %>%
      pull(cell_type)

    comp_main$cell_type <- factor(comp_main$cell_type, levels = ct_order)
    comp_main$method    <- factor(comp_main$method,
                                  levels = methods_main,
                                  labels = c(paste0("ddqc (res=", focal_resolution, ")"),
                                             "miQC", "Standard Cutoff"))

    p_bar <- ggplot(comp_main, aes(x = cell_type, y = pct, fill = method)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8),
               width = 0.75) +
      scale_fill_manual(values = c("#4e79a7", "#f28e2b", "#e15759")) +
      theme_comp +
      labs(
        title = paste0("Cell Type Composition by QC Method (", annot_label, ")"),
        x     = "Cell Type",
        y     = "% of Retained Cells",
        fill  = "QC Method"
      )

    ggsave(
      file.path(out_dir, paste0("composition_barplot_", annot_col, ".png")),
      p_bar, width = 16, height = 6, dpi = 300
    )
    message("Saved: composition_barplot_", annot_col, ".png")
  }


  # -- 5b. Heatmap: proportion of each cell type per method --------------------

  comp_all <- bind_rows(
    composition_list[grep(annot_col, names(composition_list), value = TRUE)]
  ) %>% filter(!is.na(cell_type))

  if (nrow(comp_all) > 0) {

    comp_wide <- comp_all %>%
      select(cell_type, method, pct) %>%
      pivot_wider(names_from = method, values_from = pct, values_fill = 0)

    comp_long <- comp_wide %>%
      pivot_longer(-cell_type, names_to = "method", values_to = "pct")

    # Order cell types
    ct_order <- comp_wide %>%
      mutate(mean_pct = rowMeans(select(., -cell_type))) %>%
      arrange(desc(mean_pct)) %>%
      pull(cell_type)

    comp_long$cell_type <- factor(comp_long$cell_type, levels = rev(ct_order))

    p_heat <- ggplot(comp_long, aes(x = method, y = cell_type, fill = pct)) +
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(aes(label = sprintf("%.1f", pct)), size = 3) +
      scale_fill_gradient(low = "white", high = "#4e79a7",
                          name = "% Retained") +
      theme_comp +
      theme(axis.text.x = element_text(angle = 35, hjust = 1)) +
      labs(
        title = paste0("Cell Type % per QC Method (", annot_label, ")"),
        x = "QC Method", y = "Cell Type"
      )

    ggsave(
      file.path(out_dir, paste0("composition_heatmap_", annot_col, ".png")),
      p_heat, width = max(12, length(unique(comp_long$method)) * 1.5),
      height = max(8, length(unique(comp_long$cell_type)) * 0.35),
      dpi = 300
    )
    message("Saved: composition_heatmap_", annot_col, ".png")
  }


  # -- 5c. Retention RATE per cell type ----------------------------------------
  # Shows what fraction of each cell type each method retains
  # If all methods have similar retention rates → no biological bias

  ret_main <- bind_rows(
    retention_rate_list[grep(annot_col, names(retention_rate_list), value = TRUE)]
  ) %>% filter(!is.na(cell_type))

  if (nrow(ret_main) > 0) {

    ret_main$method <- factor(ret_main$method,
                              levels = c(paste0("ddqc_res", focal_resolution),
                                         "miQC", "standard_cutoff"),
                              labels = c(paste0("ddqc (res=", focal_resolution, ")"),
                                         "miQC", "Standard Cutoff"))

    # Order by mean retention rate
    ct_order_ret <- ret_main %>%
      group_by(cell_type) %>%
      summarise(mean_ret = mean(retention_rate), .groups = "drop") %>%
      arrange(desc(mean_ret)) %>%
      pull(cell_type)

    ret_main$cell_type <- factor(ret_main$cell_type, levels = ct_order_ret)

    p_ret <- ggplot(ret_main, aes(x = cell_type, y = retention_rate * 100,
                                  color = method, group = method)) +
      geom_line(linewidth = 0.8, alpha = 0.8) +
      geom_point(size = 2.5) +
      scale_color_manual(values = c("#4e79a7", "#f28e2b", "#e15759")) +
      scale_y_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
      theme_comp +
      labs(
        title = paste0("Retention Rate per Cell Type by QC Method (", annot_label, ")"),
        x     = "Cell Type",
        y     = "% of Cell Type Retained",
        color = "QC Method"
      )

    ggsave(
      file.path(out_dir, paste0("retention_rate_by_celltype_", annot_col, ".png")),
      p_ret, width = 16, height = 6, dpi = 300
    )
    message("Saved: retention_rate_by_celltype_", annot_col, ".png")
  }


  # -- 5d. ddqc resolution comparison ------------------------------------------
  # Shows how cell type composition changes across ddqc resolutions

  keys_res <- paste(paste0("ddqc_res", resolutions), annot_col)
  comp_res <- bind_rows(composition_list[keys_res]) %>%
    filter(!is.na(cell_type))

  if (nrow(comp_res) > 0) {

    ct_order_res <- comp_res %>%
      group_by(cell_type) %>%
      summarise(mean_pct = mean(pct), .groups = "drop") %>%
      arrange(desc(mean_pct)) %>%
      pull(cell_type)

    comp_res$cell_type <- factor(comp_res$cell_type, levels = ct_order_res)
    comp_res$resolution <- gsub("ddqc_res", "", comp_res$method)

    p_res <- ggplot(comp_res, aes(x = cell_type, y = pct,
                                  fill = resolution, group = resolution)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.85),
               width = 0.8) +
      scale_fill_brewer(palette = "Blues", direction = 1) +
      theme_comp +
      labs(
        title = paste0("Cell Type Composition by ddqc Resolution (", annot_label, ")"),
        x     = "Cell Type",
        y     = "% of Retained Cells",
        fill  = "Resolution"
      )

    ggsave(
      file.path(out_dir, paste0("composition_ddqc_resolutions_", annot_col, ".png")),
      p_res, width = 16, height = 6, dpi = 300
    )
    message("Saved: composition_ddqc_resolutions_", annot_col, ".png")
  }
}


# ---- 6. SAVE COMPOSITION TABLE -----------------------------------------------

# Full table: all methods x all cell types x both annotations
all_comp <- bind_rows(composition_list) %>%
  filter(!is.na(cell_type)) %>%
  select(method, cell_type, n, pct, n_total)

write.csv(all_comp,
          file.path(out_dir, "cell_type_composition_table.csv"),
          row.names = FALSE)

message("Saved: cell_type_composition_table.csv")


# ---- SESSION INFO -------------------------------------------------------------

sink(file.path(out_dir, "session_info.txt"))
sessionInfo()
sink()

message("\n=== Done. Outputs in: ", out_dir, " ===")
message("Key plots:")
message("  composition_barplot_{annot}.png       — side by side bar, 3 methods")
message("  composition_heatmap_{annot}.png       — heatmap all methods x cell types")
message("  retention_rate_by_celltype_{annot}.png — line plot per cell type")
message("  composition_ddqc_resolutions_{annot}.png — ddqc resolution comparison")