## ============================================================
## Annotation Comparison Visualization
## Adjusted from authors' celltypist_vs_ddqc_umap_plot.R
##
## Parts:
##   A — ddqc (loop over resolutions)
##   B — miQC (no cluster column, Venn + status UMAP)
##   C — Standard cutoff
## ============================================================
install.packages("ggvenn", repos = "https://cran.r-project.org")
suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(dplyr)
  library(shadowtext)
  library(ggvenn)
  library(Seurat)
})


# ---- PATHS --------------------------------------------------------------------

BASE_DIR       <- "/projectnb/ds596/projects/Team 9/scQC_project"
annot_comp_dir <- file.path(BASE_DIR, "compare_annot")
comp_dir       <- annot_comp_dir
miqc_dir       <- file.path(BASE_DIR, "recreation/miQC")
std_dir        <- file.path(BASE_DIR, "recreation/standard_cutoff")
out_dir        <- file.path(BASE_DIR, "compare_annot/plots")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

resolutions <- c(0.4, 0.8, 1.3, 1.6, 2)


# ---- THEMES (unchanged from authors) -----------------------------------------

theme_umap <- theme(
  axis.text.x  = element_text(size = 15),
  axis.title.x = element_text(size = 15),
  axis.text.y  = element_text(size = 15),
  axis.title.y = element_text(size = 15),
  plot.title   = element_text(size = 20, face = "bold"),
  legend.title = element_text(size = 15),
  legend.text  = element_text(size = 10)
)

theme_horizontal_with_legend <- theme(
  axis.text.x  = element_text(angle = 45, size = 15, hjust = 1, face = "bold"),
  axis.text.y  = element_text(size = 15),
  axis.title.y = element_text(size = 15),
  axis.title.x = element_blank(),
  plot.title   = element_text(size = 20, face = "bold")
)

no_bkg <- theme(
  axis.line        = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border     = element_blank(),
  panel.background = element_blank()
)


# ---- PLOT FUNCTIONS (unchanged from authors) ----------------------------------

ggsave1 <- function(filename, plot, n.clusters = 30, type = "h", width.multiplier = 1) {
  if (type == "h") { height <- 10; width <- 14 / 30 * max(n.clusters, 30) }
  if (type == "v") { height <- 14 / 30 * max(n.clusters, 30); width <- 14 }
  if (type == "u") { height <- 10; width <- 10 + 2 * ceiling(n.clusters / 13) }
  ggsave(filename = filename, plot = plot,
         width = width * width.multiplier, height = height, dpi = 300)
}

DimPlotCluster <- function(obj, lbls, ttl) {
  data <- data.frame(
    UMAP1         = obj$umap1,
    UMAP2         = obj$umap2,
    cluster       = factor(obj$cluster_labels),
    other_cluster = factor(obj$other_cluster)
  )
  plot <- ggplot(data, aes(x = UMAP1, y = UMAP2, color = other_cluster)) +
    geom_point(size = 0.25) + theme_umap + ttl +
    labs(color = "Cell Type") +
    theme(
      legend.text     = element_text(size = 13),
      legend.title    = element_text(size = 14, face = "bold"),
      legend.key.size = unit(0.5, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))

  for (cl in levels(factor(obj$cluster_labels))) {
    cluster.data <- subset(data, cluster == cl)
    plot <- plot +
      geom_shadowtext(
        data = data.frame(
          x     = mean(cluster.data$UMAP1),
          y     = mean(cluster.data$UMAP2),
          label = cl
        ),
        aes(x = x, y = y, label = label),
        colour    = "black",
        bg.colour = "white",
        bg.r      = 0.18,
        size      = 7 / .pt,
        fontface  = "plain",
        inherit.aes = FALSE
      )
  }
  return(plot)
}

makeBarplot <- function(obj, lbls, ttl) {
  data <- data.frame(
    cluster       = obj$cluster_labels,
    other_cluster = obj$other_cluster
  )
  table.other_cluster <- NULL
  table.cluster       <- NULL
  table.freq          <- NULL

  for (cl in levels(as.factor(data$cluster))) {
    data.cluster <- subset(data, cluster == cl)
    table.tmp    <- as.data.frame(table(data.frame(
      other_cluster = factor(data.cluster$other_cluster),
      cluster       = as.character(data.cluster$cluster)
    )))
    table.tmp$Freq      <- table.tmp$Freq / sum(table.tmp$Freq)
    table.other_cluster <- c(table.other_cluster, as.character(table.tmp$other_cluster))
    table.cluster       <- c(table.cluster, as.integer(as.character(table.tmp$cluster)))
    table.freq          <- c(table.freq, as.character(table.tmp$Freq))
  }

  data1 <- data.frame(
    other_cluster = table.other_cluster,
    cluster       = as.factor(table.cluster),
    freq          = as.double(as.character(table.freq)) * 100
  )
  ggplot(data1, aes(x = cluster, y = freq, fill = other_cluster)) +
    geom_bar(stat = "identity") +
    theme_horizontal_with_legend + ttl +
    scale_x_discrete(labels = lbls)
}

makePlots <- function(obj, results.dir, plot.title, width.multiplier = 1) {
  message(paste("Making plots:", plot.title))
  lbls        <- seq_along(unique(obj$cluster_labels))
  names(lbls) <- seq_along(lbls)
  n.clusters  <- length(unique(obj$cluster_labels))
  ttl         <- ggtitle(plot.title)

  ggsave1(
    filename         = file.path(results.dir, paste0(plot.title, "_umap_clusters.png")),
    plot             = DimPlotCluster(obj, lbls, ttl),
    n.clusters       = n.clusters,
    type             = "u",
    width.multiplier = width.multiplier
  )
  ggsave1(
    filename   = file.path(results.dir, paste0(plot.title, "_barplot.png")),
    plot       = makeBarplot(obj, lbls, ttl),
    n.clusters = n.clusters,
    type       = "h"
  )
}


# ---- HELPERS -----------------------------------------------------------------

bool_fix <- function(cells) {
  passed_cols <- grep("_passed_qc$|^keep_", colnames(cells), value = TRUE)
  for (col in passed_cols) {
    if (is.character(cells[[col]])) cells[[col]] <- cells[[col]] == "True"
  }
  cells
}

get_umap_coords_from_rds <- function(seurat_path) {
  seu        <- readRDS(seurat_path)
  umap_embed <- Embeddings(seu, reduction = "umap")
  data.frame(
    cell  = rownames(umap_embed),
    umap1 = umap_embed[, 1],
    umap2 = umap_embed[, 2],
    stringsAsFactors = FALSE
  )
}

make_status_umap <- function(cells_umap, keep_col, compare_col,
                             compare_label, color_compare,
                             out_path) {
  status_col <- paste0("status_", compare_col)
  only_label <- paste(compare_label, "only")

  cells_umap[[status_col]] <- case_when(
    cells_umap[[keep_col]]    & cells_umap[[compare_col]]  ~ "Both retain",
    cells_umap[[keep_col]]    & !cells_umap[[compare_col]] ~ paste(gsub("_passed_qc", "", keep_col), "only"),
    !cells_umap[[keep_col]]   & cells_umap[[compare_col]]  ~ only_label,
    TRUE ~ "Both drop"
  )

  counts <- table(cells_umap[[status_col]])
  label_text <- paste(names(counts), counts, sep = ": ", collapse = "\n")

  p <- ggplot(cells_umap %>% filter(!is.na(umap1)),
              aes(x = umap1, y = umap2, color = .data[[status_col]])) +
    geom_point(size = 0.3, alpha = 0.7) +
    scale_color_manual(values = c(
      "Both retain" = "#59a14f",
      "Both drop"   = "#bab0ac"
    ) %>% c(
      setNames(color_compare,  only_label),
      setNames("#4e79a7", paste(gsub("_passed_qc", "", keep_col), "only"))
    )) +
    theme_bw(base_size = 14) +
    theme(
      plot.title      = element_text(face = "bold", hjust = 0.5),
      legend.title    = element_text(size = 13, face = "bold"),
      legend.text     = element_text(size = 12),
      legend.key.size = unit(0.5, "cm")
    ) +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    labs(
      title = paste(gsub("_passed_qc", "", keep_col), "vs", compare_label, "— Retention Status"),
      x = "UMAP1", y = "UMAP2", color = "Status"
    ) +
    annotate("text", x = Inf, y = Inf, label = label_text,
             hjust = 1.05, vjust = 1.5, size = 4, color = "black")

  ggsave(out_path, p, width = 10, height = 7, dpi = 300)
  message("Saved: ", out_path)
}


# ============================================================
# PART A — ddqc: loop over resolutions
# ============================================================

message("=== PART A: ddqc plots ===")

for (res in resolutions) {

  message(paste("\n--- Resolution:", res, "---"))

  mad_path <- file.path(annot_comp_dir,
                        paste0("ddqc_classification_mad_results_res_", res, ".csv"))
  if (!file.exists(mad_path)) {
    message(paste("  Not found:", mad_path, "— skipping"))
    next
  }

  cells <- read.csv(mad_path, row.names = 1) %>% bool_fix()

  if ("cluster" %in% colnames(cells) && !"ddqc_cluster" %in% colnames(cells))
    cells <- cells %>% rename(ddqc_cluster = cluster)

  if (!all(c("umap1", "umap2") %in% colnames(cells))) {
    message("  umap1/umap2 not found — skipping")
    next
  }

  cells$ddqc_cluster <- as.factor(cells$ddqc_cluster)
  cells$cell_typist  <- as.factor(cells$cell_typist)
  if ("single_r" %in% colnames(cells))
    cells$single_r <- as.factor(cells$single_r)

  res_out_dir <- file.path(out_dir, paste0("res_", res))
  dir.create(res_out_dir, showWarnings = FALSE)

  # -- CellTypist comparison ---------------------------------------------------
  if ("cell_typist_passed_qc" %in% colnames(cells)) {
    cells.ddqc <- cells[cells$ddqc_cluster_passed_qc == TRUE, ]
    cells.ddqc[["cluster_labels"]] <- as.factor(cells.ddqc$ddqc_cluster)
    cells.ddqc[["other_cluster"]]  <- as.factor(cells.ddqc$cell_typist)

    cells.ct <- cells[cells$cell_typist_passed_qc == TRUE, ]
    cells.ct[["cluster_labels"]] <- as.factor(cells.ct$cell_typist)
    cells.ct[["other_cluster"]]  <- as.factor(cells.ct$ddqc_cluster)

    makePlots(cells.ddqc, res_out_dir,
              plot.title = paste0("res", res, "_ddqc_colored_by_celltypist"))
    makePlots(cells.ct,   res_out_dir,
              plot.title = paste0("res", res, "_celltypist_colored_by_ddqc"))
  }

  # -- SingleR comparison ------------------------------------------------------
  if ("single_r_passed_qc" %in% colnames(cells)) {
    cells.ddqc.sr <- cells[cells$ddqc_cluster_passed_qc == TRUE, ]
    cells.ddqc.sr[["cluster_labels"]] <- as.factor(cells.ddqc.sr$ddqc_cluster)
    cells.ddqc.sr[["other_cluster"]]  <- as.factor(cells.ddqc.sr$single_r)

    cells.sr <- cells[cells$single_r_passed_qc == TRUE, ]
    cells.sr[["cluster_labels"]] <- as.factor(cells.sr$single_r)
    cells.sr[["other_cluster"]]  <- as.factor(cells.sr$ddqc_cluster)

    makePlots(cells.ddqc.sr, res_out_dir,
              plot.title = paste0("res", res, "_ddqc_colored_by_singler"))
    makePlots(cells.sr,      res_out_dir,
              plot.title = paste0("res", res, "_singler_colored_by_ddqc"))
  }
}


# ============================================================
# PART B — miQC
# ============================================================

message("\n=== PART B: miQC plots ===")

miqc_mad_path <- file.path(comp_dir, "miqc_classification_mad_results.csv")
miqc_rds_path <- file.path(miqc_dir, "miQC_filtered.rds")

if (file.exists(miqc_mad_path) && file.exists(miqc_rds_path)) {

  cells <- read.csv(miqc_mad_path, row.names = 1) %>% bool_fix()

  umap_coords           <- get_umap_coords_from_rds(miqc_rds_path)
  rownames(umap_coords) <- umap_coords$cell
  cells[["umap1"]]      <- umap_coords[rownames(cells), "umap1"]
  cells[["umap2"]]      <- umap_coords[rownames(cells), "umap2"]

  cells$cell_typist <- as.factor(cells$cell_typist)
  if ("single_r" %in% colnames(cells))
    cells$single_r <- as.factor(cells$single_r)

  miqc_out_dir <- file.path(out_dir, "miQC")
  dir.create(miqc_out_dir, showWarnings = FALSE)

  # -- Venn diagram ------------------------------------------------------------
  venn_cols <- intersect(
    c("keep_miqc", "cell_typist_passed_qc", "single_r_passed_qc"),
    colnames(cells)
  )

  if (length(venn_cols) >= 2) {
    venn_data <- lapply(
      setNames(venn_cols, gsub("_passed_qc", "", venn_cols)),
      function(col) rownames(cells)[cells[[col]] == TRUE]
    )
    p_venn <- ggvenn(venn_data,
                     fill_color = c("#4e79a7", "#f28e2b", "#e15759")[seq_along(venn_data)],
                     stroke_size = 0.5, set_name_size = 5)
    ggsave(file.path(miqc_out_dir, "miqc_venn.png"), p_venn,
           width = 8, height = 6, dpi = 300)
    message("Saved: miqc_venn.png")
  }

  cells_umap <- cells %>% filter(!is.na(umap1) & !is.na(umap2))

  # -- CellTypist status UMAP --------------------------------------------------
  if ("cell_typist_passed_qc" %in% colnames(cells_umap)) {
    make_status_umap(
      cells_umap    = cells_umap,
      keep_col      = "keep_miqc",
      compare_col   = "cell_typist_passed_qc",
      compare_label = "CellTypist MAD",
      color_compare = "#f28e2b",
      out_path      = file.path(miqc_out_dir, "miqc_vs_celltypist_retention_status.png")
    )
  }

  # -- SingleR status UMAP -----------------------------------------------------
  if ("single_r_passed_qc" %in% colnames(cells_umap)) {
    make_status_umap(
      cells_umap    = cells_umap,
      keep_col      = "keep_miqc",
      compare_col   = "single_r_passed_qc",
      compare_label = "SingleR MAD",
      color_compare = "#e15759",
      out_path      = file.path(miqc_out_dir, "miqc_vs_singler_retention_status.png")
    )
  }

  # -- CellTypist makePlots ----------------------------------------------------
  if ("cell_typist_passed_qc" %in% colnames(cells)) {
    cells.miqc <- cells[cells$keep_miqc == TRUE & !is.na(cells$umap1), ]
    cells.miqc[["cluster_labels"]] <- as.factor(cells.miqc$cell_typist)
    cells.miqc[["other_cluster"]]  <- as.factor(cells.miqc$cell_typist)

    cells.ct <- cells[cells$cell_typist_passed_qc == TRUE & !is.na(cells$umap1), ]
    cells.ct[["cluster_labels"]] <- as.factor(cells.ct$cell_typist)
    cells.ct[["other_cluster"]]  <- as.factor(cells.ct$cell_typist)

    makePlots(cells.miqc, miqc_out_dir, plot.title = "miqc_colored_by_celltypist")
    makePlots(cells.ct,   miqc_out_dir, plot.title = "celltypist_mad_filtered")
  }

  # -- SingleR makePlots -------------------------------------------------------
  if ("single_r_passed_qc" %in% colnames(cells)) {
    cells.miqc.sr <- cells[cells$keep_miqc == TRUE & !is.na(cells$umap1), ]
    cells.miqc.sr[["cluster_labels"]] <- as.factor(cells.miqc.sr$single_r)
    cells.miqc.sr[["other_cluster"]]  <- as.factor(cells.miqc.sr$single_r)

    cells.sr <- cells[cells$single_r_passed_qc == TRUE & !is.na(cells$umap1), ]
    cells.sr[["cluster_labels"]] <- as.factor(cells.sr$single_r)
    cells.sr[["other_cluster"]]  <- as.factor(cells.sr$single_r)

    makePlots(cells.miqc.sr, miqc_out_dir, plot.title = "miqc_colored_by_singler")
    makePlots(cells.sr,      miqc_out_dir, plot.title = "singler_mad_filtered")
  }

} else {
  message("miQC files not found — skipping")
}


# ============================================================
# PART C — Standard Cutoff
# ============================================================

message("\n=== PART C: Standard cutoff plots ===")

std_mad_path <- file.path(comp_dir, "std_classification_mad_results.csv")

if (file.exists(std_mad_path)) {

  cells <- read.csv(std_mad_path, row.names = 1) %>% bool_fix()

  if ("cluster" %in% colnames(cells) && !"std_cluster" %in% colnames(cells))
    cells <- cells %>% rename(std_cluster = cluster)

  cells <- cells %>% filter(!is.na(std_cluster) & !is.na(umap1))

  cells$std_cluster  <- as.factor(cells$std_cluster)
  cells$cell_typist  <- as.factor(cells$cell_typist)
  if ("single_r" %in% colnames(cells))
    cells$single_r <- as.factor(cells$single_r)

  std_out_dir <- file.path(out_dir, "standard_cutoff")
  dir.create(std_out_dir, showWarnings = FALSE)

  # -- CellTypist comparison ---------------------------------------------------
  if ("cell_typist_passed_qc" %in% colnames(cells)) {
    cells.std <- cells[cells$keep_std == TRUE, ]
    cells.std[["cluster_labels"]] <- as.factor(cells.std$std_cluster)
    cells.std[["other_cluster"]]  <- as.factor(cells.std$cell_typist)

    cells.ct <- cells[cells$cell_typist_passed_qc == TRUE, ]
    cells.ct[["cluster_labels"]] <- as.factor(cells.ct$cell_typist)
    cells.ct[["other_cluster"]]  <- as.factor(cells.ct$std_cluster)

    makePlots(cells.std, std_out_dir, plot.title = "std_colored_by_celltypist")
    makePlots(cells.ct,  std_out_dir, plot.title = "celltypist_colored_by_std")
  }

  # -- SingleR comparison ------------------------------------------------------
  if ("single_r_passed_qc" %in% colnames(cells)) {
    cells.std.sr <- cells[cells$keep_std == TRUE, ]
    cells.std.sr[["cluster_labels"]] <- as.factor(cells.std.sr$std_cluster)
    cells.std.sr[["other_cluster"]]  <- as.factor(cells.std.sr$single_r)

    cells.sr <- cells[cells$single_r_passed_qc == TRUE, ]
    cells.sr[["cluster_labels"]] <- as.factor(cells.sr$single_r)
    cells.sr[["other_cluster"]]  <- as.factor(cells.sr$std_cluster)

    makePlots(cells.std.sr, std_out_dir, plot.title = "std_colored_by_singler")
    makePlots(cells.sr,     std_out_dir, plot.title = "singler_colored_by_std")
  }

} else {
  message("Standard cutoff MAD results not found — skipping")
}


# ---- SESSION INFO -------------------------------------------------------------

sink(file.path(out_dir, "session_info.txt"))
sessionInfo()
sink()

message("\n=== Done. Plots saved in: ", out_dir, " ===")