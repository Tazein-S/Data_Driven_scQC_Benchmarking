## ============================================================
## ddqc vs. Independent Cell Annotation Comparison
## Reproduces analysis from Subramanian et al. Genome Biology (2022) 23:267
##
## Based on authors' actual scripts:
##   - single_r.Rmd          → SingleR annotation on PBMC
##   - cell_typist.py        → CellTypist annotation (Python, run separately)
##   - merge_classifications.py → joins annotation + ddqc QC stats by barcode
##   - mad_filtering.py      → MAD filter using each grouping column
##   - celltypist_vs_ddqc_umap_plot.R → UMAP + barplot visualizations
##
## NOTE: CellTypist and merge_classifications steps are Python-based.
##       Run cell_typist.py and merge_classifications.py first to produce
##       summary.csv, then run MAD filtering and visualization here in R.
##
## Workflow:
##   PART A — PBMC (SingleR + Azimuth + CellTypist)
##     1. Load raw PBMC 10x data
##     2. Minimal filter (nFeature_RNA > 10, no standard QC cutoffs)
##     3. Run SingleR using HumanPrimaryCellAtlasData reference
##     4. Save SingleR labels as CSV
##     5. Load Azimuth predictions (pre-run via web interface, loaded as TSV)
##     [CellTypist run separately in Python → cell_typist.py]
##     [merge_classifications.py joins everything into summary.csv]
##     6. MAD filtering using each grouping column
##     7. UMAP + barplot visualization
##
##   PART B — Krasnow Lung (CellTypist only, SingleR/Azimuth skipped per paper)
## ============================================================

# ---- 0. DEPENDENCIES ----------------------------------------------------------
# install.packages(c("Seurat", "ggplot2", "cowplot", "dplyr"))
# BiocManager::install(c("SingleR", "celldex"))

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(ggplot2)
  library(cowplot)
  library(dplyr)
})


# ---- PATHS — adjust to your directory structure -------------------------------

PBMC_10X_DIR             <- "data/misc/pbmc_seurat/"
PBMC_AZIMUTH_TSV         <- "data/misc/pbmc_seurat/azimuth_pred.tsv"  # from Azimuth web interface

# Output of Python scripts (run cell_typist.py + merge_classifications.py first)
PBMC_SUMMARY_CSV         <- "cell_classification/pbmc/summary.csv"
LUNG_SUMMARY_CSV         <- "cell_classification/krasnow_lung/summary.csv"
LUNG_CELLS_INITIAL_CSV   <- "cell_classification/krasnow_lung/!cells_initial.csv"

OUT_DIR <- "ddqc_comparison_output"
dir.create(OUT_DIR, showWarnings = FALSE)


# ============================================================
# PART A — PBMC DATASET
# ============================================================

# ---- 1. LOAD RAW PBMC DATA ----------------------------------------------------
# Authors use min.cells=0, min.features=0 — no filtering at object creation

message("=== Loading PBMC data ===")
data       <- Read10X(data.dir = PBMC_10X_DIR)
seurat.obj <- CreateSeuratObject(
  counts       = data,
  project      = "pbmc3k",
  min.cells    = 0,
  min.features = 0
)


# ---- 2. MINIMAL FILTER --------------------------------------------------------
# Authors only apply nFeature_RNA > 10
# No percent.mt filter, no standard QC thresholds applied here

seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > 10)
message(sprintf("After minimal filter: %d cells", ncol(seurat.obj)))


# ---- 3. SINGLER ANNOTATION ----------------------------------------------------
# Authors pass seurat.obj@assays$RNA@data directly to SingleR
# Reference: HumanPrimaryCellAtlasData with label.main

message("=== Running SingleR ===")
ref.data <- celldex::HumanPrimaryCellAtlasData()

predictions <- SingleR(
  test            = seurat.obj@assays$RNA@data,
  assay.type.test = 1,
  ref             = ref.data,
  labels          = ref.data$label.main
)

message("SingleR label distribution:")
print(table(predictions$labels))


# ---- 4. ADD SINGLER LABELS TO SEURAT OBJECT + SAVE ----------------------------

labels        <- predictions$labels
names(labels) <- rownames(predictions)

seurat.obj <- AddMetaData(
  object   = seurat.obj,
  metadata = predictions$labels,
  col.name = "predictions"
)

write.csv(labels, file.path(OUT_DIR, "pbmc_SingleR.csv"))
saveRDS(seurat.obj, file.path(OUT_DIR, "pbmc_unfiltered.rds"))
message("SingleR labels saved.")


# ---- 5. AZIMUTH ANNOTATIONS ---------------------------------------------------
# Authors ran Azimuth via https://azimuth.hubmapconsortium.org/ (web interface)
# and load the resulting TSV directly — no R-based Azimuth package needed.
# For Krasnow lung, Azimuth was skipped due to web interface problems (paper p.22).

if (file.exists(PBMC_AZIMUTH_TSV)) {
  message("=== Loading Azimuth predictions ===")
  azimuth_pred <- read.delim(PBMC_AZIMUTH_TSV, row.names = 1)
  seurat.obj   <- AddMetaData(object = seurat.obj, metadata = azimuth_pred)
  message("Azimuth predictions added.")
} else {
  message("Azimuth TSV not found — run via https://azimuth.hubmapconsortium.org/ and save output as azimuth_pred.tsv")
}

# NOTE: CellTypist for PBMC is run separately in Python (cell_typist.py)
# using the Immune_All_Low.pkl model. After cell_typist.py and
# merge_classifications.py have produced summary.csv, continue below.


# ============================================================
# MAD FILTERING FUNCTIONS
# Direct R translation of authors' mad_filtering.py
# ============================================================

INF <- 1e10

#' MAD — matches authors' Python implementation (constant = 1.4826)
mad_custom <- function(x, constant = 1.4826) {
  constant * median(abs(x - median(x, na.rm = TRUE)), na.rm = TRUE)
}

#' Per-metric MAD/outlier filter within each group
#' Matches authors' metric_filter() in mad_filtering.py exactly.
#'
#' @param data           data.frame, cells as rows
#' @param classification column name defining groups ("cell_typist" or "ddqc_cluster")
#' @param method         "mad" or "outlier"
#' @param param          MAD multiplier (2 per paper)
#' @param metric_name    QC column to filter on
#' @param do_lower_co    compute lower cutoff?
#' @param do_upper_co    compute upper cutoff?
#' @param lower_bound    floor for lower cutoff (200 for n_genes, INF otherwise)
#' @param upper_bound    ceiling for upper cutoff (10 for percent_mito, -INF otherwise)

metric_filter <- function(data, classification, method, param, metric_name,
                          do_lower_co = FALSE, do_upper_co = FALSE,
                          lower_bound = INF, upper_bound = -INF) {

  data[[paste0(metric_name, "_qc_pass")]]  <- FALSE
  data[[paste0(metric_name, "_lower_co")]] <- NA
  data[[paste0(metric_name, "_upper_co")]] <- NA

  if (method == "mad") {
    data[[paste0(metric_name, "_median")]] <- NA
    data[[paste0(metric_name, "_mad")]]    <- NA
  }

  for (ct in unique(data[[classification]])) {
    idx       <- which(data[[classification]] == ct)
    cell_vals <- data[[metric_name]][idx]

    lower_co <- -INF
    upper_co  <-  INF

    if (method == "mad") {
      med_val <- median(cell_vals, na.rm = TRUE)
      mad_val <- mad_custom(cell_vals)

      data[idx, paste0(metric_name, "_median")] <- med_val
      data[idx, paste0(metric_name, "_mad")]    <- mad_val

      if (do_lower_co) lower_co <- min(med_val - param * mad_val, lower_bound)
      if (do_upper_co) upper_co  <- max(med_val + param * mad_val, upper_bound)
    }

    if (method == "outlier") {
      q75 <- quantile(cell_vals, 0.75, na.rm = TRUE)
      q25 <- quantile(cell_vals, 0.25, na.rm = TRUE)
      iqr <- q75 - q25
      if (do_lower_co) lower_co <- min(q25 - 1.5 * iqr, lower_bound)
      if (do_upper_co) upper_co  <- max(q75 + 1.5 * iqr, upper_bound)
    }

    pass_idx <- idx[
      data[[metric_name]][idx] >= lower_co &
      data[[metric_name]][idx] <= upper_co
    ]
    data[pass_idx, paste0(metric_name, "_qc_pass")] <- TRUE

    if (do_upper_co) data[idx, paste0(metric_name, "_upper_co")] <- upper_co
    if (do_lower_co) data[idx, paste0(metric_name, "_lower_co")] <- lower_co
  }

  data
}


#' Filter cells across all QC metrics, produce cumulative <classification>_passed_qc
#' Matches authors' filter_cells() in mad_filtering.py exactly.
#' Authors call with do_counts=FALSE for the annotation comparison
#' (only n_genes and percent_mito filtered, per paper Methods p.21)

filter_cells <- function(data, classification, method = "mad", threshold = 2,
                         n_genes_lower_bound      = 200,
                         percent_mito_upper_bound = 10,
                         do_counts = TRUE,
                         do_genes  = TRUE,
                         do_mito   = TRUE,
                         do_ribo   = FALSE) {

  data_copy <- data

  if (do_counts) {
    data_copy <- metric_filter(data_copy, classification, method, threshold,
                               "n_counts", do_lower_co = TRUE)
  } else {
    data_copy[["n_counts_qc_pass"]] <- TRUE
  }

  if (do_genes) {
    data_copy <- metric_filter(data_copy, classification, method, threshold,
                               "n_genes", do_lower_co = TRUE,
                               lower_bound = n_genes_lower_bound)
  } else {
    data_copy[["n_genes_qc_pass"]] <- TRUE
  }

  if (do_mito) {
    data_copy <- metric_filter(data_copy, classification, method, threshold,
                               "percent_mito", do_upper_co = TRUE,
                               upper_bound = percent_mito_upper_bound)
  } else {
    data_copy[["percent_mito_qc_pass"]] <- TRUE
  }

  if (do_ribo) {
    data_copy <- metric_filter(data_copy, classification, method, threshold,
                               "percent_ribo", do_upper_co = TRUE)
  } else {
    data_copy[["percent_ribo_qc_pass"]] <- TRUE
  }

  data_copy[["passed_qc"]] <- (
    data_copy[["n_counts_qc_pass"]]     &
    data_copy[["n_genes_qc_pass"]]      &
    data_copy[["percent_mito_qc_pass"]] &
    data_copy[["percent_ribo_qc_pass"]]
  )

  # Transfer cumulative result to original data as <classification>_passed_qc
  data[[paste0(classification, "_passed_qc")]] <- data_copy[["passed_qc"]]
  data
}


# ============================================================
# VISUALIZATION FUNCTIONS
# Direct translation of authors' celltypist_vs_ddqc_umap_plot.R
# ============================================================

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

ggsave1 <- function(filename, plot, n.clusters = 30, type = "h", width.multiplier = 1) {
  if (type == "h") { height <- 10; width <- 14 / 30 * max(n.clusters, 30) }
  if (type == "v") { height <- 14 / 30 * max(n.clusters, 30); width <- 14 }
  if (type == "u") { height <- 10; width <- 10 + 2 * ceiling(n.clusters / 13) }
  ggsave(filename = filename, plot = plot + no_bkg,
         width = width * width.multiplier, height = height)
}

# UMAP colored by other_cluster, with ddqc cluster number labels at centroids
DimPlotCluster <- function(obj, lbls, ttl) {
  data <- data.frame(
    UMAP1         = obj$umap1,
    UMAP2         = obj$umap2,
    cluster       = factor(obj$cluster_labels),
    other_cluster = factor(obj$other_cluster)
  )
  plot <- ggplot(data, aes(x = UMAP1, y = UMAP2, color = other_cluster)) +
    geom_point(size = 0.25) + theme_umap + ttl +
    scale_fill_discrete(labels = lbls)
  for (cl in levels(factor(obj$cluster_labels))) {
    cluster.data <- subset(data, cluster == cl)
    plot <- plot + annotate("text",
                            x     = mean(cluster.data$UMAP1),
                            y     = mean(cluster.data$UMAP2),
                            label = cl, size = 7, fontface = 2)
  }
  plot
}

# Stacked barplot: per ddqc cluster, proportion of each annotation label
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
  lbls       <- seq_along(unique(obj$cluster_labels))
  names(lbls) <- seq_along(lbls)
  n.clusters  <- length(unique(obj$cluster_labels))
  ttl         <- ggtitle(plot.title)

  ggsave1(
    filename         = file.path(results.dir, paste0(plot.title, "_umap_clusters.pdf")),
    plot             = DimPlotCluster(obj, lbls, ttl),
    n.clusters       = n.clusters,
    type             = "u",
    width.multiplier = width.multiplier
  )
  ggsave1(
    filename   = file.path(results.dir, paste0(plot.title, "_barplot.pdf")),
    plot       = makeBarplot(obj, lbls, ttl),
    n.clusters = n.clusters,
    type       = "h"
  )
}


# ============================================================
# PART B — KRASNOW LUNG
# Run after cell_typist.py + merge_classifications.py → summary.csv
# SingleR skipped (no Human Lung reference in celldex at time of paper)
# Azimuth skipped (web interface had problems with this dataset)
# ============================================================

if (file.exists(LUNG_SUMMARY_CSV) && file.exists(LUNG_CELLS_INITIAL_CSV)) {

  message("=== Krasnow Lung: MAD filtering ===")

  classification_data <- read.csv(LUNG_SUMMARY_CSV, row.names = 1)

  # Run MAD filter for each grouping column — matches authors' mad_filtering.py __main__
  # do_counts=FALSE as per authors' script
  for (cl in c("cell_typist", "ddqc_cluster")) {
    classification_data <- filter_cells(
      data           = classification_data,
      classification = cl,
      method         = "mad",
      threshold      = 2,
      do_counts      = FALSE
    )
    n_pass <- sum(classification_data[[paste0(cl, "_passed_qc")]], na.rm = TRUE)
    message(sprintf("  %s: %d / %d cells pass QC (%.1f%%)",
                    cl, n_pass, nrow(classification_data),
                    100 * n_pass / nrow(classification_data)))
  }

  write.csv(classification_data,
            file.path(OUT_DIR, "krasnow_lung_classification_mad_results.csv"))

  # Concordance between the two methods
  agree <- mean(
    classification_data$ddqc_cluster_passed_qc ==
    classification_data$cell_typist_passed_qc,
    na.rm = TRUE
  )
  message(sprintf("Concordance ddqc_cluster vs CellTypist: %.2f%%", 100 * agree))

  # Cross-tabulation (Table S5)
  tab <- table(
    ddqc_cluster = classification_data$ddqc_cluster,
    cell_typist  = classification_data$cell_typist
  )
  write.csv(as.data.frame.matrix(tab),
            file.path(OUT_DIR, "krasnow_lung_table_s5_ddqc_vs_celltypist.csv"))

  # ---- VISUALIZATION ----------------------------------------------------------
  # Matches celltypist_vs_ddqc_umap_plot.R exactly

  cells.with.coords           <- read.csv(LUNG_CELLS_INITIAL_CSV, row.names = 1)
  rownames(cells.with.coords) <- paste0(rownames(cells.with.coords), "-1")

  classification_data[["umap1"]] <- cells.with.coords[rownames(classification_data), "umap1"]
  classification_data[["umap2"]] <- cells.with.coords[rownames(classification_data), "umap2"]

  classification_data$ddqc_cluster <- as.factor(classification_data$ddqc_cluster)
  classification_data$cell_typist  <- as.factor(classification_data$cell_typist)

  # ddqc-passed cells: cluster = ddqc cluster, other_cluster = CellTypist label
  cells.ddqc <- classification_data[classification_data$ddqc_cluster_passed_qc == TRUE, ]
  cells.ddqc[["cluster_labels"]] <- as.factor(cells.ddqc$ddqc_cluster)
  cells.ddqc[["other_cluster"]]  <- as.factor(cells.ddqc$cell_typist)

  # CellTypist-passed cells: cluster = CellTypist label, other_cluster = ddqc cluster
  cells.celltypist <- classification_data[classification_data$cell_typist_passed_qc == TRUE, ]
  cells.celltypist[["cluster_labels"]] <- as.factor(cells.celltypist$cell_typist)
  cells.celltypist[["other_cluster"]]  <- as.factor(cells.celltypist$ddqc_cluster)

  makePlots(cells.ddqc,       OUT_DIR, plot.title = "ddqc_colored_by_celltypist")
  makePlots(cells.celltypist, OUT_DIR, plot.title = "celltypist_colored_by_ddqc")

} else {
  message("Lung summary.csv or !cells_initial.csv not found.")
  message("Run cell_typist.py and merge_classifications.py first.")
}


# ---- SESSION INFO -------------------------------------------------------------
sink(file.path(OUT_DIR, "session_info.txt"))
sessionInfo()
sink()

message("\n=== Done. Outputs in: ", OUT_DIR, " ===")