## ============================================================
## miQC Pipeline
## Produces cells_initial.csv in same format as ddqc upstream
## for use in annotation comparison
##
## Output files:
##   cells_initial.csv
##     Columns: cell, n_counts, n_genes, percent_mito,
##               percent_ribo, keep_miqc
##   miQC_filtered.rds
##   miQC_resolution_summary.csv
##   plots: plotMetrics, plotModel, plotFiltering, UMAP
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(scater)
  library(miQC)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

# ---- PATHS --------------------------------------------------------------------

rds_path <- "/projectnb/ds596/projects/Team 9/scQC_project/colon_immune_raw.rds"
out_dir  <- "/projectnb/ds596/projects/Team 9/scQC_project/recreation/miQC"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---- 1. LOAD DATA -------------------------------------------------------------

seurat_obj <- readRDS(rds_path)
cat("Loaded object with", nrow(seurat_obj), "genes and", ncol(seurat_obj), "cells\n")

# Fix NULL dimnames (Seurat v5 Assay5 issue)
counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
rownames(counts_mat) <- rownames(seurat_obj)
colnames(counts_mat) <- colnames(seurat_obj)
seurat_obj <- SetAssayData(seurat_obj, assay = "RNA", layer = "counts",
                           new.data = counts_mat)
cat("Barcodes check:", head(colnames(seurat_obj), 3), "\n")


# ---- 2. COMPUTE QC METRICS ----------------------------------------------------

if (!"percent.mt" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
}

if (!"percent.ribo" %in% colnames(seurat_obj@meta.data)) {
  ribo_genes <- grep("^RP[SL]", rownames(seurat_obj), value = TRUE, ignore.case = TRUE)
  if (length(ribo_genes) > 0) {
    seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, features = ribo_genes)
  } else {
    seurat_obj[["percent.ribo"]] <- 0
    warning("No ribosomal genes matched; percent.ribo set to 0.")
  }
}


# ---- 3. CONVERT TO SCE -------------------------------------------------------

sce <- as.SingleCellExperiment(seurat_obj)


# ---- 4. ADD PER-CELL QC METRICS ----------------------------------------------

mt_genes      <- grepl("^MT-", rownames(sce))
feature_ctrls <- list(mito = rownames(sce)[mt_genes])
sce           <- addPerCellQC(sce, subsets = feature_ctrls)


# ---- 5. DIAGNOSTIC PLOTS -----------------------------------------------------

png(file.path(out_dir, "plot_metrics.png"), width = 800, height = 600)
print(plotMetrics(sce))
dev.off()


# ---- 6. FIT MIXTURE MODEL + FILTER -------------------------------------------

model <- mixtureModel(sce)

png(file.path(out_dir, "plot_model.png"), width = 800, height = 600)
print(plotModel(sce, model))
dev.off()

png(file.path(out_dir, "plot_filtering.png"), width = 800, height = 600)
print(plotFiltering(sce, model))
dev.off()

sce_filtered <- filterCells(sce, model)

n_before <- ncol(sce)
n_after  <- ncol(sce_filtered)
n_removed <- n_before - n_after

cat("Cells before miQC:", n_before, "\n")
cat("Cells after miQC :", n_after, "\n")
cat("Cells removed    :", n_removed, "\n")
cat("Fraction removed :", round(n_removed / n_before, 4), "\n")


# ---- 7. SUMMARY TABLE --------------------------------------------------------

passed_barcodes <- colnames(sce_filtered)

summary_df <- data.frame(
  method        = "miQC",
  cells_before  = n_before,
  cells_after   = n_after,
  cells_removed = n_removed,
  frac_removed  = round(n_removed / n_before, 4)
)

write.csv(
  summary_df,
  file      = file.path(out_dir, "miQC_summary.csv"),
  row.names = FALSE
)

cat("Saved: miQC_summary.csv\n")


# ---- 8. CONVERT BACK TO SEURAT + REPROCESS -----------------------------------
# Must happen before cells_initial.csv so UMAP coords are available

seurat_filtered <- as.Seurat(sce_filtered, counts = "counts", data = NULL)

seurat_filtered <- seurat_filtered %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30, verbose = FALSE)


# ---- 9. BUILD PER-CELL TABLE (mirrors cells_initial.csv from ddqc) -----------
# All cells from before filtering, with keep_miqc flag
# UMAP coords from seurat_filtered (now available after step 8)

umap_embed <- Embeddings(seurat_filtered, reduction = "umap")
umap_df    <- data.frame(
  cell  = rownames(umap_embed),
  umap1 = umap_embed[, 1],
  umap2 = umap_embed[, 2],
  stringsAsFactors = FALSE
)

cells_initial <- seurat_obj@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(
    n_counts     = nCount_RNA,
    n_genes      = nFeature_RNA,
    percent_mito = percent.mt,
    percent_ribo = percent.ribo,
    keep_miqc    = cell %in% passed_barcodes
  ) %>%
  select(cell, n_counts, n_genes, percent_mito, percent_ribo, keep_miqc) %>%
  left_join(umap_df, by = "cell")

write.csv(
  cells_initial,
  file      = file.path(out_dir, "cells_initial_miQC.csv"),
  row.names = FALSE
)

cat("Saved: cells_initial.csv\n")


# ---- 10. UMAP PLOT -----------------------------------------------------------

umap_plot <- DimPlot(seurat_filtered, reduction = "umap") +
  ggtitle("miQC Filtered — UMAP") +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(
  filename = file.path(out_dir, "umap.png"),
  plot     = umap_plot,
  width    = 6,
  height   = 5,
  dpi      = 300
)

cat("Saved: umap.png\n")


# ---- 11. SAVE FILTERED SEURAT OBJECT -----------------------------------------

saveRDS(seurat_filtered, file = file.path(out_dir, "miQC_filtered.rds"))
cat("Saved: miQC_filtered.rds\n")


# ---- SESSION INFO -------------------------------------------------------------

sink(file.path(out_dir, "session_info.txt"))
sessionInfo()
sink()

cat("\n=== Done. Key output: ", out_dir, "/cells_initial.csv ===\n", sep = "")
cat("Next step: join SingleR/CellTypist labels to cells_initial.csv by 'cell' barcode\n")