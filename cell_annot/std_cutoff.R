## ============================================================
## Standard Cutoff — Annotation Comparison
## Same scheme as ddqc annotation comparison:
##   1. Standard cutoff filtering
##   2. Cluster filtered cells
##   3. SingleR annotation
##   4. Save cells_initial.csv (with umap coords)
##   Then run mad_filtering_comparison.py and
##   annotation_comparison_plots.R as usual
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# ---- PATHS --------------------------------------------------------------------

infile      <- "/projectnb/ds596/projects/Team 9/scQC_project/colon_immune_raw.rds"
out_dir     <- "/projectnb/ds596/projects/Team 9/scQC_project/recreation/standard_cutoff"
singler_out <- "/projectnb/ds596/projects/Team 9/scQC_project/cell_classification/SingleR_stdcutoff.csv"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---- PARAMETERS ---------------------------------------------------------------

basic_n_genes        <- 100
basic_percent_mito   <- 80
std_n_genes_min      <- 200
std_percent_mito_max <- 10
npcs_use             <- 50
k_param              <- 20
cluster_resolution   <- 1.3
random_seed          <- 29


# ---- 1. LOAD + FIX DIMNAMES --------------------------------------------------

seurat_obj <- readRDS(infile)
cat("Loaded:", nrow(seurat_obj), "genes x", ncol(seurat_obj), "cells\n")

counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
rownames(counts_mat) <- rownames(seurat_obj)
colnames(counts_mat) <- colnames(seurat_obj)
seurat_obj <- SetAssayData(seurat_obj, assay = "RNA", layer = "counts",
                           new.data = counts_mat)
cat("Barcodes check:", head(colnames(seurat_obj), 3), "\n")


# ---- 2. QC METRICS -----------------------------------------------------------

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

ribo_genes <- grep("^RP[SL]", rownames(seurat_obj), value = TRUE, ignore.case = TRUE)
seurat_obj[["percent.ribo"]] <- if (length(ribo_genes) > 0) {
  PercentageFeatureSet(seurat_obj, features = ribo_genes)
} else 0


# ---- 3. BASIC FILTER (same as ddqc) ------------------------------------------

seu_prefilt <- subset(seurat_obj,
                      subset = nFeature_RNA >= basic_n_genes &
                               percent.mt   <  basic_percent_mito)

counts     <- GetAssayData(seu_prefilt, layer = "counts")
genes_keep <- rowSums(counts > 0) >= 3
seu_prefilt <- seu_prefilt[genes_keep, ]

cat("After basic filter:", ncol(seu_prefilt), "cells\n")


# ---- 4. STANDARD CUTOFF FILTER -----------------------------------------------

seu_std <- subset(seu_prefilt,
                  subset = nFeature_RNA >= std_n_genes_min &
                           percent.mt   <  std_percent_mito_max)

n_before  <- ncol(seu_prefilt)
n_after   <- ncol(seu_std)
n_removed <- n_before - n_after

cat("After standard cutoff:", n_after, "cells\n")
cat("Removed:", n_removed, sprintf("(%.1f%%)\n", 100 * n_removed / n_before))

passed_std_barcodes <- colnames(seu_std)


# ---- 5. NORMALIZE + CLUSTER --------------------------------------------------

set.seed(random_seed)

seu_std <- seu_std %>%
  NormalizeData(normalization.method = "LogNormalize",
                scale.factor = 10000, verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst",
                       nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(features = VariableFeatures(.), verbose = FALSE) %>%
  RunPCA(npcs = npcs_use,
         features = VariableFeatures(.), verbose = FALSE) %>%
  FindNeighbors(dims = 1:npcs_use, k.param = k_param, verbose = FALSE) %>%
  FindClusters(resolution  = cluster_resolution,
               algorithm   = 1,
               random.seed = random_seed,
               verbose     = FALSE) %>%
  RunUMAP(dims = 1:npcs_use, verbose = FALSE)

cat("Clusters found:", length(unique(seu_std$seurat_clusters)), "\n")


# ---- 7. BUILD cells_initial.csv ----------------------------------------------
# Same structure as ddqc cells_initial_res_{res}.csv:
#   cell, cluster, n_counts, n_genes, percent_mito, percent_ribo,
#   keep_std, umap1, umap2
#
# Only cells that passed standard cutoff get cluster + umap coords
# Cells filtered out get NA in those columns

umap_embed <- Embeddings(seu_std, reduction = "umap")
umap_df    <- data.frame(
  cell  = rownames(umap_embed),
  umap1 = umap_embed[, 1],
  umap2 = umap_embed[, 2],
  stringsAsFactors = FALSE
)

cluster_df <- data.frame(
  cell    = colnames(seu_std),
  cluster = as.character(seu_std$seurat_clusters),
  stringsAsFactors = FALSE
)

cells_initial <- seu_prefilt@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(
    n_counts     = nCount_RNA,
    n_genes      = nFeature_RNA,
    percent_mito = percent.mt,
    percent_ribo = percent.ribo,
    keep_std     = cell %in% passed_std_barcodes
  ) %>%
  select(cell, n_counts, n_genes, percent_mito, percent_ribo, keep_std) %>%
  left_join(cluster_df, by = "cell") %>%
  left_join(umap_df,    by = "cell")

write.csv(cells_initial,
          file.path(out_dir, "stdCutoff_cells_initial.csv"),
          row.names = FALSE)

cat("Saved: stdCutoff_cells_initial.csv\n")


# ---- 8. QC + UMAP PLOTS ------------------------------------------------------

p_vln <- VlnPlot(seu_prefilt,
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                 ncol = 3, pt.size = 0)

ggsave(file.path(out_dir, "qc_violin_pre_filter.png"),
       p_vln, width = 14, height = 5, dpi = 300)

p_umap <- DimPlot(seu_std, reduction = "umap",
                  label = TRUE, repel = TRUE, label.size = 5) +
  ggtitle(paste0("Standard Cutoff UMAP (res = ", cluster_resolution, ")")) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave(file.path(out_dir, paste0("UMAP_stdcutoff_res_", cluster_resolution, ".png")),
       p_umap, width = 10, height = 8, dpi = 300)

cat("Saved: UMAP plot\n")


# ---- 9. SUMMARY TABLE --------------------------------------------------------

summary_df <- data.frame(
  method           = "standard_cutoff",
  n_genes_min      = std_n_genes_min,
  percent_mito_max = std_percent_mito_max,
  resolution       = cluster_resolution,
  cells_before     = n_before,
  cells_after      = n_after,
  cells_removed    = n_removed,
  frac_removed     = round(n_removed / n_before, 4),
  n_clusters       = length(unique(seu_std$seurat_clusters))
)

write.csv(summary_df,
          file.path(out_dir, "stdcutoff_summary.csv"),
          row.names = FALSE)

saveRDS(seu_std, file.path(out_dir, "stdcutoff_filtered.rds"))


# ---- SESSION INFO -------------------------------------------------------------

sink(file.path(out_dir, "session_info.txt"))
sessionInfo()
sink()

cat("\n=== Done ===\n")
cat("Key outputs:\n")
cat("  cells_initial.csv     → input for mad_filtering_comparison.py\n")
cat("  SingleR_stdcutoff.csv → annotation labels\n")
cat("  UMAP_stdcutoff.png    → UMAP plot\n")