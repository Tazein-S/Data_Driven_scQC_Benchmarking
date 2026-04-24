## ============================================================
## ddqc Upstream Pipeline
## Produces cells_initial.csv for annotation comparison
##
## Output file: ddqc_analysis/cells_initial.csv
##   Columns: cell, cluster, n_counts, n_genes,
##             percent_mito, percent_ribo, keep_ddqc
##
## This file is then used as the base for the annotation
## comparison (SingleR/CellTypist vs ddqc clusters)
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# ---- PATHS --------------------------------------------------------------------

infile  <- "/projectnb/ds596/projects/Team 9/scQC_project/colon_immune_raw.rds"
out_dir <- "ddqc_analysis"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---- PARAMETERS ---------------------------------------------------------------
# Keep these consistent with whatever you use in the annotation comparison script

basic_n_counts           <- 0
basic_n_genes            <- 100
basic_percent_mito       <- 80

ddqc_threshold           <- 2.0
npcs_use                 <- 50
k_param                  <- 20
cluster_resolution       <- 1.3
random_seed              <- 29

n_genes_lower_bound      <- 200   # floor for n_genes lower cutoff
percent_mito_upper_bound <- 10    # ceiling for percent_mito upper cutoff


# ---- 1. LOAD DATA -------------------------------------------------------------

seurat_obj <- readRDS(infile)
cat("Loaded object with", nrow(seurat_obj), "genes and", ncol(seurat_obj), "cells\n")


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


# ---- 3. BASIC FILTER (empty droplet removal) ----------------------------------
# Mirrors ddqc default: <100 genes or >80% mito removed
# This produces the equivalent of !cells_initial in the authors' pipeline

seu_prefilt <- subset(
  seurat_obj,
  subset = nCount_RNA   >= basic_n_counts &
           nFeature_RNA >= basic_n_genes  &
           percent.mt   <= basic_percent_mito
)

cat("Cells after basic filtering:", ncol(seu_prefilt), "\n")


# ---- 4. NORMALIZE + CLUSTER (pre-ddqc) ----------------------------------------
# These cluster IDs are what ddqc uses as grouping variable for MAD filtering

set.seed(random_seed)

seu_cluster <- seu_prefilt %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = npcs_use, verbose = FALSE) %>%
  FindNeighbors(dims = 1:npcs_use, k.param = k_param, verbose = FALSE) %>%
  FindClusters(
    resolution   = cluster_resolution,
    algorithm    = 1,           # Louvain
    random.seed  = random_seed,
    verbose      = FALSE
  ) %>%
  RunUMAP(dims = 1:npcs_use, verbose = FALSE)

cat("Clusters found:", length(unique(seu_cluster$seurat_clusters)), "\n")


# ---- 5. BUILD PER-CELL METADATA TABLE -----------------------------------------

meta_ddqc <- seu_cluster@meta.data %>%
  rownames_to_column("cell") %>%
  mutate(
    cluster      = as.character(seurat_clusters),
    n_counts     = nCount_RNA,
    n_genes      = nFeature_RNA,
    percent_mito = percent.mt,
    percent_ribo = percent.ribo
  )


# ---- 6. COMPUTE CLUSTER-LEVEL MAD CUTOFFS ------------------------------------

ddqc_cutoffs <- meta_ddqc %>%
  group_by(cluster) %>%
  summarise(
    n_cells       = n(),

    counts_median = median(n_counts,     na.rm = TRUE),
    counts_mad    = mad(n_counts,
                        center   = median(n_counts, na.rm = TRUE),
                        constant = 1.4826,
                        na.rm    = TRUE),

    genes_median  = median(n_genes,      na.rm = TRUE),
    genes_mad     = mad(n_genes,
                        center   = median(n_genes, na.rm = TRUE),
                        constant = 1.4826,
                        na.rm    = TRUE),

    mito_median   = median(percent_mito, na.rm = TRUE),
    mito_mad      = mad(percent_mito,
                        center   = median(percent_mito, na.rm = TRUE),
                        constant = 1.4826,
                        na.rm    = TRUE),

    ribo_median   = median(percent_ribo, na.rm = TRUE),
    ribo_mad      = mad(percent_ribo,
                        center   = median(percent_ribo, na.rm = TRUE),
                        constant = 1.4826,
                        na.rm    = TRUE),

    .groups = "drop"
  ) %>%
  mutate(
    # n_counts: lower only
    n_counts_lower     = pmax(0, counts_median - ddqc_threshold * counts_mad),
    n_counts_upper     = counts_median + ddqc_threshold * counts_mad,

    # n_genes: lower with floor of n_genes_lower_bound (200)
    n_genes_lower_raw  = genes_median - ddqc_threshold * genes_mad,
    n_genes_lower      = pmin(n_genes_lower_raw, n_genes_lower_bound),
    n_genes_upper      = genes_median + ddqc_threshold * genes_mad,

    # percent_mito: upper with ceiling of percent_mito_upper_bound (10)
    percent_mito_upper_raw = mito_median + ddqc_threshold * mito_mad,
    percent_mito_upper     = pmax(percent_mito_upper_raw, percent_mito_upper_bound),

    # percent_ribo: upper only
    percent_ribo_upper = ribo_median + ddqc_threshold * ribo_mad
  ) %>%
  select(
    cluster,
    n_cells,
    n_counts_lower, n_counts_upper,
    n_genes_lower,  n_genes_upper,
    percent_mito_upper,
    percent_ribo_upper
  )


# ---- 7. APPLY CUTOFFS PER CELL ------------------------------------------------

meta_ddqc_labeled <- meta_ddqc %>%
  left_join(ddqc_cutoffs, by = "cluster") %>%
  mutate(
    keep_ddqc =
      n_counts     >= n_counts_lower     &
      n_counts     <= n_counts_upper     &
      n_genes      >= n_genes_lower      &
      n_genes      <= n_genes_upper      &
      percent_mito <= percent_mito_upper &
      percent_ribo <= percent_ribo_upper
  )

cells_keep_ddqc <- meta_ddqc_labeled %>%
  filter(keep_ddqc) %>%
  pull(cell)

cat("Cells kept after ddqc :", length(cells_keep_ddqc), "\n")
cat("Cells removed by ddqc :", ncol(seu_prefilt) - length(cells_keep_ddqc), "\n")
cat("Fraction removed       :",
    round((ncol(seu_prefilt) - length(cells_keep_ddqc)) / ncol(seu_prefilt), 4), "\n")


# ---- 8. SAVE OUTPUT FILES -----------------------------------------------------

# -- 8a. cells_initial.csv ------------------------------------------------------
# This is the KEY file for the annotation comparison.
# Contains one row per cell with:
#   cell         → barcode (join key for SingleR/CellTypist labels)
#   cluster      → ddqc cluster ID (grouping variable for MAD filtering)
#   n_counts, n_genes, percent_mito, percent_ribo → QC metrics
#   keep_ddqc    → pass/fail from ddqc MAD filtering (comparison baseline)

cells_initial <- meta_ddqc_labeled %>%
  select(cell, cluster, n_counts, n_genes, percent_mito, percent_ribo, keep_ddqc)

write.csv(
  cells_initial,
  file      = file.path(out_dir, "cells_initial.csv"),
  row.names = FALSE
)

cat("Saved: cells_initial.csv\n")

# -- 8b. Cluster-level cutoffs (for reporting) ----------------------------------

cluster_removal_summary <- meta_ddqc_labeled %>%
  group_by(cluster) %>%
  summarise(
    n_before     = n(),
    n_after      = sum(keep_ddqc),
    n_removed    = sum(!keep_ddqc),
    frac_removed = n_removed / n_before,
    .groups      = "drop"
  )

ddqc_cutoffs_full <- ddqc_cutoffs %>%
  left_join(cluster_removal_summary, by = "cluster")

write.csv(
  ddqc_cutoffs_full,
  file      = file.path(out_dir, "ddqc_cluster_cutoffs_with_removal_for_annot.csv"),
  row.names = FALSE
)

cat("Saved: ddqc_cluster_cutoffs_with_removal_for_annot.csv\n")


# ---- 9. QC PLOTS --------------------------------------------------------------

p_box_mt <- ggplot(meta_ddqc_labeled, aes(x = factor(cluster), y = percent_mito)) +
  geom_boxplot(outlier.size = 0.1) +
  geom_hline(yintercept = percent_mito_upper_bound, color = "red", linewidth = 0.8) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = "Percent mitochondrial reads by cluster", x = "Cluster", y = "percent_mito")

p_box_genes <- ggplot(meta_ddqc_labeled, aes(x = factor(cluster), y = n_genes)) +
  geom_boxplot(outlier.size = 0.1) +
  geom_hline(yintercept = n_genes_lower_bound, color = "red", linewidth = 0.8) +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  labs(title = "Number of genes by cluster", x = "Cluster", y = "n_genes")

ggsave(
  filename = file.path(out_dir, "ddqc_cluster_boxplots.pdf"),
  plot     = p_box_mt / p_box_genes,
  width    = 12, height = 10
)

cat("Saved: ddqc_cluster_boxplots.pdf\n")


# ---- SESSION INFO -------------------------------------------------------------

sink(file.path(out_dir, "session_info.txt"))
sessionInfo()
sink()

cat("\n=== Done. Key output: ddqc_analysis/cells_initial.csv ===\n")
cat("Next step: join SingleR/CellTypist labels to cells_initial.csv by 'cell' barcode\n")