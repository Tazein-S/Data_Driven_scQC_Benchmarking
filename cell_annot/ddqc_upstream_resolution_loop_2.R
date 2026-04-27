## ============================================================
## ddqc Upstream Pipeline — Resolution Loop
## Produces cells_initial_{res}.csv for each resolution,
## for use in annotation comparison (SingleR/CellTypist vs ddqc)
##
## Output files per resolution:
##   cells_initial_{res}.csv
##     Columns: cell, cluster, n_counts, n_genes,
##               percent_mito, percent_ribo, keep_ddqc
##   ddqc_cluster_cutoffs_with_removal_{res}.csv
##   ddqc_cluster_boxplots_{res}.pdf
##
## Final summary:
##   ddqc_resolution_summary.csv
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
out_dir <- "/projectnb/ds596/projects/Team 9/scQC_project/ddqc_annot/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---- PARAMETERS ---------------------------------------------------------------

basic_n_counts           <- 0
basic_n_genes            <- 100
basic_percent_mito       <- 80

ddqc_threshold           <- 2.0
npcs_use                 <- 50
k_param                  <- 20
random_seed              <- 29

n_genes_lower_bound      <- 200
percent_mito_upper_bound <- 10

resolutions              <- c(0.4, 0.8, 1.3, 1.6, 2)


# ---- 1. LOAD DATA (once, outside loop) ----------------------------------------

seurat_obj <- readRDS(infile)
cat("Loaded object with", nrow(seurat_obj), "genes and", ncol(seurat_obj), "cells\n")

# Fix NULL dimnames in counts matrix (Seurat v5 Assay5 issue)
# Barcodes exist in the cells slot but are not propagated to matrix dimnames
counts_mat <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
rownames(counts_mat) <- rownames(seurat_obj)
colnames(counts_mat) <- colnames(seurat_obj)
seurat_obj <- SetAssayData(seurat_obj, assay = "RNA", layer = "counts",
                           new.data = counts_mat)
cat("Barcodes check:", head(colnames(seurat_obj), 3), "\n")


# ---- 2. COMPUTE QC METRICS (once) ---------------------------------------------

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


# ---- 3. BASIC FILTER (once) ---------------------------------------------------

seu_prefilt <- subset(
  seurat_obj,
  subset = nCount_RNA   >= basic_n_counts &
           nFeature_RNA >= basic_n_genes  &
           percent.mt   <= basic_percent_mito
)

cat("Cells after basic filtering:", ncol(seu_prefilt), "\n")


# ---- 4. SHARED PREPROCESSING (once, before resolution loop) -------------------
# Normalize, HVG, scale, PCA are resolution-independent — only clustering changes

set.seed(random_seed)

seu_base <- seu_prefilt %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = npcs_use, verbose = FALSE) %>%
  FindNeighbors(dims = 1:npcs_use, k.param = k_param, verbose = FALSE)

cat("Shared preprocessing done. Starting resolution loop...\n")


# ---- 5. RESOLUTION LOOP -------------------------------------------------------

results_list <- list()

for (res in resolutions) {

  cat("\n==============================\n")
  cat("Resolution:", res, "\n")
  cat("==============================\n")

  # -- 5a. Cluster at this resolution ------------------------------------------

  seu_cluster <- FindClusters(
    seu_base,
    resolution  = res,
    algorithm   = 1,          # Louvain
    random.seed = random_seed,
    verbose     = FALSE
  )

  seu_cluster <- RunUMAP(seu_cluster, dims = 1:npcs_use, verbose = FALSE)

  n_clusters <- length(unique(seu_cluster$seurat_clusters))
  cat("Clusters found:", n_clusters, "\n")

  # -- 5b. Build per-cell metadata table ---------------------------------------

  meta_ddqc <- seu_cluster@meta.data %>%
    rownames_to_column("cell") %>%
    mutate(
      cluster      = as.character(seurat_clusters),
      n_counts     = nCount_RNA,
      n_genes      = nFeature_RNA,
      percent_mito = percent.mt,
      percent_ribo = percent.ribo
    )

  # -- 5c. Compute cluster-level MAD cutoffs -----------------------------------

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
      n_counts_lower         = pmax(0, counts_median - ddqc_threshold * counts_mad),
      n_counts_upper         = counts_median + ddqc_threshold * counts_mad,

      # n_genes: lower with floor of n_genes_lower_bound (200)
      n_genes_lower_raw      = genes_median - ddqc_threshold * genes_mad,
      n_genes_lower          = pmin(n_genes_lower_raw, n_genes_lower_bound),
      n_genes_upper          = genes_median + ddqc_threshold * genes_mad,

      # percent_mito: upper with ceiling of percent_mito_upper_bound (10)
      percent_mito_upper_raw = mito_median + ddqc_threshold * mito_mad,
      percent_mito_upper     = pmax(percent_mito_upper_raw, percent_mito_upper_bound),

      # percent_ribo: upper only
      percent_ribo_upper     = ribo_median + ddqc_threshold * ribo_mad
    ) %>%
    select(
      cluster,
      n_cells,
      n_counts_lower, n_counts_upper,
      n_genes_lower,  n_genes_upper,
      percent_mito_upper,
      percent_ribo_upper
    )

  # -- 5d. Apply cutoffs per cell ----------------------------------------------

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

  n_kept    <- length(cells_keep_ddqc)
  n_removed <- ncol(seu_prefilt) - n_kept

  cat("Cells kept   :", n_kept, "\n")
  cat("Cells removed:", n_removed, "\n")
  cat("Fraction removed:", round(n_removed / ncol(seu_prefilt), 4), "\n")

  # -- 5e. Save cells_initial_{res}.csv ----------------------------------------
  # KEY output for annotation comparison — one file per resolution
  # Includes umap1/umap2 to match Pegasus !cells_initial.csv format

  umap_embed <- Embeddings(seu_cluster, reduction = "umap")
  umap_df    <- data.frame(
    cell  = rownames(umap_embed),
    umap1 = umap_embed[, 1],
    umap2 = umap_embed[, 2],
    stringsAsFactors = FALSE
  )

  cells_initial <- meta_ddqc_labeled %>%
    select(cell, cluster, n_counts, n_genes, percent_mito, percent_ribo, keep_ddqc) %>%
    left_join(umap_df, by = "cell")

  write.csv(
    cells_initial,
    file      = file.path(out_dir, paste0("cells_initial_res_", res, ".csv")),
    row.names = FALSE
  )

  cat("Saved: cells_initial_res_", res, ".csv\n", sep = "")

  # -- 5f. Save cluster cutoffs with removal stats -----------------------------

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
    file      = file.path(out_dir, paste0("ddqc_cluster_cutoffs_res_", res, ".csv")),
    row.names = FALSE
  )

  # -- 5g. QC boxplots ---------------------------------------------------------

  p_box_mt <- ggplot(meta_ddqc_labeled,
                     aes(x = factor(cluster, levels = as.character(sort(as.numeric(unique(cluster))))),
                         y = percent_mito)) +
    geom_boxplot(outlier.size = 0.1) +
    geom_hline(yintercept = percent_mito_upper_bound, color = "red", linewidth = 0.8) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = paste0("% Mito by cluster (res = ", res, ")"),
         x = "Cluster", y = "percent_mito")

  p_box_genes <- ggplot(meta_ddqc_labeled,
                        aes(x = factor(cluster, levels = as.character(sort(as.numeric(unique(cluster))))),
                            y = n_genes)) +
    geom_boxplot(outlier.size = 0.1) +
    geom_hline(yintercept = n_genes_lower_bound, color = "red", linewidth = 0.8) +
    theme_bw(base_size = 14) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    labs(title = paste0("n_genes by cluster (res = ", res, ")"),
         x = "Cluster", y = "n_genes")

  ggsave(
    filename = file.path(out_dir, paste0("ddqc_boxplots_res_", res, ".pdf")),
    plot     = p_box_mt / p_box_genes,
    width    = 12, height = 10
  )

  # -- 5h. Collect summary for across-resolution comparison --------------------

  results_list[[as.character(res)]] <- data.frame(
    resolution       = res,
    n_clusters       = n_clusters,
    cells_before_ddqc = ncol(seu_prefilt),
    cells_after_ddqc  = n_kept,
    cells_removed     = n_removed,
    frac_removed      = round(n_removed / ncol(seu_prefilt), 4)
  )
}


# ---- 6. SAVE CROSS-RESOLUTION SUMMARY ----------------------------------------

results_df <- do.call(rbind, results_list)
print(results_df)

write.csv(
  results_df,
  file      = file.path(out_dir, "ddqc_resolution_summary.csv"),
  row.names = FALSE
)

cat("\n=== Done ===\n")
cat("Output files per resolution in:", out_dir, "\n")
cat("  cells_initial_res_{res}.csv       → input for annotation comparison\n")
cat("  ddqc_cluster_cutoffs_res_{res}.csv → cluster-level cutoffs + removal stats\n")
cat("  ddqc_boxplots_res_{res}.pdf        → QC boxplots\n")
cat("  ddqc_resolution_summary.csv        → cross-resolution summary\n")


# ---- SESSION INFO -------------------------------------------------------------

sink(file.path(out_dir, "session_info.txt"))
sessionInfo()
sink()
