## ============================================================
## Panel E Equivalent: Canonical Marker Gene Validation
##
## Validates that ddqc-retained clusters express canonical
## immune cell marker genes, confirming biological coherence
## of the clustering.
##
## Approach:
##   1. Signature score per cluster (mean log2(TP10K+1))
##      — same method as Smillie et al. Cell 2019
##   2. Marker gene dotplot per cluster
##
## Canonical markers sourced from:
##   - Smillie et al. Cell 178:714-730 (2019) — colon immune atlas
##   - Gut et al. Science 371:eabb5793 (2021) — human gut cell atlas
##   - Dominguez Conde et al. Science 376:eabl5197 (2022) — cross-tissue immune
##   - Jardine et al. Nature 600:285-291 (2021) — fetal immune reference
##
## Input:
##   ddqc_analysis/seu_ddqc_final_processed.rds  (res=1.3)
##   OR provide path to any processed Seurat object with
##   seurat_clusters in meta.data and normalized RNA data
## ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
})


# ---- PATHS --------------------------------------------------------------------

BASE_DIR   <- "/projectnb/ds596/projects/Team 9/scQC_project"
seurat_rds <- file.path(BASE_DIR, "ddqc_analysis/seu_ddqc_final_processed.rds")
out_dir    <- file.path(BASE_DIR, "compare_annot/marker_validation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---- CANONICAL IMMUNE MARKERS ------------------------------------------------
#
# Sources:
#   T cells (CD3D, CD3E, CD3G):
#     Smillie et al. Cell 2019; Dominguez Conde et al. Science 2022
#   CD4 T cells (CD4, IL7R):
#     Smillie et al. Cell 2019
#   CD8 T cells (CD8A, CD8B):
#     Smillie et al. Cell 2019; Gut et al. Science 2021
#   Regulatory T cells (FOXP3, IL2RA):
#     Smillie et al. Cell 2019
#   NK cells (NCAM1, NKG7, KLRD1, GNLY):
#     Dominguez Conde et al. Science 2022; Jardine et al. Nature 2021
#   ILC (RORC, IL1R1):
#     Smillie et al. Cell 2019
#   B cells (CD79A, CD79B, MS4A1, CD19):
#     Smillie et al. Cell 2019; Dominguez Conde et al. Science 2022
#   Plasma cells (MZB1, IGHG1, JCHAIN, SDC1):
#     Smillie et al. Cell 2019; Gut et al. Science 2021
#   Monocytes/Macrophages (CD14, LYZ, CSF1R, CD68):
#     Smillie et al. Cell 2019; Dominguez Conde et al. Science 2022
#   Dendritic cells (ITGAX, CLEC9A, CD1C, FCER1A):
#     Smillie et al. Cell 2019
#   Mast cells (TPSAB1, CPA3, KIT):
#     Smillie et al. Cell 2019; Gut et al. Science 2021

MARKERS <- list(

  # Pan-T cell
  "T cells"              = c("CD3D", "CD3E", "CD3G"),

  # CD4 T subtypes
  "CD4+ T cells"         = c("CD4", "IL7R"),
  "Regulatory T cells"   = c("FOXP3", "IL2RA", "CTLA4"),

  # CD8 T subtypes
  "CD8+ T cells"         = c("CD8A", "CD8B"),
  "Cytotoxic T cells"    = c("GZMB", "GZMK", "PRF1"),

  # NK / ILC
  "NK cells"             = c("NCAM1", "NKG7", "KLRD1", "GNLY"),
  "ILC"                  = c("RORC", "IL1R1", "KIT"),

  # B cells
  "B cells"              = c("CD79A", "CD79B", "MS4A1", "CD19"),
  "Plasma cells"         = c("MZB1", "IGHG1", "JCHAIN", "SDC1"),
  "Germinal center B"    = c("AICDA", "BCL6", "CXCR4"),

  # Myeloid
  "Monocytes/Macrophages" = c("CD14", "LYZ", "CSF1R", "CD68"),
  "DC1"                  = c("CLEC9A", "XCR1", "CADM1"),
  "DC2"                  = c("CD1C", "FCER1A", "CLEC10A"),
  "Mast cells"           = c("TPSAB1", "CPA3", "KIT")
)

# Flat list of all markers for dotplot
ALL_MARKERS <- unique(unlist(MARKERS))


# ---- 1. LOAD SEURAT OBJECT ---------------------------------------------------

message("Loading Seurat object...")
seu <- readRDS(seurat_rds)
cat("Loaded:", ncol(seu), "cells,", nrow(seu), "genes\n")
cat("Clusters:", length(unique(seu$seurat_clusters)), "\n")

# Ensure cluster identity is set
Idents(seu) <- "seurat_clusters"

# Filter to markers present in the dataset
markers_present <- intersect(ALL_MARKERS, rownames(seu))
markers_missing <- setdiff(ALL_MARKERS, rownames(seu))
if (length(markers_missing) > 0)
  message("Markers not found in dataset: ",
          paste(markers_missing, collapse = ", "))


# ---- 2. SIGNATURE SCORES -----------------------------------------------------
# Method: mean log2(TP10K+1) per cluster per signature
# Matches Smillie et al. Cell 2019 Methods

message("Computing signature scores...")

# Get normalized data matrix
norm_mat <- GetAssayData(seu, assay = "RNA", layer = "data")

# Compute per-cell signature scores
for (sig_name in names(MARKERS)) {
  genes <- intersect(MARKERS[[sig_name]], rownames(norm_mat))
  if (length(genes) == 0) next
  score_col <- paste0("score_", gsub("[^A-Za-z0-9]", "_", sig_name))
  seu[[score_col]] <- colMeans(norm_mat[genes, , drop = FALSE])
}

# Compute per-cluster mean signature scores
score_cols <- grep("^score_", colnames(seu@meta.data), value = TRUE)

cluster_scores <- seu@meta.data %>%
  select(seurat_clusters, all_of(score_cols)) %>%
  group_by(seurat_clusters) %>%
  summarise(across(everything(), mean), .groups = "drop") %>%
  rename(cluster = seurat_clusters)

# Rename columns back to signature names
sig_clean <- gsub("[^A-Za-z0-9]", "_", names(MARKERS))
colnames(cluster_scores) <- c("cluster",
                               names(MARKERS)[match(
                                 gsub("^score_", "", score_cols),
                                 sig_clean
                               )])

write.csv(cluster_scores,
          file.path(out_dir, "cluster_signature_scores.csv"),
          row.names = FALSE)


# ---- 3. SIGNATURE SCORE HEATMAP ----------------------------------------------
# Shows which cell type signature each cluster expresses most strongly
# Analogous to Fig 3E in Subramanian et al. Genome Biology 2022

message("Plotting signature score heatmap...")

score_long <- cluster_scores %>%
  pivot_longer(-cluster,
               names_to  = "signature",
               values_to = "score") %>%
  mutate(cluster = factor(cluster,
                           levels = as.character(
                             sort(as.numeric(as.character(
                               unique(cluster_scores$cluster)
                             )))
                           )))

# Identify dominant signature per cluster
dominant <- score_long %>%
  group_by(cluster) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(cluster, dominant_sig = signature)

# Order signatures by cell type family
sig_order <- c(
  "T cells", "CD4+ T cells", "Regulatory T cells",
  "CD8+ T cells", "Cytotoxic T cells",
  "NK cells", "ILC",
  "B cells", "Plasma cells", "Germinal center B",
  "Monocytes/Macrophages", "DC1", "DC2", "Mast cells"
)
sig_order <- intersect(sig_order, unique(score_long$signature))
score_long$signature <- factor(score_long$signature, levels = sig_order)

p_sig_heat <- ggplot(score_long,
                     aes(x = signature, y = cluster, fill = score)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradient2(
    low      = "#f7f7f7",
    mid      = "#fee08b",
    high     = "#d73027",
    midpoint = median(score_long$score),
    name     = "Mean\nlog2(TP10K+1)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.caption = element_text(size = 8, color = "grey50", hjust = 0),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y  = element_text(size = 8),
    axis.title   = element_blank(),
    panel.grid   = element_blank(),
    legend.title = element_text(size = 9)
  ) +
  labs(
    title   = "Cluster Signature Scores — Canonical Immune Markers",
    caption = paste0(
      "Signature score = mean log2(TP10K+1) across marker genes per cluster\n",
      "Markers sourced from: Smillie et al. Cell 2019; ",
      "Dominguez Conde et al. Science 2022;\n",
      "Gut et al. Science 2021; Jardine et al. Nature 2021"
    )
  )

n_clusters <- nrow(cluster_scores)
n_sigs     <- length(sig_order)

ggsave(
  file.path(out_dir, "signature_score_heatmap.png"),
  p_sig_heat,
  width  = max(7, n_sigs * 0.55),
  height = max(5, n_clusters * 0.28),
  dpi    = 450,
  limitsize = FALSE
)
message("Saved: signature_score_heatmap.png")


# ---- 4. SIGNATURE SCORE VIOLIN PLOT ------------------------------------------
# Per-cluster, per-signature — matches Fig 3E style from paper
# Shows distribution across cells within each cluster

message("Plotting signature score violins...")

# Add cluster labels to metadata for plotting
meta_scores <- seu@meta.data %>%
  select(seurat_clusters, all_of(score_cols)) %>%
  rename(cluster = seurat_clusters)

colnames(meta_scores)[-1] <- names(MARKERS)[match(
  gsub("^score_", "", score_cols), sig_clean
)]

meta_long <- meta_scores %>%
  pivot_longer(-cluster,
               names_to  = "signature",
               values_to = "score") %>%
  mutate(
    cluster   = factor(cluster,
                        levels = as.character(
                          sort(as.numeric(as.character(unique(meta_scores$cluster))))
                        )),
    signature = factor(signature, levels = sig_order)
  )

p_violin <- ggplot(meta_long,
                   aes(x = cluster, y = score, fill = signature)) +
  geom_violin(scale = "width", linewidth = 0.2, color = NA) +
  facet_wrap(~signature, ncol = 2, scales = "free_y") +
  theme_bw(base_size = 10) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.caption = element_text(size = 7, color = "grey50", hjust = 0),
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                 size = 6),
    axis.text.y  = element_text(size = 7),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 9),
    strip.text   = element_text(face = "bold", size = 9),
    legend.position = "none",
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title   = "Signature Score Distribution per Cluster",
    y       = "Mean log2(TP10K+1)",
    caption = paste0(
      "Markers: Smillie et al. Cell 2019; Dominguez Conde et al. Science 2022; ",
      "Gut et al. Science 2021; Jardine et al. Nature 2021"
    )
  )

ggsave(
  file.path(out_dir, "signature_score_violin.png"),
  p_violin,
  width  = 14,
  height = max(10, ceiling(n_sigs / 2) * 2.5),
  dpi    = 450,
  limitsize = FALSE
)
message("Saved: signature_score_violin.png")


# ---- 5. MARKER GENE DOTPLOT --------------------------------------------------
# Standard Seurat dotplot showing % expressing and mean expression
# per cluster for all canonical markers

message("Plotting marker gene dotplot...")

# Only use markers present in dataset
markers_for_dot <- intersect(ALL_MARKERS, rownames(seu))

p_dot <- DotPlot(
  seu,
  features = markers_for_dot,
  group.by = "seurat_clusters",
  dot.scale = 4
) +
  scale_color_gradient2(
    low      = "#2166ac",
    mid      = "#f7f7f7",
    high     = "#d73027",
    midpoint = 0,
    name     = "Avg\nExpression"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5, size = 12),
    plot.caption = element_text(size = 7, color = "grey50", hjust = 0),
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y  = element_text(size = 8),
    axis.title   = element_blank(),
    legend.title = element_text(size = 9),
    panel.grid   = element_blank()
  ) +
  labs(
    title   = "Canonical Immune Marker Expression per ddqc Cluster",
    caption = paste0(
      "Dot size = % cells expressing marker (> 0); ",
      "Color = scaled average expression\n",
      "Markers: Smillie et al. Cell 2019; Dominguez Conde et al. Science 2022; ",
      "Gut et al. Science 2021; Jardine et al. Nature 2021"
    )
  )

# Add vertical separators between marker groups
marker_group_sizes <- sapply(MARKERS, function(g) length(intersect(g, rownames(seu))))
marker_group_sizes <- marker_group_sizes[marker_group_sizes > 0]
sep_positions <- cumsum(marker_group_sizes)[-length(marker_group_sizes)] + 0.5

for (xpos in sep_positions) {
  p_dot <- p_dot + geom_vline(xintercept = xpos,
                               linetype = "dashed",
                               color = "grey60",
                               linewidth = 0.4)
}

ggsave(
  file.path(out_dir, "marker_dotplot.png"),
  p_dot,
  width  = max(10, length(markers_for_dot) * 0.45),
  height = max(6,  n_clusters * 0.28),
  dpi    = 450,
  limitsize = FALSE
)
message("Saved: marker_dotplot.png")


# ---- 6. DOMINANT SIGNATURE UMAP ----------------------------------------------
# UMAP colored by dominant cell type signature per cluster
# Provides biological labels derived from marker genes independently

message("Plotting dominant signature UMAP...")

seu <- RunUMAP(seu, dims = 1:50, verbose = FALSE)

# Add dominant signature per cell (via cluster)
seu$dominant_signature <- dominant$dominant_sig[
  match(as.character(seu$seurat_clusters), as.character(dominant$cluster))
]

p_umap_sig <- DimPlot(
  seu,
  group.by   = "dominant_signature",
  label      = TRUE,
  repel      = TRUE,
  label.size = 3,
  pt.size    = 0.15,
  raster     = TRUE
) +
  ggtitle("ddqc Clusters — Dominant Canonical Signature") +
  theme_bw(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold", hjust = 0.5),
    legend.title    = element_text(size = 10, face = "bold"),
    legend.text     = element_text(size = 9),
    legend.key.size = unit(0.4, "cm")
  ) +
  labs(
    caption = paste0(
      "Cell type assigned by dominant canonical marker gene signature\n",
      "Markers: Smillie et al. Cell 2019; Dominguez Conde et al. Science 2022"
    )
  )

ggsave(
  file.path(out_dir, "umap_dominant_signature.png"),
  p_umap_sig,
  width = 12, height = 8, dpi = 450
)
message("Saved: umap_dominant_signature.png")


# ---- SESSION INFO -------------------------------------------------------------

sink(file.path(out_dir, "session_info.txt"))
sessionInfo()
sink()

message("\n=== Done. Outputs in: ", out_dir, " ===")
message("Outputs:")
message("  signature_score_heatmap.png  — cluster x signature score matrix")
message("  signature_score_violin.png   — per-signature score distributions")
message("  marker_dotplot.png           — canonical marker dotplot")
message("  umap_dominant_signature.png  — UMAP colored by dominant signature")
message("  cluster_signature_scores.csv — raw scores table")
message("")
message("Citation for markers:")
message("  Smillie et al. Cell 178:714-730 (2019)")
message("  Dominguez Conde et al. Science 376:eabl5197 (2022)")
message("  Gut et al. Science 371:eabb5793 (2021)")
message("  Jardine et al. Nature 600:285-291 (2021)")