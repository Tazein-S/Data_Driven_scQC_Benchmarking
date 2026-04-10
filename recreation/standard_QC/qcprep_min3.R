#1. Load the MTX + barcodes + genes into Seurat with ReadMtx
#2. Add metadata from all.meta2.txt, subset to immune barcodes
#3. Compute QC metrics (nGene, nUMI, percent mito)

# Steps 1-3: Load, metadata, QC metrics
# Output: smillie_colon_immune_qc.rds
# =============================================================================

library(Seurat)
library(Matrix)

# -- Paths --------------------------------------------------------------------
# Update these to match your data directory on the SCC
mtx_file      <- "/projectnb/ds596/projects/Team 9/scQC_project/data/recreate_colon/gene_sorted-Imm.matrix.mtx"
barcodes_file <- "/projectnb/ds596/projects/Team 9/scQC_project/data/recreate_colon/Imm.barcodes2.tsv"
genes_file    <- "/projectnb/ds596/projects/Team 9/scQC_project/data/recreate_colon/Imm.genes.tsv"
meta_file     <- "/projectnb/ds596/projects/Team 9/scQC_project/data/recreate_colon/all.meta2.txt"
out_rds       <- "colon_qc_min3.rds"

# -----------------------------------------------------------------------------
# Step 1: Load MTX into Seurat
# -----------------------------------------------------------------------------
message("Loading count matrix...")

counts <- ReadMtx(
  mtx      = mtx_file,
  cells    = barcodes_file,
  features = genes_file,
  feature.column = 1  # genes file has one column (gene symbol)
)

seurat_obj <- CreateSeuratObject(
  counts  = counts,
  project = "smillie_colon_immune",
  min.cells = 3 
)

message(sprintf("Loaded %d genes x %d cells", nrow(seurat_obj), ncol(seurat_obj)))

# -----------------------------------------------------------------------------
# Step 2: Add metadata, subset to immune barcodes
# -----------------------------------------------------------------------------
message("Adding metadata...")

meta <- read.table(
  meta_file,
  header    = TRUE,
  sep       = "\t",
  row.names = 1,
  check.names = FALSE
)

# The Seurat object was built from the immune MTX, so all barcodes are already
# immune. This step pulls only the matching rows from the full metadata file.
immune_barcodes <- colnames(seurat_obj)
meta_immune     <- meta[rownames(meta) %in% immune_barcodes, ]

message(sprintf(
  "Matched %d / %d barcodes in metadata",
  nrow(meta_immune), length(immune_barcodes)
))

# Warn if any barcodes are unmatched (e.g. suffix mismatch like -1)
unmatched <- sum(!immune_barcodes %in% rownames(meta))
if (unmatched > 0) {
  warning(sprintf(
    "%d barcodes not found in metadata. Check for suffix mismatches (e.g. '-1').",
    unmatched
  ))
}

seurat_obj <- AddMetaData(seurat_obj, metadata = meta_immune)

# -----------------------------------------------------------------------------
# Step 3: Compute QC metrics
# -----------------------------------------------------------------------------
message("Computing QC metrics...")

# nFeature_RNA (nGene) and nCount_RNA (nUMI) are added automatically
# by CreateSeuratObject. Adding percent mito and ribo here.

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

# Quick sanity check
message("QC metric summary:")
print(summary(seurat_obj@meta.data[, c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb")]))

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
message(sprintf("Saving to %s ...", out_rds))
saveRDS(seurat_obj, file = out_rds)
message("Done.")