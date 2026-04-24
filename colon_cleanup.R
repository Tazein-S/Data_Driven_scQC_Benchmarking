library(Seurat)
library(Matrix)

# -- Paths --------------------------------------------------------------------
mtx_file      <- "data/recreate_colon/gene_sorted-Imm.matrix.mtx"
barcodes_file <- "data/recreate_colon/Imm.barcodes2.tsv"
genes_file    <- "data/recreate_colon/Imm.genes.tsv"
meta_file     <- "data/recreate_colon/all.meta2.txt"
out_rds       <- "colon_immune_raw.rds"

# -- Load ---------------------------------------------------------------------
message("Loading count matrix...")
counts <- ReadMtx(
  mtx            = mtx_file,
  cells          = barcodes_file,
  features       = genes_file,
  feature.column = 1
)

seurat_obj <- CreateSeuratObject(
  counts  = counts,
  project = "smillie_colon_immune",
  min.cells    = 0,
  min.features = 0
)
message(sprintf("Loaded %d genes x %d cells", nrow(seurat_obj), ncol(seurat_obj)))

# -- Metadata -----------------------------------------------------------------
message("Adding metadata...")
meta <- read.table(meta_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)
meta_immune <- meta[rownames(meta) %in% colnames(seurat_obj), ]
seurat_obj <- AddMetaData(seurat_obj, metadata = meta_immune)

# -- Save ---------------------------------------------------------------------
message(sprintf("Saving to %s ...", out_rds))
saveRDS(seurat_obj, file = out_rds)
message("Done.")