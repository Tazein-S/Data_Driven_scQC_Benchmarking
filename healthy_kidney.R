library(Seurat)

mat <- read.table(
  gzfile("data/extension_kidney/GSM5222734_s200707N5.expression_matrix.txt.gz"),
  header = TRUE, row.names = 1
)

healthy <- CreateSeuratObject(counts = mat, project = "GSM5222734")
saveRDS(healthy, "data/extension_kidney/GSE171314_scRNA_healthy.rds")