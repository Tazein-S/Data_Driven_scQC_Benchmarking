library(Seurat)

files <- c(
  "data/extension_kidney/GSM5222730_s200525a4.expression_matrix.txt.gz",
  "data/extension_kidney/GSM5222731_s200702a5.expression_matrix.txt.gz",
  "data/extension_kidney/GSM5222732_s200707a7.expression_matrix.txt.gz",
  "data/extension_kidney/GSM5222733_s200709a8.expression_matrix.txt.gz"
)

sample_ids <- c("GSM5222730", "GSM5222731", "GSM5222732", "GSM5222733")

samples <- lapply(seq_along(files), function(i) {
  mat <- read.table(gzfile(files[i]), header = TRUE, row.names = 1)
  obj <- CreateSeuratObject(counts = mat, project = sample_ids[i])
  obj
})

disease <- merge(samples[[1]], y = samples[2:4])
disease <- JoinLayers(disease)

saveRDS(disease, "data/extension_kidney/GSE171314_scRNA_disease.rds")