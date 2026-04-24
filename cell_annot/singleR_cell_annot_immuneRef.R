# Load packages
library(Seurat)
library(SingleR)
library(celldex)

# Load Seurat object
seurat.obj <- readRDS("/projectnb/ds596/projects/Team 9/scQC_project/colon_immune_raw.rds")

#Matching count matrix with cell barcode
#Extract the counts matrix
counts_mat <- GetAssayData(seurat.obj, assay = "RNA", layer = "counts")

# Assign dimnames from the assay's feature and cell slots
rownames(counts_mat) <- rownames(seurat.obj)
colnames(counts_mat) <- colnames(seurat.obj)

# Put it back
seurat.obj <- SetAssayData(seurat.obj, assay = "RNA", layer = "counts", new.data = counts_mat)

# 1. Filter cells FIRST (still a Seurat object)
seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > 10)
message(sprintf("After minimal filter: %d cells", ncol(seurat.obj)))
seurat.obj <- NormalizeData(seurat.obj)
# 2. Extract expression matrix (Seurat v5 style)
expr <- GetAssayData(seurat.obj, assay = "RNA", layer = "data")
head(colnames(expr))

# 3. Run SingleR
message("=== Running SingleR ===")

ref.data <- celldex::MonacoImmuneData()

predictions <- SingleR(
  test   = expr,
  ref    = ref.data,
  labels = ref.data$label.main
)

message("SingleR label distribution:")
print(table(predictions$labels))

# 4. Add labels back to Seurat object
seurat.obj <- AddMetaData(
  object   = seurat.obj,
  metadata = predictions$labels,
  col.name = "SingleR_label"
)

# 5. Save labels
labels <- predictions$labels
names(labels) <- rownames(predictions)
write.csv(labels, "SingleR_colon_immuneRef.csv")
message("SingleR labels saved.")