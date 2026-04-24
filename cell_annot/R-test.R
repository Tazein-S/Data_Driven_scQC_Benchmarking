library(dplyr)

cells <- read.csv("/projectnb/ds596/projects/Team 9/scQC_project/compare_annot/ddqc_classification_mad_results_res_0.4.csv", row.names = 1)

if ("cluster" %in% colnames(cells)) cells <- cells %>% rename(ddqc_cluster = cluster)

# Check 1: how many pass each filter
cat("ddqc_cluster_passed_qc TRUE:", sum(cells$ddqc_cluster_passed_qc == TRUE, na.rm = TRUE), "\n")
cat("cell_typist_passed_qc TRUE:", sum(cells$cell_typist_passed_qc == TRUE, na.rm = TRUE), "\n")
cat("single_r_passed_qc TRUE:", sum(cells$single_r_passed_qc == TRUE, na.rm = TRUE), "\n")

# Check 2: NAs in key columns
cat("NA umap1:", sum(is.na(cells$umap1)), "\n")
cat("NA umap2:", sum(is.na(cells$umap2)), "\n")
cat("NA cell_typist:", sum(is.na(cells$cell_typist)), "\n")
cat("NA single_r:", sum(is.na(cells$single_r)), "\n")
cat("NA ddqc_cluster:", sum(is.na(cells$ddqc_cluster)), "\n")

# Check 3: sample values
cat("\nSample ddqc_cluster values:\n")
print(head(cells$ddqc_cluster))
cat("\nSample cell_typist values:\n")
print(head(cells$cell_typist))
cat("\nSample umap1 values:\n")
print(head(cells$umap1))
class(cells$ddqc_cluster_passed_qc)
head(cells$ddqc_cluster_passed_qc)
unique(cells$ddqc_cluster_passed_qc)