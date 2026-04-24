import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models

# ---- PATHS --------------------------------------------------------------------

mtx_dir      = "/projectnb/ds596/projects/Team 9/scQC_project/data/recreate_colon"
matrix_file  = "gene_sorted-Imm.matrix.mtx"
barcodes_file = "Imm.barcodes2.tsv"
genes_file   = "Imm.genes.tsv"
out_csv      = "/projectnb/ds596/projects/Team 9/scQC_project/cell_annot/colon_immune_celltypist.csv"

# ---- DOWNLOAD MODELS (skip if already cached) ---------------------------------

models.download_models(force_update=False)

# List available models — run once to pick the right one
print(models.models_description())

# ---- LOAD DATA ----------------------------------------------------------------
# Read MTX directly — matches your three-file structure

adata = sc.read_mtx(f"{mtx_dir}/{matrix_file}").T  # transpose: cells x genes

# Assign barcodes and gene names
barcodes = pd.read_csv(f"{mtx_dir}/{barcodes_file}", header=None)[0].values
genes    = pd.read_csv(f"{mtx_dir}/{genes_file}",    header=None)[0].values

adata.obs_names = barcodes   # cells
adata.var_names = genes      # genes

print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")
print("Barcode check:", adata.obs_names[:3].tolist())

# ---- BASIC FILTER (mirrors R script) ------------------------------------------
# Same thresholds as ddqc_upstream: min_genes=100, percent_mito < 80

sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)

adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)
adata = adata[adata.obs.pct_counts_mt < 80, :]

print(f"After basic filter: {adata.n_obs} cells")

# ---- NORMALIZE (required by CellTypist) ---------------------------------------

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# ---- CELLTYPIST ANNOTATION ----------------------------------------------------
# Choose model appropriate for colon immune cells
# Run models.models_description() above to confirm model name

model_name  = "Immune_All_Low.pkl"   # adjust if needed after checking models_description()
model       = models.Model.load(model=model_name)

predictions = celltypist.annotate(adata, model=model, majority_voting=True)
labels      = predictions.predicted_labels

print("CellTypist label distribution:")
print(labels["majority_voting"].value_counts())

# ---- SAVE ---------------------------------------------------------------------
# Save with barcodes as index — matches cells_initial.csv 'cell' column format

labels.to_csv(out_csv)
print(f"Saved: {out_csv}")
