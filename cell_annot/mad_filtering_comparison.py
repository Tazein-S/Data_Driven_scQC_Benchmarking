import numpy as np
import pandas as pd
import os

# ---- PATHS --------------------------------------------------------------------

base_dir        = "/projectnb/ds596/projects/Team 9/scQC_project"
ddqc_dir        = f"{base_dir}/ddqc_annot"
miqc_dir        = f"{base_dir}/recreation/miQC"
annot_dir       = f"{base_dir}/cell_annot"
out_dir         = f"{base_dir}/compare_annot"

os.makedirs(out_dir, exist_ok=True)

celltypist_csv  = f"{annot_dir}/colon_immune_celltypist.csv"
singler_csv     = f"{annot_dir}/SingleR_colon_immuneRef.csv"

resolutions     = [0.4, 0.8, 1.3, 1.6, 2.0]

INF = 1e10


# ---- MAD FILTERING FUNCTIONS (unchanged from authors' mad_filtering.py) ------

def mad(data, axis=0, constant=1.4826):
    return constant * np.median(np.absolute(data - np.median(data, axis)), axis)


def metric_filter(data, classification, method, param, metric_name,
                  do_lower_co=False, do_upper_co=False,
                  lower_bound=INF, upper_bound=-INF):

    data[metric_name + "_qc_pass"]  = False
    data[metric_name + "_lower_co"] = None
    data[metric_name + "_upper_co"] = None

    if method == "mad":
        data[metric_name + "_median"] = None
        data[metric_name + "_mad"]    = None

    for ct in data[classification].unique():
        lower_co = -INF
        upper_co =  INF
        cell_type = data[data[classification] == ct]

        if method == "mad":
            data.loc[data[classification] == ct, metric_name + "_median"] = np.median(cell_type[metric_name])
            data.loc[data[classification] == ct, metric_name + "_mad"]    = mad(cell_type[metric_name])
            if do_lower_co:
                lower_co = min(np.median(cell_type[metric_name]) - param * mad(cell_type[metric_name]), lower_bound)
            if do_upper_co:
                upper_co = max(np.median(cell_type[metric_name]) + param * mad(cell_type[metric_name]), upper_bound)

        if method == "outlier":
            q75, q25 = np.percentile(cell_type[metric_name], [75, 25])
            if do_lower_co:
                lower_co = min(q25 - 1.5 * (q75 - q25), lower_bound)
            if do_upper_co:
                upper_co = max(q75 + 1.5 * (q75 - q25), upper_bound)

        filters = [
            data[classification] == ct,
            data[metric_name] >= lower_co,
            data[metric_name] <= upper_co
        ]
        if do_upper_co:
            data.loc[data[classification] == ct, metric_name + "_upper_co"] = upper_co
        if do_lower_co:
            data.loc[data[classification] == ct, metric_name + "_lower_co"] = lower_co

        data.loc[np.logical_and.reduce(filters), metric_name + "_qc_pass"] = True
    return data


def filter_cells(data, classification, method, threshold,
                 n_genes_lower_bound=200, percent_mito_upper_bound=10,
                 do_counts=True, do_genes=True, do_mito=True, do_ribo=False,
                 return_full=False):

    data_copy = data.copy()

    if do_counts:
        data_copy = metric_filter(data_copy, classification, method, threshold,
                                  "n_counts", do_lower_co=True)
    else:
        data_copy["n_counts_qc_pass"] = True

    if do_genes:
        data_copy = metric_filter(data_copy, classification, method, threshold,
                                  "n_genes", do_lower_co=True,
                                  lower_bound=n_genes_lower_bound)
    else:
        data_copy["n_genes_qc_pass"] = True

    if do_mito:
        data_copy = metric_filter(data_copy, classification, method, threshold,
                                  "percent_mito", do_upper_co=True,
                                  upper_bound=percent_mito_upper_bound)
    else:
        data_copy["percent_mito_qc_pass"] = True

    if do_ribo:
        data_copy = metric_filter(data_copy, classification, method, threshold,
                                  "percent_ribo", do_upper_co=True)
    else:
        data_copy["percent_ribo_qc_pass"] = True

    filters = [
        data_copy["n_counts_qc_pass"],
        data_copy["n_genes_qc_pass"],
        data_copy["percent_mito_qc_pass"],
        data_copy["percent_ribo_qc_pass"],
    ]
    data_copy["passed_qc"] = False
    data_copy.loc[np.logical_and.reduce(filters), "passed_qc"] = True

    if not return_full:
        data[classification + "_passed_qc"] = data_copy["passed_qc"]
        return data
    else:
        data_copy[classification + "_passed_qc"] = data_copy["passed_qc"]
        return data_copy


# ---- LOAD ANNOTATION LABELS (once — resolution-independent) ------------------

print("Loading annotation labels...")

# CellTypist — index is barcode, use majority_voting column
celltypist_df = pd.read_csv(celltypist_csv, index_col=0)
celltypist_df.index.name = "cell"
print(f"CellTypist: {len(celltypist_df)} cells")
print(f"CellTypist columns: {celltypist_df.columns.tolist()}")

# SingleR — index is barcode, label column (named 'x' from write.csv in R)
singler_df = pd.read_csv(singler_csv, index_col=0)
singler_df.index.name = "cell"
singler_df.columns = ["single_r"]   # rename to standard column name
print(f"SingleR: {len(singler_df)} cells")


# ---- HELPER: BUILD SUMMARY (join ddqc/miQC cells with annotation labels) -----

def build_summary(cells_df, celltypist_df, singler_df, has_cluster=True):
    """
    Join QC cell table with annotation labels by barcode.

    cells_df   : per-cell dataframe with 'cell' column as barcode
                 columns: cell, [cluster,] n_counts, n_genes,
                          percent_mito, percent_ribo, keep_*
    has_cluster: True for ddqc (has cluster column), False for miQC
    """
    df = cells_df.set_index("cell")

    # Join CellTypist
    ct_col = "majority_voting" if "majority_voting" in celltypist_df.columns else celltypist_df.columns[0]
    df["cell_typist"] = celltypist_df[ct_col]

    # Join SingleR
    df["single_r"] = singler_df["single_r"]

    # Report join success
    n_ct_matched = df["cell_typist"].notna().sum()
    n_sr_matched = df["single_r"].notna().sum()
    print(f"  CellTypist matched: {n_ct_matched} / {len(df)} cells")
    print(f"  SingleR matched   : {n_sr_matched} / {len(df)} cells")

    return df


# ---- CONCORDANCE HELPER ------------------------------------------------------

def concordance(df, col_a, col_b):
    valid = df[col_a].notna() & df[col_b].notna()
    return (df.loc[valid, col_a] == df.loc[valid, col_b]).mean()


# ============================================================
# PART A — ddqc: loop over resolutions
# ============================================================

print("\n" + "=" * 50)
print("PART A — ddqc annotation comparison")
print("=" * 50)

ddqc_concordance_rows = []

for res in resolutions:

    print(f"\n--- Resolution: {res} ---")

    cells_path = os.path.join(ddqc_dir, f"cells_initial_res_{res}.csv")
    if not os.path.exists(cells_path):
        print(f"  File not found: {cells_path} — skipping")
        continue

    cells_df = pd.read_csv(cells_path)
    print(f"  Loaded {len(cells_df)} cells")

    # Build summary table
    summary = build_summary(cells_df, celltypist_df, singler_df, has_cluster=True)

    # Drop cells with missing annotation labels
    summary_ct = summary.dropna(subset=["cell_typist"])
    summary_sr = summary.dropna(subset=["single_r"])

    # MAD filtering using each grouping column
    # do_counts=False — matches authors' mad_filtering.py exactly
    for classification in ["ddqc_cluster", "cell_typist", "single_r"]:
        sub = summary if classification == "ddqc_cluster" else \
              summary_ct if classification == "cell_typist" else summary_sr

        # Rename cluster col if needed
        if classification == "ddqc_cluster" and "cluster" in sub.columns:
            sub = sub.rename(columns={"cluster": "ddqc_cluster"})

        sub = filter_cells(
            sub.copy(),
            classification = classification,
            method         = "mad",
            threshold      = 2,
            do_counts      = False
        )

        # Write back passed_qc column to summary
        summary[classification + "_passed_qc"] = sub[classification + "_passed_qc"]

    # Save full results for this resolution
    out_path = os.path.join(out_dir, f"ddqc_classification_mad_results_res_{res}.csv")
    summary.to_csv(out_path)
    print(f"  Saved: {out_path}")

    # Concordance vs keep_ddqc
    for col in ["cell_typist_passed_qc", "single_r_passed_qc"]:
        if col in summary.columns:
            c = concordance(summary, "keep_ddqc", col)
            print(f"  Concordance keep_ddqc vs {col}: {c:.4f} ({c*100:.2f}%)")
            ddqc_concordance_rows.append({
                "resolution"  : res,
                "comparison"  : f"keep_ddqc vs {col}",
                "concordance" : round(c, 4)
            })

    # Cross-tabulation: ddqc cluster vs cell_typist
    if "cell_typist" in summary.columns and "cluster" in cells_df.columns:
        tab = pd.crosstab(
            summary["ddqc_cluster"] if "ddqc_cluster" in summary.columns
            else summary.rename(columns={"cluster": "ddqc_cluster"})["ddqc_cluster"],
            summary["cell_typist"]
        )
        tab.to_csv(os.path.join(out_dir, f"crosstab_ddqc_vs_celltypist_res_{res}.csv"))

    if "single_r" in summary.columns and "cluster" in cells_df.columns:
        tab = pd.crosstab(
            summary["ddqc_cluster"] if "ddqc_cluster" in summary.columns
            else summary.rename(columns={"cluster": "ddqc_cluster"})["ddqc_cluster"],
            summary["single_r"]
        )
        tab.to_csv(os.path.join(out_dir, f"crosstab_ddqc_vs_singler_res_{res}.csv"))

# Save ddqc concordance summary
if ddqc_concordance_rows:
    pd.DataFrame(ddqc_concordance_rows).to_csv(
        os.path.join(out_dir, "ddqc_concordance_summary.csv"), index=False
    )
    print(f"\nSaved: ddqc_concordance_summary.csv")


# ============================================================
# PART B — miQC: single file, no cluster column
# ============================================================

print("\n" + "=" * 50)
print("PART B — miQC annotation comparison")
print("=" * 50)

miqc_cells_path = os.path.join(miqc_dir, "cells_initial_miQC.csv")

if os.path.exists(miqc_cells_path):

    cells_df = pd.read_csv(miqc_cells_path)
    print(f"Loaded {len(cells_df)} cells")

    # Build summary — no cluster column for miQC
    summary = build_summary(cells_df, celltypist_df, singler_df, has_cluster=False)

    summary_ct = summary.dropna(subset=["cell_typist"])
    summary_sr = summary.dropna(subset=["single_r"])

    # MAD filtering using annotation labels as grouping
    for classification in ["cell_typist", "single_r"]:
        sub = summary_ct if classification == "cell_typist" else summary_sr

        sub = filter_cells(
            sub.copy(),
            classification = classification,
            method         = "mad",
            threshold      = 2,
            do_counts      = False
        )

        summary[classification + "_passed_qc"] = sub[classification + "_passed_qc"]
    
    # In mad_filtering_comparison.py, after the miQC MAD filtering block
    retention = pd.DataFrame({
    "method": ["miQC", "CellTypist MAD", "SingleR MAD"],
    "cells_retained": [
        summary["keep_miqc"].sum(),
        summary["cell_typist_passed_qc"].sum(),
        summary["single_r_passed_qc"].sum()
    ],
    "total_cells": len(summary),
    })
    retention["frac_retained"] = retention["cells_retained"] / retention["total_cells"]
    retention.to_csv(os.path.join(out_dir, "miqc_retention_summary.csv"), index=False)

    # Save full results
    out_path = os.path.join(out_dir, "miqc_classification_mad_results.csv")
    summary.to_csv(out_path)
    print(f"Saved: {out_path}")

    # Concordance vs keep_miqc
    miqc_concordance_rows = []
    for col in ["cell_typist_passed_qc", "single_r_passed_qc"]:
        if col in summary.columns:
            c = concordance(summary, "keep_miqc", col)
            print(f"Concordance keep_miqc vs {col}: {c:.4f} ({c*100:.2f}%)")
            miqc_concordance_rows.append({
                "comparison"  : f"keep_miqc vs {col}",
                "concordance" : round(c, 4)
            })

    if miqc_concordance_rows:
        pd.DataFrame(miqc_concordance_rows).to_csv(
            os.path.join(out_dir, "miqc_concordance_summary.csv"), index=False
        )
        print("Saved: miqc_concordance_summary.csv")

else:
    print(f"miQC cells_initial.csv not found at {miqc_cells_path} — skipping")


print("\n=== Done. Outputs in:", out_dir, "===")
