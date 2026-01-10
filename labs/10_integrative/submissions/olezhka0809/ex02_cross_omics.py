"""
Exercise 10.2 — Identify top SNP–Gene correlations

- load integrated multi-omics matrix (samples x features)
- split columns into SNPs vs Genes (heuristics on column names)
- compute correlations r(SNP, Gene)
- filter |r| > threshold (default 0.5)
- export snp_gene_pairs_<handle>.csv

Input expected:
  labs/10_integrative/submissions/<handle>/multiomics_concat_<handle>.csv

Output:
  labs/10_integrative/submissions/<handle>/snp_gene_pairs_<handle>.csv
"""

from pathlib import Path
import argparse
import numpy as np
import pandas as pd


def is_snp_col(col: str) -> bool:
    c = str(col).strip()
    c_low = c.lower()
    return (
        c_low.startswith("rs") or
        c_low.startswith("snp_") or
        c_low.startswith("snp-")
    )


def is_gene_col(col: str) -> bool:
    c = str(col).strip()
    c_low = c.lower()
    # Ensembl gene IDs or common gene symbols or explicit prefixes
    return (
        c.startswith("ENSG") or
        c_low.startswith("gene_") or
        c_low.startswith("gene-") or
        (c.isalpha() and c.isupper() and len(c) >= 2)  # TP53, APOE etc.
    )


def zscore_matrix(X: np.ndarray) -> np.ndarray:
    """Z-score per column; safe for constant columns."""
    mean = np.nanmean(X, axis=0)
    std = np.nanstd(X, axis=0)
    std[std == 0] = 1.0
    return (X - mean) / std


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--handle", required=True)
    parser.add_argument(
        "--threshold", type=float, default=0.5,
        help="Filter pairs by |r| > threshold. Default 0.5"
    )
    parser.add_argument(
        "--chunk_genes", type=int, default=500,
        help="Process genes in chunks to limit memory. Default 500."
    )
    args = parser.parse_args()

    handle = args.handle
    JOINT_CSV = Path(f"/workspaces/bioinf-y4-lab/labs/10_integrative/submissions/olezhka0809/multiomics_concat_olezhka0809.csv")
    OUT_CSV = Path(f"labs/10_integrative/submissions/olezhka0809/snp_gene_pairs_olezhka0809.csv")
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)

    print("[INFO] Loading joint multi-omics matrix...")
    df = pd.read_csv(JOINT_CSV, index_col=0)
    print("  Shape:", df.shape)

    # Keep numeric only (critical)
    df = df.apply(pd.to_numeric, errors="coerce")
    # drop all-NaN columns
    df = df.dropna(axis=1, how="all")

    # Identify SNP vs Gene columns
    snp_cols = [c for c in df.columns if is_snp_col(c)]
    gene_cols = [c for c in df.columns if is_gene_col(c) and c not in snp_cols]

    print(f"[INFO] Detected SNP cols:  {len(snp_cols)}")
    print(f"[INFO] Detected Gene cols: {len(gene_cols)}")

    if len(snp_cols) == 0 or len(gene_cols) == 0:
        raise ValueError(
            "Nu am putut separa SNP vs Gene pe baza numelor de coloane.\n"
            "Verifică naming-ul coloanelor din multiomics_concat și ajustează "
            "funcțiile is_snp_col / is_gene_col."
        )

    # Extract matrices (samples x snps) and (samples x genes)
    Xs = df[snp_cols].to_numpy(dtype=float)
    Xg = df[gene_cols].to_numpy(dtype=float)

    # Drop rows with too many NaNs? For simplicity: fill NaN with column mean.
    # (Alternative: drop rows with any NaN, but that may discard many samples.)
    print("[INFO] Imputare NaN cu media pe coloană...")
    # SNP
    snp_means = np.nanmean(Xs, axis=0)
    inds = np.where(np.isnan(Xs))
    Xs[inds] = np.take(snp_means, inds[1])
    # Gene
    gene_means = np.nanmean(Xg, axis=0)
    inds = np.where(np.isnan(Xg))
    Xg[inds] = np.take(gene_means, inds[1])

    # Z-score each matrix per column (so corr becomes dot product / (n-1))
    print("[INFO] Z-score normalizare pentru SNP și Gene...")
    Xs_z = zscore_matrix(Xs)
    Xg_z = zscore_matrix(Xg)

    n = Xs_z.shape[0]
    denom = max(n - 1, 1)

    print("[INFO] Calcul corelații SNP–Gene (chunked)...")
    results = []

    # Correlation matrix would be (num_snps x num_genes) -> can be huge.
    # We chunk genes to keep memory bounded.
    for start in range(0, Xg_z.shape[1], args.chunk_genes):
        end = min(start + args.chunk_genes, Xg_z.shape[1])
        chunk = Xg_z[:, start:end]  # (samples x chunk_genes)

        # corr(snps, genes_chunk) = (Xs_z.T @ chunk) / (n-1)
        corr_block = (Xs_z.T @ chunk) / denom  # shape: (num_snps x chunk_genes)

        # filter by threshold
        mask = np.abs(corr_block) > args.threshold
        if not mask.any():
            continue

        snp_idx, gene_idx = np.where(mask)
        for i, j in zip(snp_idx.tolist(), gene_idx.tolist()):
            results.append((snp_cols[i], gene_cols[start + j], float(corr_block[i, j])))

        print(f"  Chunk genes {start}:{end} -> pairs found: {mask.sum()}")

    out_df = pd.DataFrame(results, columns=["SNP", "Gene", "r"])
    out_df["abs_r"] = out_df["r"].abs()
    out_df = out_df.sort_values("abs_r", ascending=False).drop(columns=["abs_r"])

    print(f"[INFO] Total pairs |r|>{args.threshold}: {len(out_df)}")
    out_df.to_csv(OUT_CSV, index=False)
    print(f"[OK] Saved: {OUT_CSV}")


if __name__ == "__main__":
    main()
