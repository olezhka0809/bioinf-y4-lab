"""
Exercise 10 — PCA Single-Omics vs Joint

Steps:
- load SNP and Expression matrices
- align samples
- z-score normalize each layer
- run PCA on:
    1) SNP
    2) Expression
    3) Joint (concatenation)
- export 3 PNG figures
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# =========================
# CONFIG
# =========================
HANDLE = "olezhka0809"

SNP_CSV = Path(f"/workspaces/bioinf-y4-lab/data/work/olezhka0809/lab10/snp_matrix_olezhka0809.csv")
EXP_CSV = Path(f"/workspaces/bioinf-y4-lab/data/work/olezhka0809/lab10/expression_matrix_olezhka0809.csv")

OUT_DIR = Path(f"/workspaces/bioinf-y4-lab/labs/10_integrative/submissions/olezhka0809")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# =========================
# LOAD DATA
# =========================
print("[INFO] Loading data...")

snp = pd.read_csv(SNP_CSV, index_col=0)
exp = pd.read_csv(EXP_CSV, index_col=0)

print("SNP shape:", snp.shape)
print("Expression shape:", exp.shape)

# =========================
# ALIGN SAMPLES
# =========================
common_samples = snp.index.intersection(exp.index)

snp = snp.loc[common_samples]
exp = exp.loc[common_samples]

print("[INFO] Aligned samples:", len(common_samples))

# =========================
# NORMALIZATION (z-score)
# =========================
print("[INFO] Z-score normalization...")

scaler_snp = StandardScaler()
scaler_exp = StandardScaler()

snp_z = scaler_snp.fit_transform(snp)
exp_z = scaler_exp.fit_transform(exp)

# =========================
# PCA FUNCTION
# =========================
def run_pca(X, title, out_png):
    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(X)

    plt.figure(figsize=(6, 5))
    plt.scatter(coords[:, 0], coords[:, 1], s=30, alpha=0.7)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title(
        f"{title}\nExplained variance: "
        f"{pca.explained_variance_ratio_[0]:.2f}, "
        f"{pca.explained_variance_ratio_[1]:.2f}"
    )
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

    print(f"[OK] Saved: {out_png}")
    print("Explained variance:", pca.explained_variance_ratio_)

# =========================
# PCA 1 — SNP
# =========================
print("[INFO] PCA on SNP layer...")
run_pca(
    snp_z,
    "PCA — SNP layer",
    OUT_DIR / f"pca_snp_{HANDLE}.png"
)

# =========================
# PCA 2 — Expression
# =========================
print("[INFO] PCA on Expression layer...")
run_pca(
    exp_z,
    "PCA — Expression layer",
    OUT_DIR / f"pca_expression_{HANDLE}.png"
)

# =========================
# PCA 3 — Joint (concat)
# =========================
print("[INFO] PCA on Joint (SNP + Expression)...")

joint = pd.concat(
    [
        pd.DataFrame(snp_z, index=common_samples),
        pd.DataFrame(exp_z, index=common_samples)
    ],
    axis=1
)

run_pca(
    joint.values,
    "PCA — Joint SNP + Expression",
    OUT_DIR / f"pca_joint_{HANDLE}.png"
)

print("\n✔ Exercise 10 completed successfully.")
