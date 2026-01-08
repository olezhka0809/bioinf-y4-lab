"""
Exercise 9.1 — Drug–Gene Bipartite Network & Drug Similarity Network (SOLUTION)

Scop:
- Construiți o rețea bipartită drug–gene
- Proiectați layer-ul de medicamente folosind Jaccard similarity
- Exportați tabelul cu muchiile de similaritate

Output files:
  1. drug_summary_{HANDLE}.csv          — drug, num_targets
  2. drug_similarity_{HANDLE}.csv       — drug1, drug2, similarity
  3. network_drug_gene_{HANDLE}.gpickle — bipartite graph
  4. network_drug_similarity_{HANDLE}.gpickle — drug similarity graph
  5. network_drug_gene_{HANDLE}.png     — visualization
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Set, Tuple, List

import itertools
import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import numpy as np

# =====================================================
# CONFIG
# =====================================================
HANDLE = "olezhka0809"

# ✅ Use absolute paths or resolve from current file location
# Support both relative and absolute paths
_script_dir = Path(__file__).parent.resolve()
_project_root = _script_dir.parent.parent.parent  # Go up to project root

DRUG_GENE_CSV = Path(f"/workspaces/bioinf-y4-lab/data/work/{HANDLE}/lab09/drug_gene_{HANDLE}.csv")
if not DRUG_GENE_CSV.exists():
    # Fallback: try relative path
    DRUG_GENE_CSV = Path(f"data/work/{HANDLE}/lab09/drug_gene_{HANDLE}.csv")

OUT_DIR = Path(f"/workspaces/bioinf-y4-lab/labs/09_repurposing/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_DRUG_SUMMARY = OUT_DIR / f"drug_summary_{HANDLE}.csv"
OUT_DRUG_SIMILARITY = OUT_DIR / f"drug_similarity_{HANDLE}.csv"
OUT_GRAPH_DRUG_GENE = OUT_DIR / f"network_drug_gene_{HANDLE}.gpickle"
OUT_GRAPH_DRUG_SIM = OUT_DIR / f"network_drug_similarity_{HANDLE}.gpickle"
OUT_VIZ_BIPARTITE = OUT_DIR / f"network_drug_gene_{HANDLE}.png"
OUT_VIZ_DRUG_SIM = OUT_DIR / f"network_drug_similarity_{HANDLE}.png"

# =====================================================
# UTILITY FUNCTIONS
# =====================================================
def ensure_exists(path: Path) -> None:
    """Verify file exists, raise FileNotFoundError if not."""
    if not path.is_file():
        raise FileNotFoundError(f"File not found: {path}\nMake sure to run generate_dataset_alzheimer.py first!")
    print(f"✓ Found: {path}")

def load_drug_gene_table(path: Path) -> pd.DataFrame:
    """Load and validate drug-gene CSV."""
    df = pd.read_csv(path)
    required_cols = {"drug", "gene"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"CSV must contain columns: {required_cols}. Got: {df.columns.tolist()}")
    return df

def build_drug2genes(df: pd.DataFrame) -> Dict[str, Set[str]]:
    """Create mapping: drug -> set of target genes."""
    drug2genes = {}
    for drug, group in df.groupby("drug"):
        genes = set(group["gene"].unique())
        drug2genes[drug] = genes
    return drug2genes

def build_bipartite_graph(drug2genes: Dict[str, Set[str]]) -> nx.Graph:
    """
    Build bipartite drug-gene network.
    
    Nodes:
      - drug nodes with attribute bipartite="drug"
      - gene nodes with attribute bipartite="gene"
    Edges:
      - drug-gene edges (unweighted)
    """
    B = nx.Graph()
    
    # Add drug nodes
    for drug in drug2genes.keys():
        B.add_node(drug, bipartite="drug")
    
    # Add gene nodes and edges
    for drug, genes in drug2genes.items():
        for gene in genes:
            B.add_node(gene, bipartite="gene")
            B.add_edge(drug, gene)
    
    return B

def summarize_drugs(drug2genes: Dict[str, Set[str]]) -> pd.DataFrame:
    """Create summary DataFrame: drug, num_targets."""
    data = [
        {"drug": drug, "num_targets": len(genes)}
        for drug, genes in drug2genes.items()
    ]
    return pd.DataFrame(data).sort_values("num_targets", ascending=False)

def jaccard_similarity(s1: Set[str], s2: Set[str]) -> float:
    """
    Jaccard similarity: J(A, B) = |A ∩ B| / |A ∪ B|
    """
    if not s1 and not s2:
        return 0.0
    inter = len(s1 & s2)
    union = len(s1 | s2)
    return inter / union if union > 0 else 0.0

def compute_drug_similarity_edges(
    drug2genes: Dict[str, Set[str]],
    min_sim: float = 0.1,
) -> List[Tuple[str, str, float]]:
    """
    Compute Jaccard similarity between all drug pairs.
    
    Returns:
      List of (drug1, drug2, similarity) tuples with similarity >= min_sim
    """
    drugs = list(drug2genes.keys())
    edges = []
    
    for drug1, drug2 in itertools.combinations(drugs, 2):
        sim = jaccard_similarity(drug2genes[drug1], drug2genes[drug2])
        if sim >= min_sim:
            edges.append((drug1, drug2, sim))
    
    return edges

def edges_to_dataframe(edges: List[Tuple[str, str, float]]) -> pd.DataFrame:
    """Convert edge list to DataFrame."""
    return pd.DataFrame(edges, columns=["drug1", "drug2", "similarity"])

def build_drug_similarity_graph(edges: List[Tuple[str, str, float]]) -> nx.Graph:
    """Build weighted drug-drug similarity network."""
    G = nx.Graph()
    for drug1, drug2, weight in edges:
        G.add_edge(drug1, drug2, weight=weight)
    return G

# =====================================================
# VISUALIZATION FUNCTIONS
# =====================================================
def visualize_bipartite_network(B: nx.Graph, drug2genes: Dict[str, Set[str]], output_path: Path) -> None:
    """
    Visualize bipartite network with drugs and genes.
    Color drugs by number of targets.
    """
    print(f"\n[VIZ] Drawing bipartite network...")
    
    fig, ax = plt.subplots(figsize=(14, 10), dpi=150)
    
    # Bipartite layout: drugs on left, genes on right
    pos = {}
    drugs = [n for n, d in B.nodes(data=True) if d.get("bipartite") == "drug"]
    genes = [n for n, d in B.nodes(data=True) if d.get("bipartite") == "gene"]
    
    # Position drugs on left
    for i, drug in enumerate(drugs):
        pos[drug] = (0, i - len(drugs) / 2)
    
    # Position genes on right
    for i, gene in enumerate(genes):
        pos[gene] = (10, i - len(genes) / 2)
    
    # Color drugs by num_targets
    drug_colors = [len(drug2genes[d]) for d in drugs]
    
    # Draw edges
    nx.draw_networkx_edges(B, pos, alpha=0.2, width=0.5, ax=ax)
    
    # Draw drug nodes
    nx.draw_networkx_nodes(
        B, pos,
        nodelist=drugs,
        node_color=drug_colors,
        node_size=500,
        node_shape='o',
        cmap='YlOrRd',
        ax=ax,
        label='Drugs'
    )
    
    # Draw gene nodes
    nx.draw_networkx_nodes(
        B, pos,
        nodelist=genes,
        node_color='lightblue',
        node_size=300,
        node_shape='s',
        ax=ax,
        label='Genes'
    )
    
    # Draw labels (only drugs for clarity)
    drug_labels = {d: d for d in drugs}
    nx.draw_networkx_labels(B, pos, drug_labels, font_size=8, ax=ax)
    
    ax.set_title(f"Bipartite Drug-Gene Network (Alzheimer's Disease)\nDrugs: {len(drugs)}, Genes: {len(genes)}", fontsize=14, fontweight='bold')
    ax.legend(scatterpoints=1, loc='upper left')
    ax.axis('off')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"✓ Saved visualization: {output_path}")
    plt.close()

def visualize_drug_similarity_network(G: nx.Graph, top_n: int = 50, output_path: Path = None) -> None:
    """
    Visualize drug-drug similarity network (projected graph).
    Show only top edges by similarity.
    """
    print(f"\n[VIZ] Drawing drug similarity network (top {top_n} edges)...")
    
    # Keep only top edges
    edges_sorted = sorted(
        [(u, v, d['weight']) for u, v, d in G.edges(data=True)],
        key=lambda x: x[2],
        reverse=True
    )
    G_top = nx.Graph()
    for u, v, w in edges_sorted[:top_n]:
        G_top.add_edge(u, v, weight=w)
    
    fig, ax = plt.subplots(figsize=(14, 10), dpi=150)
    
    # Spring layout
    pos = nx.spring_layout(G_top, k=0.5, iterations=50, seed=42)
    
    # Get edge weights for colors
    edges = G_top.edges()
    weights = [G_top[u][v]['weight'] for u, v in edges]
    
    # Draw edges (width and color by weight)
    nx.draw_networkx_edges(
        G_top, pos,
        width=[w * 5 for w in weights],
        edge_color=weights,
        edge_cmap=plt.cm.RdYlGn,
        ax=ax,
        alpha=0.7
    )
    
    # Draw nodes (size by degree)
    node_sizes = [G_top.degree(n) * 200 for n in G_top.nodes()]
    nx.draw_networkx_nodes(
        G_top, pos,
        node_size=node_sizes,
        node_color='lightblue',
        ax=ax
    )
    
    # Draw labels
    nx.draw_networkx_labels(G_top, pos, font_size=9, ax=ax)
    
    ax.set_title(f"Drug-Drug Similarity Network (Top {top_n} edges)\nJaccard Similarity", fontsize=14, fontweight='bold')
    ax.axis('off')
    
    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"✓ Saved visualization: {output_path}")
    plt.close()

# =====================================================
# MAIN EXECUTION
# =====================================================
def main():
    print("=" * 70)
    print("EXERCISE 9.1 — DRUG-GENE BIPARTITE & DRUG SIMILARITY NETWORK")
    print("=" * 70)
    
    # Step 1: Load data
    print("\n[1] Loading drug-gene interactions...")
    print(f"  Looking for: {DRUG_GENE_CSV}")
    ensure_exists(DRUG_GENE_CSV)
    df = load_drug_gene_table(DRUG_GENE_CSV)
    print(f"  Rows: {len(df)}")
    print(f"  Unique drugs: {df['drug'].nunique()}")
    print(f"  Unique genes: {df['gene'].nunique()}")
    
    # Step 2: Build drug->genes mapping
    print("\n[2] Building drug → gene mapping...")
    drug2genes = build_drug2genes(df)
    print(f"  Drugs in mapping: {len(drug2genes)}")
    
    # Step 3: Build bipartite graph
    print("\n[3] Building bipartite network...")
    B = build_bipartite_graph(drug2genes)
    print(f"  Nodes: {B.number_of_nodes()}")
    print(f"  Edges: {B.number_of_edges()}")
    
    # Step 4: Summarize drugs
    print("\n[4] Generating drug summary...")
    drug_summary = summarize_drugs(drug2genes)
    drug_summary.to_csv(OUT_DRUG_SUMMARY, index=False)
    print(f"✓ Saved: {OUT_DRUG_SUMMARY}")
    print(f"\nTop 10 drugs by target count:")
    print(drug_summary.head(10).to_string(index=False))
    
    # Step 5: Compute drug similarity
    print("\n[5] Computing drug-drug similarities (Jaccard)...")
    edges = compute_drug_similarity_edges(drug2genes, min_sim=0.05)
    print(f"  Edges (sim >= 0.05): {len(edges)}")
    
    # Step 6: Export similarity table
    sim_df = edges_to_dataframe(edges)
    sim_df.to_csv(OUT_DRUG_SIMILARITY, index=False)
    print(f"✓ Saved: {OUT_DRUG_SIMILARITY}")
    print(f"\nTop 15 similar drug pairs:")
    print(sim_df.nlargest(15, 'similarity').to_string(index=False))
    
    # Step 7: Build drug similarity graph
    print("\n[6] Building drug similarity graph...")
    G = build_drug_similarity_graph(edges)
    print(f"  Nodes: {G.number_of_nodes()}")
    print(f"  Edges: {G.number_of_edges()}")
    
    # Step 8: Save graphs
    print("\n[7] Saving graphs...")
    # NetworkX 3.0+ uses write_graphml instead of write_gpickle
    import pickle
    with open(OUT_GRAPH_DRUG_GENE, 'wb') as f:
        pickle.dump(B, f)
    print(f"✓ Saved bipartite graph: {OUT_GRAPH_DRUG_GENE}")
    
    with open(OUT_GRAPH_DRUG_SIM, 'wb') as f:
        pickle.dump(G, f)
    print(f"✓ Saved similarity graph: {OUT_GRAPH_DRUG_SIM}")
    
    # Step 9: Visualizations
    print("\n[8] Creating visualizations...")
    visualize_bipartite_network(B, drug2genes, OUT_VIZ_BIPARTITE)
    visualize_drug_similarity_network(G, top_n=40, output_path=OUT_VIZ_DRUG_SIM)
    
    print("\n" + "=" * 70)
    print("✓ EXERCISE 9.1 COMPLETED")
    print("=" * 70)
    print(f"\nOutput files:")
    print(f"  1. {OUT_DRUG_SUMMARY}")
    print(f"  2. {OUT_DRUG_SIMILARITY}")
    print(f"  3. {OUT_GRAPH_DRUG_GENE}")
    print(f"  4. {OUT_GRAPH_DRUG_SIM}")
    print(f"  5. {OUT_VIZ_BIPARTITE}")
    print(f"  6. {OUT_VIZ_DRUG_SIM}")

if __name__ == "__main__":
    main()