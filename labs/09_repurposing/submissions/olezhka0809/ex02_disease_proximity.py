"""
Exercise 9.2 — Disease Proximity and Drug Ranking (SOLUTION)

Scop:
- Calculați distanța medie dintre fiecare medicament și genele bolii
- Prioritizați medicamentele pentru repurposing pe baza network proximity
- Exportați ranking-ul medicamentelor

Output files:
  1. drug_priority_{HANDLE}.csv      — drug, distance, num_disease_genes_hit, rank
  2. disease_proximity_{HANDLE}.txt  — statistici și interpretări
  3. drug_proximity_heatmap_{HANDLE}.png — vizualizare
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Set, List, Tuple
import warnings

import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')

# =====================================================
# CONFIG
# =====================================================
HANDLE = "olezhka0809"

# ✅ Use absolute paths
GRAPH_DRUG_GENE = Path(f"/workspaces/bioinf-y4-lab/labs/09_repurposing/submissions/{HANDLE}/network_drug_gene_{HANDLE}.gpickle")
DRUG_GENE_CSV = Path(f"/workspaces/bioinf-y4-lab/data/work/{HANDLE}/lab09/drug_gene_{HANDLE}.csv")
DISEASE_GENES_TXT = Path(f"/workspaces/bioinf-y4-lab/data/work/{HANDLE}/lab09/disease_genes_{HANDLE}.txt")

# Output paths
OUT_DIR = Path(f"/workspaces/bioinf-y4-lab/labs/09_repurposing/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_DRUG_PRIORITY = OUT_DIR / f"drug_priority_{HANDLE}.csv"
OUT_PROXIMITY_STATS = OUT_DIR / f"disease_proximity_{HANDLE}.txt"
OUT_HEATMAP = OUT_DIR / f"drug_proximity_heatmap_{HANDLE}.png"

# Parameters
MAX_DIST = 10  # Penalize disconnected drug-gene pairs
MIN_PROXIMITY_DRUGS = 20  # Show top N drugs in report

# =====================================================
# UTILITY FUNCTIONS
# =====================================================
def ensure_exists(path: Path) -> None:
    """Verify file exists."""
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")
    print(f"✓ Found: {path}")

def load_bipartite_graph_or_build() -> nx.Graph:
    """
    Load bipartite graph from pickle, or rebuild from CSV if not available.
    """
    if GRAPH_DRUG_GENE.exists():
        print(f"  Loading from pickle...")
        import pickle
        with open(GRAPH_DRUG_GENE, 'rb') as f:
            B = pickle.load(f)
        return B
    
    print(f"  Rebuilding from CSV...")
    df = pd.read_csv(DRUG_GENE_CSV)
    B = nx.Graph()
    
    for _, row in df.iterrows():
        drug = row['drug']
        gene = row['gene']
        B.add_node(drug, bipartite="drug")
        B.add_node(gene, bipartite="gene")
        B.add_edge(drug, gene)
    
    return B

def load_disease_genes(path: Path) -> Set[str]:
    """
    Load disease genes from TXT file (one per line, skip comments).
    """
    genes = set()
    with open(path, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                genes.add(line)
    return genes

def get_drug_nodes(B: nx.Graph) -> List[str]:
    """Extract all drug nodes from bipartite graph."""
    drugs = [n for n, d in B.nodes(data=True) if d.get("bipartite") == "drug"]
    return sorted(drugs)

def get_gene_nodes(B: nx.Graph) -> List[str]:
    """Extract all gene nodes from bipartite graph."""
    genes = [n for n, d in B.nodes(data=True) if d.get("bipartite") == "gene"]
    return sorted(genes)

def compute_shortest_path_distances(
    B: nx.Graph,
    source: str,
    targets: Set[str],
    max_dist: int = 10
) -> Dict[str, float]:
    """
    Compute shortest path distance from source node to each target.
    
    Returns:
      Dict[target] -> distance or max_dist+1 if unreachable
    """
    distances = {}
    for target in targets:
        if target not in B:
            distances[target] = max_dist + 1
            continue
        
        try:
            dist = nx.shortest_path_length(B, source, target)
            distances[target] = dist
        except nx.NetworkXNoPath:
            distances[target] = max_dist + 1
    
    return distances

def compute_drug_disease_distance(
    B: nx.Graph,
    drug: str,
    disease_genes: Set[str],
    mode: str = "mean",
    max_dist: int = 10,
) -> Tuple[float, int]:
    """
    Compute network proximity of a drug to disease genes.
    
    Returns:
      (proximity_score, num_disease_genes_hit)
    
    Proximity score is the average or minimum distance to disease genes.
    Lower is better (closer to disease genes).
    """
    if drug not in B:
        return max_dist + 1, 0
    
    # Filter disease genes that exist in graph
    disease_genes_in_graph = disease_genes & set(B.nodes())
    if not disease_genes_in_graph:
        return max_dist + 1, 0
    
    # Compute distances
    distances = compute_shortest_path_distances(B, drug, disease_genes_in_graph, max_dist)
    dist_values = list(distances.values())
    
    # Count how many disease genes we hit (distance <= max_dist)
    num_hits = sum(1 for d in dist_values if d <= max_dist)
    
    # Compute proximity metric
    if mode == "mean":
        proximity = np.mean(dist_values)
    elif mode == "min":
        proximity = np.min(dist_values)
    elif mode == "median":
        proximity = np.median(dist_values)
    else:
        raise ValueError(f"Unknown mode: {mode}")
    
    return proximity, num_hits

def rank_drugs_by_proximity(
    B: nx.Graph,
    disease_genes: Set[str],
    mode: str = "mean",
    max_dist: int = 10,
) -> pd.DataFrame:
    """
    Rank all drugs by network proximity to disease genes.
    
    Returns:
      DataFrame with columns: drug, proximity_distance, num_disease_genes_hit, rank
    """
    drugs = get_drug_nodes(B)
    results = []
    
    for drug in drugs:
        proximity, num_hits = compute_drug_disease_distance(
            B, drug, disease_genes, mode=mode, max_dist=max_dist
        )
        results.append({
            "drug": drug,
            "proximity_distance": proximity,
            "num_disease_genes_hit": num_hits,
        })
    
    df = pd.DataFrame(results)
    
    # Rank: lower distance = better (rank 1)
    df = df.sort_values("proximity_distance", ascending=True).reset_index(drop=True)
    df['rank'] = df.index + 1
    
    return df

# =====================================================
# ANALYSIS & REPORTING
# =====================================================
def generate_report(
    df_priority: pd.DataFrame,
    disease_genes: Set[str],
    B: nx.Graph,
    output_path: Path
) -> None:
    """Generate text report with analysis."""
    
    genes_in_graph = disease_genes & set(B.nodes())
    drugs = get_drug_nodes(B)
    
    with open(output_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("NETWORK-BASED DRUG REPURPOSING ANALYSIS — ALZHEIMER'S DISEASE\n")
        f.write("=" * 80 + "\n\n")
        
        # Dataset statistics
        f.write("1. DATASET STATISTICS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Total drugs: {len(drugs)}\n")
        f.write(f"Total genes: {len(get_gene_nodes(B))}\n")
        f.write(f"Disease genes: {len(disease_genes)}\n")
        f.write(f"Disease genes in network: {len(genes_in_graph)}\n")
        f.write(f"Network edges: {B.number_of_edges()}\n\n")
        
        # Proximity metrics
        f.write("2. PROXIMITY METRICS\n")
        f.write("-" * 80 + "\n")
        f.write(f"Mean proximity distance: {df_priority['proximity_distance'].mean():.3f}\n")
        f.write(f"Median proximity distance: {df_priority['proximity_distance'].median():.3f}\n")
        f.write(f"Std dev proximity distance: {df_priority['proximity_distance'].std():.3f}\n")
        f.write(f"Min proximity distance: {df_priority['proximity_distance'].min():.3f}\n")
        f.write(f"Max proximity distance: {df_priority['proximity_distance'].max():.3f}\n\n")
        
        # Disease genes hit
        f.write("3. DISEASE GENE TARGETING\n")
        f.write("-" * 80 + "\n")
        f.write(f"Drugs with ≥1 disease gene hit: {(df_priority['num_disease_genes_hit'] > 0).sum()}\n")
        f.write(f"Mean disease genes hit per drug: {df_priority['num_disease_genes_hit'].mean():.2f}\n")
        f.write(f"Max disease genes hit: {df_priority['num_disease_genes_hit'].max()}\n\n")
        
        # Top repurposing candidates
        f.write("4. TOP DRUG REPURPOSING CANDIDATES (by network proximity)\n")
        f.write("-" * 80 + "\n")
        top_drugs = df_priority.head(MIN_PROXIMITY_DRUGS)
        for idx, row in top_drugs.iterrows():
            f.write(f"  Rank {row['rank']}: {row['drug']:<40} ")
            f.write(f"Distance: {row['proximity_distance']:.3f}  ")
            f.write(f"Hits: {int(row['num_disease_genes_hit'])}\n")
        
        f.write("\n5. INTERPRETATION\n")
        f.write("-" * 80 + "\n")
        f.write("""
Network proximity scoring:
  - Lower distance = drug targets genes closer to disease genes in the network
  - Smaller hops = potentially more direct mechanism of action
  
Disease gene hits:
  - Count of disease genes directly reachable from the drug
  - Higher = more direct targeting of AD pathology genes
  
Top candidates are promising for experimental validation because:
  - They're topologically close to known AD risk genes
  - They may work through mechanisms related to AD pathology
  - They can be tested using cell/animal models of AD
        """)

# =====================================================
# VISUALIZATION
# =====================================================
def create_proximity_heatmap(df_priority: pd.DataFrame, output_path: Path, top_n: int = 30) -> None:
    """
    Create heatmap showing top drugs and their characteristics.
    """
    print(f"\n[VIZ] Creating proximity heatmap...")
    
    top_drugs_df = df_priority.head(top_n).copy()
    
    # Normalize metrics for heatmap
    top_drugs_df['proximity_normalized'] = (
        1 - (top_drugs_df['proximity_distance'] / top_drugs_df['proximity_distance'].max())
    )
    top_drugs_df['hits_normalized'] = (
        top_drugs_df['num_disease_genes_hit'] / top_drugs_df['num_disease_genes_hit'].max()
    )
    
    # Create heatmap data
    heatmap_data = top_drugs_df[['drug', 'proximity_normalized', 'hits_normalized']].copy()
    heatmap_data = heatmap_data.set_index('drug')
    
    fig, ax = plt.subplots(figsize=(10, 12), dpi=150)
    
    sns.heatmap(
        heatmap_data,
        annot=True,
        fmt='.2f',
        cmap='RdYlGn',
        cbar_kws={'label': 'Score (normalized)'},
        ax=ax,
        vmin=0,
        vmax=1
    )
    
    ax.set_title(f"Top {top_n} Drugs by Network Proximity (Alzheimer's Disease)\nGreen=Better", 
                 fontsize=12, fontweight='bold')
    ax.set_xlabel("Metric", fontsize=11)
    ax.set_ylabel("Drug", fontsize=11)
    ax.set_xticklabels(['Proximity Score', 'Disease Gene Hits'], rotation=45)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"✓ Saved visualization: {output_path}")
    plt.close()

# =====================================================
# MAIN
# =====================================================
def main():
    print("=" * 70)
    print("EXERCISE 9.2 — NETWORK PROXIMITY & DRUG RANKING")
    print("=" * 70)
    
    # Step 1: Verify inputs
    print("\n[1] Verifying input files...")
    ensure_exists(DISEASE_GENES_TXT)
    
    # Step 2: Load graph
    print("\n[2] Loading bipartite graph...")
    B = load_bipartite_graph_or_build()
    print(f"  Nodes: {B.number_of_nodes()}")
    print(f"  Edges: {B.number_of_edges()}")
    
    # Step 3: Load disease genes
    print("\n[3] Loading disease genes...")
    disease_genes = load_disease_genes(DISEASE_GENES_TXT)
    genes_in_graph = disease_genes & set(B.nodes())
    print(f"  Total disease genes: {len(disease_genes)}")
    print(f"  Disease genes in network: {len(genes_in_graph)}")
    print(f"  Coverage: {100 * len(genes_in_graph) / len(disease_genes):.1f}%")
    
    # Step 4: Compute drug ranking
    print("\n[4] Computing drug-disease proximity...")
    df_priority = rank_drugs_by_proximity(B, disease_genes, mode="mean", max_dist=MAX_DIST)
    print(f"  Drugs ranked: {len(df_priority)}")
    
    # Step 5: Export results
    print("\n[5] Exporting results...")
    df_priority.to_csv(OUT_DRUG_PRIORITY, index=False)
    print(f"✓ Saved: {OUT_DRUG_PRIORITY}")
    
    print(f"\nTop 15 repurposing candidates:")
    print(df_priority[['rank', 'drug', 'proximity_distance', 'num_disease_genes_hit']].head(15).to_string(index=False))
    
    # Step 6: Generate report
    print("\n[6] Generating analysis report...")
    generate_report(df_priority, disease_genes, B, OUT_PROXIMITY_STATS)
    print(f"✓ Saved: {OUT_PROXIMITY_STATS}")
    
    # Step 7: Visualization
    print("\n[7] Creating visualizations...")
    create_proximity_heatmap(df_priority, OUT_HEATMAP, top_n=30)
    
    print("\n" + "=" * 70)
    print("✓ EXERCISE 9.2 COMPLETED")
    print("=" * 70)
    print(f"\nOutput files:")
    print(f"  1. {OUT_DRUG_PRIORITY}")
    print(f"  2. {OUT_PROXIMITY_STATS}")
    print(f"  3. {OUT_HEATMAP}")

if __name__ == "__main__":
    main()