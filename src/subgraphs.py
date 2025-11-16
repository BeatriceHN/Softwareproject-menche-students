import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import re
from pathlib import Path

# read PPI 
ppi_path = "SpatialPPI/data/consensus_ppi_bioplex_biogrid_intact_huri_edgelist.tsv"

edges = pd.read_csv(
    ppi_path,
    sep="\t",
    header=None,
    usecols=[0, 1],
    names=["u", "v"],
    dtype=str,
    low_memory=False
)
# remove self-loops
edges = edges[edges["u"] != edges["v"]]

# read Localization Data
loc_df = pd.read_csv("location_with_uniprot.tsv", sep="\t", dtype=str)

# get all unique proteins
proteins = loc_df["uniprot_id"].unique().tolist()
locations = loc_df["location"].unique().tolist()

#print(proteins)
#print(locations)

#Mapping dictionaries
loc_to_proteins = (
    loc_df
    .groupby("location")["uniprot_id"]
    .apply(set)
    .to_dict()
)

protein_to_locations = (
    loc_df
    .groupby("uniprot_id")["location"]
    .apply(lambda s: sorted(set(s)))
    .to_dict()
)

#print(loc_to_proteins["Intermediate filaments"])

#NetworkX Graph

G = nx.Graph()
G = nx.from_pandas_edgelist(edges, "u", "v")
N = G.number_of_nodes()
M = G.number_of_edges()
all_nodes = set(G.nodes()) #all nodes in the PPI network

# PPI nodes without localisation
nodes_without_loc = all_nodes.difference(loc_df["uniprot_id"])
#print(nodes_without_loc)
Path("nodes_without_loc.tsv").write_text("\n".join(sorted(nodes_without_loc)))


location_subgraphs = {}

for loc, prots in loc_to_proteins.items():
    #only proteins which are in the PPI network AND localization data
    nodes_here = list(all_nodes.intersection(prots))

    if not nodes_here:
        print(f"no nodes for {loc} location in the PPI")
        continue

    H_loc = G.subgraph(nodes_here).copy()   #extract subgraph via intersecting nodes
    location_subgraphs[loc] = H_loc         #save subgraph in dictionary with location as key

    print(f"{loc} Subnetwork: {H_loc.number_of_nodes()} nodes, {H_loc.number_of_edges()} edges")

#Check if node is in Subnetwork
print("P05981" in location_subgraphs["Plasma membrane"].nodes())
print("P05981" in location_subgraphs["Endoplasmic reticulum"].nodes())  
#should both be true as P05981 is multilocalized

#print(G.edges())

#Search Protein Edges
protein = "Q9UJQ7"
edges_with_protein = list(G.edges(protein))

print(f"Protein {protein} Edges in Network:")
for e in edges_with_protein:
    print(e)

'''
#Plot Subgraphs
MAX_NODES = 50 

for loc, H in location_subgraphs.items():
    if H.number_of_nodes() > MAX_NODES:
        deg = dict(H.degree())
        top_nodes = [
            n for n, _ in sorted(deg.items(), key=lambda kv: kv[1], reverse=True)[:MAX_NODES]
        ]
        H_plot = H.subgraph(top_nodes).copy()
    else:
        H_plot = H

    if H_plot.number_of_nodes() == 0:
        continue

    plt.figure(figsize=(10, 8))
    pos = nx.spring_layout(H_plot, seed=42)

    degH = dict(H_plot.degree())
    sizes = np.array([max(d, 1) for d in degH.values()], dtype=float)
    sizes = 80 * np.sqrt(sizes)

    with_labels = H_plot.number_of_nodes() <= 60

    nx.draw_networkx(
        H_plot,
        pos=pos,
        with_labels=with_labels,
        node_size=sizes,
        font_size=8 if with_labels else 0,
        width=0.6,
    )

    plt.axis("off")
    plt.title(
        f"{loc} subgraph: {H_plot.number_of_nodes()} nodes, "
        f"{H_plot.number_of_edges()} edges"
    )

    out_png = Path(f"/Plots/ppi_{loc.replace(' ', '_')}.png")
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved plot for {loc} to {out_png}")

'''