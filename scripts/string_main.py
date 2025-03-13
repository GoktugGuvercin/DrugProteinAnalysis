
from utils import network
from networks import PPI
from protein import ProteinDB
from embedder import ProtTransEmbedder

ppi_graphs = "/Users/goktug/Desktop/Cancer-Research/ppi_graphs"
prot_db_dir = "./data/human_proteome_reviewed.tsv"

# creating PPI and ProteinDB objects
ppi = PPI()
protein_db = ProteinDB(prot_db_dir)

# ProtTrans Embeddings
# ====================
mras = protein_db.search_gene("MRAS")
shoc2 = protein_db.search_gene("SHOC2")

mras_seq = mras["Sequence"].values[0]
shoc2_seq = mras["Sequence"].values[0]

prot_trans = ProtTransEmbedder(device="cpu")
embeds = prot_trans.compute_res_embeds([mras_seq, shoc2_seq])
print(embeds.shape)

# PPI Graph of Phosphatase Protein Complex
# ========================================
nodes, edges = network(
    method_name="network",
    query_genes=["MRAS", "SHOC2", "PP1C"],
    species=9606,
    network_type="functional",
    confidence=350.0,
    add_color_nodes=5,
)

ppi.add_nodes(nodes)
ppi.add_edges(edges)
ppi.draw_d3_graph("phosphatase complex", "#FFA500", ppi_graphs)
