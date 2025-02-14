from argparse import Namespace
import numpy as np
from pfam import Pfam
from utils import load_embeds_pickle

pfams_dir = "/Users/goktug/Desktop/Cancer-Research/data/genes_pfams.tsv"
embeds_dir = "/Users/goktug/Desktop/Cancer-Research/data"

protbert_embeds = load_embeds_pickle(embeds_dir, "prot_bert_embeds.p")
pfamily = Pfam(pfams_dir, "Gene")
pfam_groups = pfamily.detect_large_pfams(120)

print("Large Pfam Groups:")
print(pfam_groups)

query_genes = list(protbert_embeds.keys())
query_embeds = np.array(list(protbert_embeds.values()))

config = Namespace()
config.n_components = 2
config.perplexity = 50
config.init = "random"
config.learning_rate = "auto"

target_pfams = ["PF00096;PF01352", "PF00069"]

pfamily.apply_tsne(query_genes, query_embeds, config, target_pfams)
