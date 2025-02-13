from argparse import Namespace
import numpy as np
from scripts.pfam import Pfam
from scripts.utils import load_embeds_pickle

pfams_dir = "/Users/goktug/Desktop/Cancer-Research/data/genes_pfams.tsv"
embeds_dir = "/Users/goktug/Desktop/Cancer-Research/data"

protbert_embeds = load_embeds_pickle(embeds_dir, "prot_bert_embeds.p")
pfamily = Pfam(pfams_dir, "Gene")

query_genes = list(protbert_embeds.keys())
query_embeds = np.array(list(protbert_embeds.values()))

config = Namespace()
config.n_components = 2
config.perplexity = 50
config.init = "random"
config.learning_rate = "auto"

target_pfams = ["PF00096", "PF00069"]

pfamily.apply_tsne(query_genes, query_embeds, config, target_pfams)
