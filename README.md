# Cancer Research

This repository is dedicated to cancer research. It will higlight theoretical knowledge-base together with related deep learning projects. 

1. [Cancer Drugs](https://github.com/GoktugGuvercin/Cancer-Research/blob/main/Drugs.md)
2. [Protein Database](https://github.com/GoktugGuvercin/Cancer-Research/blob/main/utils.py)


## ProteinDB

To construct a protein database, we need to have protein entries in a `.tsv` file format. I specifically opted for human proteome entries, which is provided by [UniProt](https://www.uniprot.org/proteomes/UP000005640) database. The given code block below is standard `main.py` to realize a protein database and search for any protein id or gene name.
Our `ProteinDB` also allows for saving specific columns of its database entries; passing requested column names is enough for this.

```python
from scripts.utils import ProteinDB

data_dir = "./data/human_proteome_reviewed.tsv"
protein_db = ProteinDB(data_dir)

print(len(protein_db))

# search by uniprot id
uniprot_id = "Q7LG56"
entry = protein_db.search_uniprot_id(uniprot_id)
print("\nUniprot Id: ", uniprot_id)
print(entry)
print()

# search by gene name
proteins = protein_db.search_gene("AKAP7")
print(proteins)
print()

# saving genes and pfams in a separate tsv file
save_dir = "./data/genes_pfams.tsv"
protein_db.save_columns(["Entry Name", "Gene", "Pfam"], save_dir, True)

```

## Pfam

Protein families are categorization groups for the structure and functionality of proteins. Our `Pfam` provides an information retrieval system between genes or proteins and related pfam entries. To instantiate `Pfam` database, we can use `genes_pfams.tsv` file, which will construct gene $\rightarrow$ pfam and pfam $\rightarrow$ gene dictionaries, as well as assign color codes to distinct pfam groups. 

By using `ProtTrans` model, which is not added to this repository, we created embedding vectors for thousands of proteins and matched these embeddings with gene names of those proteins in pickle file `prot_bert_embeds.p`. By using these embeddings, we can generate TSNE visualization of the proteins, and we can assign color codes to these proteins in TSNE depending on protein family entries in our `Pfam` database. To achieve that, we can access pfam entries of target proteins by their gene names in `prot_bert_embeds.p` file. 

```python
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

```



## Proteins 

The primary structure of a protein is essentially a sequence of amino acids linked by peptide bonds, forming a polypeptide chain. This sequence is not random but is precisely determined by the corresponding gene via messenger RNA. The diversity of polypeptide chains arises from various combinations of 20 distinct amino acids, though selenocysteine, the 21st amino acid, also appears in 25 human selenoproteins, some of which are involved in antioxidant defense and thyroid hormone regulation.

Within a polypeptide chain, amino acids can form hydrogen bonds, which lay the foundation for the protein's secondary structure. This structure arises when the chain undergoes local folding, typically forming one of two patterns: the α-helix or β-sheet. In an α-helix, hydrogen bonds occur between the carbonyl group (C=O) of one amino acid and the amino group (N-H) of another, four residues further along the chain, creating a helical shape when these bonds follow a regular pattern. In contrast, β-sheets are stabilized by hydrogen bonds between the carbonyl group of one polypeptide segment and the amino group of a parallel segment, aligning multiple chain segments side by side.

Each amino acid has a unique side chain, or R group, which defines its properties, such as polarity, charge, and size. These R groups play a crucial role in folding the polypeptide into its tertiary structure through various interactions. For instance, ionic bonds form between oppositely charged R groups, stabilizing the structure through salt bridges [MB06; DKD11; Kuo+13]. Similarly, hydrogen bonds between polar side chains contribute to stabilization, and disulfide bonds between cysteine residues facilitate folding and structural integrity. Additionally, hydrophobic R groups tend to cluster inside the protein to avoid water, while hydrophilic R groups are typically found on the surface, interacting with water molecules.

The tertiary structure results from the integrity of the primary and secondary structures, producing a single monomeric protein. When multiple monomeric chains interact, they form a quaternary structure, an organized assembly of subunits.

# Protein Existence (PE) in Uniprot

PE stands for “Protein Existence in UniProt. It assigns a value from 1 to 5 to reflects the type of evidence supporting the existence of the protein:

PE=1: Evidence at protein level (direct protein sequencing or mass spectrometry).
PE=2: Evidence at transcript level (mRNA data but no direct protein-level evidence).
PE=3: Inferred from homology (similarity to proteins in other organisms).
PE=4: Predicted (no experimental evidence, no strong homology evidence).
PE=5: Uncertain (very little or contradictory evidence).

# Proteins in Pfam Domain

Pfam is a comprehensive database of protein families. Each protein is classified into one family or multiple families at the same time,
based on the presence of conserved domains. In general each domain defines a family, and qualify the proteins to have some functionalities.

For example, PF00069 refers to protein-kinase domain. It is highly conserved region found in a vast array of proteins. These proteins are involved in various cellular processes. A protein kinase is an enzyme responsible for phosphorylation, that is the process of adding  phosphates to proteins to regulate their functionality and creating a signal through the cell. At this point, structural conformation of
target protein can be changed, and its activity would be modified. In this process,  phosphates are sourced by ATP molecules. 

Protein kinases are categorized in same family group, but this does not mean that their primary amino-acid chains are completely same.
In fact, they have quite different ordering of amino-acids in polypeptide chain with some subtle exceptions of motifs. To phosphorylate other molecules, all kinases have

- ATP-binding motifs: It faciliates ATP binding and positioning for the phosphoryl transfer reaction, known as Walker-A motif or P-loop. It is capable of interacting with the phosphate groups of ATP. These motifs are generaly formed in G-x(4)-GK-\[TS\] pattern

- Catalytic loop motif: It is HRD motif (Histidine-Arginine-Aspartic acid)

- Activation loop motif: It is DFG motif (Aspartic acid-Phenylalanine-Glycine) to regulate kinase activity. 

* https://en.wikipedia.org/wiki/Walker_motifs
* https://en.wikipedia.org/wiki/Protein_kinase
* https://en.wikipedia.org/wiki/Protein_phosphorylation#cite_note-1