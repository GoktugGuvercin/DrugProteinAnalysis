from utils import ProteinDB

fasta_dir = "./human_proteome_reviewed.fasta"
protein_db = ProteinDB("FASTA", "", fasta_dir)
protein_db.query_to_uniprot("A1L4Q6")

# search by uniprot id
uniprot_id = "Q7LG56"
entry = protein_db.search_uniprot_id(uniprot_id)
print("\nUniprot Id: ", uniprot_id)
print(entry[:-1])
print(entry[-1])
print()

# search by gene name
proteins = protein_db.search_gene("AKAP7")
print(proteins)
print()

uniprot_id = "Q95330"  # rabbit tp53
entry = protein_db.query_to_uniprot(uniprot_id)
print(entry)
print()

entry = protein_db.search_uniprot_id(uniprot_id)
print(entry)
print()
