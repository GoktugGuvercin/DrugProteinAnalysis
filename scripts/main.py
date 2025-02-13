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
