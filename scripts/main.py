from protein import ProteinDB
from utils import default_fields

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

# search for mouse p53 in our database
uniprot_id = "P02340"
entry = protein_db.query_to_uniprot(uniprot_id, default_fields)
print(entry)
print()

"""# saving genes and pfams in a separate tsv file
save_dir = "./data/genes_pfams.tsv"
protein_db.save_columns(["Entry Name", "Gene", "Pfam"], save_dir, True)"""
