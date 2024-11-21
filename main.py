from utils import ProteinDB

fasta_dir = "./human_proteome_reviewed.fasta"
protein_db = ProteinDB(fasta_dir)

# search by uniprot id
seq, gene, exp = protein_db.search_uniprot_id("Q7LG56")
print("Protein Q7LG56: ")
print("Sequence: ", seq)
print("Gene name: ", gene)
print("Description: ", exp)
print()

# search by gene name
proteins = protein_db.search_gene("AKAP7")
for protein in proteins:
    print("Uni-prot id: ", protein[0])
    print("Sequence: ", protein[1])
    print("Description: ", protein[2])
    print()
