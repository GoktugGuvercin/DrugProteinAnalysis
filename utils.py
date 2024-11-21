from Bio import SeqIO


def is_string(name: str) -> bool:
    """ checks if input argument is a string object or not

    Args:
      name: uniprot id or gene name in str format
    """
    if isinstance(name, str):
        return True
    return False


def parse_fasta_header(title: str) -> tuple:
    """ Derive protein description and gene name from fasta title.

    * Fasta titles contain uniprot entry, protein description, organism
    species (OS), organism taxonomoy (OX), protein existence (PE), and
    sequence version (SV). This function parses them and fetch protein
    explanation with respective gene name.

    Args:
      title: title of a fasta file in string format

    Returns: a tuple of protein description and gene name
    """
    content = " ".join(title.split(" ")[1:])
    os_index = content.find("OS=")  # complexity = O(3N)
    gn_index = content.find("GN=")  # complexity = O(3N)

    gene = content[gn_index + 3:].split(" ")[0]
    gene = "##" if len(gene) == 1 else gene
    description = content[: (os_index - 1)]

    return description, gene


def read_fasta(input_file: str) -> tuple:

    """ Returns all prot-ids, sequences, genes, prot-descriptions in fasta.

    * A fasta file can contain more than one protein entry. Each entry
    contains specific fasta header, and protein sequence. This function
    reads all the context and decomposes it into uniprot ids, sequences,
    genes, and descriptions.

    Args:
      input_file: the directory of fasta file in str format

    Returns: A tuple of uniprot ids, sequences, genes and descriptions in
      list format
    """

    uniprot_ids, fasta_seqs = [], []
    explanations, genes = [], []
    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')

    for fasta in fasta_sequences:
        uniprot_ids.append(fasta.id.split("|")[1])
        fasta_seqs.append(str(fasta.seq))

        exp, gene = parse_fasta_header(fasta.description)
        explanations.append(exp)
        genes.append(gene)

    return uniprot_ids, fasta_seqs, genes, explanations,


class ProteinDB:

    def __init__(
            self,
            huproteome_reviewed_dir: str,
    ) -> None:

        """ Defines protein database.

        Args:
          huproteome_reviewed_dir: the directory of human proteins fasta file [1]

        [1]: https://www.uniprot.org/proteomes/UP000005640
        """

        contents = read_fasta(huproteome_reviewed_dir)
        prot_ids, seqs, genes, explanations = contents

        self.pid_db = {}
        self.gene_db: dict[str, list] = {gene: [] for gene in genes}

        for pid, seq, gene, exp in zip(*contents):
            self.pid_db[pid] = (seq, gene, exp)
            self.gene_db[gene].append((pid, seq, exp))

    def search_uniprot_id(self, uniprot_id: str) -> tuple:

        """It returns sequence, gene, and description of target protein.

        * The function search for uniprot id given as input argument in
        human proteome dictionary. If it exists, it returns respective
        sequence, gene name, and description context.

        Args:
          uniprot_id: human uniprot id in string format.

        Returns: a tuple of sequence, gene, protein description.
        """

        is_string(uniprot_id)
        uniprot_id = uniprot_id.upper()

        if uniprot_id not in self.pid_db:
            raise ValueError("Invalid uniprot id for reviewed human proteome")
        return self.pid_db[uniprot_id]

    def search_gene(self, gene: str) -> list[tuple]:

        """It returns uniprot id, sequence, and description of target protein.

        * The function search for gene name given as input argument in
        human proteome dictionary. If it exists, it returns respective
        uniprot ids, sequences, and description context. One gene name can
        refer to multiple uniprot ids, because multiple protein isoforms can
        be synthesized from one gene by alternative splicing such as AKAP7.

        Args:
          gene: gene name in string format.

        Returns: a tuple of three lists (sequences, genes, protein descriptions).
        """

        is_string(gene)
        gene = gene.upper()
        if gene not in self.gene_db:
            raise ValueError("Invalid gene name for reviewed human proteome")
        return self.gene_db[gene]
