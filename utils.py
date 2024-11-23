import requests as r
import pandas as pd
from Bio import SeqIO
from io import StringIO

base_url = "http://www.uniprot.org/uniprot/"


def is_string(name: str) -> bool:
    """ checks if input argument is a string object or not

    Args:
      name: uniprot id or gene name in str format
    """
    if isinstance(name, str):
        return True
    return False


def parse_fasta_header(title: str) -> list:
    """ Derive entries (protein name, gene, species etc.) from fasta header.

    * Fasta headers contain uniprot entry, protein name, organism
    species (OS), organism taxonomoy (OX), protein existence (PE), and
    sequence version (SV). This function at first parses them and then
    return in a list.

    Args:
      title: title of a fasta file in string format

    Returns: a list of fasta-header entries
    """

    prot_name = title[:title.find("_")].split("|")[-1]
    entries = {"Protein": prot_name, "OS": "##", "OX": "##",
               "GN": "##", "PE": "##", "SV": "##"}
    parts = title.split()

    for part in parts:
        if "=" in part:
            key, value = part.split("=", 1)
            entries[key] = value

    return list(entries.values())


def read_fasta(input_file: str) -> tuple:

    """ Returns prot-ids, sequences, and headers retrieved from fasta.

    * A fasta file can contain more than one protein entry. Each entry
    contains specific fasta header, uniprot id, protein sequence. This
    function reads all the context and decomposes it into uniprot ids,
    sequences, and headers in list format.

    Args:
      input_file: the directory of fasta file in str format

    Returns: A tuple of uniprot ids, sequences, and headers in
      list format
    """

    uniprot_ids, seqs, headers = [], [], []
    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')

    for fasta in fasta_sequences:
        uniprot_ids.append(fasta.id.split("|")[1])
        seqs.append(str(fasta.seq))
        header = parse_fasta_header(fasta.description)
        headers.append(header)

    return uniprot_ids, seqs, headers


class ProteinDB:

    def __init__(
            self,
            db_choice: str,
            database_csv_dir: str,
            proteome_reviewed_dir: str,
    ) -> None:

        """ Defines protein database.

        Args:
          db_choice: which dir is chosen for db construction (CSV or FASTA)
          database_csv_dir: the directory of csv file containing protein entries
          huproteome_reviewed_dir: the directory of human proteins fasta file [1]

        [1]: https://www.uniprot.org/proteomes/UP000005640
        """

        if db_choice.upper() == "FASTA":
            contents = read_fasta(proteome_reviewed_dir)
            prot_ids, seqs, headers = contents
            prots, os, ox, gn, pe, sv = zip(*headers)

            pid_db = {"ID": prot_ids, "Protein Name": prots,
                      "Gene": gn, "Species": os, "Taxonomy": ox,
                      "PE": pe, "SV": sv, "Sequence": seqs}

            self.database = pd.DataFrame(pid_db).set_index("ID")

        else:
            self.database = pd.read_csv(database_csv_dir)

    def search_uniprot_id(self, uniprot_id: str) -> list:
        """It returns sequence, gene, and description of target protein.

        * The function searches for uniprot id given as input argument in
        the DataFrame. If it exists, the features of target protein is returned
        [protein name, gene, species, taxonomy, PE, SV, sequence].

        Args:
          uniprot_id: uniprot id in string format.

        Returns: a list of corresponding protein features in str.
        """

        is_string(uniprot_id)
        uniprot_id = uniprot_id.upper()

        if uniprot_id not in self.database.index:
            raise ValueError("Invalid uniprot id for proteome database")
        return self.database.loc[uniprot_id].tolist()

    def search_gene(self, gene: str) -> pd.DataFrame:
        """It returns uniprot id, sequence, and description of target protein.

        * The function searches for gene name given as input argument in
        the DataFrame. One gene can refer to more then one uniprot id, because
        multiple protein isoforms can be synthesized from one gene by
        alternative splicing such as AKAP7. Hence, all proteins matched with
        given gene name are returned as a DataFrame object.

        Args:
          gene: gene name in string format.

        Returns: a DataFrame object that can contain multiple protein entries.
        """

        is_string(gene)
        gene = gene.upper()
        if gene not in self.database["Gene"].to_list():
            raise ValueError("Invalid gene name for proteome database")

        matches = (self.database["Gene"] == gene)
        return self.database[matches]

    def query_to_uniprot(self, uniprot_id: str) -> pd.DataFrame:
        """ Fetches a protein entry from uniprot and save it to database.

        * The function sends an API request to UniProt to search for
        particular uniprot id. If it exists, the response of the request
        is converted to fasta and parsed into a new protein entry, which
        is saved into protein database.

        Args:
          uniprot_id: uniprot id in string format.
        """

        # sending query to UniProt
        fasta_url = f"{base_url}{uniprot_id}.fasta"
        response = r.post(fasta_url)

        if response.status_code == 404:
            raise ValueError("Invalid uniprot id is passed as input argument")

        # converting query response to fasta
        fasta_str = StringIO(response.text)
        fasta = list(SeqIO.parse(fasta_str, 'fasta'))[0]

        # parsing fasta content
        seq = str(fasta.seq)
        uniprot_id = fasta.id.split("|")[1]
        header = parse_fasta_header(fasta.description)
        prot, os, ox, gn, pe, sv = header

        # adding new entry
        entry = {"Protein Name": prot, "Gene": gn, "Species": os,
                 "Taxonomy": ox, "PE": pe, "SV": sv, "Sequence": seq}
        df_entry = pd.DataFrame(entry, index=[uniprot_id])
        self.database = pd.concat([self.database, df_entry])

        return df_entry

    def save_to_csv(self, save_dir: str) -> None:
        """Saves protein database into csv file."""
        self.database.to_csv(save_dir)
