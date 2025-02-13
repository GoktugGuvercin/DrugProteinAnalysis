import os
import pickle
import requests as r
import pandas as pd
import numpy as np

default_fields = ["accession", "id", "protein_name", "gene_primary",
                  "organism_name", "organism_id", "protein_existence",
                  "sequence_version", "xref_pfam", "protein_families",
                  "sequence"]


def load_embeds_pickle(source_dir: str, file_name: str) -> dict:
    pickle_dir = os.path.join(source_dir, file_name)
    return pickle.load(open(pickle_dir, "rb"))


def is_string(name: str) -> bool:
    """ checks if input argument is a string object or not

    Args:
      name: uniprot id or gene name in str format
    """
    if isinstance(name, str):
        return True
    return False


def has_pyrrolysine(seq: str) -> bool:
    """Checks the sequence contains pyrrolsine amino acid.

    * Pyrrolysine, encoded by the 'amber' stop codon UAG, is
    an a-amino acid involved in protein biosynthesis in certain
    methanogenic archaea and bacteria. Notably, it is absent in
    humans.
    * Example: Uniprot id: O30642, Prot. Name: MTMB1_METBA, Gene: mtmB1

    Args:
      seq: primary structure of a protein in str format
    """
    if not is_string(seq):
        raise TypeError(f"Invalid input type; should be str, got {type(seq)}")
    return "O" in seq


def has_selenocysteine(seq: str) -> bool:
    """Checks the sequence contains selenocysteine amino acid.

    * Most of human proteins comprise 20 major amino acids.
    Selenocysteine (U), as 21th amino acid, appears in 25 human
    selenoproteins, some of which are involved in antioxidant
    defense and thyroid hormone regulation. In protein decomposition,
    selenocysteine is not considered.

    Args:
      seq: primary structure of a protein in str format
    """
    if not is_string(seq):
        raise TypeError(f"Invalid input type; should be str, got {type(seq)}")
    return "U" in seq


def index_pfams(pfams: list, index_type: str = "int") -> dict:

    if index_type == "color":
        indices = generate_colors()
    else:
        indices = np.arange(0, len(pfams), 1)

    index = 1
    pfam_dict = {"null": indices[0]}

    for pfam in pfams:
        pfam_items = sorted(pfam.split(";")[:-1])
        pfam_key = ";".join(pfam_items)
        if pfam_key not in pfam_dict:
            pfam_dict[pfam_key] = indices[index]
            index += 1

    return pfam_dict


def generate_colors() -> np.ndarray:

    colors = np.indices((256, 256, 256), dtype=np.uint8)
    colors = colors.reshape(3, -1).T
    np.random.shuffle(colors)
    return colors


class ProteinDB:

    def __init__(
            self,
            database_tsv_dir: str,
    ) -> None:

        """ Defines protein database.

        Args:
          database_tsv_dir: the directory of tsv file containing protein entries [1]

        [1]: https://www.uniprot.org/proteomes/UP000005640
        """

        database = pd.read_csv(
            database_tsv_dir,
            sep="\t",
            dtype=str,
            na_filter=False,
            keep_default_na=False
        )

        # database["Gene"] = [gene.split(" ")[0] for gene in database.iloc[:]["Gene"]]
        self.database = pd.DataFrame(database).set_index("ID")

    def search_uniprot_id(self, uniprot_id: str) -> dict:
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
        return self.database.loc[uniprot_id].to_dict()

    def search_gene(self, gene_name: str) -> pd.DataFrame:
        """It returns uniprot id, sequence, and description of target protein.

        * The function searches for gene name given as input argument in
        the DataFrame. One gene can refer to more then one uniprot id, because
        multiple protein isoforms can be synthesized from one gene by
        alternative splicing such as AKAP7. Hence, all proteins matched with
        given gene name are returned as a DataFrame object.

        Args:
          gene_name: gene name in string format.

        Returns: a DataFrame object that can contain multiple protein entries.
        """

        is_string(gene_name)
        gene_name = gene_name.upper()
        genes_list = self.database["Gene"].to_list()

        indices = []
        for i, gene in enumerate(genes_list):
            if gene_name in gene:
                indices.append(i)
        if len(indices) == 0:
            raise ValueError("Invalid gene name for proteome database")

        return self.database.iloc[indices]

    def save_columns(self, columns: list, save_dir: str, index: bool) -> None:

        if all([name in self.database.columns for name in columns]):
            sub_df = self.database[columns].replace("", "null")
            sub_df.to_csv(save_dir, index=index, index_label="ID")
        else:
            raise ValueError("Invalid column names exist")

    def save_to_tsv(self, save_dir: str) -> None:
        """Saves protein database into tsv file."""
        self.database.to_csv(save_dir)

    def __len__(self) -> int:
        return len(self.database)
