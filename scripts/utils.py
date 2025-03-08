import os
import pickle
import numpy as np
import requests

string_api_url = "https://version-12-0.string-db.org/api"
base_pfam_url = "https://www.ebi.ac.uk/interpro/api/entry/pfam"
base_query_url = "https://rest.uniprot.org/uniprotkb/search?"

prot_headers = {"Entry": "ID",
                "Entry Name": "Entry Name",
                "Protein names": "Protein Name",
                "Gene Names (primary)": "Gene",
                "Organism": "Organism",
                "Organism (ID)": "Taxonomy",
                "Protein existence": "PE",
                "Sequence version": "SV",
                "Pfam": "Pfam",
                "Sequence": "Sequence"}

default_fields = ["accession", "id", "protein_name", "gene_primary",
                  "organism_name", "organism_id", "protein_existence",
                  "sequence_version", "xref_pfam", "sequence"]


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


def network(
        method_name: str,
        genes: list,
        species: int,
        network_type: str,
        confidence: float,
        add_color_nodes: int,
        unk_prots_dir: str = "",
        output_format: str = "tsv-no-header"
) -> tuple:

    params = {
        "identifiers": "%0d".join(genes),
        "species": species,
        "network_type": network_type,
        "add_nodes": add_color_nodes,
        "required_score": confidence,
        "caller_identity": "www.awesome_app.org"
    }

    # sending a request to STRING API
    request_url = "/".join([string_api_url, output_format, method_name])
    response = requests.post(request_url, data=params)

    # reading unknown proteins into a set
    if unk_prots_dir != "":
        with open(unk_prots_dir, "r") as file:
            unknown_prots = set([line.rstrip() for line in file])
    else:
        unknown_prots = set()

    edges: list[list] = []
    nodes: set[str] = set()

    for line in response.text.strip().split("\n"):
        line = line.strip().split("\t")

        if len(line) == 1 and line[0] == "":
            continue

        node1, node2, edge_score = line[2], line[3], line[5]
        if (node1 not in unknown_prots) and (node2 not in unknown_prots):
            nodes.update([node1, node2])
            edges.append([node1, node2, edge_score])

    return nodes, edges
