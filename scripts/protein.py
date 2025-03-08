import requests as r
import pandas as pd

from utils import (
    is_string,
    prot_headers,
    base_query_url
)


class ProteinDB:

    def __init__(
            self,
            database_tsv_dir: str,
    ) -> None:

        """ It defines protein database.
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

    def search_uniprot_id(self, uniprot_id: str) -> pd.DataFrame:
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
        return self.database.loc[uniprot_id].to_frame().T

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

    def query_to_uniprot(self, uniprot_id: str, fields: list) -> pd.DataFrame:
        """ Fetches a protein entry from uniprot and save it to database.

        * The function sends a query to UniProt to search for particular
        uniprot id. If it exists, the response of the request
        is converted to fasta and parsed into a new protein entry, which
        is saved into protein database.

        Args:
          uniprot_id: uniprot id in string format.

        Fields: https://www.uniprot.org/help/return_fields
        Crosslink Fields: https://www.uniprot.org/help/return_fields_databases

        """

        # creating query url
        format = "tsv"
        query = f"accession:{uniprot_id}"
        fields_tag = ",".join(fields)
        query_url = f"{base_query_url}query={query}&format={format}&fields={fields_tag}"
        print(query_url)

        # sending query for response
        response = r.post(query_url)
        if response.status_code == 404:
            raise ValueError("Invalid uniprot id is passed as input argument")

        # In query response,
        # headers and entries are located in different lines.
        # there is tab char between entries and headers.
        # we fetch headers and entries separately.
        headers = response.text.split("\n")[0].split("\t")
        headers = [prot_headers[h] for h in headers]  # query headers to our headers
        entries = response.text.split("\n")[1].split("\t")

        # Adding new entry to protein database
        outcome = dict(zip(headers, entries))
        df_entry = pd.DataFrame(outcome, index=[uniprot_id])
        self.database = pd.concat([self.database, df_entry])

        return df_entry

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
