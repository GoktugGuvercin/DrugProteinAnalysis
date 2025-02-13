from argparse import Namespace
import numpy as np
import pandas as pd

from sklearn.manifold import TSNE
import plotly.express as px


class Pfam:

    def __init__(self, database_tsv_dir: str, gene_or_prot: str = "Gene") -> None:
        """ Defines protein family database.
        This database is maintained as a DataFrame of
            - protein ids
            - gene names
            - protein families

        Args:
          * database_tsv_dir: the directory of tsv file for pfam db entries.
          * gene_or_prot: a string of "Gene" or "ID" to choose database type:
            - gene-pfam
            - protein-pfam """

        self.database = pd.read_csv(
            database_tsv_dir,
            sep=",",
            dtype=str,
            na_filter=False,
            keep_default_na=False
        )

        print(self.database)

        self.colors = Pfam.yield_color_codes()
        self.gp_pfam, self.pfam_color = self.gp_to_pfam(gene_or_prot)
        self.pfam_gp = Pfam.reverse_gp_pfam(self.gp_pfam)

    def detect_large_pfams(self, size: int = 30) -> list:

        """Detect large protein family groups.

        * Multiple genes or proteins can share same pfam entries.
        In that case, these genes/proteins are collected under that
        pfam group. This function helps to choose pfam groups larger
        than the size given as input.

        Args:
          * size: the size of pfam group as threshold"""

        large_pfams = []
        for pfam, gp_list in self.pfam_gp.items():
            if len(gp_list) >= size:
                large_pfams.append(pfam)
        return large_pfams

    def gp_to_pfam(self, choice: str = "Gene") -> tuple:

        """ Creates a dict of gp-pfam pairs.
        * Either gene names or protein names are used as keys
        to access related protein family ids.

        * Some proteins do not have a certain protein family, which
        is denoted as null and annotated with white color.

        [1] There can be multiple copies (paralogs) of a gene.
        For example, human beta-defensin genes (DEFB104) reside
        in a region of chromosome 8 (8p23.1), and undergoes
        multiple dublications. As a result, multiple prologs of
        this gene (DEFB104A, DEFB104B) occur. They are closely
        related copies, and encode almost same protein sequence.
        These gene names are generally separated by ; sign.

        https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=DetailsSearch&Term=503618
        """

        index = 0
        white = np.array(3 * [255], dtype=np.uint8)

        gene_pfam_dict = {}
        pfam_color_dict = {"null": white}

        gps = self.database[choice].values
        pfams = self.database["Pfam"].values

        for gp, pfam in zip(gps, pfams):
            if gp == "null":
                continue

            # one gene can have multiple synonym names
            pfam = Pfam.sort_pfam(pfam)
            for synon in gp.split(" "):
                if synon.endswith(";"):  # [1]
                    gene_pfam_dict[synon[:-1]] = pfam
                else:
                    gene_pfam_dict[synon] = pfam

            if pfam not in pfam_color_dict:
                pfam_color_dict[pfam] = self.colors[index]
                index += 1

        return gene_pfam_dict, pfam_color_dict

    def apply_tsne(
            self,
            query_genes: list,
            query_embeds: np.ndarray,
            config: Namespace,
            target_pfams: list = [],
    ) -> tuple:

        colors: list[np.ndarray] = []
        awhite = np.array(3 * [245], dtype=np.uint8)

        # getting colors for query genes
        # some query genes may not exist in Pfam database
        for qgene in query_genes:
            if qgene in self.gp_pfam:
                pfam = self.gp_pfam[qgene]
                color = self.pfam_color[pfam]
                colors.append(color)
            else:
                colors.append(awhite)

        # embeddings of query genes are projected onto 2/3D coords
        print("Fitting TSNE ...")
        tsne = TSNE(n_components=config.n_components,
                    perplexity=config.perplexity,
                    init=config.init,
                    learning_rate=config.learning_rate)
        tsne_embeds = tsne.fit_transform(X=query_embeds).T

        # creating a dict of query genes across tsne coords and colors
        xs, ys = tsne_embeds[0], tsne_embeds[1]
        tsne_dict = {g: (x, y) for g, x, y in zip(query_genes, xs, ys)}

        # visualizing entire pfam scatter plot
        print("Visualizing entire pfam scatter")
        hex_colors = Pfam.rgb_to_hex(np.array(colors))
        df1 = pd.DataFrame({"x": xs, "y": ys})
        fig = px.scatter(df1, x="x", y="y")

        fig.update_traces(marker=dict(size=4, symbol="circle", color=hex_colors))
        fig.update_xaxes(range=[-120, 120], dtick=20)
        fig.update_yaxes(range=[-120, 120], dtick=20)
        fig.show()

        # visualizing some part of pfam scatter plot
        print("Visualizing some part of pfam scatter")
        filt_xs, filt_ys, filt_legends = [], [], []
        for tpfam in target_pfams:
            if tpfam in self.pfam_gp:

                tgene = self.pfam_gp[tpfam]  # a pfam: gene or genes
                txyl = [(*tsne_dict[g], self.gp_pfam[g])
                        for g in tgene if g in tsne_dict]
                txs, tys, tlegends = map(list, zip(*txyl))

                filt_xs.extend(txs)
                filt_ys.extend(tys)
                filt_legends.extend(tlegends)
            else:
                print("Given target pfam does not exist")

        color_map = {legend: Pfam.rgb_to_hex(self.pfam_color[legend])[0]
                     for legend in filt_legends}

        df2 = pd.DataFrame({
            "x": filt_xs,
            "y": filt_ys,
            "family": filt_legends}
        )

        fig = px.scatter(df2, x="x", y="y", color="family",
                         color_discrete_map=color_map)

        fig.update_traces(marker=dict(size=4, symbol="circle"))
        fig.update_xaxes(range=[-120, 120], dtick=20)
        fig.update_yaxes(range=[-120, 120], dtick=20)
        fig.show()

        return df1, df2

    @staticmethod
    def rgb_to_hex(rgb_colors: np.ndarray) -> list:

        """ Maps rgb colors to hex string
        Args:
          * rgb_colors: a 2D numpy array of integers in [N, 3] shape"""

        if len(rgb_colors.shape) == 1:
            rgb_colors = np.array([rgb_colors])

        hex_colors = []
        for r, g, b in rgb_colors:
            hex = "#" + ('{:02X}' * 3).format(r, g, b)
            hex_colors.append(hex)
        return hex_colors

    @staticmethod
    def yield_color_codes(max_val: int = 235) -> np.ndarray:
        colors = np.indices(3 * [max_val], dtype=np.uint8)
        colors = colors.reshape(3, -1).T[1:]
        np.random.shuffle(colors)
        return colors

    @staticmethod
    def sort_pfam(pfams: str) -> str:
        """ Sorts multiple pfam entries

        * While creating a dict of pfam-gp pairs, we need to consider
        two details.
            1) Some genes/proteins (gp) can have more than one pfam entry.
            2) Pfam entries of some genes/proteins can be same.

        Example:
        --------
        ESYT2 gene (A0FGR8 protein) --> PF00168;PF17047;
        ESYT3 gene (A0FGR9 protein) --> PF17047;PF00168;

        These two genes/proteins have multiple and same pfam entries. In that case,
        to create a unique pfam id as a key in the dict, we need to sort them.
        In that way, these two genes/proteins would be located in same pfam group:
        <"PF00168;PF17047", [ESYT2, ESYT3]>

        Args:
          * pfams: a string of single/multiple pfam ids"""

        pfams_list = sorted(pfams.split(";")[:-1])
        pfams_key = ";".join(pfams_list)
        return pfams_key

    @staticmethod
    def reverse_gp_pfam(gene_pfam: dict) -> dict:

        """ Reverses gene to pfam dict.

        We might need to choose or cluster the genes/proteins in
        the same pfam group, where we need pfam to gene dict."""

        pfam_gp = {}
        for gene, pfam in gene_pfam.items():
            if pfam not in pfam_gp:
                pfam_gp[pfam] = [gene]
            else:
                pfam_gp[pfam].append(gene)

        return pfam_gp
