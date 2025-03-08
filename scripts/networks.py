import os
import networkx as nx
import matplotlib.pyplot as plt
from d3graph import d3graph, vec2adjmat


class PPI:
    def __init__(self) -> None:

        # Data structures to allocate PPI graph are defined.
        self.node_to_id: dict[str, int] = {}
        self.id_to_node: dict[int, str] = {}
        self.node_edges: set[tuple] = set()
        self.id_edges: set[tuple] = set()
        self.count = 0

    def add_nodes(self, nodes: list[str]) -> None:
        """ It adds the nodes to two dictionaries.
        * Each node refers to a protein.
        * These proteins are represented by names and id numbers.
        * Protein nodes are stored in "name to id" and "id to name" dicts. """

        for node in nodes:
            if node not in self.node_to_id:
                self.node_to_id[node] = self.count
                self.id_to_node[self.count] = node
                self.count += 1

    def add_edge(self, node1: str, node2: str, weight: float) -> None:
        """ It adds the edge to two sets.
        * Each edge is represented by a triplet of two nodes and weight.
        * Two node names and indices are aggregated with edge weight to
        contruct a triplet. These triplets are added to two different sets. 

        Args:
          - node1: the name of node 1
          - node2: the name of node 2
          - weight: connectivity score of two nodes as weight """

        if node1 not in self.node_to_id or node2 not in self.node_to_id:
            raise ValueError("The nodes of that edge are not in the graph")

        edge = tuple(sorted((node1, node2)) + [weight])
        self.node_edges.add(edge)

        id1 = self.node_to_id[edge[0]]
        id2 = self.node_to_id[edge[1]]
        self.id_edges.add((id1, id2, weight))

    def add_edges(self, edges: list) -> None:
        """ It adds the edges to two sets()."""
        for edge in edges:
            n1, n2, score = edge
            self.add_edge(n1, n2, score)

    def list_graph_nodes(self) -> None:
        """It prints all nodes by indices and protein names."""
        for key, value in self.node_to_id.items():
            print(f"Node {key}: {value}")

    def build_nx_graph(self) -> nx.Graph:
        """It builds a networkx graph for visualization."""
        graph = nx.Graph()
        for node1, node2, score in self.node_edges:
            graph.add_edge(node1, node2, weight=score)
        return graph

    def draw_nx_graph(self) -> None:

        """It draws a PPI graph.
        * It builds a networkx graph, and visualize them by protein names. """

        graph = self.build_nx_graph()
        pos = nx.spring_layout(graph, seed=7)

        nx.draw(graph, pos, with_labels=False, node_size=20)
        nx.draw_networkx_labels(
            G=graph,
            pos=pos,
            font_size=7,
            verticalalignment='bottom',
            horizontalalignment='right')

        plt.show()

    def draw_d3_graph(self, name: str, color: str, save_dir: str) -> None:

        """It draws a PPI graph.
        * It builds a d3 graph, and visualize them by protein names.

        Args:
          - name: the name of plot
          - color: the color of nodes, only hex codes
          - save_dir: the directory where ppi plot will be saved."""

        source, target = [], []
        for node1, node2, _ in self.node_edges:
            source.append(node1)
            target.append(node2)

        weight = [1 for i in range(len(source))]
        adjmat = vec2adjmat(source, target, weight=weight)

        colors = [color for _ in range(self.count)]
        file_name = f"{name} proteins"

        d3 = d3graph()
        d3.graph(adjmat)
        d3.set_node_properties(color=colors)
        d3.show(filepath=os.path.join(save_dir, f"{file_name}.html"))
