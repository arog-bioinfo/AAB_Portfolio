from Algorithms.graph import Graph

#=================================================================
#-----------------------Assembly Algorythm------------------------
#=================================================================

def composition(k, seq):
    return sorted([seq[i:i+k] for i in range(len(seq)-k+1)])

def prefix(seq): return seq[:-1]
def suffix(seq): return seq[1:]

class DeBruijnGraph(Graph):
    """
    Class representing a De Bruijn graph constructed from k-mers (fragments).
    Inherits from the base Graph class.
    """

    def __init__(self, frags):
        """
        Initializes a De Bruijn graph from a list of k-mers.

        Each k-mer contributes a directed edge from its prefix (first k-1 chars)
        to its suffix (last k-1 chars).

        :param frags: List of k-mer strings.
        """
        super().__init__()
        for seq in frags:
            self.add_edge(prefix(seq), suffix(seq))

    def check_nearly_balanced_graph(self):
        """
        Checks if the graph is nearly balanced (Eulerian path exists).

        A graph is nearly balanced if:
        - One vertex has out-degree = in-degree + 1 (start node)
        - One vertex has in-degree = out-degree + 1 (end node)
        - All other nodes have equal in-degree and out-degree

        :return: Tuple (start_node, end_node) if nearly balanced, else (None, None)
        """
        res = None, None  # (start, end)
        for n in self.graph:
            indeg = self.in_degree(n)
            outdeg = self.out_degree(n)
            if indeg - outdeg == 1:
                res = res[0], n  # candidate for end node
            elif outdeg - indeg == 1:
                res = n, res[1]  # candidate for start node
            elif indeg != outdeg:
                return None, None
        return res

    def eulerian_path(self):
        """
        Constructs a Eulerian path using Hierholzer's algorithm.

        Requires the graph to be nearly balanced. Adds a temporary edge to
        make the graph Eulerian, finds the cycle, then removes the temporary edge.

        :return: List of nodes in Eulerian path, or None if not found.
        """
        start, end = self.check_nearly_balanced_graph()
        if not start or not end:
            return None

        self.add_edge(end, start)  # temporary edge

        path = []
        stack = [start]
        edges = {u: list(vs) for u, vs in self.graph.items()}  # deep copy

        while stack:
            u = stack[-1]
            if edges.get(u):
                stack.append(edges[u].pop())
            else:
                path.append(stack.pop())

        path.reverse()

        # Identify and remove the temporary edge
        for i in range(len(path) - 1):
            if path[i] == end and path[i + 1] == start:
                return path[i + 1:] + path[1:i + 1]

        return None

    def seq_from_path(self, path):
        """
        Reconstructs the original sequence from the Eulerian path.

        Each node represents a (k-1)-mer, so the full sequence is constructed
        by appending the last character of each successive node.

        :param path: List of node strings forming a Eulerian path.
        :return: Reconstructed string sequence, or None if path is invalid.
        """
        if not path:
            return None
        return path[0] + ''.join(n[-1] for n in path[1:])


class OverlapGraph(Graph):
    """
    Class representing an overlap graph built from a set of DNA fragments.

    An edge exists from fragment A to fragment B if suffix(A) == prefix(B).
    Handles identical sequences by uniquely labeling them.
    """

    def __init__(self, frags):
        """
        Initializes the OverlapGraph from a list of fragments.

        :param frags: List of sequence fragments (strings).
        """
        super().__init__()
        self.frags = [f"{f}-{i}" for i, f in enumerate(frags, 1)]

        for f in self.frags:
            self.add_vertex(f)

        for f1 in self.frags:
            s1 = suffix(f1.split('-')[0])
            for f2 in self.frags:
                if prefix(f2.split('-')[0]) == s1:
                    self.add_edge(f1, f2)

    def search_hamiltonian_path(self):
        """
        Attempts to find a Hamiltonian path using backtracking.

        A Hamiltonian path visits each vertex exactly once.

        :return: List of vertex labels in path, or None if not found.
        """
        def bt(path):
            if len(path) == len(self.graph):
                return path
            for neighbor in self.graph[path[-1]]:
                if neighbor not in path:
                    res = bt(path + [neighbor])
                    if res:
                        return res
            return None

        for start in self.graph:
            res = bt([start])
            if res:
                return res
        return None

    def get_seq(self, node):
        """
        Extracts the original sequence from a node label.

        :param node: Node label (e.g., 'ATG-3').
        :return: Sequence part before the suffix (e.g., 'ATG').
        """
        return node.split('-')[0]

    def seq_from_path(self, path):
        """
        Reconstructs a DNA sequence from a Hamiltonian path in the overlap graph.

        The full sequence is built by appending the last character of each successive node.

        :param path: List of vertex labels forming a Hamiltonian path.
        :return: Reconstructed sequence, or None if path is invalid.
        """
        if not path:
            return None
        return self.get_seq(path[0]) + ''.join(self.get_seq(n)[-1] for n in path[1:])
