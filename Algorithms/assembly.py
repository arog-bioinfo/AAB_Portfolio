from Algorithms.graph import Graph

#=================================================================
#-----------------------Assembly Algorythm------------------------
#=================================================================

def composition(k, seq):
    return sorted([seq[i:i+k] for i in range(len(seq)-k+1)])

def prefix(seq): return seq[:-1]
def suffix(seq): return seq[1:]

    
class DeBruijnGraph(Graph):
    def __init__(self, frags):
        super().__init__()
        # Add edges where each k-mer contributes an edge from its prefix to its suffix
        for seq in frags:
            self.add_edge(prefix(seq), suffix(seq))

    def check_nearly_balanced_graph(self):
        # Identify if the graph is nearly balanced:
        # One node with out-degree = in-degree + 1 (start)
        # One node with in-degree = out-degree + 1 (end)
        res = None, None  # (start, end)
        for n in self.graph:
            indeg = self.in_degree(n)
            outdeg = self.out_degree(n)
            if indeg - outdeg == 1:
                res = res[0], n  # candidate to be end node
            elif outdeg - indeg == 1:
                res = n, res[1]  # candidate to be start node
            elif indeg != outdeg:
                return None, None  # not balanced or nearly balanced
        return res

    def eulerian_path(self):
        # Find a Eulerian path using Hierholzer's algorithm
        start, end = self.check_nearly_balanced_graph()
        if not start or not end:
            return None

        # Add a temporary edge to make the graph Eulerian
        self.add_edge(end, start)

        path = []
        stack = [start]
        # Copy of graph edges to allow mutation during traversal
        edges = {u: list(vs) for u, vs in self.graph.items()}

        while stack:
            u = stack[-1]
            if edges.get(u):
                stack.append(edges[u].pop())
            else:
                path.append(stack.pop())

        path.reverse()

        # Remove the temporary edge to recover the original path
        for i in range(len(path) - 1):
            if path[i] == end and path[i + 1] == start:
                return path[i + 1:] + path[1:i + 1]

        return None

    def seq_from_path(self, path):
        # Reconstruct the original sequence from a Eulerian path
        if not path:
            return None
        return path[0] + ''.join(n[-1] for n in path[1:])

class OverlapGraph(Graph):
    def __init__(self, frags):
        super().__init__()

        # Add unique suffix to each fragment to handle duplicates
        self.frags = [f"{f}-{i}" for i, f in enumerate(frags, 1)]

        # Add all vertices to the graph
        for f in self.frags:
            self.add_vertex(f)

        # Add edges based on overlap: suffix of f1 matches prefix of f2
        for f1 in self.frags:
            s1 = suffix(f1.split('-')[0])  # suffix of the sequence
            for f2 in self.frags:
                if prefix(f2.split('-')[0]) == s1:
                    self.add_edge(f1, f2)

    def search_hamiltonian_path(self):
        # Try to find a Hamiltonian path using backtracking
        def bt(path):
            if len(path) == len(self.graph):
                return path
            for neighbor in self.graph[path[-1]]:
                if neighbor not in path:
                    res = bt(path + [neighbor])
                    if res:
                        return res
            return None

        # Attempt to start from every vertex
        for start in self.graph:
            res = bt([start])
            if res:
                return res
        return None

    def get_seq(self, node):
        # Extract the original sequence from node label (e.g., 'ATG-3' -> 'ATG')
        return node.split('-')[0]

    def seq_from_path(self, path):
        # Reconstruct sequence from Hamiltonian path
        if not path:
            return None
        return self.get_seq(path[0]) + ''.join(self.get_seq(n)[-1] for n in path[1:])
