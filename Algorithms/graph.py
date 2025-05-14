from collections import deque
import heapq
import graphviz

#=================================================================
#---------------Oriented Graph and Weightned Graph----------------
#=================================================================

def is_in_tuple_list(tl, val):
    """
    Checks whether a value exists as the first element in any tuple in a list.

    :param tl: List of tuples.
    :param val: Value to search for.
    :return: True if val is found as the first element in any tuple, False otherwise.
    """
    for (x, _) in tl:
        if val == x:
            return True
    return False


class Graph:
    """
    Class representing a directed, unweighted graph using adjacency lists.
    """

    def __init__(self, g=None):
        """
        Initializes the graph. If a dictionary `g` is provided, adds its vertices and edges.

        :param g: Optional dictionary representing the graph (e.g., {'A': ['B', 'C']}).
        """
        self.graph = {}
        if g:
            for v, neighbors in g.items():
                self.add_vertex(v)
                for d in neighbors:
                    self.add_edge(v, d)

    def add_vertex(self, v):
        """
        Adds a vertex to the graph if not already present.

        :param v: Vertex label.
        """
        if v not in self.graph:
            self.graph[v] = []

    def add_edge(self, o, d):
        """
        Adds a directed edge from vertex `o` to vertex `d`.

        :param o: Origin vertex.
        :param d: Destination vertex.
        """
        self.add_vertex(o)
        self.add_vertex(d)
        if d not in self.graph[o]:
            self.graph[o].append(d)

    def print_graph(self):
        """
        Prints the adjacency list representation of the graph.
        """
        for v in self.graph:
            print(v, "->", self.graph[v])

    def size(self):
        """
        Returns the number of vertices and edges in the graph.

        :return: Tuple (number of nodes, number of edges).
        """
        return len(self.get_nodes()), len(self.get_edges())

    def get_nodes(self):
        """
        Returns a list of all vertices in the graph.

        :return: List of vertex labels.
        """
        return list(self.graph.keys())

    def get_edges(self):
        """
        Returns a list of all edges in the graph as (origin, destination) tuples.

        :return: List of edge tuples.
        """
        edges = []
        for v in self.graph:
            for d in self.graph[v]:
                edges.append((v, d))
        return edges

    def get_successors(self, v):
        """
        Returns the outgoing neighbors (successors) of a vertex.

        :param v: Vertex label.
        :return: List of successors.
        """
        return self.graph.get(v, [])

    def get_predecessors(self, v):
        """
        Returns the incoming neighbors (predecessors) of a vertex.

        :param v: Vertex label.
        :return: List of predecessors.
        """
        return [u for u in self.graph if v in self.graph[u]]

    def get_adjacents(self, v):
        """
        Returns all adjacent vertices (union of predecessors and successors).

        :param v: Vertex label.
        :return: List of adjacent vertices.
        """
        return list(set(self.get_successors(v)) | set(self.get_predecessors(v)))

    def in_degree(self, v):
        """
        Returns the number of incoming edges for a vertex.

        :param v: Vertex label.
        :return: In-degree.
        """
        return len(self.get_predecessors(v))

    def out_degree(self, v):
        """
        Returns the number of outgoing edges for a vertex.

        :param v: Vertex label.
        :return: Out-degree.
        """
        return len(self.get_successors(v))

    def degree(self, v):
        """
        Returns the total degree (in + out) of a vertex.

        :param v: Vertex label.
        :return: Total degree.
        """
        return self.in_degree(v) + self.out_degree(v)

    def reachable_bfs(self, v):
        """
        Returns all vertices reachable from `v` using Breadth-First Search.

        :param v: Starting vertex.
        :return: List of reachable vertices.
        """
        queue = [v]
        visited = set()
        while queue:
            node = queue.pop(0)
            for neighbor in self.get_successors(node):
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
        return list(visited)

    def reachable_dfs(self, v):
        """
        Returns all vertices reachable from `v` using Depth-First Search.

        :param v: Starting vertex.
        :return: List of reachable vertices.
        """
        stack = [v]
        visited = set()
        while stack:
            node = stack.pop()
            for neighbor in reversed(self.get_successors(node)):
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)
        return list(visited)

    def distance(self, s, d):
        """
        Computes the shortest distance between `s` and `d` using BFS.

        :param s: Source vertex.
        :param d: Destination vertex.
        :return: Distance (int) or None if unreachable.
        """
        if s == d:
            return 0
        queue = [(s, 0)]
        visited = set([s])
        while queue:
            node, dist = queue.pop(0)
            for neighbor in self.get_successors(node):
                if neighbor == d:
                    return dist + 1
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append((neighbor, dist + 1))
        return None

    def shortest_path(self, s, d):
        """
        Finds the shortest path from `s` to `d` using BFS.

        :param s: Source vertex.
        :param d: Destination vertex.
        :return: List of vertices forming the path, or None if no path.
        """
        if s == d:
            return [s]
        queue = [(s, [])]
        visited = set([s])
        while queue:
            node, path = queue.pop(0)
            for neighbor in self.get_successors(node):
                if neighbor == d:
                    return path + [node, neighbor]
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append((neighbor, path + [node]))
        return None

    def reachable_with_dist(self, s):
        """
        Returns all vertices reachable from `s` with their respective distances.

        :param s: Starting vertex.
        :return: List of tuples (vertex, distance).
        """
        res = []
        l = [(s, 0)]
        while l:
            node, dist = l.pop(0)
            if node != s:
                res.append((node, dist))
            for elem in self.get_successors(node):
                if not is_in_tuple_list(l, elem) and not is_in_tuple_list(res, elem):
                    l.append((elem, dist + 1))
        return res

    def node_has_cycle(self, v):
        """
        Checks if there's a cycle starting and ending at vertex `v`.

        :param v: Vertex label.
        :return: True if cycle found, else False.
        """
        queue = [v]
        visited = set([v])
        while queue:
            node = queue.pop(0)
            for neighbor in self.get_successors(node):
                if neighbor == v:
                    return True
                if neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
        return False

    def has_cycle(self):
        """
        Checks if the graph has any cycles.

        :return: True if the graph contains a cycle, else False.
        """
        return any(self.node_has_cycle(v) for v in self.graph)

    def visualize(self):
        """
        Generates a PNG visualization of the graph using Graphviz.
        The output image is saved and opened.
        """
        dot = graphviz.Digraph(comment='Graph', format='png')
        for v in self.graph:
            dot.node(v)
        for v in self.graph:
            for d in self.graph[v]:
                dot.edge(v, d)
        dot.render('output/unweighted_graph', view=True)


class WeightedGraph(Graph):
    """
    Class representing a directed, weighted graph. Inherits from Graph.
    """

    def __init__(self, g=None):
        """
        Initializes the weighted graph. Accepts a dictionary in the form:
        {'A': [('B', weight), ('C', weight)]}.

        :param g: Optional weighted graph representation.
        """
        self.graph = {}
        if g:
            for v, neighbors in g.items():
                self.add_vertex(v)
                for d, w in neighbors:
                    self.add_edge(v, d, w)

    def add_edge(self, o, d, w):
        """
        Adds a weighted edge from vertex `o` to vertex `d`.

        :param o: Origin vertex.
        :param d: Destination vertex.
        :param w: Edge weight (numeric).
        """
        self.add_vertex(o)
        self.add_vertex(d)
        self.graph[o].append((d, w))

    def get_edges(self):
        """
        Returns a list of all edges with weights.

        :return: List of tuples (origin, destination, weight).
        """
        edges = []
        for v in self.graph:
            for d, w in self.graph[v]:
                edges.append((v, d, w))
        return edges

    def get_successors(self, v):
        """
        Returns the outgoing neighbors of a vertex, ignoring weights.

        :param v: Vertex label.
        :return: List of successor vertices.
        """
        return [d for d, _ in self.graph.get(v, [])]

    def dijkstra(self, start):
        """
        Computes shortest paths from `start` to all other vertices using Dijkstra's algorithm.

        :param start: Starting vertex.
        :return: Dictionary mapping each vertex to its shortest distance from `start`.
        """
        distances = {v: float('inf') for v in self.graph}
        distances[start] = 0
        heap = [(0, start)]
        visited = set()

        while heap:
            dist_u, u = heapq.heappop(heap)
            if u in visited:
                continue
            visited.add(u)
            for v, weight in self.graph[u]:
                alt = dist_u + weight
                if alt < distances[v]:
                    distances[v] = alt
                    heapq.heappush(heap, (alt, v))
        return distances

    def visualize(self):
        """
        Generates a PNG visualization of the weighted graph using Graphviz.
        Each edge is labeled with its weight.
        """
        dot = graphviz.Digraph(comment='Weighted Graph', format='png')
        for v in self.graph:
            dot.node(v)
        for v in self.graph:
            for d, w in self.graph[v]:
                dot.edge(v, d, label=str(w))
        dot.render('output/weighted_graph', view=True)
