from collections import deque
import heapq
import graphviz

#=================================================================
#---------------Oriented Graph and Weightned Graph----------------
#=================================================================

class Graph:
    def __init__(self, g=None):
        self.graph = {}
        if g:
            for v, neighbors in g.items():
                self.add_vertex(v)
                for d in neighbors:
                    self.add_edge(v, d)

    def add_vertex(self, v):
        if v not in self.graph:
            self.graph[v] = []

    def add_edge(self, o, d):
        self.add_vertex(o)
        self.add_vertex(d)
        if d not in self.graph[o]:
            self.graph[o].append(d)

    def print_graph(self):
        for v in self.graph:
            print(v, "->", self.graph[v])

    def size(self):
        return len(self.get_nodes()), len(self.get_edges())

    def get_nodes(self):
        return list(self.graph.keys())

    def get_edges(self):
        edges = []
        for v in self.graph:
            for d in self.graph[v]:
                edges.append((v, d))
        return edges

    def get_successors(self, v):
        return self.graph.get(v, [])

    def get_predecessors(self, v):
        return [u for u in self.graph if v in self.graph[u]]

    def get_adjacents(self, v):
        return list(set(self.get_successors(v)) | set(self.get_predecessors(v)))

    def in_degree(self, v):
        return len(self.get_predecessors(v))

    def out_degree(self, v):
        return len(self.get_successors(v))

    def degree(self, v):
        return self.in_degree(v) + self.out_degree(v)

    def reachable_bfs(self, v):
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

    def node_has_cycle(self, v):
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
        return any(self.node_has_cycle(v) for v in self.graph)

    def visualize(self):
        dot = graphviz.Digraph(comment='Graph', format='png')
        for v in self.graph:
            dot.node(v)
        for v in self.graph:
            for d in self.graph[v]:
                dot.edge(v, d)
        dot.render('output/unweighted_graph', view=True)


class WeightedGraph(Graph):
    def __init__(self, g=None):
        self.graph = {}
        if g:
            for v, neighbors in g.items():
                self.add_vertex(v)
                for d, w in neighbors:
                    self.add_edge(v, d, w)

    def add_edge(self, o, d, w):
        self.add_vertex(o)
        self.add_vertex(d)
        self.graph[o].append((d, w))

    def get_edges(self):
        edges = []
        for v in self.graph:
            for d, w in self.graph[v]:
                edges.append((v, d, w))
        return edges

    def get_successors(self, v):
        return [d for d, _ in self.graph.get(v, [])]

    def dijkstra(self, start):
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
        dot = graphviz.Digraph(comment='Weighted Graph', format='png')
        for v in self.graph:
            dot.node(v)
        for v in self.graph:
            for d, w in self.graph[v]:
                dot.edge(v, d, label=str(w))
        dot.render('output/weighted_graph', view=True)




