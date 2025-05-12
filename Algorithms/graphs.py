import heapq
import graphviz

class Graph:
    def __init__(self, g=None):
        self.graph = {}
        if g:
            for vertex, neighbors in g.items():
                self.add_vertex(vertex)
                for neighbor in neighbors:
                    if isinstance(neighbor, tuple):
                        self.add_edge(vertex, neighbor[0], neighbor[1])
                    else:
                        self.add_edge(vertex, neighbor)
            self.id = self.id_graph()
        else:
            self.id = None

    def id_graph(self):
        for v in self.graph:
            for neighbor in self.graph[v]:
                if isinstance(neighbor, tuple):
                    return 'grw'
        return 'gr'

    def add_vertex(self, v):
        if v not in self.graph:
            self.graph[v] = []

    def add_edge(self, o, d, w='bin'):
        self.add_vertex(o)
        self.add_vertex(d)
        if w == 'bin':
            if d not in self.graph[o]:
                self.graph[o].append(d)
        else:
            self.graph[o].append((d, w))

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
            for dest in self.graph[v]:
                if isinstance(dest, tuple):
                    edges.append((v, dest[0], dest[1]))
                else:
                    edges.append((v, dest))
        return edges

    def get_successors(self, v):
        if self.id == 'gr':
            return self.graph[v]
        return [dest[0] for dest in self.graph[v]]

    def get_predecessors(self, v):
        result = []
        for node in self.graph:
            for dest in self.graph[node]:
                if (dest == v) or (isinstance(dest, tuple) and dest[0] == v):
                    result.append(node)
        return result

    def get_adjacents(self, v):
        return list(set(self.get_successors(v)) | set(self.get_predecessors(v)))

    def in_degree(self, v):
        return len(self.get_predecessors(v))

    def out_degree(self, v):
        return len(self.graph[v])

    def degree(self, v):
        return len(self.get_adjacents(v))

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
        for node in self.graph:
            if self.node_has_cycle(node):
                return True
        return False

    def dijkstra(self, start):
        min_heap = []
        visited = set()
        distances = {vertex: float('inf') for vertex in self.graph}
        distances[start] = 0

        heapq.heappush(min_heap, (0, start))

        while min_heap:
            current_distance, current_vertex = heapq.heappop(min_heap)
            if current_vertex in visited:
                continue
            visited.add(current_vertex)
            for neighbor in self.graph[current_vertex]:
                if isinstance(neighbor, tuple):
                    dest, weight = neighbor
                else:
                    dest, weight = neighbor, 1  # unweighted edge default
                distance = current_distance + weight
                if distance < distances[dest]:
                    distances[dest] = distance
                    heapq.heappush(min_heap, (distance, dest))

        return distances

    def visualize(self):
        dot = graphviz.Digraph(comment='Graph', format='png')
        for v in self.graph:
            dot.node(v, v)
        for v in self.graph:
            for neighbor in self.graph[v]:
                if isinstance(neighbor, tuple):
                    d, w = neighbor
                    dot.edge(v, d, label=str(w))
                else:
                    dot.edge(v, neighbor)
        dot.render('output/graph', view=True)


if __name__ == "__main__":
    print("=== Unweighted Graph ===")
    g1 = {
        'A': ['B', 'C'],
        'B': ['D'],
        'C': ['D'],
        'D': []
    }

    G1 = Graph(g1)
    G1.print_graph()
    print("Nodes:", G1.get_nodes())
    print("Edges:", G1.get_edges())
    print("Size:", G1.size())
    print("Successors of A:", G1.get_successors("A"))
    print("Predecessors of B:", G1.get_predecessors("B"))
    print("Adjacents of C:", G1.get_adjacents("C"))
    print("Out-degree of B:", G1.out_degree("B"))
    print("In-degree of C:", G1.in_degree("C"))
    print("Degree of A:", G1.degree("A"))
    print("Shortest path from A to D:", G1.shortest_path('A', 'D'))
    print("Distance from A to D:", G1.distance("A", "D"))
    print("Reachable from A (BFS):", G1.reachable_bfs('A'))
    print("Reachable from A (DFS):", G1.reachable_dfs('A'))
    print("Graph has cycle:", G1.has_cycle())
    print()

    print("=== Weighted Graph ===")
    g2 = {
        'A': [('B', 1), ('C', 4)],
        'B': [('C', 2), ('D', 5)],
        'C': [('D', 1)],
        'D': []
    }

    G2 = Graph(g2)
    G2.print_graph()
    print("Nodes:", G2.get_nodes())
    print("Edges:", G2.get_edges())
    print("Size:", G2.size())
    print("Successors of A:", G2.get_successors("A"))
    print("Predecessors of B:", G2.get_predecessors("B"))
    print("Adjacents of C:", G2.get_adjacents("C"))
    print("Out-degree of B:", G2.out_degree("B"))
    print("In-degree of C:", G2.in_degree("C"))
    print("Degree of A:", G2.degree("A"))
    print("Dijkstra distances from A:", G2.dijkstra('A'))
    print("Reachable from A (DFS):", G2.reachable_dfs('A'))
    print("Graph has cycle:", G2.has_cycle())
