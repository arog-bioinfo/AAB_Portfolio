
def is_in_tuple_list(tl, val):
    for (x, _) in tl:
        if val == x:
            return True
    return False

class MyGraph:
    def __init__(self, g={}):
        self.graph = g

    def print_graph(self):
        for v in self.graph.keys():
            print(v, " -> ", self.graph[v])

    def get_nodes(self):
        return list(self.graph.keys())

    def get_edges(self):
        edges = []
        for v in self.graph.keys():
            for d in self.graph[v]:
                edges.append((v, d))
        return edges

    def size(self):
        return len(self.get_nodes()), len(self.get_edges())

    def add_vertex(self, v):
        if v not in self.graph:
            self.graph[v] = []

    def add_edge(self, o, d):
        if o not in self.graph:
            self.add_vertex(o)
        if d not in self.graph:
            self.add_vertex(d)
        if d not in self.graph[o]:
            self.graph[o].append(d)

    def get_successors(self, v):
        return list(self.graph[v])

    def get_predecessors(self, v):
        res = []
        for k in self.graph:
            if v in self.graph[k]:
                res.append(k)
        return res

    def get_adjacents(self, v):
        suc = self.get_successors(v)
        pred = self.get_predecessors(v)
        res = pred[:]
        for p in suc:
            if p not in res:
                res.append(p)
        return res

    def out_degree(self, v):
        return len(self.graph[v])

    def in_degree(self, v):
        return len(self.get_predecessors(v))

    def degree(self, v):
        return len(self.get_adjacents(v))

    def all_degrees(self, deg_type="inout"):
        degs = {}
        for v in self.graph:
            if deg_type == "out" or deg_type == "inout":
                degs[v] = len(self.graph[v])
            else:
                degs[v] = 0
        if deg_type == "in" or deg_type == "inout":
            for v in self.graph:
                for d in self.graph[v]:
                    if deg_type == "in" or v not in self.graph[d]:
                        degs[d] = degs.get(d, 0) + 1
        return degs

    def highest_degrees(self, all_deg=None, deg_type="inout", top=10):
        if all_deg is None:
            all_deg = self.all_degrees(deg_type)
        ord_deg = sorted(all_deg.items(), key=lambda x: x[1], reverse=True)
        return [x[0] for x in ord_deg[:top]]

    def mean_degree(self, deg_type="inout"):
        degs = self.all_degrees(deg_type)
        return sum(degs.values()) / float(len(degs))

    def prob_degree(self, deg_type="inout"):
        degs = self.all_degrees(deg_type)
        res = {}
        for k in degs:
            res[degs[k]] = res.get(degs[k], 0) + 1
        for k in res:
            res[k] /= float(len(degs))
        return res

    def reachable_bfs(self, v):
        l = [v]
        res = []
        while l:
            node = l.pop(0)
            if node != v:
                res.append(node)
            for elem in self.graph[node]:
                if elem not in res and elem not in l and elem != node:
                    l.append(elem)
        return res

    def reachable_dfs(self, v):
        l = [v]
        res = []
        while l:
            node = l.pop(0)
            if node != v:
                res.append(node)
            s = 0
            for elem in self.graph[node]:
                if elem not in res and elem not in l:
                    l.insert(s, elem)
                    s += 1
        return res

    def distance(self, s, d):
        if s == d:
            return 0
        l = [(s, 0)]
        visited = [s]
        while l:
            node, dist = l.pop(0)
            for elem in self.graph[node]:
                if elem == d:
                    return dist + 1
                elif elem not in visited:
                    l.append((elem, dist + 1))
                    visited.append(elem)
        return None

    def shortest_path(self, s, d):
        if s == d:
            return 0
        l = [(s, [])]
        visited = [s]
        while l:
            node, preds = l.pop(0)
            for elem in self.graph[node]:
                if elem == d:
                    return preds + [node, elem]
                elif elem not in visited:
                    l.append((elem, preds + [node]))
                    visited.append(elem)
        return None

    def reachable_with_dist(self, s):
        res = []
        l = [(s, 0)]
        while l:
            node, dist = l.pop(0)
            if node != s:
                res.append((node, dist))
            for elem in self.graph[node]:
                if not is_in_tuple_list(l, elem) and not is_in_tuple_list(res, elem):
                    l.append((elem, dist + 1))
        return res

    def mean_distances(self):
        tot = 0
        num_reachable = 0
        for k in self.graph:
            distsk = self.reachable_with_dist(k)
            for _, dist in distsk:
                tot += dist
            num_reachable += len(distsk)
        meandist = float(tot) / num_reachable
        n = len(self.get_nodes())
        return meandist, float(num_reachable) / ((n - 1) * n)

    def closeness_centrality(self, node):
        dist = self.reachable_with_dist(node)
        if len(dist) == 0:
            return 0.0
        s = sum(d[1] for d in dist)
        return len(dist) / s

    def highest_closeness(self, top=10):
        cc = {k: self.closeness_centrality(k) for k in self.graph}
        ord_cl = sorted(cc.items(), key=lambda x: x[1], reverse=True)
        return [x[0] for x in ord_cl[:top]]

    def betweenness_centrality(self, node):
        total_sp = 0
        sps_with_node = 0
        for s in self.graph:
            for t in self.graph:
                if s != t and s != node and t != node:
                    sp = self.shortest_path(s, t)
                    if sp:
                        total_sp += 1
                        if node in sp:
                            sps_with_node += 1
        return sps_with_node / total_sp if total_sp > 0 else 0.0

    def node_has_cycle(self, v):
        l = [v]
        visited = [v]
        while l:
            node = l.pop(0)
            for elem in self.graph[node]:
                if elem == v:
                    return True
                elif elem not in visited:
                    l.append(elem)
                    visited.append(elem)
        return False

    def has_cycle(self):
        return any(self.node_has_cycle(v) for v in self.graph)

    def clustering_coef(self, v):
        adjs = self.get_adjacents(v)
        if len(adjs) <= 1:
            return 0.0
        ligs = 0
        for i in adjs:
            for j in adjs:
                if i != j and (j in self.graph[i] or i in self.graph[j]):
                    ligs += 1
        return float(ligs) / (len(adjs) * (len(adjs) - 1))

    def all_clustering_coefs(self):
        return {k: self.clustering_coef(k) for k in self.graph}

    def mean_clustering_coef(self):
        ccs = self.all_clustering_coefs()
        return sum(ccs.values()) / float(len(ccs)) if ccs else 0.0

    def mean_clustering_perdegree(self, deg_type="inout"):
        degs = self.all_degrees(deg_type)
        ccs = self.all_clustering_coefs()
        degs_k = {}
        for k in degs:
            degs_k.setdefault(degs[k], []).append(k)
        ck = {}
        for k in degs_k:
            tot = sum(ccs[v] for v in degs_k[k])
            ck[k] = tot / len(degs_k[k])
        return ck