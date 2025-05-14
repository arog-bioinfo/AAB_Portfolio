from collections import deque
from Algorithms.graph import Graph
import regex as re
import heapq

#=================================================================
#-----------------------Metabolic Network-------------------------
#=================================================================

class MN_Graph(Graph):
    """
    Metabolite Network Graph extending a base Graph class with additional methods 
    for degree analysis, centrality, clustering, and network metrics.
    """

    def all_degrees(self, deg_type="inout"):
        """
        Compute the degree of each node based on the specified type.

        Parameters:
            deg_type (str): "in", "out", or "inout" (default) indicating degree type.

        Returns:
            dict: Mapping of node to degree count.
        """
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
        """
        Return nodes with the highest degrees.

        Parameters:
            all_deg (dict, optional): Precomputed degree dict.
            deg_type (str): Degree type.
            top (int): Number of top nodes to return.

        Returns:
            list: Top nodes sorted by degree.
        """
        if all_deg is None:
            all_deg = self.all_degrees(deg_type)
        ord_deg = sorted(all_deg.items(), key=lambda x: x[1], reverse=True)
        return [x[0] for x in ord_deg[:top]]

    def mean_degree(self, deg_type="inout"):
        """
        Compute the mean degree of all nodes.

        Parameters:
            deg_type (str): Degree type.

        Returns:
            float: Mean degree.
        """
        degs = self.all_degrees(deg_type)
        return sum(degs.values()) / float(len(degs))

    def prob_degree(self, deg_type="inout"):
        """
        Calculate the degree distribution as a probability.

        Parameters:
            deg_type (str): Degree type.

        Returns:
            dict: Mapping of degree value to probability.
        """
        degs = self.all_degrees(deg_type)
        res = {}
        for k in degs:
            res[degs[k]] = res.get(degs[k], 0) + 1
        for k in res:
            res[k] /= float(len(degs))
        return res

    def mean_distances(self):
        """
        Compute mean shortest path length and density.

        Returns:
            tuple: (mean distance, network density)
        """
        tot = 0
        num_reachable = 0
        for k in self.get_nodes():
            distsk = self.reachable_with_dist(k)
            for _, dist in distsk:
                tot += dist
            num_reachable += len(distsk)
        meandist = float(tot) / num_reachable if num_reachable else 0
        n = len(self.get_nodes())
        density = float(num_reachable) / ((n - 1) * n) if n > 1 else 0
        return meandist, density

    def closeness_centrality(self, node):
        """
        Compute closeness centrality of a node.

        Parameters:
            node: Target node.

        Returns:
            float: Closeness centrality.
        """
        dist = self.reachable_with_dist(node)
        if len(dist) == 0:
            return 0.0
        s = sum(d[1] for d in dist)
        return len(dist) / s

    def highest_closeness(self, top=10):
        """
        Return nodes with the highest closeness centrality.

        Parameters:
            top (int): Number of top nodes.

        Returns:
            list: Top nodes.
        """
        cc = {k: self.closeness_centrality(k) for k in self.get_nodes()}
        ord_cl = sorted(cc.items(), key=lambda x: x[1], reverse=True)
        return [x[0] for x in ord_cl[:top]]

    def betweenness_centrality(self, node):
        """
        Estimate betweenness centrality using path enumeration.

        Parameters:
            node: Target node.

        Returns:
            float: Betweenness score.
        """
        total_sp = 0
        sps_with_node = 0
        for s in self.get_nodes():
            for t in self.get_nodes():
                if s != t and s != node and t != node:
                    sp = self.shortest_path(s, t)
                    if sp:
                        total_sp += 1
                        if node in sp:
                            sps_with_node += 1
        return sps_with_node / total_sp if total_sp > 0 else 0.0

    def clustering_coef(self, v):
        """
        Calculate local clustering coefficient for a node.

        Parameters:
            v: Node ID.

        Returns:
            float: Clustering coefficient.
        """
        adjs = self.get_adjacents(v)
        if len(adjs) <= 1:
            return 0.0
        ligs = 0
        for i in adjs:
            for j in adjs:
                if i != j and (j in self.get_successors(i) or i in self.get_successors(j)):
                    ligs += 1
        return float(ligs) / (len(adjs) * (len(adjs) - 1))

    def all_clustering_coefs(self):
        """
        Compute clustering coefficients for all nodes.

        Returns:
            dict: Node to coefficient mapping.
        """
        return {k: self.clustering_coef(k) for k in self.get_nodes()}

    def mean_clustering_coef(self):
        """
        Compute mean clustering coefficient.

        Returns:
            float: Mean coefficient.
        """
        ccs = self.all_clustering_coefs()
        return sum(ccs.values()) / float(len(ccs)) if ccs else 0.0

    def mean_clustering_perdegree(self, deg_type="inout"):
        """
        Average clustering coefficient grouped by node degree.

        Parameters:
            deg_type (str): Degree type.

        Returns:
            dict: Degree to average coefficient.
        """
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


class CentralityAnalyzer:
    """
    Analyzer class for computing centrality metrics on a graph.
    """

    def __init__(self, graph):
        """
        Initialize with a graph instance.
        """
        self.graph = graph

    def degree_centrality(self):
        """
        Calculate degree centrality for all nodes.

        Returns:
            dict: Node to centrality value.
        """
        return {node: len(self.graph.get_successors(node)) for node in self.graph.get_nodes()}

    def closeness_centrality(self):
        """
        Compute closeness centrality for all nodes.

        Returns:
            dict: Node to centrality value.
        """
        centrality = {}
        for node in self.graph.get_nodes():
            total_dist, reachable_count = self._bfs_total_distance_and_reach_count(node)
            centrality[node] = (reachable_count / total_dist) if total_dist > 0 else 0.0
        return centrality

    def _bfs_total_distance_and_reach_count(self, start):
        """
        Helper to perform BFS and compute distance metrics.

        Returns:
            tuple: (total distance, reachable node count)
        """
        visited = set()
        queue = deque([(start, 0)])
        total = 0
        reachable_count = 0

        while queue:
            node, dist = queue.popleft()
            if node not in visited:
                visited.add(node)
                if node != start:
                    total += dist
                    reachable_count += 1
                for neighbor in self.graph.get_successors(node):
                    if neighbor not in visited:
                        queue.append((neighbor, dist + 1))
        return total, reachable_count

    def betweenness_centrality(self):
        """
        Calculate betweenness centrality for all nodes using Brandes' algorithm.

        Returns:
            dict: Node to centrality value.
        """
        centrality = dict.fromkeys(self.graph.get_nodes(), 0.0)
        for s in self.graph.get_nodes():
            stack = []
            pred = {w: [] for w in self.graph.get_nodes()}
            sigma = dict.fromkeys(self.graph.get_nodes(), 0)
            dist = dict.fromkeys(self.graph.get_nodes(), -1)
            sigma[s] = 1
            dist[s] = 0
            queue = deque([s])

            while queue:
                v = queue.popleft()
                stack.append(v)
                for w in self.graph.get_successors(v):
                    if dist[w] < 0:
                        dist[w] = dist[v] + 1
                        queue.append(w)
                    if dist[w] == dist[v] + 1:
                        sigma[w] += sigma[v]
                        pred[w].append(v)

            delta = dict.fromkeys(self.graph.get_nodes(), 0)
            while stack:
                w = stack.pop()
                for v in pred[w]:
                    delta[v] += (sigma[v] / sigma[w]) * (1 + delta[w])
                if w != s:
                    centrality[w] += delta[w]
        return centrality

    def top_nodes(self, centrality_dict, top_n=5):
        """
        Return top N nodes by centrality score.

        Parameters:
            centrality_dict (dict): Centrality mapping.
            top_n (int): Number of top nodes.

        Returns:
            list: Top node-score pairs.
        """
        return heapq.nlargest(top_n, centrality_dict.items(), key=lambda x: x[1])


# REACTION PARSER
def parse_reactions(file_path):
    """
    Parses reactions from file into structured dicts.

    Parameters:
        file_path (str): Path to the reaction file.

    Returns:
        list: List of reaction dicts with substrates and products.
    """
    reactions = []
    with open(file_path, 'r') as f:
        for line in f:
            if ':' not in line:
                continue
            parts = re.split(r':\s*', line.strip(), maxsplit=1)
            if len(parts) != 2:
                continue
            reaction_id, formula = parts
            match = re.search(r"^(.*?)\s*(<=>|=>)\s*(.*?)$", formula)
            if not match:
                continue
            substrates = [m.strip() for m in match.group(1).split('+')]
            products = [m.strip() for m in match.group(3).split('+')]
            reactions.append({
                'id': reaction_id,
                'substrates': substrates,
                'products': products
            })
    return reactions


def build_metabolite_graph(reactions):
    """
    Builds a graph where edges represent shared participation in a reaction.

    Parameters:
        reactions (list): List of parsed reactions.

    Returns:
        MN_Graph: Constructed metabolite graph.
    """
    g = MN_Graph()
    for r in reactions:
        metabolites = r['substrates'] + r['products']
        for i in range(len(metabolites)):
            for j in range(i + 1, len(metabolites)):
                g.add_edge(metabolites[i], metabolites[j])
    return g


# ASSESSMENT FUNCTIONS
def get_active_reactions(metabolites_set, reactions):
    """
    Get reactions whose substrates are all in the known metabolite set.

    Parameters:
        metabolites_set (set): Set of available metabolites.
        reactions (list): List of all reactions.

    Returns:
        list: Active reactions.
    """
    return [r for r in reactions if all(sub in metabolites_set for sub in r['substrates'])]

def get_produced_metabolites(active_reactions):
    """
    Get products produced by active reactions.

    Parameters:
        active_reactions (list): List of active reactions.

    Returns:
        set: Produced metabolites.
    """
    produced = set()
    for r in active_reactions:
        produced.update(r['products'])
    return produced

def compute_final_metabolites(initial_metabolites, reactions):
    """
    Iteratively compute the full set of metabolites that can be produced.

    Parameters:
        initial_metabolites (set): Starting set of metabolites.
        reactions (list): All reactions.

    Returns:
        set: Complete reachable metabolite set.
    """
    known_metabolites = set(initial_metabolites)
    while True:
        active = get_active_reactions(known_metabolites, reactions)
        new = get_produced_metabolites(active)
        if new.issubset(known_metabolites):
            break
        known_metabolites.update(new)
    return known_metabolites