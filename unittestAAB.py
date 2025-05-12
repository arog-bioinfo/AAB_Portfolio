import random
import unittest
from Algorithms.BB import *
from Algorithms.exaustivo import *
from Algorithms.gibbs import *
from Algorithms.BWT import *
from Algorithms.auto_finito import *
from Algorithms.tries import *
from Algorithms.graph import * 
from Algorithms.metabolic_networks import *
from Algorithms.assembly import *




def gerar_seq(n_seqs, tam, alphabet = "ACTG" ):

    seqs = [''.join(random.choices(alphabet, k=tam)) for _ in range(n_seqs)]
    return seqs

def inserir_motif(seqs, motif):
    tam_motif = len(motif)
    posicoes = []

    for i in range(len(seqs)):
        tam_seq = len(seqs[i])
        pos = random.randint(0, tam_seq - tam_motif)  # Escolhe uma posição aleatória válida
        seqs[i] = seqs[i][:pos] + motif + seqs[i][pos + tam_motif:]  # Substitui parte da sequência pelo motif
        posicoes.append(pos)

    return seqs, posicoes

#Checkpoint 1

class TestMotifs(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.m = str(*gerar_seq(1, 4, "XYZV"))
        cls.tam_motif = len(cls.m)
        cls.interations = (10 - cls.tam_motif + 1) ** 3  # (seq_len - motif_len + 1) ** n_seqs

        score = 12
        while score == 12:
            seqs = gerar_seq(3, 10)
            _, score = motif(seqs, 3, 10, cls.tam_motif)

        cls.seqs, cls.posicoes = inserir_motif(seqs, cls.m)

    def test_exhaustive_search(self):
        self.assertEqual(
            motif(self.seqs, 3, 10, self.tam_motif, counting=True),
            (tuple(self.posicoes), 12, self.interations), self.seqs
        )

    def test_bb(self):
        self.assertEqual(
            motif_bb(self.seqs, 3, 10, self.tam_motif),
            (self.posicoes, 12), self.seqs
        )

    def test_ES_versus_BB(self):
        *_, a = motif(self.seqs, 3, 10, self.tam_motif)
        *_, b = motif_bb(self.seqs, 3, 10, self.tam_motif)
        self.assertGreaterEqual(a, b)

class TestGibbs(unittest.TestCase):
    def test_Gibbs_flaky_with_threshold(self):
        # Parameters
        n_trials = 50  # Number of test runs
        success_threshold = 0.5  # Require 50% success rate 
        
        # Generate different sequences with inserted motif
        seqs = ["AGCTTACGGA", "TGCAGTCTAC", "CCGTAGGAAT"]
        while seqs == None or len(seqs) != len(set(seqs)):
            seqs = gerar_seq(3, 10)

        m = str(*gerar_seq(1, 4))
        seqs, true_positions = inserir_motif(seqs, m)
        tam_motif = len(m)
        
        # Run Gibbs sampling multiple times
        successes = 0
        all_false_predicted_positions = []
        for _ in range(n_trials):
            gibbs = GibbsSampling(seqs, tam_motif)
            predicted_positions, v = gibbs.gibbs_sampling()

            # Check if predicted positions match true positions
            if list(predicted_positions.values()) == true_positions:
                successes += 1
            else:
                all_false_predicted_positions.append([list(predicted_positions.values()), v])
        
        success_rate = successes / n_trials
        self.assertGreaterEqual(
            success_rate,
            success_threshold,
            f"Gibbs succeeded in {successes}/{n_trials} trials ({success_rate*100:.1f}%), "
            f"but needed {success_threshold*100}%. \n Seqs:{seqs}"
            f"\nPosition: {true_positions}\nFalse_predictions: {all_false_predicted_positions}"
        )

#Checkpoint 2

class TestAutomata(unittest.TestCase):

    def test_no_ocurrence(self):
        seq = str(*gerar_seq(1,20))
        alphabet = {"A","C","T","G","N"}
        pattern = "NNNNN"

        automaton = AutomatosFinitos(alphabet, pattern)
        count, positions = automaton.find_occurrences(seq)
        
        self.assertTrue((count,positions), (0,[]))

    def test_one_ocurrence(self):
        seq = gerar_seq(1,20)
        alphabet = {"A","C","T","G","H","E","R"}
        pattern ="HERE"
        seq, p = inserir_motif(seq, pattern)

        automaton = AutomatosFinitos(alphabet, pattern)
        count, positions = automaton.find_occurrences(*seq)

        self.assertEqual((count,positions), (1, p), f"{seq}")
    
    def test_more_ocurrences(self):
        seq = ""
        pos = []
        for _ in range(4):
            s = gerar_seq(1,8)
            pattern ="HERE"
            s, p = inserir_motif(s, pattern)

            seq += s[0]
            pos.append(p[0] + (8*_))


        alphabet = {"A","C","T","G","H","E","R"}
        automaton = AutomatosFinitos(alphabet, pattern)
        count, positions = automaton.find_occurrences(seq)

        self.assertEqual((count,positions), (4, pos), f"{seq}")

class TestTrees(unittest.TestCase):

    def test_prefix_trie_findpattern(self):
        seq = gerar_seq(1, 10)
        pattern = gerar_seq(1,4, alphabet="XYZW")
        seq, p = inserir_motif(seq, *pattern)

        t = Trie()
        t.trie_from_patterns(pattern)
        self.assertEqual(t.trie_matches(seq[0]), [(p[0], pattern[0])])

    def test_prefix_find2patterns(self):
        seq = gerar_seq(1, 10)
        pattern = gerar_seq(1,2, alphabet="XYZW")
        m = pattern[0] + "A" + pattern[0]
        seq, p = inserir_motif(seq, m)

        t = Trie()
        t.trie_from_patterns(pattern)
        self.assertEqual(t.trie_matches(seq[0]), [(p[0], pattern[0]),(p[0] +len(pattern[0]) + 1, pattern[0])])

    def test_prefix_match(self):
        seq = gerar_seq(1, 10)

        real_pattern = seq[0][0:5]
        patterns = gerar_seq(4, 4, alphabet="XYZW")
        patterns.append(real_pattern)

        t = Trie()
        t.trie_from_patterns(patterns)
        self.assertEqual(t.prefix_trie_match(seq[0]), real_pattern)


    def test_suffix_trie(self):
        seq = gerar_seq(1, 10)
        pattern = gerar_seq(1,2, alphabet="XYZW")
        m = str(*pattern) + "A" + str(*pattern)
        seq, p = inserir_motif(seq, m)

        st = SuffixTree()
        st.build(*seq)
        self.assertEqual(st.find_pattern(pattern[0]), [p[0],p[0]+3])

class TestBwt(unittest.TestCase):

    def test_BWT_backtracking(self):
        seq = str(*gerar_seq(1,10))
        bwt_instance = BWT(seq)
        resultado = bwt_instance.bwt
        recuperado = bwt_reverse(resultado)
        
        self.assertEqual(recuperado, seq)

    def test_BWT_subwords(self):
        seq = gerar_seq(1, 10)
        pattern = gerar_seq(1,2, alphabet="XYZW")
        m = str(*pattern) + "A" + str(*pattern)
        seq, p = inserir_motif(seq, m)

        bwt_instance = BWT(*seq, buildsufarray=True)

        self.assertEqual(bwt_instance.bw_matching_pos(pattern[0]), [p[0], p[0]+3])

#Checkpoint 3

class TestGraph(unittest.TestCase):

    def setUp(self):
        self.Graph = Graph

    def test_add_vertex_and_edge(self):
        g = self.Graph()
        g.add_vertex("A")
        g.add_edge("A", "B")
        self.assertIn("A", g.graph)
        self.assertIn("B", g.graph)
        self.assertIn("B", g.graph["A"])
        self.assertNotIn("A", g.graph["B"])

    def test_get_nodes_and_edges(self):
        g = self.Graph({"A": ["B", "C"], "B": ["C"]})
        self.assertCountEqual(g.get_nodes(), ["A", "B", "C"])
        self.assertCountEqual(g.get_edges(), [("A", "B"), ("A", "C"), ("B", "C")])

    def test_degrees(self):
        g = self.Graph({"A": ["B", "C"], "C": ["A"]})
        self.assertEqual(g.in_degree("A"), 1)
        self.assertEqual(g.out_degree("A"), 2)
        self.assertEqual(g.degree("A"), 3)  # in + out without double-count

    def test_successors_and_predecessors(self):
        g = self.Graph({"A": ["B"], "B": ["C"], "C": []})
        self.assertEqual(g.get_successors("A"), ["B"])
        self.assertEqual(g.get_predecessors("C"), ["B"])

    def test_bfs_reachability(self):
        g = self.Graph({"A": ["B"], "B": ["C"], "C": []})
        self.assertCountEqual(g.reachable_bfs("A"), ["B", "C"])
        self.assertCountEqual(g.reachable_bfs("B"), ["C"])

    def test_dfs_reachability(self):
        g = self.Graph({"A": ["B"], "B": ["C"], "C": []})
        self.assertCountEqual(g.reachable_dfs("A"), ["B", "C"])

    def test_distance_and_shortest_path(self):
        g = self.Graph({"A": ["B"], "B": ["C"], "C": []})
        self.assertEqual(g.distance("A", "C"), 2)
        self.assertEqual(g.shortest_path("A", "C"), ["A", "B", "C"])
        self.assertIsNone(g.distance("C", "A"))

    def test_cycle_detection(self):
        g1 = self.Graph({"A": ["B"], "B": ["C"], "C": ["A"]})  # cycle
        g2 = self.Graph({"A": ["B"], "B": ["C"], "C": []})     # no cycle
        self.assertTrue(g1.has_cycle())
        self.assertFalse(g2.has_cycle())

    def test_size(self):
        g = self.Graph({"A": ["B"], "B": ["C"]})
        nodes, edges = g.size()
        self.assertEqual(nodes, 3)
        self.assertEqual(edges, 2)

class TestWeightnedGraph(unittest.TestCase):
    def setUp(self):
        self.WeightedGraph = WeightedGraph

    def test_add_edge_and_get_edges(self):
        g = self.WeightedGraph()
        g.add_edge('A', 'B', 5)
        g.add_edge('A', 'C', 2)
        self.assertIn(('A', 'B', 5), g.get_edges())
        self.assertIn(('A', 'C', 2), g.get_edges())

    def test_graph_construction_from_dict(self):
        g = self.WeightedGraph({
            'A': [('B', 3), ('C', 1)],
            'B': [('C', 7)]
        })
        self.assertCountEqual(g.get_edges(), [
            ('A', 'B', 3),
            ('A', 'C', 1),
            ('B', 'C', 7)
        ])

    def test_get_successors(self):
        g = self.WeightedGraph({
            'X': [('Y', 2), ('Z', 4)],
            'Y': [('Z', 1)]
        })
        self.assertEqual(g.get_successors('X'), ['Y', 'Z'])
        self.assertEqual(g.get_successors('Y'), ['Z'])
        self.assertEqual(g.get_successors('Z'), [])

    def test_dijkstra(self):
        g = self.WeightedGraph({
            'A': [('B', 1), ('C', 4)],
            'B': [('C', 2), ('D', 5)],
            'C': [('D', 1)],
            'D': []
        })
        distances = g.dijkstra('A')
        expected = {'A': 0, 'B': 1, 'C': 3, 'D': 4}
        self.assertEqual(distances, expected)

    def test_dijkstra_disconnected(self):
        g = self.WeightedGraph({
            'A': [('B', 2)],
            'B': [],
            'C': []  # disconnected node
        })
        distances = g.dijkstra('A')
        self.assertEqual(distances['A'], 0)
        self.assertEqual(distances['B'], 2)
        self.assertEqual(distances['C'], float('inf'))

class TestBiologicalNetwork(unittest.TestCase):
    def setUp(self):
        self.graph = MN_Graph({
            'A': ['B', 'C'],
            'B': ['C', 'D'],
            'C': ['D'],
            'D': ['E'],
            'E': []
        })

    def test_reachable_with_dist(self):
        reachable = self.graph.reachable_with_dist('A')
        self.assertIn(('B', 1), reachable)
        self.assertIn(('C', 1), reachable)
        self.assertIn(('D', 2), reachable)
        self.assertIn(('E', 3), reachable)

    def test_mean_distances(self):
        mean_dist, density = self.graph.mean_distances()
        self.assertGreater(mean_dist, 0)
        self.assertGreaterEqual(density, 0)
        self.assertLessEqual(density, 1)

    def test_closeness_centrality(self):
        cc = self.graph.closeness_centrality('A')
        self.assertGreater(cc, 0)
        self.assertLessEqual(cc, 1)

    def test_highest_closeness(self):
        top = self.graph.highest_closeness(top=3)
        self.assertTrue(all(n in self.graph.get_nodes() for n in top))
        self.assertEqual(len(top), 3)

    def test_betweenness_centrality(self):
        bc = self.graph.betweenness_centrality('C')
        self.assertGreaterEqual(bc, 0)
        self.assertLessEqual(bc, 1)

    def test_clustering_coef(self):
        coef = self.graph.clustering_coef('B')
        self.assertGreaterEqual(coef, 0)
        self.assertLessEqual(coef, 1)

    def test_mean_clustering_coef(self):
        mean_cc = self.graph.mean_clustering_coef()
        self.assertGreaterEqual(mean_cc, 0)
        self.assertLessEqual(mean_cc, 1)

    def test_mean_clustering_perdegree(self):
        ck = self.graph.mean_clustering_perdegree()
        for k, v in ck.items():
            self.assertGreaterEqual(v, 0)
            self.assertLessEqual(v, 1)

class TestAssembly(unittest.TestCase):
    pass

if __name__ == "__main__":
    unittest.main()