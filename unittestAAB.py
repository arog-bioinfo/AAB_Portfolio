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

#=================================================================
#---------------Oriented Graph and Weightned Graph----------------
#=================================================================

class TestGraph(unittest.TestCase):

    def setUp(self):
        self.Graph = Graph

    def test_add_vertex_and_edge_complex(self):
        g = self.Graph()
        g.add_vertex("A")
        g.add_edge("A", "B")
        g.add_edge("B", "A")  # making it bidirectional
        g.add_edge("C", "C")  # self-loop
        g.add_edge("A", "C")
        self.assertIn("A", g.graph)
        self.assertIn("B", g.graph)
        self.assertIn("C", g.graph)
        self.assertIn("B", g.graph["A"])
        self.assertIn("A", g.graph["B"])
        self.assertIn("C", g.graph["C"])  # self-loop
        self.assertNotIn("D", g.graph)

    def test_get_nodes_and_edges_with_self_loop(self):
        g = self.Graph({"A": ["B", "A"], "B": ["C"], "C": []})
        self.assertCountEqual(g.get_nodes(), ["A", "B", "C"])
        self.assertCountEqual(g.get_edges(), [("A", "B"), ("A", "A"), ("B", "C")])

    def test_degrees_with_self_and_multiple_edges(self):
        g = self.Graph({"A": ["B", "C", "A"], "C": ["A"], "B": []})
        self.assertEqual(g.in_degree("A"), 2)
        self.assertEqual(g.out_degree("A"), 3)
        self.assertEqual(g.degree("A"), 5)

    def test_successors_and_predecessors_with_disconnected_node(self):
        g = self.Graph({"A": ["B"], "B": ["C"], "C": [], "D": []})
        self.assertEqual(g.get_successors("D"), [])
        self.assertEqual(g.get_predecessors("C"), ["B"])

    def test_bfs_and_dfs_on_disconnected_graph(self):
        g = self.Graph({"A": ["B"], "B": ["C"], "C": [], "X": ["Y"], "Y": []})
        self.assertCountEqual(g.reachable_bfs("A"), ["B", "C"])
        self.assertCountEqual(g.reachable_bfs("X"), ["Y"])
        self.assertCountEqual(g.reachable_dfs("A"), ["B", "C"])

    def test_distance_and_shortest_path_extended(self):
        g = self.Graph({
            "A": ["B", "C"],
            "B": ["D"],
            "C": ["D"],
            "D": [],
            "E": []
        })
        self.assertEqual(g.distance("A", "D"), 2)
        self.assertIn(g.shortest_path("A", "D"), [["A", "B", "D"], ["A", "C", "D"]])
        self.assertIsNone(g.shortest_path("D", "A"))
        self.assertIsNone(g.distance("E", "A"))

    def test_reachable_with_dist_and_isolated(self):
        g = self.Graph({
            "A": ["B", "C"],
            "B": ["C"],
            "C": ["D"],
            "D": [],
            "E": [],
            "F": ["F"]  # self-loop node
        })
        reachable = g.reachable_with_dist("A")
        expected = [('B', 1), ('C', 1), ('D', 2)]
        for node in expected:
            self.assertIn(node, reachable)
        self.assertNotIn(('E', 1), reachable)
        self.assertNotIn(('F', 1), reachable)

    def test_cycle_detection_extended(self):
        g1 = self.Graph({"A": ["B"], "B": ["C"], "C": ["A"], "D": []})  # cycle + isolated
        g2 = self.Graph({"A": ["B"], "B": ["C"], "C": [], "D": []})    # no cycle
        g3 = self.Graph({"X": ["X"]})  # self-loop = cycle
        self.assertTrue(g1.has_cycle())
        self.assertFalse(g2.has_cycle())
        self.assertTrue(g3.has_cycle())

    def test_size_with_disconnected_and_self_loops(self):
        g = self.Graph({
            "A": ["B", "A"],
            "B": ["C"],
            "C": [],
            "D": ["D"],
            "E": []
        })
        nodes, edges = g.size()
        self.assertEqual(nodes, 5)
        self.assertEqual(edges, 4)  # A→B, A→A, B→C, D→D

    def test_random_large_graph_consistency(self):
        g = self.Graph()
        nodes = [chr(i) for i in range(65, 75)]  # A-J
        for node in nodes:
            g.add_vertex(node)
        for _ in range(20):
            u, v = random.sample(nodes, 2)
            g.add_edge(u, v)

        all_nodes = g.get_nodes()
        self.assertCountEqual(all_nodes, nodes)
        size_nodes, size_edges = g.size()
        self.assertEqual(size_nodes, len(nodes))
        self.assertTrue(size_edges <= 20)

    def test_invalid_node_operations(self):
        g = self.Graph({"A": ["B"], "B": []})
        with self.assertRaises(KeyError):
            g.in_degree("Z")
        with self.assertRaises(KeyError):
            g.out_degree("Z")
        with self.assertRaises(KeyError):
            g.get_successors("Z")
        with self.assertRaises(KeyError):
            g.get_predecessors("Z")

class TestWeightnedGraph(unittest.TestCase):

    def setUp(self):
        self.WeightedGraph = WeightedGraph

    def test_add_edge_and_get_edges_complex(self):
        g = self.WeightedGraph()
        g.add_edge('A', 'B', 5)
        g.add_edge('A', 'C', 2)
        g.add_edge('B', 'A', 1)  # reverse edge
        g.add_edge('C', 'C', 7)  # self-loop
        edges = g.get_edges()
        self.assertIn(('A', 'B', 5), edges)
        self.assertIn(('A', 'C', 2), edges)
        self.assertIn(('B', 'A', 1), edges)
        self.assertIn(('C', 'C', 7), edges)

    def test_graph_construction_from_dict_complex(self):
        g = self.WeightedGraph({
            'A': [('B', 3), ('C', 1)],
            'B': [('C', 7)],
            'D': []  # disconnected node
        })
        self.assertCountEqual(g.get_edges(), [
            ('A', 'B', 3),
            ('A', 'C', 1),
            ('B', 'C', 7)
        ])
        self.assertIn('D', g.get_nodes())
        self.assertEqual(g.get_successors('D'), [])

    def test_get_successors_with_weights(self):
        g = self.WeightedGraph({
            'X': [('Y', 2), ('Z', 4)],
            'Y': [('Z', 1)],
            'Z': []
        })
        self.assertEqual(set(g.get_successors('X')), {'Y', 'Z'})
        self.assertEqual(g.get_successors('Y'), ['Z'])
        self.assertEqual(g.get_successors('Z'), [])

    def test_dijkstra_with_cycle_and_self_loop(self):
        g = self.WeightedGraph({
            'A': [('B', 1)],
            'B': [('C', 2)],
            'C': [('A', 3), ('C', 0)],  # cycle and self-loop
            'D': []  # disconnected
        })
        distances = g.dijkstra('A')
        self.assertEqual(distances['A'], 0)
        self.assertEqual(distances['B'], 1)
        self.assertEqual(distances['C'], 3)
        self.assertEqual(distances['D'], float('inf'))

    def test_dijkstra_disconnected_extended(self):
        g = self.WeightedGraph({
            'A': [('B', 2)],
            'B': [],
            'C': [],
            'D': [('E', 1)],
            'E': []
        })
        distances = g.dijkstra('A')
        self.assertEqual(distances['A'], 0)
        self.assertEqual(distances['B'], 2)
        self.assertEqual(distances['C'], float('inf'))
        self.assertEqual(distances['D'], float('inf'))
        self.assertEqual(distances['E'], float('inf'))

    def test_update_existing_edge(self):
        g = self.WeightedGraph()
        g.add_edge('X', 'Y', 10)
        g.add_edge('X', 'Y', 5)  # update
        self.assertIn(('X', 'Y', 5), g.get_edges())
        self.assertNotIn(('X', 'Y', 10), g.get_edges())

    def test_add_edge_invalid_weight(self):
        g = self.WeightedGraph()
        with self.assertRaises(TypeError):
            g.add_edge('A', 'B', 'high')  # invalid weight

    def test_nonexistent_node_handling(self):
        g = self.WeightedGraph({'A': [('B', 1)], 'B': []})
        with self.assertRaises(KeyError):
            g.get_successors('Z')

        with self.assertRaises(KeyError):
            g.dijkstra('Z')

    def test_nodes_and_edges_consistency(self):
        g = self.WeightedGraph({
            'M': [('N', 2), ('O', 4)],
            'N': [('P', 1)],
            'O': [],
            'P': []
        })
        self.assertCountEqual(g.get_nodes(), ['M', 'N', 'O', 'P'])
        edges = g.get_edges()
        self.assertIn(('M', 'N', 2), edges)
        self.assertIn(('M', 'O', 4), edges)
        self.assertIn(('N', 'P', 1), edges)
        self.assertEqual(len(edges), 3)

#=================================================================
#-----------------------Metabolic Network-------------------------
#=================================================================

class TestMNGraph(unittest.TestCase):

    def setUp(self):
        # Simulate a directed metabolic network with cycles, hubs, and isolates
        self.graph = MN_Graph({
            'A': ['B', 'C', 'D'],
            'B': ['C'],
            'C': ['D', 'E'],
            'D': ['A', 'F'],
            'E': ['F'],
            'F': [],
            'G': [],         # isolated node
            'H': ['I', 'J'],
            'I': ['J'],
            'J': []
        })

    def test_mean_distances_debug(self):
        g = MN_Graph({
            'A': ['B'],
            'B': ['C'],
            'C': ['D'],
            'D': []
        })
        
        mean_dist, density = g.mean_distances()
        self.assertAlmostEqual(mean_dist, 1.6667, delta=0.0001)
        self.assertAlmostEqual(density, 0.5, delta=0.0001) #Manually calculated

    def test_closeness_centrality_values(self):
        cc_A = self.graph.closeness_centrality('A')
        cc_G = self.graph.closeness_centrality('G')  # isolated node
        cc_D = self.graph.closeness_centrality('D')

        self.assertGreater(cc_A, 0)
        self.assertEqual(cc_G, 0)  # no reachable nodes
        self.assertGreater(cc_D, cc_G)

    def test_highest_closeness_top3(self):
        top_nodes = self.graph.highest_closeness(top=3)
        self.assertEqual(len(top_nodes), 3)
        scores = [self.graph.closeness_centrality(n) for n in top_nodes]
        self.assertTrue(all(x >= y for x, y in zip(scores, scores[1:])))

    def test_betweenness_centrality_distribution(self):
        bc_C = self.graph.betweenness_centrality('C')
        bc_G = self.graph.betweenness_centrality('G')
        bc_D = self.graph.betweenness_centrality('D')

        self.assertGreaterEqual(bc_C, 0)
        self.assertLessEqual(bc_C, 1)
        self.assertEqual(bc_G, 0)
        self.assertGreater(bc_C, bc_G)
        self.assertNotEqual(bc_D, 0)  # D is in multiple paths

    def test_clustering_coef_boundaries(self):
        coef_B = self.graph.clustering_coef('B')
        coef_A = self.graph.clustering_coef('A')
        coef_G = self.graph.clustering_coef('G')  # isolated node

        self.assertGreaterEqual(coef_B, 0)
        self.assertLessEqual(coef_B, 1)
        self.assertGreaterEqual(coef_A, 0)
        self.assertLessEqual(coef_A, 1)
        self.assertEqual(coef_G, 0)

    def test_mean_clustering_coef_accuracy(self):
        mean_cc = self.graph.mean_clustering_coef()
        self.assertGreaterEqual(mean_cc, 0)
        self.assertLessEqual(mean_cc, 1)
        # Optional: check against manual average if expected

    def test_mean_clustering_perdegree_structure(self):
        ck = self.graph.mean_clustering_perdegree()
        for degree, coef in ck.items():
            self.assertIsInstance(degree, int)
            self.assertGreaterEqual(coef, 0)
            self.assertLessEqual(coef, 1)

        unique_degrees = set(self.graph.all_degrees().values())
        self.assertEqual(set(ck.keys()), unique_degrees)


class TestCentralityAnalyzer(unittest.TestCase):
    def setUp(self):
        self.graph = MN_Graph({
            'A': ['B', 'C'],
            'B': ['C', 'D'],
            'C': ['D'],
            'D': ['E'],
            'E': []
        })
        self.analyzer = CentralityAnalyzer(self.graph)

    def test_degree_centrality(self):
        degree_c = self.analyzer.degree_centrality()
        expected = {
            'A': 2,
            'B': 2,
            'C': 1,
            'D': 1,
            'E': 0
        }
        self.assertEqual(degree_c, expected)

    def test_closeness_centrality_values(self):
        closeness = self.analyzer.closeness_centrality()
        for val in closeness.values():
            self.assertGreaterEqual(val, 0.0)
            self.assertLessEqual(val, 1.0)

    def test_betweenness_centrality_values(self):
        betweenness = self.analyzer.betweenness_centrality()
        for val in betweenness.values():
            self.assertGreaterEqual(val, 0.0)

    def test_top_nodes_degree(self):
        degree = self.analyzer.degree_centrality()
        top = self.analyzer.top_nodes(degree, top_n=2)
        self.assertEqual(len(top), 2)
        # Should contain nodes with out-degree 2: A and B
        self.assertIn('A', [n for n, _ in top])
        self.assertIn('B', [n for n, _ in top])

    def test_top_nodes_betweenness(self):
        bc = self.analyzer.betweenness_centrality()
        top = self.analyzer.top_nodes(bc, top_n=2)
        self.assertEqual(len(top), 2)
        self.assertTrue(all(n in self.graph.get_nodes() for n, _ in top))

#=================================================================
#-----------------------Assembly Algorythm------------------------
#=================================================================

class TestDeBruijnAssembly(unittest.TestCase):
    def test_composition_basic(self):
        self.assertCountEqual(
            composition(3, "ATGCG"),
            ["ATG", "TGC", "GCG"]
        )

    def test_composition_with_duplicates(self):
        self.assertCountEqual(
            composition(3, "ATGCGATG"),
            ["ATG", "TGC", "GCG", "CGA", "GAT", "ATG"]
        )

    def test_debruijn_graph_construction_complex(self):
        kmers = ["ATG", "TGC", "GCG", "CGA", "GAT"]
        dbg = DeBruijnGraph(kmers)
        expected_nodes = {"AT", "TG", "GC", "CG", "GA"}
        self.assertTrue(expected_nodes.issubset(set(dbg.graph.keys())))
        self.assertIn("GC", dbg.graph["TG"])
        self.assertIn("GA", dbg.graph["CG"])

    def test_debruijn_eulerian_path_validity(self):
        kmers = ["CTTA", "ACCA", "TACC", "GGCT", "GCTT", "TTAC"]
        dbg = DeBruijnGraph(kmers)
        path = dbg.eulerian_path()
        self.assertIsNotNone(path)
        seq = dbg.seq_from_path(path)
        for kmer in kmers:
            self.assertIn(kmer, seq)

    def test_debruijn_no_path_for_disconnected_graph(self):
        kmers = ["AAA", "TTT", "GGG"]
        dbg = DeBruijnGraph(kmers)
        path = dbg.eulerian_path()
        self.assertIsNone(path)

    def test_debruijn_seq_from_path_none(self):
        dbg = DeBruijnGraph([])
        self.assertIsNone(dbg.seq_from_path(None))

class TestOverlapAssembly(unittest.TestCase):

    def test_overlap_graph_construction_simple(self):
        frags = ["ATT", "TTA", "TAC"]
        og = OverlapGraph(frags)

        def clean(fragment_id):
            return fragment_id.split('-')[0]

        for src in og.graph:
            for dst in og.graph[src]:
                src_clean = clean(src)
                dst_clean = clean(dst)
                self.assertTrue(
                    src_clean[-2:] == dst_clean[:2],
                    msg=f"Failed overlap: {src_clean} → {dst_clean}"
                )

    def test_overlap_graph_with_multiple_overlaps(self):
        frags = ["ATTA", "TTAC", "TACC", "ACCG"]
        og = OverlapGraph(frags)
        path = og.search_hamiltonian_path()
        self.assertIsNotNone(path)
        seq = og.seq_from_path(path)
        for frag in frags:
            self.assertIn(frag, seq)

    def test_overlap_graph_hamiltonian_fails_for_disconnected(self):
        frags = ["AAA", "CCC", "GGG"]
        og = OverlapGraph(frags)
        self.assertIsNone(og.search_hamiltonian_path())

    def test_overlap_seq_from_path_none(self):
        og = OverlapGraph([])
        self.assertIsNone(og.seq_from_path(None))

if __name__ == "__main__":
    unittest.main()