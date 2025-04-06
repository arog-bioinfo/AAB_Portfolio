import random
import unittest
from Algorithms.BB import *
from Algorithms.exaustivo import *
from Algorithms.gibbs import *
from Algorithms.BWT import *
from Algorithms.auto_finito import *
# from Algorithms.tsTrees import *

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
            (tuple(self.posicoes), 12, self.interations)
        )

    def test_bb(self):
        self.assertEqual(
            motif_bb(self.seqs, 3, 10, self.tam_motif),
            (self.posicoes, 12)
        )

    def test_ES_versus_BB(self):
        *_, a = motif(self.seqs, 3, 10, self.tam_motif)
        *_, b = motif_bb(self.seqs, 3, 10, self.tam_motif)
        self.assertGreaterEqual(a, b)

class TestGibbs(unittest.TestCase):
    def test_Gibbs_flaky_with_threshold(self):
        # Parameters
        n_trials = 10  # Number of test runs
        success_threshold = 0.7  # Require 70% success rate (adjust as needed)
        
        # Generate different sequences with inserted motif
        seqs = None
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
    pass

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
        pass

class TestTrees(unittest.TestCase):
    def test_prefix_trie(self):
        pass

    def test_suffix_trie(self):
        pass

    def test_suffixcom_trie(self):
        pass

class TestBwt(unittest.TestCase):

    def test_BWT_backtracking(self):

        seq = str(gerar_seq(1,10))
        resultado = bwt_transform(seq)
        recuperado = bwt_reverse(resultado[1])
        
        self.assertEqual(recuperado, seq)

if __name__ == "__main__":
    unittest.main()