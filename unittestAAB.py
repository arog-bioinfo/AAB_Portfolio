import random
import unittest
from Algorithms.BB import *
from Algorithms.exaustivo import *
from Algorithms.gibbs import *
# from Algorithms.BWT import *
# from Algorithms.automato import *
# from Algorithms.tsTrees import *

def gerar_seq(n_seqs, tam):
    nucleotidos = "ACTG"
    seqs = [''.join(random.choices(nucleotidos, k=tam)) for _ in range(n_seqs)]
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

class TestExhaustiveSearch(unittest.TestCase):
    def test_exaustiive_search(self):
        seqs = gerar_seq(3, 10)  
        m = str(*gerar_seq(1, 4))
        seqs, posicoes = inserir_motif(seqs, m)
        tam_motif = len(m)
        
        interations = (10 - tam_motif + 1) ** 3 #(seq_len - motig_len + 1) ** n_seqs

        self.assertEqual(motif(seqs, 3, 10, tam_motif, counting = True), (tuple(posicoes), 12, interations))

class TestBB(unittest.TestCase):
    def test_bb(self):
        seqs = gerar_seq(3, 10)  
        m = str(*gerar_seq(1, 4))
        seqs, posicoes = inserir_motif(seqs, m)
        tam_motif = len(m)

        self.assertEqual(motif_bb(seqs, 3, 10, tam_motif), (posicoes,12))
    
    def test_ES_versus_BB(self):
        seqs = gerar_seq(3, 10)  
        m = "ACGT"
        seqs, posicoes = inserir_motif(seqs, m)
        tam_motif = len(m)

        *_, a = motif(seqs, 3, 10, tam_motif)
        *_, b = motif_bb(seqs, 3, 10, tam_motif)

        self.assertGreaterEqual(a,b)

class TestGibbs(unittest.TestCase):
    def test_Gibbs(self):
        seqs = gerar_seq(3, 10)  
        m = "ACGT"
        seqs, posicoes = inserir_motif(seqs, m)
        tam_motif = len(m)
        gibbs = GibbsSampling(seqs, 4)
        p , s = gibbs.gibbs_sampling()
        self.assertEqual( (list(p.values())), (posicoes))

#Checkpoint 2

class TestAutomata(unittest.TestCase):
    pass

    def test_no_ocurrence(self):
        seq = str(*gerar_seq(1,20))
        padrao = "NNNNN"

        self.assertFalse(auto_fin(seq,padrao))

    def test_one_ocurrence(self):
        seq = gerar_seq(1,20)
        padrao ="AQUI"
        seq, posicoes = inserir_motif(seq, padrao)

        self.assertEqual(auto_fin(seq, padrao), (1,posicoes))
    
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
    def test_BWT(self):
        pass

if __name__ == "__main__":
    unittest.main()