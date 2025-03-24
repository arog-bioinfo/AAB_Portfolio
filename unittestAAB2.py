import random
import unittest
from Algorithms.BB import *
from Algorithms.exaustivo import *
from Algorithms.gibbs import *

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

class counter():
    def __init__(self):
        self.count = 0
        
        pass
    pass
class TestExhaustiveSearch(unittest.TestCase):
    def test_exaustiive_search(self):
        seqs = gerar_seq(3, 10)  
        m = "ACGT"
        seqs, posicoes = inserir_motif(seqs, m)
        tam_motif = len(m)

        self.assertEqual(motif(seqs, 3, 10, tam_motif), (tuple(posicoes),12))

class TestBB(unittest.TestCase):
    def test_bb(self):
        seqs = gerar_seq(3, 10)  
        m = "ACGT"
        seqs, posicoes = inserir_motif(seqs, m)
        tam_motif = len(m)

        self.assertEqual(motif_bb(seqs, 3, 10, tam_motif), (posicoes,12))

class TestGibbs(unittest.TestCase):
    def test_Gibbs(self):
        seqs = gerar_seq(3, 10)  
        m = "ACGT"
        seqs, posicoes = inserir_motif(seqs, m)
        tam_motif = len(m)
        gibbs = GibbsSampling(seqs, 4)
        p , s = gibbs.gibbs_sampling()
        self.assertEqual( (list(p.values())), (posicoes))

if __name__ == "__main__":
    unittest.main()