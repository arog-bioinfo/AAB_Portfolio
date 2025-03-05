# UNITEST
import random
import unittest

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

# (score_bb, motif_bb, proximo, enumerar, score, motif, gibbs_sampling) devem ser definidas 


class TestBranchAndBound(unittest.TestCase):
    
    def test_score_bb(self):
        seqs = gerar_seq(3, 10)  
        motif = "ACGT"
        seqs, posicoes = inserir_motif(seqs, motif)
        tam_motif = len(motif)
        
        # Verifica se o score é consistente com o motif inserido
        self.assertEqual(score_bb(seqs, posicoes, tam_motif), 4)
    
    def test_motif_bb(self):
        seqs = gerar_seq(3, 10)  
        motif = "ACGT"
        seqs, posicoes = inserir_motif(seqs, motif)
        num_seqs = 3
        tam_seq = 10
        tam_motif = len(motif)
        
        melhor_pos, maior_s = motif_bb(seqs, num_seqs, tam_seq, tam_motif)
        
        # Verifica se o score máximo foi alcançado
        self.assertEqual(maior_s, 12)
        
        # Verifica se as posições encontradas são válidas
        for pos in melhor_pos:
            self.assertTrue(0 <= pos <= tam_seq - tam_motif)

class TestExhaustiveSearch(unittest.TestCase):
    
    def test_proximo(self):
        offset = [0, 0, 0]
        num_seqs = 3
        limite = 4
        
        novo_offset = proximo(offset, num_seqs, limite)
        self.assertEqual(novo_offset, [1, 0, 0])
        
        novo_offset = proximo(novo_offset, num_seqs, limite)
        self.assertEqual(novo_offset, [2, 0, 0])
        
        offset = [3, 0, 0]
        novo_offset = proximo(offset, num_seqs, limite)
        self.assertEqual(novo_offset, [0, 1, 0])
        
        offset = [3, 3, 3]
        novo_offset = proximo(offset, num_seqs, limite)
        self.assertIsNone(novo_offset)
    
    def test_enumerar(self):
        num_seqs = 2
        tam_seqs = 4
        tam_motifs = 2
        limite = tam_seqs - tam_motifs + 1
        
        estados_esperados = [
            (0, 0),
            (1, 0),
            (2, 0),
            (0, 1),
            (1, 1),
            (2, 1),
            (0, 2),
            (1, 2),
            (2, 2),
        ]
        
        estados_gerados = list(enumerar(num_seqs, tam_seqs, tam_motifs))
        self.assertEqual(estados_gerados, estados_esperados)
    
    def test_score(self):
        seqs = gerar_seq(3, 10)  
        motif = "ACGT"
        seqs, posicoes = inserir_motif(seqs, motif)
        tam_motif = len(motif)
        
        # Verifica se o score é consistente com o motif inserido
        self.assertEqual(score(seqs, posicoes, tam_motif), 4)
    
    def test_motif(self):
        seqs = gerar_seq(3, 10)  
        motif = "ACGT"
        seqs, posicoes = inserir_motif(seqs, motif)
        num_seqs = 3
        tam_seq = 10
        tam_motif = len(motif)
        
        melhor_pos, maior_s = motif(seqs, num_seqs, tam_seq, tam_motif)
        
        # Verifica se o score máximo foi alcançado
        self.assertEqual(maior_s, 12)
        
        # Verifica se as posições encontradas são válidas
        for pos in melhor_pos:
            self.assertTrue(0 <= pos <= tam_seq - tam_motif)

class TestGibbsSampling(unittest.TestCase):
    def test_gibbs_sampling(self):
        seqs = gerar_seq(3, 10)  
        motif = "ACGT"
        seqs, posicoes = inserir_motif(seqs, motif)
        num_seqs = 3
        tam_seq = 10
        tam_motif = len(motif)
        num_iterations = 10
        
        offsets, score_final = gibbs_sampling(seqs, num_seqs, tam_seq, tam_motif, num_iterations)
        
        # Verifica se o score é consistente com o motif inserido
        self.assertEqual(score_final, 12)
        
        # Verifica se as posições encontradas são válidas
        for offset in offsets:
            self.assertTrue(0 <= offset <= tam_seq - tam_motif)

if __name__ == "__main__":
    unittest.main()
