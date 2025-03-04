#UNITEST
import random

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
#BB 
import unittest

class TestBranchAndBound(unittest.TestCase):
    
    def test_score_bb(self):
        seqs = ["ACGT", "ACGT", "ACGT"]
        offset = [0, 0, 0]
        tam_motif = 4
        self.assertEqual(score_bb(seqs, offset, tam_motif), 4)
        
        seqs = ["ACGT", "TGCA", "ACGT"]
        offset = [0, 0, 0]
        tam_motif = 4
        self.assertEqual(score_bb(seqs, offset, tam_motif), 2)
        
        seqs = ["ACGT", "ACGT", "ACGT"]
        offset = [1, 1, 1]
        tam_motif = 3
        self.assertEqual(score_bb(seqs, offset, tam_motif), 3)
    
    def test_motif_bb(self):
        seqs = ["ACGT", "ACGT", "ACGT"]
        num_seqs = 3
        tam_seq = 4
        tam_motif = 2
        melhor_pos, maior_s = motif_bb(seqs, num_seqs, tam_seq, tam_motif)
        self.assertEqual(melhor_pos, [0, 0, 0])
        self.assertEqual(maior_s, 6)
        
        seqs = ["ACGT", "TGCA", "ACGT"]
        num_seqs = 3
        tam_seq = 4
        tam_motif = 2
        melhor_pos, maior_s = motif_bb(seqs, num_seqs, tam_seq, tam_motif)
        self.assertEqual(melhor_pos, [0, 0, 0])
        self.assertEqual(maior_s, 4)
        
        seqs = ["ACGT", "ACGT", "ACGT"]
        num_seqs = 3
        tam_seq = 4
        tam_motif = 3
        melhor_pos, maior_s = motif_bb(seqs, num_seqs, tam_seq, tam_motif)
        self.assertEqual(melhor_pos, [0, 0, 0])
        self.assertEqual(maior_s, 9)
        
        seqs = ["ACGT", "TGCA", "CGTA"]
        num_seqs = 3
        tam_seq = 4
        tam_motif = 2
        melhor_pos, maior_s = motif_bb(seqs, num_seqs, tam_seq, tam_motif)
        self.assertEqual(melhor_pos, [0, 0, 0])
        self.assertEqual(maior_s, 4)

if __name__ == "__main__":
    unittest.main()
    
    
#Procura exaustiva
import unittest

class TestExhaustiveSearch(unittest.TestCase):
    
    def test_proximo(self):
        # Testa a função `proximo` para verificar se ela avança corretamente os offsets
        offset = [0, 0, 0]
        num_seqs = 3
        limite = 4
        
        # Primeira iteração
        novo_offset = proximo(offset, num_seqs, limite)
        self.assertEqual(novo_offset, [1, 0, 0])
        
        # Segunda iteração
        novo_offset = proximo(novo_offset, num_seqs, limite)
        self.assertEqual(novo_offset, [2, 0, 0])
        
        # Testa o "carry-over" (quando um offset atinge o limite)
        offset = [3, 0, 0]
        novo_offset = proximo(offset, num_seqs, limite)
        self.assertEqual(novo_offset, [0, 1, 0])
        
        # Testa o caso em que todos os offsets atingem o limite
        offset = [3, 3, 3]
        novo_offset = proximo(offset, num_seqs, limite)
        self.assertIsNone(novo_offset)
    
    def test_enumerar(self):
        # Testa a função `enumerar` para verificar se gera todos os estados possíveis
        num_seqs = 2
        tam_seqs = 4
        tam_motifs = 2
        limite = tam_seqs - tam_motifs + 1  # limite = 3
        
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
        # Testa a função `score` para verificar se calcula corretamente o score
        seqs = ["ACGT", "ACGT", "ACGT"]
        offset = [0, 0, 0]
        tam_motif = 4
        self.assertEqual(score(seqs, offset, tam_motif), 4)
        
        seqs = ["ACGT", "TGCA", "ACGT"]
        offset = [0, 0, 0]
        tam_motif = 4
        self.assertEqual(score(seqs, offset, tam_motif), 2)
        
        seqs = ["ACGT", "ACGT", "ACGT"]
        offset = [1, 1, 1]
        tam_motif = 3
        self.assertEqual(score(seqs, offset, tam_motif), 3)
    
    def test_motif(self):
        # Testa a função `motif` para verificar se encontra o melhor 
        seqs = ["ACGT", "ACGT", "ACGT"]
        num_seqs = 3
        tam_seq = 4
        tam_motif = 2
        melhor_pos, maior_s = motif(seqs, num_seqs, tam_seq, tam_motif)
        self.assertEqual(melhor_pos, (0, 0, 0))
        self.assertEqual(maior_s, 6)
        
        seqs = ["ACGT", "TGCA", "ACGT"]
        num_seqs = 3
        tam_seq = 4
        tam_motif = 2
        melhor_pos, maior_s = motif(seqs, num_seqs, tam_seq, tam_motif)
        self.assertEqual(melhor_pos, (0, 0, 0))
        self.assertEqual(maior_s, 4)
        
        seqs = ["ACGT", "ACGT", "ACGT"]
        num_seqs = 3
        tam_seq = 4
        tam_motif = 3
        melhor_pos, maior_s = motif(seqs, num_seqs, tam_seq, tam_motif)
        self.assertEqual(melhor_pos, (0, 0, 0))
        self.assertEqual(maior_s, 9)
        
        seqs = ["ACGT", "TGCA", "CGTA"]
        num_seqs = 3
        tam_seq = 4
        tam_motif = 2
        melhor_pos, maior_s = motif(seqs, num_seqs, tam_seq, tam_motif)
        self.assertEqual(melhor_pos, (0, 0, 0))
        self.assertEqual(maior_s, 4)

if __name__ == "__main__":
    unittest.main()
    
    
#Gibs

class TestGibbsSampling(unittest.TestCase):
    def test_gibbs_sampling(self):
        # Todas as sequências são iguais
        seqs = ["ACGT", "ACGT", "ACGT"]
        num_seqs = 3
        tam_seq = 4
        tam_motif = 2
        num_iterations = 10
        
        offsets, score_final = gibbs_sampling(seqs, num_seqs, tam_seq, tam_motif, num_iterations)
        self.assertEqual(score_final, 6)  # Score máximo esperado
        self.assertEqual(offsets, [0, 0, 0])  # Melhor posição esperada
        
        # Sequências diferentes
        seqs = ["ACGT", "TGCA", "CGTA"]
        num_seqs = 3
        tam_seq = 4
        tam_motif = 2
        num_iterations = 10
        
        offsets, score_final = gibbs_sampling(seqs, num_seqs, tam_seq, tam_motif, num_iterations)
        self.assertTrue(score_final >= 4)  # Score mínimo esperado
        self.assertTrue(all(0 <= o <= tam_seq - tam_motif for o in offsets))  # Offsets válidos
        
        # Motif maior
        seqs = ["ACGTACGT", "ACGTACGT", "ACGTACGT"]
        num_seqs = 3
        tam_seq = 8
        tam_motif = 4
        num_iterations = 10
        
        offsets, score_final = gibbs_sampling(seqs, num_seqs, tam_seq, tam_motif, num_iterations)
        self.assertEqual(score_final, 12)  # Score máximo esperado
        self.assertEqual(offsets, [0, 0, 0])  # Melhor posição esperada

if __name__ == "__main__":
    unittest.main()
