import random

class GibbsSampling:
    def __init__(self, sequences, motif_length, pseudo=1):
        '''
        Inicializa a classe GibbsSampling.
        
        Parâmetros:
            sequences (list): Lista de sequências de DNA.
            motif_length (int): Tamanho do motif a ser descoberto.
            pseudo (float): Pequena constante para evitar probabilidades zero na PWM (padrão=1).
        '''
        self.seqs = sequences
        self.w = motif_length
        self.pseudo = pseudo
        self.n = len(sequences)
        self.t = len(sequences[0])

    def random_init_positions(self):
        '''
        Seleciona posições iniciais aleatórias para os motifs em cada sequência.
        
        Retorna:
            dict: Dicionário com sequências como chaves e posições iniciais como valores.
        '''
        return {seq: random.randint(0, self.t - self.w) for seq in self.seqs}

    def create_motifs(self, positions):
        '''
        Extrai motifs das sequências com base nas posições fornecidas.
        
        Parâmetros:
            positions (dict): Dicionário com posições iniciais para cada sequência.
        
        Retorna:
            list: Lista de motifs extraídos.
        '''
        return [seq[positions[seq]:positions[seq] + self.w] for seq in positions]

    def pwm(self, motifs):
        '''
        Calcula a Matriz de Pesos de Posição (PWM) a partir de um conjunto de motifs.
        
        Parâmetros:
            motifs (list): Lista de motifs.
        
        Retorna:
            list: PWM representada como uma lista de dicionários.
        '''
        bases = 'ATCG'
        pwm_matrix = []
        for pos in zip(*motifs):
            counts = {base: pos.count(base) + self.pseudo for base in bases}
            total = sum(counts.values())
            pwm_matrix.append({base: counts[base] / total for base in bases})
        return pwm_matrix

    def prob_seq(self, seq, pwm):
        '''
        Calcula a probabilidade de uma sequência dada a PWM.
        
        Parâmetros:
            seq (str): Sequência de DNA.
            pwm (list): Matriz de Pesos de Posição.
        
        Retorna:
            float: Probabilidade da sequência dada a PWM.
        '''
        prob = 1.0
        for i, base in enumerate(seq):
            prob *= pwm[i][base]
        return prob

    def prob_positions(self, seq, pwm):
        '''
        Calcula as probabilidades de cada posição da sequência conter o motif.
        
        Parâmetros:
            seq (str): Sequência de DNA.
            pwm (list): Matriz de Pesos de Posição.
        
        Retorna:
            list: Lista de probabilidades para cada posição possível do motif.
        '''
        probabilities = [self.prob_seq(seq[i:i+self.w], pwm) for i in range(self.t - self.w + 1)]
        total = sum(probabilities)
        return [p / total for p in probabilities]

    def roulette_wheel(self, probabilities):
        '''
        Realiza uma seleção estocástica usando o método da roleta.
        
        Parâmetros:
            probabilities (list): Lista de probabilidades.
        
        Retorna:
            int: Índice selecionado com base nas probabilidades.
        '''
        r = random.uniform(0, sum(probabilities))
        s = 0
        for i, p in enumerate(probabilities):
            s += p
            if s >= r:
                return i
        return len(probabilities) - 1

    def gibbs_sampling(self, max_iter=100, threshold=50):
        '''
        Executa o algoritmo de Amostragem de Gibbs para encontrar motifs.
        
        Parâmetros:
            max_iter (int): Número máximo de iterações (padrão=100).
            threshold (int): Critério de parada baseado na ausência de melhoria (padrão=50).
        
        Retorna:
            tuple: Melhores posições dos motifs e maior pontuação encontrada.
        '''
        positions = self.random_init_positions()
        best_positions = positions.copy()
        best_score = 0
        count = 0

        for _ in range(max_iter):
            if count < threshold:
                seq_to_remove = random.choice(self.seqs)
                temp_positions = positions.copy()
                temp_positions.pop(seq_to_remove)
                motifs = self.create_motifs(temp_positions)
                pwm = self.pwm(motifs)
                probabilities = self.prob_positions(seq_to_remove, pwm)
                new_pos = self.roulette_wheel(probabilities)
                positions[seq_to_remove] = new_pos

                current_motifs = self.create_motifs(positions)
                pwm_all = self.pwm(current_motifs)
                score = sum(max(col.values()) for col in pwm_all)

                if score > best_score:
                    best_score = score
                    best_positions = positions.copy()
                    count = 0
                else:
                    count += 1

        return best_positions, best_score

    def print_motif(self, positions):
        '''
        Imprime o alinhamento do motif descoberto.
        
        Parâmetros:
            positions (dict): Dicionário das melhores posições dos motifs.
        '''
        for seq in self.seqs:
            if seq in positions:
                start = positions[seq]
                print(seq[start:start+self.w])


def motif_gibbs(seqs, num_seqs, tam_seq, tam_motif, max_iter=100, threshold=50):
    assert len(seqs) == num_seqs
    assert all(len(s) == tam_seq for s in seqs)

    gibbs = GibbsSampling(seqs, tam_motif)
    positions, score = gibbs.gibbs_sampling(max_iter=max_iter, threshold=threshold)
    return positions, score
