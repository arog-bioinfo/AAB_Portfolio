import random

class GibbsSampling:
    def __init__(self, sequences, motif_length, pseudo=1):
        self.seqs = sequences
        self.w = motif_length
        self.pseudo = pseudo
        self.n = len(sequences)
        self.t = len(sequences[0])

    def random_init_positions(self):
        return {seq: random.randint(0, self.t - self.w) for seq in self.seqs}

    def create_motifs(self, positions):
        return [seq[positions[seq]:positions[seq] + self.w] for seq in positions]

    def pwm(self, motifs):
        bases = 'ATCG'
        pwm_matrix = []
        for pos in zip(*motifs):
            counts = {base: pos.count(base) + self.pseudo for base in bases}
            total = sum(counts.values())
            pwm_matrix.append({base: counts[base] / total for base in bases})
        return pwm_matrix

    def prob_seq(self, seq, pwm):
        prob = 1.0
        for i, base in enumerate(seq):
            prob *= pwm[i][base]
        return prob

    def prob_positions(self, seq, pwm):
        probabilities = [self.prob_seq(seq[i:i+self.w], pwm) for i in range(self.t - self.w + 1)]
        total = sum(probabilities)
        return [p / total for p in probabilities]

    def roulette_wheel(self, probabilities):
        r = random.uniform(0, sum(probabilities))
        s = 0
        for i, p in enumerate(probabilities):
            s += p
            if s >= r:
                return i
        return len(probabilities) - 1

    def gibbs_sampling(self, max_iter=100, threshold=50):
        positions = self.random_init_positions()
        best_positions, best_score, count = positions.copy(), 0, 0
        
        for _ in range(max_iter):
            seq_to_remove = random.choice(self.seqs)
            temp_positions = positions.copy()
            temp_positions.pop(seq_to_remove)
            motifs = self.create_motifs(temp_positions)
            pwm = self.pwm(motifs)
            probabilities = self.prob_positions(seq_to_remove, pwm)
            new_pos = self.roulette_wheel(probabilities)
            positions[seq_to_remove] = new_pos
            
            score = sum(max(col.values()) for col in self.pwm(self.create_motifs(positions)))
            if score > best_score:
                best_positions, best_score, count = positions.copy(), score, 0
            else:
                count += 1
                if count >= threshold:
                    break
        
        return best_positions, best_score

    def print_motif(self, positions):
        for seq in self.seqs:
            if seq in positions:
                start = positions[seq]
                print(seq[start:start+self.w])

if __name__ == "__main__":
    n = int(input('Introduza o tamanho do motif: '))
    no = int(input('Introduza o número de sequências: '))
    seqs = [input(f'Introduza a sequência {i+1}: ').strip() for i in range(no)]
    
    gibbs = GibbsSampling(seqs, n)
    positions, score = gibbs.gibbs_sampling()
    print(f"Posições: {positions}\nPontuação: {score}\n")
    print("Motif encontrado:")
    gibbs.print_motif(positions)
