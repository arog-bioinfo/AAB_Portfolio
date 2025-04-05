"""
Gibbs Sampling implementation for motif discovery in biological sequences.

This module provides a GibbsSampling class and helper functions to identify 
conserved motifs in DNA sequences using a probabilistic sampling approach.
"""

import random

class GibbsSampling:
    """Implements Gibbs Sampling algorithm for motif discovery.
    
    Attributes:
        seqs (list): List of DNA sequences to analyze
        w (int): Motif length to search for
        pseudo (int): Pseudocount value to prevent zero probabilities
        n (int): Number of sequences
        t (int): Length of each sequence (assumes equal-length sequences)
    """
    
    def __init__(self, sequences, motif_length, pseudo=1):
        """Initialize GibbsSampling parameters.
        
        Args:
            sequences (list): List of DNA sequences (must be equal length)
            motif_length (int): Length of motif to search for
            pseudo (int, optional): Pseudocount for Laplace smoothing. Default=1
        """
        self.seqs = sequences
        self.w = motif_length
        self.pseudo = pseudo
        self.n = len(sequences)
        self.t = len(sequences[0])

    def random_init_positions(self):
        """Generate random starting positions for motif search.
        
        Returns:
            dict: Initial positions mapping {sequence: start_index}
        """
        return {seq: random.randint(0, self.t - self.w) for seq in self.seqs}

    def create_motifs(self, positions):
        """Extract motifs based on current positions.
        
        Args:
            positions (dict): Current motif positions {seq: start_index}
            
        Returns:
            list: List of motif sequences corresponding to current positions
        """
        return [seq[positions[seq]:positions[seq] + self.w] for seq in positions]

    def pwm(self, motifs):
        """Build Position Weight Matrix (PWM) from motifs.
        
        Args:
            motifs (list): List of motif sequences to analyze
            
        Returns:
            list: PWM matrix as list of dictionaries {base: probability}
        """
        bases = 'ATCG'
        pwm_matrix = []
        for pos in zip(*motifs):
            counts = {base: pos.count(base) + self.pseudo for base in bases}
            total = sum(counts.values())
            pwm_matrix.append({base: counts[base] / total for base in bases})
        return pwm_matrix

    def prob_seq(self, seq, pwm):
        """Calculate probability of sequence given PWM.
        
        Args:
            seq (str): DNA sequence to evaluate
            pwm (list): Position Weight Matrix
            
        Returns:
            float: Probability score for the sequence
        """
        prob = 1.0
        for i, base in enumerate(seq):
            prob *= pwm[i][base]
        return prob

    def prob_positions(self, seq, pwm):
        """Calculate normalized probabilities for all motif positions in sequence.
        
        Args:
            seq (str): Sequence to analyze
            pwm (list): Current Position Weight Matrix
            
        Returns:
            list: Normalized probabilities for each possible starting position
        """
        probabilities = [self.prob_seq(seq[i:i+self.w], pwm) for i in range(self.t - self.w + 1)]
        total = sum(probabilities)
        return [p / total for p in probabilities]

    def roulette_wheel(self, probabilities):
        """Select position using probabilistic roulette wheel selection.
        
        Args:
            probabilities (list): List of position probabilities
            
        Returns:
            int: Selected position index
        """
        r = random.uniform(0, sum(probabilities))
        s = 0
        for i, p in enumerate(probabilities):
            s += p
            if s >= r:
                return i
        return len(probabilities) - 1

    def gibbs_sampling(self, max_iter=100, threshold=50):
        """Execute Gibbs Sampling algorithm.
        
        Args:
            max_iter (int): Maximum iterations. Default=100
            threshold (int): Convergence threshold (iterations without improvement). Default=50
            
        Returns:
            tuple: (best_positions dict, best_score float)
        """
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
        """Display motifs based on final positions.
        
        Args:
            positions (dict): Final motif positions {seq: start_index}
        """
        for seq in self.seqs:
            if seq in positions:
                start = positions[seq]
                print(seq[start:start+self.w])


def motif_gibbs(seqs, num_seqs, tam_seq, tam_motif, max_iter=100, threshold=50):
    """Wrapper function for Gibbs Sampling motif discovery.
    
    Args:
        seqs (list): List of DNA sequences
        num_seqs (int): Expected number of sequences (validation)
        tam_seq (int): Expected sequence length (validation)
        tam_motif (int): Motif length to search for
        max_iter (int): Maximum iterations. Default=100
        threshold (int): Convergence threshold. Default=50
        
    Returns:
        tuple: (positions dict, score float)
        
    Raises:
        AssertionError: If input validation fails
    """
    assert len(seqs) == num_seqs
    assert all(len(s) == tam_seq for s in seqs)

    gibbs = GibbsSampling(seqs, tam_motif)
    positions, score = gibbs.gibbs_sampling(max_iter=max_iter, threshold=threshold)
    return positions, score
