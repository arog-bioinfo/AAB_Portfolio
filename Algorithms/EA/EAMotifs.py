from EvolAlgorithm import EvolAlgorithm
from Popul import PopulInt, PopulReal
from MotifFinding import MotifFinding
from MyMotifs import MyMotifs


def createMatZeros(nl, nc):
    res = []
    for _ in range(0, nl):
        res.append([0]*nc)
    return res


def printMat(mat):
    for i in range(0, len(mat)):
        for j in range(len(mat[i])):
            print(f"{mat[i][j]:.3f}", end=' ')
        print()


class EAMotifsInt (EvolAlgorithm):
    def __init__(self, popsize, numits, noffspring, filename):
        self.motifs = MotifFinding()
        self.motifs.readFile(filename, "dna")
        indsize = len(self.motifs)
        EvolAlgorithm.__init__(self, popsize, numits, noffspring, indsize)

    def initPopul(self, indsize):
        maxvalue = self.motifs.seqSize(0) - self.motifs.motifSize
        self.popul = PopulInt(self.popsize, indsize,
                              maxvalue, [])

    def evaluate(self, indivs):
        for i in range(len(indivs)):
            ind = indivs[i]
            sol = ind.getGenes()
            fit = self.motifs.score(sol)
            ind.setFitness(fit)


class EAMotifsReal (EvolAlgorithm):
    def __init__(self, popsize, numits, noffspring, filename):
        self.motifs = MotifFinding()
        self.motifs.readFile(filename, "dna")
        
        L = self.motifs.motifSize
        A = 4
        indsize = L * A
        
        super().__init__(popsize, numits, noffspring, indsize)

    def initPopul(self, indsize):

        self.popul = PopulReal(self.popsize, indsize, lb=0.0, ub=1.0)

    def evaluate(self, indivs):
        """Evaluates each individual by constructing a PWM, finding best motif positions, and computing fitness."""
        for ind in indivs:
            pwm = self.reshape_to_pwm(ind.getGenes(), 4, self.motifs.motifSize)
            pwm = self.normalizePWM(pwm)

            # Use MyMotifs to extract motif from best positions
            motif_positions = self.findBestMotifPositions(pwm)
            motif = self.motifs.createMotifFromIndexes(motif_positions)
            motif.createPWM()  

            # Compute fitness score using MotifFinding's scoring function
            fit = self.motifs.score(motif_positions)
            ind.setFitness(fit)

    def reshape_to_pwm(self, genes, rows, cols):
        """Reshapes a flat list into a PWM matrix of dimensions rows × cols."""
        return [genes[i * cols:(i + 1) * cols] for i in range(rows)]

    def normalizePWM(self, pwm):
        """Uses MyMotifs to normalize the PWM."""
        motif = MyMotifs(pwm=pwm, alphabet="ACGT")  
        motif.createPWM()  
        return motif.pwm  

    def findBestMotifPositions(self, pwm):
        """Finds the most probable motif positions in each sequence using the PWM."""
        return [self.mostProbableMotif(seq, pwm) for seq in self.motifs.seqs]

    def mostProbableMotif(self, seq, pwm):
        """Uses MyMotifs to find the best start position for the motif."""
        motif = MyMotifs(pwm=pwm, alphabet="ACGT")
        return motif.mostProbableSeq(seq)


def test1():
    ea = EAMotifsInt(100, 500, 50, "exemploMotifs.txt")
    ea.run()
    ea.printBestSolution()


def test2():
    ea = EAMotifsReal(100, 500, 50, "exemploMotifs.txt")
    ea.run()
    ea.printBestSolution()


# test1()
test2()
