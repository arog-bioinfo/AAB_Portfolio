class Counter():
    def __init__(self):
        self.count = 0
        
    def increment(self):
        self.count += 1
    
    def __int__(self):
        return self.count

#Procura Exaustiva

def proximo(offset, num_seqs, limite):
    pos = 0

    while pos < num_seqs:
        offset[pos] +=1

        if offset[pos] < limite:
            return offset

        offset[pos] = 0
        pos += 1

def enumerar(num_seqs,tam_seqs,tam_motifs):
    estado = [0] * num_seqs
    limite = tam_seqs - tam_motifs + 1

    while estado is not None:
        yield tuple(estado)
        estado = proximo(estado, num_seqs, limite)

def score(seqs, offset, tam_motif):
    snips = [s[p: p + tam_motif] for p,s in zip(offset, seqs)]
    return sum(max([col.count(x) for x in col]) for col in zip(*snips))

def motif(seqs, num_seqs, tam_seq, tam_motif, counting = False):

    assert all(len(s) == tam_seq for s in seqs)
    
    C = Counter() if counting else None

    maior_s = 0
    for estado in enumerar(num_seqs, tam_seq, tam_motif):
        if C:
            C.increment()
        atual = score(seqs, estado, tam_motif)
        if atual > maior_s:
            melhor_pos = estado
            maior_s = atual

    return (melhor_pos, maior_s, int(C)) if counting else (melhor_pos, maior_s)
