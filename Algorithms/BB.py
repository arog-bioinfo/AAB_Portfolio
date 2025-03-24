#Branch and Bond


def score_bb(seqs, offset, tam_motif):
    snips = [s[p: p + tam_motif] for p, s in zip(offset, seqs)]
    return sum(max(col.count(x) for x in set(col)) for col in zip(*snips))

def branch_and_bound(offset, num_seqs, limite, tam_motifs, estado_global, seqs, level=0, counting = False):
    """ Implementação de Branch and Bound para busca do melhor motivo """
    
    if level == num_seqs:  # Todos os offsets foram preenchidos
        atual = score_bb(seqs, offset, tam_motifs)
        
        # Se for melhor que o maior score encontrado, atualizar
        if atual > estado_global["maior_s"]:
            estado_global["maior_s"] = atual
            estado_global["melhor_pos"] = offset[:]
        return

    # Gerar próximos candidatos (Branch)
    for pos in range(limite):
        offset[level] = pos
        atual = score_bb(seqs, offset[:level+1] + [0] * (num_seqs - level - 1), tam_motifs)

        # Poda (Bound): Se o score máximo possível já for pior, descartar o ramo
        limite_superior = atual + (num_seqs - level - 1) * tam_motifs
        if limite_superior > estado_global["maior_s"]:
            counter.increment()
            branch_and_bound(offset, num_seqs, limite, tam_motifs, estado_global, seqs, level+1)

def motif_bb(seqs, num_seqs, tam_seq, tam_motif, counting = False):
    assert all(len(s) == tam_seq for s in seqs)
    
    estado_global = {"maior_s": 0, "melhor_pos": None}
    offset = [0] * num_seqs
    limite = tam_seq - tam_motif + 1
    if counting == True:
        counter = Counter()

    branch_and_bound(offset, num_seqs, limite, tam_motif, estado_global, seqs, counting)

    if counting == False:
        return estado_global["melhor_pos"], estado_global["maior_s"]
    
    return estado_global["melhor_pos"], estado_global["maior_s"], counter