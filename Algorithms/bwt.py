from typing import List, Union

def bwt_transform(seq: str) -> List[Union[str, int]]:
    """
    Executa a transformação de Burrows-Wheeler (BWT) numa sequência.

    Parâmetros:
        seq (str): A sequência de entrada (sem o terminador).

    Retorna:
        List[Union[str, int]]: Uma lista com:
            - a sequência original;
            - a sequência transformada BWT;
            - o índice da sequência original na lista de rotações ordenadas.
    """
    original = seq 
    if '$' not in seq:
        seq += '$'

    rotacoes = [seq[i:] + seq[:i] for i in range(len(seq))]
    rotacoes.sort()
    bwt = ''.join(rot[-1] for rot in rotacoes)
    pos_original = rotacoes.index(seq)

    return [original, bwt, pos_original] 


def bwt_reverse(bwt: str, original_index: int) -> str:
    """
    Reverte a transformação BWT e recupera a sequência original.

    Parâmetros:
        bwt (str): A sequência transformada pela BWT.
        original_index (int): Índice da sequência original na matriz ordenada.

    Retorna:
        str: A sequência original (sem o terminador).
    """
    n = len(bwt)
    table = [''] * n

    for _ in range(n):
        # Adiciona bwt como prefixo a cada linha da tabela
        table = sorted([bwt[i] + table[i] for i in range(n)])

    # A linha original é aquela no índice original_index
    original = table[original_index]

    # Remover o símbolo de fim '$' antes de devolver
    return original.rstrip('$')

if __name__ == "__main__":
    resultado = bwt_transform("banana")
    print(resultado)
    original_recuperado = bwt_reverse(resultado[1], resultado[2])
    print("Sequência recuperada:", original_recuperado)


