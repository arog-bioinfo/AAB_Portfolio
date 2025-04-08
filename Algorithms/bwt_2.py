from typing import List, Union

class BWT:
    def __init__(self, seq: str, buildsufarray: bool = False):
        """
        Inicializa a classe BWT e constrói a transformada BWT e o Suffix Array (opcional).
        """
        self.seq = seq + "$" if "$" not in seq else seq
        self.bwt = self.build_bwt(self.seq)
        self.sa = self.build_suffix_array(self.seq) if buildsufarray else []

    def build_bwt(self, text: str) -> str:
        """
        Constrói a transformada de Burrows-Wheeler (BWT) para uma sequência.
        """
        rotations = [text[i:] + text[:i] for i in range(len(text))]
        rotations.sort()
        return ''.join(rot[-1] for rot in rotations)

    def build_suffix_array(self, text: str) -> List[int]:
        """
        Constrói o Suffix Array para uma sequência.
        """
        suffixes = [(text[i:], i) for i in range(len(text))]
        suffixes.sort()  # Ordena os sufixos lexicograficamente
        return [suffix[1] for suffix in suffixes]

    def last_to_first(self) -> List[int]:
        """
        Cria o mapeamento da última coluna para a primeira coluna na matriz BWT.
        """
        first_col = sorted(self.bwt)
        last_to_first_mapping = []
        
        count = {}
        for i, char in enumerate(self.bwt):
            count[char] = count.get(char, 0) + 1
            occurrence_index = count[char]
            pos_in_first_col = find_ith_occurrence(first_col, char, occurrence_index)
            last_to_first_mapping.append(pos_in_first_col)
        
        return last_to_first_mapping

    def bw_matching(self, pattern: str) -> List[int]:
        """
        Procura todas as ocorrências de um padrão na sequência original usando BWT.
        
        Retorna os índices das linhas da matriz M onde o padrão ocorre.
        """
        lf_mapping = self.last_to_first()
        
        top = 0
        bottom = len(self.bwt) - 1
        
        while top <= bottom and pattern:
            symbol = pattern[-1]
            pattern = pattern[:-1]
            
            sub_bwt = self.bwt[top:bottom + 1]
            
            if symbol in sub_bwt:
                top_index = sub_bwt.index(symbol) + top
                bottom_index = bottom - sub_bwt[::-1].index(symbol)
                
                top = lf_mapping[top_index]
                bottom = lf_mapping[bottom_index]
            else:
                return []  # Padrão não encontrado
        
        return list(range(top, bottom + 1))

    def bw_matching_pos(self, pattern: str) -> List[int]:
        """
        Procura as posições iniciais do padrão na sequência original usando o Suffix Array.
        
        Retorna uma lista com as posições iniciais ordenadas.
        """
        matches = self.bw_matching(pattern)
        
        if not matches:
            return []  # Nenhuma ocorrência encontrada
        
        positions = [self.sa[m] for m in matches]
        
        positions.sort()
        
        return positions
    
    #  Função auxiliar para encontrar a n-ésima ocorrência de um símbolo em uma lista
    def find_ith_occurrence(lst: List[str], elem: str, index: int) -> int:
        count = 0
        for i, c in enumerate(lst):
            if c == elem:
                count += 1
                if count == index:
                    return i
        return -1


# Teste completo
if __name__ == "__main__":
    seq = "BANANA"
    bwt_instance = BWT(seq, buildsufarray=True)

    print(f'[ {seq} , {bwt_instance.bwt} ]')
    print("Suffix Array:", bwt_instance.sa)

    pattern = "ANA"
    positions = bwt_instance.bw_matching_pos(pattern)
    count = len(positions)

    
    print(f"'{pattern}', {count}", positions)
