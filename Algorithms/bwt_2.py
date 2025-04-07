from typing import List, Union

class BWT:
    def __init__(self, seq: str, buildsufarray: bool = False):
        """
        Initializes the BWT class and constructs the Burrows-Wheeler Transform (BWT) 
        and optionally the Suffix Array.
        
        :param seq: The input sequence.
        :param buildsufarray: Whether to construct the Suffix Array (default: False).
        """
        self.seq = seq + "$" if "$" not in seq else seq
        self.bwt = self.build_bwt(self.seq)
        self.sa = self.build_suffix_array(self.seq) if buildsufarray else []

    def build_bwt(self, text: str) -> str:
        """
        Constructs the Burrows-Wheeler Transform (BWT) for a given sequence.
        
        :param text: The input sequence.
        :return: The BWT transformed sequence.
        """
        rotations = [text[i:] + text[:i] for i in range(len(text))]
        rotations.sort()
        return ''.join(rot[-1] for rot in rotations)

    def build_suffix_array(self, text: str) -> List[int]:
        """
        Constructs the Suffix Array for a given sequence.
        
        :param text: The input sequence.
        :return: The suffix array as a list of starting indices.
        """
        suffixes = [(text[i:], i) for i in range(len(text))]
        suffixes.sort()  # Sort suffixes lexicographically
        return [suffix[1] for suffix in suffixes]

    def last_to_first(self) -> List[int]:
        """
        Creates the mapping from the last column to the first column in the BWT matrix.
        
        :return: A list mapping each index in the last column to its corresponding index in the first column.
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
        Searches for all occurrences of a pattern in the original sequence using the BWT.
        
        :param pattern: The pattern to search for.
        :return: The indices of the rows in the BWT matrix where the pattern occurs.
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
                return []  # Pattern not found
        
        return list(range(top, bottom + 1))

    def bw_matching_pos(self, pattern: str) -> List[int]:
        """
        Searches for the starting positions of a pattern in the original sequence using the Suffix Array.
        
        :param pattern: The pattern to search for.
        :return: A sorted list of starting positions.
        """
        matches = self.bw_matching(pattern)
        
        if not matches:
            return []  # No occurrences found
        
        positions = [self.sa[m] for m in matches]
        
        positions.sort()
        
        return positions

# Helper function to find the n-th occurrence of an element in a list
def find_ith_occurrence(lst: List[str], elem: str, index: int) -> int:
    """
    Finds the n-th occurrence of a given element in a list.
    
    :param lst: The list to search in.
    :param elem: The element to find.
    :param index: The occurrence index (1-based).
    :return: The position of the n-th occurrence, or -1 if not found.
    """
    count = 0
    for i, c in enumerate(lst):
        if c == elem:
            count += 1
            if count == index:
                return i
    return -1

# Full test
if __name__ == "__main__":
    seq = "BANANA"
    bwt_instance = BWT(seq, buildsufarray=True)

    print(f'[ {seq} , {bwt_instance.bwt} ]')
    print("Suffix Array:", bwt_instance.sa)

    pattern = "ANA"
    positions = bwt_instance.bw_matching_pos(pattern)
    occurrences = bw_pattern_search(bwt, pattern)
    
    print(f"'{pattern}', {occurrences}", positions)
