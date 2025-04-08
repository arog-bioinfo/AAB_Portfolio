from typing import List

def find_ith_occurrence(lst: List[str], char: str, i: int) -> int:
    """
    Finds the index of the i-th occurrence of a character in a list.

    Args:
        lst (List[str]): The list to search in.
        char (str): The character to search for.
        i (int): The occurrence number (1-based).

    Returns:
        int: The index of the i-th occurrence.

    Raises:
        ValueError: If the character does not occur i times.
    """
    count = 0
    for idx, c in enumerate(lst):
        if c == char:
            count += 1
            if count == i:
                return idx
    raise ValueError(f"{char} does not occur {i} times.")

class BWT:
    """
    A class to perform Burrows-Wheeler Transform (BWT) operations, including:
    - BWT construction
    - Suffix Array generation
    - Pattern matching using the BWT
    """

    def __init__(self, seq: str, buildsufarray: bool = False):
        """
        Initializes the BWT class with the given sequence.

        Args:
            seq (str): The input sequence.
            buildsufarray (bool): If True, builds the suffix array.
        """
        self.seq = seq + "$" if "$" not in seq else seq
        self.bwt = self.build_bwt(self.seq)
        self.sa = self.build_suffix_array(self.seq) if buildsufarray else []

    def build_bwt(self, text: str) -> str:
        """
        Constructs the Burrows-Wheeler Transform of the sequence.

        Args:
            text (str): The input string with a terminal character.

        Returns:
            str: The BWT string.
        """
        rotations = [text[i:] + text[:i] for i in range(len(text))]
        rotations.sort()
        return ''.join(rot[-1] for rot in rotations)

    def build_suffix_array(self, text: str) -> List[int]:
        """
        Constructs the suffix array of the sequence.

        Args:
            text (str): The input string.

        Returns:
            List[int]: The suffix array as a list of starting positions.
        """
        suffixes = [(text[i:], i) for i in range(len(text))]
        suffixes.sort()
        return [suffix[1] for suffix in suffixes]

    def last_to_first(self) -> List[int]:
        """
        Builds the Last-to-First mapping used in BWT matching.

        Returns:
            List[int]: A list mapping each index in the BWT (last column)
                       to its corresponding index in the first column.
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
        Performs backward search using BWT to find the pattern.

        Args:
            pattern (str): The pattern to search for.

        Returns:
            List[int]: A list of BWT matrix row indices where the pattern occurs.
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
                return []

        return list(range(top, bottom + 1))

    def bw_matching_pos(self, pattern: str) -> List[int]:
        """
        Returns the original positions where the pattern occurs,
        using the suffix array and BWT matching.

        Args:
            pattern (str): The pattern to find.

        Returns:
            List[int]: Sorted list of starting positions of the pattern.
        """
        matches = self.bw_matching(pattern)

        if not matches:
            return []

        positions = [self.sa[m] for m in matches]
        positions.sort()
        return positions

def bwt_reverse(bwt: str) -> str:
    """
    Reverses the Burrows-Wheeler Transform and recovers the original sequence.

    Args:
        bwt (str): The BWT-transformed string.

    Returns:
        str: The original sequence (without the terminal character).
    """
    n = len(bwt)
    table = [''] * n
    for _ in range(n):
        table = sorted([bwt[i] + table[i] for i in range(n)])
    for row in table:
        if row.endswith('$'):
            return row.rstrip('$')
    return ""

# ================== TEST ====================
if __name__ == "__main__":
    seq = "BANANA"
    bwt_instance = BWT(seq, buildsufarray=True)

    res = seq, bwt_instance.bwt
    suffix = bwt_instance.sa

    pattern = "ANA"
    positions = bwt_instance.bw_matching_pos(pattern)
    padrao = pattern, positions

    print(res)                       
    print(suffix)                    
    print(padrao)                    

    original_recuperado = bwt_reverse(bwt_instance.bwt)
    print("Recovered sequence:", original_recuperado)  
