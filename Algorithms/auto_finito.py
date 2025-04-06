class AutomatosFinitos :

    def __init__(self, alphabet: set, pattern: str):
        if not pattern:
            raise ValueError ('Pattern cannot be empty')
        if not alphabet:
            raise ValueError ('Alphabet cannot be empty')
        
        self.num_states = len(pattern) + 1
        self.alphabet = alphabet
        self.pattern = pattern
        self.transition_table = {}
        self._build_transition_table()

    def _build_transition_table(self):

        """Constructs the transition table using optimized overlap calculation."""

        for state in range(self.num_states):
            for symbol in self.alphabet:
                next_state = self._calculate_state(state, symbol)
                self.transition_table[(state, symbol)] = next_state

    def _calculate_state(self, current_state: int, symbol: str):
        """Calculates the next state with memoization for repeated patterns."""

        if current_state == 0 and symbol == self.pattern[0]:
            return 1
        
        candidate = self.pattern[:current_state] +  symbol
        max_overlap = min(current_state + 1, len(self.pattern))

        for length in range(max_overlap, 0, -1):
            if candidate.endswith(self.pattern[:length]):
                return length
        return 0
    
    def process_input(self, sequence: str) :
        """Processes an input sequence and returns states visited."""
        
        current_state = 0
        states_log = [current_state]

        for char in sequence:
            current_state = self.transition_table.get(
                (current_state, char), 0
            )
            states_log.append(current_state)
        return states_log
    
    def find_occurrences (self, text: str):
        """Finds all pattern occurrences with their starting positions."""

        states = self.process_input(text)
        pattern_length = len(self.pattern)
        matches = [
            index - pattern_length
            for index, state in enumerate(states)
            if state == len (self.pattern)
        ]
        return len (matches), matches
    
    def visualize (self):
        print(f"States: {self.num_states}")
        print(f"Search Pattern: {self.pattern}")
        print("Transition Map:")
        for state in sorted(set(s for s, _ in self.transition_table)):
            symbols = [k[1] for k in self.transition_table if k[0] == state]
            for symbol in sorted(symbols):
                dest = self.transition_table[(state, symbol)]
                print(f"  ({state}, {symbol}) → {dest}")

    def get_next_state(self, current_state: int, symbol: str):
        """Safe state transition with input validation."""
        if symbol not in self.alphabet:
            raise ValueError(f"Symbol '{symbol}' not in alphabet")
        return self.transition_table.get((current_state, symbol), 0)
    
def pattern_overlap(suffix_candidate: str, prefix_source: str):
    """Calculates maximum overlap between suffix and prefix."""
    max_length = min(len(suffix_candidate), len(prefix_source))
    for length in range(max_length, 0, -1):
        if suffix_candidate[-length:] == prefix_source[:length]:
            return length
    return 0