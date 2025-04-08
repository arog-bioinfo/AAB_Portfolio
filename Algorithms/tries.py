class Trie:
    def __init__(self) -> None:
        """
        Initializes an empty Trie with a root node (node 0).
        """
        self.nodes = {0: {}}  # root node
        self.num = 0          # node counter

    def print_trie(self) -> None:
        """
        Prints the structure of the Trie for debugging.
        """
        for k in self.nodes.keys():
            print(k, "->", self.nodes[k])

    def add_node(self, origin: int, symbol: str) -> None:
        """
        Adds a new node to the Trie.

        :param origin: ID of the node where the edge starts
        :param symbol: character for the edge label
        """
        self.num += 1
        self.nodes[origin][symbol] = self.num
        self.nodes[self.num] = {}

    def add_pattern(self, p: str) -> None:
        """
        Inserts a single pattern into the Trie.

        :param p: pattern to insert
        """
        pos = 0
        node = 0
        while pos < len(p):
            if p[pos] not in self.nodes[node]:
                self.add_node(node, p[pos])
            node = self.nodes[node][p[pos]]
            pos += 1

    def trie_from_patterns(self, pats: list[str]) -> None:
        """
        Inserts a list of patterns into the Trie.

        :param pats: list of patterns
        """
        for p in pats:
            self.add_pattern(p)

    def prefix_trie_match(self, text: str) -> str | None:
        """
        Searches for the longest prefix in the Trie that matches the beginning of text.

        :param text: input string to match
        :return: matched prefix string or None
        """
        pos = 0
        match = ""
        node = 0
        while pos < len(text):
            char = text[pos]
            if char in self.nodes[node]:
                node = self.nodes[node][char]
                match += char
                if self.nodes[node] == {}:  # terminal node
                    return match
                pos += 1
            else:
                return None
        return None

    def trie_matches(self, text: str) -> list[tuple[int, str]]:
        """
        Searches for all patterns in the Trie that match substrings in text.

        :param text: the input string
        :return: list of tuples (start_position, matched_string)
        """
        res = []
        for i in range(len(text)):
            match = self.prefix_trie_match(text[i:])
            if match is not None:
                res.append((i, match))
        return res


class SuffixTree:
    def __init__(self):
        """
        Initializes the suffix tree with a root node.
        Each node is represented as a tuple: (leaf_index, children_dict).
        """
        self.nodes = {0: (-1, {})}  # node_id: (leaf_index or -1, {symbol: child_node_id})
        self.node_count = 0
        self.text = ""

    def _add_node(self, parent, symbol, leaf_index=-1):
        """
        Adds a new node as a child of the given parent node.

        :param parent: index of the parent node
        :param symbol: character representing the edge to the new node
        :param leaf_index: suffix index if it's a leaf node, otherwise -1
        """
        self.node_count += 1
        self.nodes[parent][1][symbol] = self.node_count
        self.nodes[self.node_count] = (leaf_index, {})

    def _add_suffix(self, suffix, index):
        """
        Adds a suffix to the tree.

        :param suffix: suffix string to be inserted
        :param index: starting index of the suffix in the original text
        """
        node = 0
        for i, char in enumerate(suffix):
            if char not in self.nodes[node][1]:
                is_leaf = (i == len(suffix) - 1)
                self._add_node(node, char, index if is_leaf else -1)
            node = self.nodes[node][1][char]

    def build(self, text):
        """
        Builds the suffix tree from the given input string.

        :param text: input string to build the suffix tree from
        """
        self.text = text + "$"
        for i in range(len(self.text)):
            self._add_suffix(self.text[i:], i)

    def _get_leaves_below(self, node_id):
        """
        Recursively collects all leaf indices under the given node.

        :param node_id: ID of the starting node
        :return: list of leaf indices
        """
        leaf_index, children = self.nodes[node_id]
        if leaf_index >= 0:
            return [leaf_index]

        results = []
        for child in children.values():
            results.extend(self._get_leaves_below(child))
        return results

    def find_pattern(self, pattern):
        """
        Searches for a pattern in the suffix tree.

        :param pattern: string pattern to search for
        :return: list of starting positions in the text where the pattern occurs, or None
        """
        node = 0
        for char in pattern:
            if char in self.nodes[node][1]:
                node = self.nodes[node][1][char]
            else:
                return None
        return self._get_leaves_below(node)

    def print_tree(self):
        """
        Prints the structure of the suffix tree.
        """
        for node_id, (leaf_index, children) in self.nodes.items():
            if leaf_index >= 0:
                print(f"{node_id} : leaf index {leaf_index}")
            else:
                print(f"{node_id} -> children: {children}")