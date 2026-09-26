
from itertools import combinations

from .alignment_group import AlignmentGroup
from .constants import ALPHABET
from .matrix import Matrix
from .needleman_wunsch import NeedleWunsch


class MultipleAlign:
    """Drives progressive alignment of a set of input sequences."""

    def __init__(self, input_sequences):
        self.alignment_orders = None
        self.input_sequences = input_sequences
        self.num_combinations = (len(self.input_sequences) * (len(self.input_sequences) - 1)) / 2
        self.similarity_matrix = Matrix(input_sequences, input_sequences, "SM")
        self.alphabet = ALPHABET

        self.alignments = []

    def fill_matrix(self):
        """Score every pairwise combination of input sequences via NW."""
        for seq1ind, seq2ind in combinations(range(len(self.input_sequences)), 2):
            seq1 = self.input_sequences[seq1ind]
            seq2 = self.input_sequences[seq2ind]

            alignment = NeedleWunsch(seq1, seq2)
            alignment.run(False)

            self.similarity_matrix.matrix[seq2ind][seq1ind] = alignment.max_score

    def guide_tree(self):
        """Rank sequence pairs by similarity, most similar first.

        Sets self.alignment_orders to a list of (index, index) tuples that
        can be used to look up sequences in self.input_sequences.
        """
        self.alignment_orders = []

        while len(self.alignment_orders) < self.num_combinations:
            self.max_score = -1000
            val1 = 0
            val2 = 0

            for i, row in enumerate(self.similarity_matrix.matrix):
                for j, value in enumerate(row):
                    if value == "-":
                        continue

                    if value > self.max_score:
                        self.max_score = value
                        val1 = i
                        val2 = j

            self.alignment_orders.append((val2, val1))
            self.similarity_matrix.matrix[val1][val2] = "-"

    def generate_matrices(self, seq1, seq2):
        """Align two raw input sequences and start a new AlignmentGroup."""
        seq1name = self.input_sequences.index(seq1)
        seq2name = self.input_sequences.index(seq2)

        nw = NeedleWunsch(seq1, seq2)
        nw.run(False)
        s1 = nw.alignment1
        s2 = nw.alignment2

        ag = AlignmentGroup([seq1name, seq2name], [s1, s2])

        for char in self.alphabet:
            ag.count_matrix[char] = [0 for i in range(len(s1))]

        counter = 0
        for char1, char2 in zip(s1, s2):
            if char1 == char2:
                ag.count_matrix[char1][counter] += 2
            else:
                ag.count_matrix[char1][counter] += 1
                ag.count_matrix[char2][counter] += 1

            counter += 1

        ag.generate_profile()

        self.alignments.append(ag)

    def update_alignment_groups(self):
        """Drop alignment groups that have been merged into another group.

        for alignment in self.alignments:
            if alignment.merged:
                self.alignments.remove(alignment)
