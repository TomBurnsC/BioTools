"""AlignmentGroup: a profile of already-aligned sequences."""

import math

from .constants import ALPHABET


class AlignmentGroup:
    """A set of pre-aligned sequences tracked as count/log-odds profiles"""

    def __init__(self, seq_names, sequences):
        self.seq_names = seq_names
        self.alphabet = ALPHABET
        self.sequences = sequences
        self.count_matrix = {}
        self.log_odds_matrix = None
        self.merged = False

    def generate_profile(self):
        """Build log-odds profile from current count matrix."""

        self.log_odds_matrix = {
            char: values.copy()
            for char, values in self.count_matrix.items()
        }

        for col in range(len(self.count_matrix["A"])):
            col_total = 0

            for char in self.alphabet:
                col_total += self.count_matrix[char][col]

            for char in self.alphabet:
                if self.count_matrix[char][col] == 0:
                    self.log_odds_matrix[char][col] = 0
                    continue

                p = self.count_matrix[char][col] / col_total
                self.log_odds_matrix[char][col] = round(math.log(p / 0.25), 3)

    def update_matrices(self):
        """Recompute the count matrix from self.sequences, then the profile."""
        count_mat = {
            char: [0] * len(self.sequences[0])
            for char in ALPHABET
        }

        for sequence in self.sequences:

            for col, char in enumerate(sequence):

                count_mat[char][col] += 1

        self.count_matrix = count_mat

        self.generate_profile()

    def insert_gap(self, index=None):
        """Insert a gap column into every sequence in the group.

        With no index, appends a gap to the end of each sequence. With an
        index, inserts at that position — reversed, since NW traceback runs
        backwards.
        """
        if index is None:
            for seq in self.sequences:
                seq.append("-")

        else:
            rev_ind = len(self.sequences[0]) - (index + 1)

            for seq in self.sequences:
                seq.insert(rev_ind, "-")
