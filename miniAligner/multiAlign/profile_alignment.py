"""NW-based alignment of a sequence or profile against an existing profile.

NWProfileAndSeq aligns a raw sequence against an AlignmentGroup; NWProfiles
aligns two AlignmentGroups against each other. Both fill a DP matrix and
trace back through it the same way NeedleWunsch does for two raw sequences,
just scored from a log-odds profile instead of a simple substitution matrix.
"""

from .constants import ALPHABET
from .matrix import Matrix


class NWProfileAndSeq:
    """Aligns a raw sequence against an existing profile (AlignmentGroup)."""

    def __init__(self, alignment_group, seq, seq_ind):
        self.alignment_group = alignment_group
        self.lo_mat = alignment_group.log_odds_matrix
        self.seq = seq
        self.seq_ind = seq_ind

        self.alphabet = ALPHABET

        self.alignment = []

        self.max_score = 0
        self.Matrix = Matrix(self.alignment_group.sequences[0], self.seq, "NW")

    def fill_matrix(self):
        """Fill the DP matrix using the profile's log-odds scores."""
        for index, row in enumerate(self.Matrix.matrix):
            if index > 0:
                for i in range(1, len(row)):
                    upper_diag = self.Matrix.matrix[index - 1][i - 1]
                    upper = self.Matrix.matrix[index - 1][i]
                    left = self.Matrix.matrix[index][i - 1]

                    diag = round(upper_diag + self.lo_mat[self.seq[i - 1]][index - 1], 3)
                    down = round(upper - 2, 3)
                    across = round(left - 2, 3)

                    self.Matrix.matrix[index][i] += max(diag, down, across)

    def traceback(self):
        """Trace back the DP matrix and fold the new sequence into the group."""
        i = len(self.lo_mat["A"])
        j = len(self.seq)

        position = self.Matrix.matrix[i][j]
        self.max_score = position

        while i != 0 or j != 0:
            upper_diag = self.Matrix.matrix[i - 1][j - 1]
            upper = self.Matrix.matrix[i - 1][j]
            left = self.Matrix.matrix[i][j - 1]

            # Diagonal: match or mismatch
            if position == round(upper_diag + self.lo_mat[self.seq[j - 1]][i - 1], 3):
                self.alignment.append(self.seq[j - 1])

                position = upper_diag
                i -= 1
                j -= 1

            # Up: gap in sequence 2
            elif position == round(upper - 2, 3):
                self.alignment.append("-")

                position = upper
                i -= 1

            # Left: gap in sequence 1
            elif position == round(left - 2, 3):
                self.alignment_group.insert_gap(i)
                self.alignment.append(self.seq[j - 1])

                position = left
                j -= 1

            else:
                print("TRACEBACK STUCK")
                print("Position:", position)
                print(f"{position} = {upper_diag} + {self.lo_mat[self.seq[j - 1]][i - 1]}")
                print("i:", i, "j:", j)
                break

        # If sequence 1 still has characters remaining
        while i > 0:
            self.alignment.append("-")
            i -= 1

        # If sequence 2 still has characters remaining
        while j > 0:
            self.alignment_group.insert_gap()
            j -= 1

        self.alignment.reverse()
        self.alignment_group.seq_names.append(self.seq_ind)
        self.alignment_group.sequences.append(self.alignment)

        self.alignment_group.update_matrices()


class NWProfiles:
    """Aligns two existing profiles (AlignmentGroups) against each other."""

    def __init__(self, alignment_group1, alignment_group2):
        self.alignment_group1 = alignment_group1
        self.alignment_group2 = alignment_group2
        self.lo_mat1 = alignment_group1.log_odds_matrix
        self.lo_mat2 = alignment_group2.log_odds_matrix

        self.alphabet = ALPHABET

        self.max_score = 0
        self.Matrix = Matrix(self.alignment_group1.sequences[0], self.alignment_group2.sequences[0], "NW")

        self.col_scores = []

    def col_scorer(self, col1, col2):
        """Sum-of-pairs style score between a column from each profile."""
        column_score = 0

        for char_a in self.alphabet:
            for char_b in self.alphabet:
                score = self.lo_mat1[char_a][col1 - 1] * self.lo_mat2[char_b][col2 - 1]

                if char_a != char_b:
                    score *= -1

                column_score += score

        return column_score

    def fill_matrix(self):
        """Fill the DP matrix using profile-profile column scores."""
        for index, row in enumerate(self.Matrix.matrix):
            if index > 0:
                for i in range(1, len(row)):
                    upper_diag = self.Matrix.matrix[index - 1][i - 1]
                    upper = self.Matrix.matrix[index - 1][i]
                    left = self.Matrix.matrix[index][i - 1]

                    col_score = self.col_scorer(col1=index, col2=i)

                    diag = round(upper_diag + col_score, 3)
                    down = round(upper - 2, 3)
                    across = round(left - 2, 3)

                    self.Matrix.matrix[index][i] += max(diag, down, across)

    def traceback(self):
        """Trace back the DP matrix and merge the two profiles into one group."""
        i = len(self.lo_mat1["A"])
        j = len(self.lo_mat2["A"])

        position = self.Matrix.matrix[i][j]
        self.max_score = position

        while i != 0 or j != 0:
            upper_diag = self.Matrix.matrix[i - 1][j - 1]
            upper = self.Matrix.matrix[i - 1][j]
            left = self.Matrix.matrix[i][j - 1]

            # Diagonal: match or mismatch
            if position == round(upper_diag + self.col_scorer(i, j), 3):
                position = upper_diag
                i -= 1
                j -= 1

            # Up: gap in sequence 2
            elif position == round(upper - 2, 3):
                self.alignment_group2.insert_gap(j)

                position = upper
                i -= 1

            # Left: gap in sequence 1
            elif position == round(left - 2, 3):
                self.alignment_group1.insert_gap(i)

                position = left
                j -= 1

            else:
                print("TRACEBACK STUCK")
                print("Position:", position)
                print("i:", i, "j:", j)
                break

        # If sequence 1 still has characters remaining
        while i > 0:
            self.alignment_group2.insert_gap(j)
            i -= 1

        # If sequence 2 still has characters remaining
        while j > 0:
            self.alignment_group1.insert_gap(i)
            j -= 1

        for num in range(len(self.alignment_group2.seq_names)):
            self.alignment_group1.seq_names.append(self.alignment_group2.seq_names[num])
            self.alignment_group1.sequences.append(self.alignment_group2.sequences[num])

        self.alignment_group1.update_matrices()
        self.alignment_group2.merged = True

        return self.alignment_group1
