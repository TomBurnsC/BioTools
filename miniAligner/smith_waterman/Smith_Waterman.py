import matrix

class NeedleWunsch:
    """
    Performs local sequence alignment using the Smith-Waterman algorithm.

    Scoring:
        Match    = +1
        Mismatch = -1
        Gap      = -2
    """

    MATCH = 1
    MISMATCH = -1
    GAP = -2

    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2

        self.method = "SW"

        self.Matrix = matrix.Matrix(
            self.seq1,
            self.seq2,
            self.method
        )

        self.max_score = 0
        self.max_positions = []

    # ---------------------------------------------------------
    # Scoring
    # ---------------------------------------------------------

    def match_check(self, a, b):
        """
        Returns the score for aligning two bases.
        """
        if a == b:
            return self.MATCH
        return self.MISMATCH

    # ---------------------------------------------------------
    # Matrix filling
    # ---------------------------------------------------------

    def fill_matrix(self):
        """
        Fills the Smith-Waterman scoring matrix.

        Each cell is calculated as:

        H(i,j) = max(
            0,
            H(i-1,j-1) + match/mismatch,
            H(i-1,j) + gap,
            H(i,j-1) + gap
        )
        """

        for i in range(1, len(self.Matrix.matrix)):

            for j in range(1, len(self.Matrix.matrix[i])):

                diagonal = (
                    self.Matrix.matrix[i - 1][j - 1]
                    + self.match_check(
                        self.seq1[i - 1],
                        self.seq2[j - 1]
                    )
                )

                up = (
                    self.Matrix.matrix[i - 1][j]
                    + self.GAP
                )

                left = (
                    self.Matrix.matrix[i][j - 1]
                    + self.GAP
                )

                self.Matrix.matrix[i][j] = max(
                    0,
                    diagonal,
                    up,
                    left
                )

    # ---------------------------------------------------------
    # Find maximum score
    # ---------------------------------------------------------

    def find_maximum(self):
        """
        Finds the maximum score in the matrix.

        All positions containing the maximum score are stored,
        allowing multiple optimal local alignments to be found.
        """

        self.max_score = 0
        self.max_positions = []

        for i, row in enumerate(self.Matrix.matrix):

            for j, value in enumerate(row):

                # Ignore the first row/column
                if i == 0 or j == 0:
                    continue

                if value > self.max_score:

                    self.max_score = value
                    self.max_positions = [(i, j)]

                elif value == self.max_score:

                    self.max_positions.append((i, j))

    # ---------------------------------------------------------
    # Traceback
    # ---------------------------------------------------------

    def traceback(self, position):
        """
        Performs traceback from a specified maximum-scoring cell.

        Returns:

            alignment1
            alignment2
            start_position
            end_position
            score

        """

        i, j = position

        end_position = (i, j)

        score = self.Matrix.matrix[i][j]

        alignment1 = []
        alignment2 = []

        while self.Matrix.matrix[i][j] > 0:

            current = self.Matrix.matrix[i][j]

            diagonal = self.Matrix.matrix[i - 1][j - 1]
            up = self.Matrix.matrix[i - 1][j]
            left = self.Matrix.matrix[i][j - 1]

            # -------------------------------------------------
            # Diagonal move
            # -------------------------------------------------

            diagonal_score = (
                diagonal
                + self.match_check(
                    self.seq1[i - 1],
                    self.seq2[j - 1]
                )
            )

            if current == diagonal_score:

                alignment1.append(self.seq1[i - 1])
                alignment2.append(self.seq2[j - 1])

                i -= 1
                j -= 1

            # -------------------------------------------------
            # Up move
            # -------------------------------------------------

            elif current == up + self.GAP:

                alignment1.append(self.seq1[i - 1])
                alignment2.append("-")

                i -= 1

            # -------------------------------------------------
            # Left move
            # -------------------------------------------------

            elif current == left + self.GAP:

                alignment1.append("-")
                alignment2.append(self.seq2[j - 1])

                j -= 1

            else:
                raise RuntimeError(
                    "Traceback Error - No Valid Previous Values"
                )

        start_position = (i, j)

        alignment1 = "".join(reversed(alignment1))
        alignment2 = "".join(reversed(alignment2))

        return {
            "alignment1": alignment1,
            "alignment2": alignment2,
            "start": start_position,
            "end": end_position,
            "score": score
        }


    def format_alignment(self, result):
        """
        Prints an alignment in easily interpretable format.
        """

        alignment1 = result["alignment1"]
        alignment2 = result["alignment2"]

        match_line = ""

        for a, b in zip(alignment1, alignment2):

            if a == b:
                match_line += "|"

            elif a == "-" or b == "-":
                match_line += " "

            else:
                match_line += "."

        start1, start2 = result["start"]
        end1, end2 = result["end"]

        # Convert matrix coordinates into 1-based sequence
        # coordinates.

        sequence_start1 = start1 + 1
        sequence_end1 = end1

        sequence_start2 = start2 + 1
        sequence_end2 = end2

        print()
        print("Alignment")
        print("=" * 60)

        print(
            f"Seq1 : "
            f"{alignment1} [{sequence_start1}-{sequence_end1}]"
        )

        print(
            f"       "
            f"{match_line}"
        )

        print(
            f"Seq2 : "
            f"{alignment2} [{sequence_start2}-{sequence_end2}]"
        )

        print()
        print(f"Score: {result['score']}")

    # ---------------------------------------------------------
    # Run all optimal alignments
    # ---------------------------------------------------------

    def generate_alignments(self):
        """
        prints all optimal local alignments.
        """

        self.find_maximum()

        print("\nSmith-Waterman Local Alignment")
        print("=" * 60)

        print(f"Maximum score: {self.max_score}")
        print(f"Number of maximum-scoring positions: "
              f"{len(self.max_positions)}")

        for number, position in enumerate(
            self.max_positions,
            start=1
        ):

            print()
            print(f"Optimal alignment {number}")
            print("-" * 60)

            result = self.traceback(position)

            self.format_alignment(result)


    def run(self):
        """
        Runs the Smith-Waterman algorithm.
        """

        self.fill_matrix()

        self.generate_alignments()

        matrix_query = input("Generate matrix? y or n: ")

        if matrix_query == "y":
            self.Matrix.generate_matrix()

