import matrix

seq1 = "AACGATAGCTAGCTAG"
seq2 = "ACCGAAGCGAG"


class NeedleWunsch(object):

    '''
    Fills matrix and performs traceback before building optimal alignment
    '''

    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2

        self.alignment1 = []
        self.alignment2 = []

        self.max_score = 0
        self.Matrix = matrix.Matrix(self.seq1, self.seq2, "NW")

    def match_check(self, i, j):
        if i == j:
            return 1
        else:
            return -1

    def fill_matrix(self):

        for index, row in enumerate(self.Matrix.matrix):

            if index > 0:

                for i in range(1, len(row)):

                    upper_diag = self.Matrix.matrix[index - 1][i - 1]
                    upper = self.Matrix.matrix[index - 1][i]
                    left = self.Matrix.matrix[index][i - 1]

                    if self.seq1[index - 1] == self.seq2[i - 1]:
                        diag = upper_diag + 1
                    else:
                        diag = upper_diag - 1

                    down = upper - 2
                    across = left - 2

                    self.Matrix.matrix[index][i] += max(
                        diag,
                        down,
                        across
                    )

    def traceback(self):

        i = len(self.seq1)
        j = len(self.seq2)

        position = self.Matrix.matrix[i][j]
        self.max_score = position

        while i != 0 or j != 0:

            upper_diag = self.Matrix.matrix[i - 1][j - 1]
            upper = self.Matrix.matrix[i - 1][j]
            left = self.Matrix.matrix[i][j - 1]

            # Diagonal: match or mismatch
            if position == upper_diag + self.match_check(
                    self.seq1[i - 1],
                    self.seq2[j - 1]
            ):

                self.alignment1.append(self.seq1[i - 1])
                self.alignment2.append(self.seq2[j - 1])

                position = upper_diag
                i -= 1
                j -= 1

            # Up: gap in sequence 2
            elif position == upper - 2:

                self.alignment1.append(self.seq1[i - 1])
                self.alignment2.append("-")

                position = upper
                i -= 1

            # Left: gap in sequence 1
            elif position == left - 2:

                self.alignment1.append("-")
                self.alignment2.append(self.seq2[j - 1])

                position = left
                j -= 1

            else:
                print("TRACEBACK STUCK")
                print("Position:", position)
                print("i:", i, "j:", j)
                break

        # If sequence 2 still has characters remaining
        while i > 0:

            self.alignment1.append(self.seq1[i - 1])
            self.alignment2.append("-")

            i -= 1

        # If sequence 2 still has characters remaining
        while j > 0:

            self.alignment1.append("-")
            self.alignment2.append(self.seq2[j - 1])

            j -= 1

        self.alignment1.reverse()
        self.alignment2.reverse()

    def generate_alignment(self):
        """
        Prints the optimal global alignment.
        """

        print("\nNeedleman-Wunsch Global Alignment")
        print("=" * 60)

        print(
            f"Alignment score: "
            f"{self.Matrix.matrix[len(self.seq1)][len(self.seq2)]}"
        )

        print("\nOptimal Alignment")
        print("-" * 60)

        alignment1 = "".join(self.alignment1)
        alignment2 = "".join(self.alignment2)

        match_line = ""

        for base1, base2 in zip(alignment1, alignment2):

            if base1 == base2:
                match_line += "|"
            elif base1 == "-" or base2 == "-":
                match_line += " "
            else:
                match_line += "."

        print(alignment1)
        print(match_line)
        print(alignment2)


    def run(self):
        """
        Runs the Needleman-Wunsch algorithm - Included for convenience.
        """

        self.fill_matrix()
        self.traceback()
        self.generate_alignment()

        matrix_query = input("Generate matrix? y or n: ")

        if matrix_query == "y":
            self.Matrix.generate_matrix()




