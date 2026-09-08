class Matrix(object):

    '''
    Constructs matrix according to sequence lengths
    '''

    def __init__(self, seq1, seq2, method = "NW"):
        self.seq1 = seq1
        self.seq2 = seq2
        self.invalid = False

        # construct raw matrix
        row = len(seq1) + 1
        col = len(seq2) + 1
        self.matrix = [[0 for x in range(col)] for j in range(row)]

        # fill the boundary rows

        if method == "NW":
            for index, row in enumerate(self.matrix):
                if index == 0:
                    counter = 0
                    for j in range(len(self.matrix[index])):
                        self.matrix[index][j] += counter * -2
                        counter += 1
                else:
                    self.matrix[index][0] += index * -2

        elif method != "SW":
            print("!Invalid method!")
            self.invalid = True

    def generate_matrix(self):
        if not self.invalid:
            for i in self.matrix:
                print(i)
        else:
            print("!Invalid Method! Cannot generate matrix!")