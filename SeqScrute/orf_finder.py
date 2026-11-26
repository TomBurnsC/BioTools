from itertools import combinations
import math

start_codon = "ATG"

stop_codons = {
    "TAA",
    "TAG",
    "TGA"
}

class ORF_finder:

    def __init__(self, sequence):
        self.sequence = list(sequence)

        self.reverse_complement = self.sequence[::-1]

        self.frames = {}

        self.frame_orfs = {}

        self.frames = {

            "+1" : self.sequence,
            "+2" : self.sequence[1:],
            "+3" : self.sequence[2:],

            "-1": self.reverse_complement[0:],
            "-2": self.reverse_complement[1:],
            "-3": self.reverse_complement[2:]

        }

        self.orf_indexes = {
        "+1": [],
        "+2": [],
        "+3": [],

        "-1": [],
        "-2": [],
        "-3": []
        }

    def orf_logger(self, frame, seq, start_codon_index):

        self.orf_indexes[frame].append(start_codon_index)

        if frame not in self.frame_orfs:
            self.frame_orfs[frame] = [start_codon]
            orf_index = 0

        else:
            self.frame_orfs[frame].append(start_codon)

            orf_index = len(self.frame_orfs[frame])-1
            self.frame_orfs[frame][orf_index] = start_codon


        for codon_index in range(start_codon_index+3, len(seq), 3):

            if codon_index < len(seq) - 2:

                codon = "".join(seq[codon_index:codon_index + 3])

                if codon == start_codon:
                    self.orf_logger(frame, seq, codon_index)

                self.frame_orfs[frame][orf_index] += codon

                if codon in stop_codons:
                    break

    def start_codon_screener(self, frame, seq):

        for codon_start_ind in range(0, len(seq), 3):

            if codon_start_ind < len(seq) - 2:

                codon = seq[codon_start_ind] + seq[codon_start_ind + 1] + seq[codon_start_ind + 2]

                if codon == start_codon:

                    self.orf_logger(frame, seq, codon_start_ind)

    def format_output(self):
        width = 60

        sequence = "".join(self.sequence)

        print("\n" + "=" * width)
        print("ORF FINDER RESULTS".center(width))
        print("=" * width)

        print(f"\nOriginal Sequence ({len(sequence)} bp):")
        print(sequence)

        for frame, alignments in self.frame_orfs.items():

            print(f"\n{'-' * width}")
            print(f"Reading Frame: {frame}")
            print(f"ORFs identified: {len(alignments)}")
            print("-" * width)

            if not alignments:
                print("No ORFs found.")
                continue

            # Select the appropriate template
            offset = int(frame[1]) - 1

            if frame[0] == "-":
                template = sequence[::-1]
            else:
                template = sequence

            for i, orf_sequence in enumerate(alignments):

                start = self.orf_indexes[frame][i] + offset

                # Convert reverse-strand coordinates to original sequence coordinates
                if frame[0] == "-":
                    position = len(sequence) - start - len(orf_sequence) + 1
                else:
                    position = start + 1

                print(
                    f"\nORF {i + 1} | Position: {position} | "
                    f"Length: {len(orf_sequence)} bp"
                )

                print(f"Template: {template}")
                print(f"          {' ' * start}{orf_sequence}")

        print("\n" + "=" * width)

    def run_orf_scan(self):


        for frame, seq in self.frames.items():
            self.start_codon_screener(frame, seq)


        self.format_output()




