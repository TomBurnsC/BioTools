from ORF_finder import ORF_finder


def main():

    sequence = "ATGAAATAGCCCATGGTAGGGATGCCCTAA"

    orf_finder = ORF_finder(sequence)
    orf_finder.run_orf_scan()


if __name__ == "__main__":
    main()
