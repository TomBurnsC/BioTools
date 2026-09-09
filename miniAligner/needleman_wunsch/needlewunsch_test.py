from Needleman_Wunsch import NeedleWunsch

def test_identical_sequences():
    nw = NeedleWunsch("TTAAGGGG", "TTAAGGGG")
    nw.run()

    assert nw.max_score == 8, "Needleman-Wunsch score should be 8"

test_identical_sequences()

def test_negative_alignment():
    nw = NeedleWunsch("CCCCCC", "GGGGGG")
    nw.run()

    assert nw.max_score == -6, "Smith-Waterman score should be -6"

test_negative_alignment()

seq1 = ""
seq2 = ""

sw = NeedleWunsch(seq1, seq2)
sw.run()
