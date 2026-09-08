from Smith_Waterman import NeedleWunsch

def test_identical_sequences():
    sw = NeedleWunsch("AAAA", "AAAA")
    sw.run()

    assert sw.max_score == 4, "Smith-Waterman score should be 4"

test_identical_sequences()

def test_no_similarity():
    sw = NeedleWunsch("AAAA", "TTTT")
    sw.fill_matrix()
    sw.find_maximum()

    assert sw.max_score == 0, "Smith-Waterman score should be 0"

test_no_similarity()

def test_local_alignment():
    sw = NeedleWunsch("XXABCYY", "ZZABCWW")
    sw.fill_matrix()
    sw.find_maximum()

    result = sw.traceback(sw.max_positions[0])

    assert result["alignment1"] == "ABC", "Alignment1 should be ABC"
    assert result["alignment2"] == "ABC", "Alignment2 should be ABC"
    assert result["score"] == 3, "Scoring error"

#
# seq1 = ""
# seq2 = ""
#
# NW = NeedleWunsch(seq1, seq2)
#
# NW.run()
#
# SW = SmithWaterman(seq1, seq2)
#
# SW.run()
