"""Run progressive multiple sequence alignment on an example sequence set."""

from progressive_msa import MultipleAlign, NWProfileAndSeq, NWProfiles

if __name__ == "__main__":

    # Input sequences
    sequences = [
        "ATGCGCTTGA",
        "ATGCTTGA",
        "CTCCCG",
        "CTGCTCTG",
        "TCTT",
    ]

    # Initialise alignment
    multi = MultipleAlign(sequences)

    # 1. Calculate pairwise sequence similarities
    multi.fill_matrix()
    multi.similarity_matrix.generate_matrix()

    # 2. Generate guide order
    #
    # Sequences are ranked from most to least similar. This determines the
    # order in which sequences/alignment groups are progressively combined.
    multi.guide_tree()

    print("Guide order:")
    for i, j in multi.alignment_orders:
        print(f"  {i} + {j}")

    # 3. Progressively combine sequences and alignment groups
    for order_num, (i, j) in enumerate(multi.alignment_orders):

        # First pair: create the initial alignment group
        if order_num == 0:
            multi.generate_matrices(multi.input_sequences[i], multi.input_sequences[j])
            print(f"Initial alignment: {i} + {j}")
            continue

        # Determine whether each sequence is already part of a profile
        profile1 = None
        profile2 = None
        merge_needed = True

        for alignment_group in multi.alignments:

            i_in_group = i in alignment_group.seq_names
            j_in_group = j in alignment_group.seq_names

            # Both sequences are already aligned together
            if i_in_group and j_in_group:
                print(f"{i} and {j} skipped")
                merge_needed = False
                break

            if i_in_group:
                profile1 = alignment_group

            if j_in_group:
                profile2 = alignment_group

        # Determine the type of alignment required
        if merge_needed:

            # Neither sequence belongs to an existing group
            if profile1 is None and profile2 is None:
                multi.generate_matrices(multi.input_sequences[i], multi.input_sequences[j])
                print(f"Sequence + sequence: {i} + {j}")

            # Sequence i belongs to a group
            elif profile1 is not None and profile2 is None:
                alignment = NWProfileAndSeq(profile1, multi.input_sequences[j], j)
                alignment.fill_matrix()
                alignment.traceback()
                print(f"Profile + sequence: {i} + {j}")

            # Sequence j belongs to a group
            elif profile1 is None and profile2 is not None:
                alignment = NWProfileAndSeq(profile2, multi.input_sequences[i], i)
                alignment.fill_matrix()
                alignment.traceback()
                print(f"Sequence + profile: {i} + {j}")

            # Both sequences belong to different groups
            elif profile1 is not None and profile2 is not None:
                alignment = NWProfiles(profile1, profile2)
                alignment.fill_matrix()
                alignment.traceback()
                print(f"Profile + profile: {i} + {j}")

                multi.update_alignment_groups()

    # 4. Display final alignment
    #
    # Groups are printed in merge order (most to least similar).
    print("\nFinal alignment:")
    print("-" * 40)

    for alignment_group in multi.alignments:
        print("Sequences:", alignment_group.seq_names)

        for name, sequence in zip(alignment_group.seq_names, alignment_group.sequences):
            print(f"{name}: {''.join(sequence)}")

        print()
