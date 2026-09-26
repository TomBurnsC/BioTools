"""
Progressive multiple sequence alignment.

Pairwise sequences are first scored and ranked by similarity (guide tree),
then combined in order of decreasing similarity using Needleman-Wunsch,
either sequence-to-sequence, profile-to-sequence, or profile-to-profile,
until all sequences are merged into alignment groups.
"""

from .alignment_group import AlignmentGroup
from .multiple_align import MultipleAlign
from .profile_alignment import NWProfileAndSeq, NWProfiles

__all__ = [
    "AlignmentGroup",
    "MultipleAlign",
    "NWProfileAndSeq",
    "NWProfiles",
]
