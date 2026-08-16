#!/usr/bin/env python
"""Unit tests for transcript-relative loci."""

import pytest

from GenomeUtils.Genome import TranscriptLocus


class TestTranscriptLocus:
    @pytest.mark.parametrize(
        ("transcript_id", "start", "end", "message"),
        [
            ("", 0, 1, "cannot be empty"),
            ("TRANS001", -1, 1, "less than zero"),
            ("TRANS001", 5, 5, "greater than start"),
            ("TRANS001", 5, 4, "greater than start"),
        ],
    )
    def test_coordinate_validation(self, transcript_id, start, end, message):
        with pytest.raises(ValueError, match=message):
            TranscriptLocus(transcript_id, start, end)

    def test_length(self):
        locus = TranscriptLocus("TRANS001", 3, 11)

        assert len(locus) == 8

    def test_equality_and_hashing(self):
        locus = TranscriptLocus("TRANS001", 3, 11)
        equivalent = TranscriptLocus("TRANS001", 3, 11)
        different = TranscriptLocus("TRANS001", 4, 11)

        assert locus == equivalent
        assert locus != different
        assert hash(locus) == hash(equivalent)
        assert len({locus, equivalent, different}) == 2
