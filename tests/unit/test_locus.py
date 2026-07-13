#!/usr/bin/env python
"""
Filename: tests/unit/test_locus.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.3
Description: Unit tests for the Locus class.
License: LGPL-3.0-or-later
"""

import pytest

from GenomeUtils.Genome import Locus


class TestLocus:
    """Test cases for the Locus class."""

    def test_locus_creation_valid(self):
        """Test valid locus creation."""
        locus = Locus("chr1", 100, 200, "+")
        assert locus.chr == "chr1"
        assert locus.start == 100
        assert locus.end == 200
        assert locus.strand == "+"

    def test_locus_creation_default_strand(self):
        """Test locus creation with default strand."""
        locus = Locus("chr2", 50, 150)
        assert locus.strand == "+"

    def test_locus_creation_negative_strand(self):
        """Test locus creation with negative strand."""
        locus = Locus("chr3", 300, 400, "-")
        assert locus.strand == "-"

    def test_locus_invalid_coordinates_start_greater_than_end(self):
        """Test that start > end raises ValueError."""
        with pytest.raises(ValueError, match="Start coordinate cannot be greater than end coordinate"):
            Locus("chr1", 200, 100)

    def test_locus_invalid_coordinates_start_less_than_one(self):
        """Test that start < 1 raises ValueError."""
        with pytest.raises(ValueError, match="Start coordinate cannot be less than 1"):
            Locus("chr1", 0, 100)

    def test_locus_equal_start_end(self):
        """Test that start == end is valid (single base)."""
        locus = Locus("chr1", 100, 100)
        assert len(locus) == 1

    def test_locus_length(self):
        """Test locus length calculation."""
        locus = Locus("chr1", 100, 200)
        assert len(locus) == 101  # 1-based inclusive

        locus_single = Locus("chr1", 50, 50)
        assert len(locus_single) == 1

    def test_locus_repr(self):
        """Test string representation of locus."""
        locus = Locus("chr1", 100, 200, "+")
        expected = "Locus(chr1:100-200, strand=+)"
        assert repr(locus) == expected

        locus_neg = Locus("chr2", 300, 400, "-")
        expected_neg = "Locus(chr2:300-400, strand=-)"
        assert repr(locus_neg) == expected_neg

    def test_locus_overlaps_same_chromosome(self):
        """Test overlap detection on same chromosome."""
        locus1 = Locus("chr1", 100, 200)
        locus2 = Locus("chr1", 150, 250)
        locus3 = Locus("chr1", 201, 300)
        locus4 = Locus("chr1", 50, 99)

        assert locus1.overlaps(locus2)
        assert locus2.overlaps(locus1)
        assert not locus1.overlaps(locus3)
        assert not locus1.overlaps(locus4)

    def test_locus_overlaps_different_chromosome(self):
        """Test that loci on different chromosomes don't overlap."""
        locus1 = Locus("chr1", 100, 200)
        locus2 = Locus("chr2", 150, 250)

        assert not locus1.overlaps(locus2)
        assert not locus2.overlaps(locus1)

    def test_locus_overlaps_edge_cases(self):
        """Test overlap edge cases."""
        locus1 = Locus("chr1", 100, 200)
        locus_touching_start = Locus("chr1", 200, 300)  # Touching at end/start
        locus_touching_end = Locus("chr1", 50, 100)    # Touching at start/end

        assert locus1.overlaps(locus_touching_start)
        assert locus1.overlaps(locus_touching_end)

    def test_locus_contains_same_chromosome(self):
        """Test containment on same chromosome."""
        locus_large = Locus("chr1", 100, 300)
        locus_small = Locus("chr1", 150, 250)
        locus_outside = Locus("chr1", 50, 99)
        locus_partial = Locus("chr1", 150, 350)

        assert locus_large.contains(locus_small)
        assert not locus_small.contains(locus_large)
        assert not locus_large.contains(locus_outside)
        assert not locus_large.contains(locus_partial)

    def test_locus_contains_different_chromosome(self):
        """Test that loci on different chromosomes don't contain each other."""
        locus1 = Locus("chr1", 100, 300)
        locus2 = Locus("chr2", 150, 250)

        assert not locus1.contains(locus2)
        assert not locus2.contains(locus1)

    def test_locus_contains_exact_same(self):
        """Test that a locus contains itself."""
        locus = Locus("chr1", 100, 200)
        locus_same = Locus("chr1", 100, 200)

        assert locus.contains(locus_same)
        assert locus_same.contains(locus)

    def test_locus_contains_edge_cases(self):
        """Test containment edge cases."""
        locus = Locus("chr1", 100, 200)
        locus_start_boundary = Locus("chr1", 100, 150)
        locus_end_boundary = Locus("chr1", 150, 200)
        locus_full_boundary = Locus("chr1", 100, 200)

        assert locus.contains(locus_start_boundary)
        assert locus.contains(locus_end_boundary)
        assert locus.contains(locus_full_boundary)

    def test_locus_dataclass_frozen(self):
        """Test that Locus is immutable (frozen dataclass)."""
        locus = Locus("chr1", 100, 200)
        
        with pytest.raises(AttributeError):
            locus.chr = "chr2"
        
        with pytest.raises(AttributeError):
            locus.start = 150

    def test_locus_dataclass_ordering(self):
        """Test that Locus objects can be ordered."""
        locus1 = Locus("chr1", 100, 200)
        locus2 = Locus("chr1", 150, 250)
        locus3 = Locus("chr2", 100, 200)

        # Test ordering works
        loci = [locus2, locus1, locus3]
        sorted_loci = sorted(loci)
        
        # Should be sorted by chromosome first, then by coordinates
        assert sorted_loci[0] == locus1  # chr1:100-200
        assert sorted_loci[1] == locus2  # chr1:150-250  
        assert sorted_loci[2] == locus3  # chr2:100-200

    def test_locus_equality(self):
        """Test locus equality comparison."""
        locus1 = Locus("chr1", 100, 200, "+")
        locus2 = Locus("chr1", 100, 200, "+")
        locus3 = Locus("chr1", 100, 200, "-")
        locus4 = Locus("chr1", 100, 201, "+")

        assert locus1 == locus2
        assert locus1 != locus3  # Different strand
        assert locus1 != locus4  # Different end coordinate

    def test_locus_hash(self):
        """Test that Locus objects are hashable (can be used in sets/dicts)."""
        locus1 = Locus("chr1", 100, 200)
        locus2 = Locus("chr1", 100, 200)
        locus3 = Locus("chr1", 150, 250)

        # Should be able to create sets
        locus_set = {locus1, locus2, locus3}
        assert len(locus_set) == 2  # locus1 and locus2 are equal

        # Should be able to use as dict keys
        locus_dict = {locus1: "value1", locus3: "value2"}
        assert len(locus_dict) == 2

    def test_locus_from_string_valid(self):
        """Test from_string method with a valid locus string."""
        locus_str = "12:25205246-25250936,-"
        locus = Locus.from_string(locus_str)

        assert locus.chr == "12"
        assert locus.start == 25205246
        assert locus.end == 25250936
        assert locus.strand == "-"

    def test_locus_from_string_invalid_format(self):
        """Test from_string method with invalid string formats."""
        invalid_strings = [
            "12:25205246",  # Missing end and strand
            "12:25205246-", # Missing end and strand
            "12:25205246-25250936", # Missing strand
            "12-25205246-25250936,-", # Missing colon
            "invalid_string",
            "12::,-",
            "12:25-,"
        ]

        for s in invalid_strings:
            with pytest.raises(ValueError, match="Invalid locus string"):
                Locus.from_string(s)

    def test_locus_from_string_positive_strand(self):
        """Test from_string method with positive strand."""
        locus_str = "chr1:1000-2000,+"
        locus = Locus.from_string(locus_str)

        assert locus.chr == "chr1"
        assert locus.start == 1000
        assert locus.end == 2000
        assert locus.strand == "+"

    def test_locus_from_string_various_chromosomes(self):
        """Test from_string method with various chromosome formats."""
        test_cases = [
            ("chr1:100-200,+", "chr1", 100, 200, "+"),
            ("1:100-200,+", "1", 100, 200, "+"),
            ("X:500-600,-", "X", 500, 600, "-"),
            ("chrX:500-600,-", "chrX", 500, 600, "-"),
            ("MT:1-100,+", "MT", 1, 100, "+"),
            ("chrM:1-100,-", "chrM", 1, 100, "-"),
        ]

        for locus_str, expected_chr, expected_start, expected_end, expected_strand in test_cases:
            locus = Locus.from_string(locus_str)
            assert locus.chr == expected_chr
            assert locus.start == expected_start
            assert locus.end == expected_end
            assert locus.strand == expected_strand

    def test_locus_from_string_single_base(self):
        """Test from_string method with single base locus."""
        locus_str = "chr1:1000-1000,+"
        locus = Locus.from_string(locus_str)

        assert locus.chr == "chr1"
        assert locus.start == 1000
        assert locus.end == 1000
        assert len(locus) == 1

    def test_locus_from_string_large_coordinates(self):
        """Test from_string method with large coordinates."""
        locus_str = "chr1:123456789-987654321,+"
        locus = Locus.from_string(locus_str)

        assert locus.chr == "chr1"
        assert locus.start == 123456789
        assert locus.end == 987654321

    def test_locus_from_string_roundtrip(self):
        """Test that from_string and __str__ are inverse operations."""
        test_loci = [
            Locus("chr1", 100, 200, "+"),
            Locus("12", 25205246, 25250936, "-"),
            Locus("X", 1000, 2000, "+"),
            Locus("MT", 1, 100, "-"),
        ]

        for original_locus in test_loci:
            locus_str = str(original_locus)
            reconstructed_locus = Locus.from_string(locus_str)
            assert reconstructed_locus == original_locus
            assert reconstructed_locus.chr == original_locus.chr
            assert reconstructed_locus.start == original_locus.start
            assert reconstructed_locus.end == original_locus.end
            assert reconstructed_locus.strand == original_locus.strand

    def test_locus_from_string_non_numeric_coordinates(self):
        """Test from_string method with non-numeric coordinates."""
        invalid_strings = [
            "chr1:abc-200,+",
            "chr1:100-xyz,+",
            "chr1:10.5-200,+",  # Float instead of int
            "chr1:100-200.5,+",
        ]

        for s in invalid_strings:
            with pytest.raises(ValueError, match="Invalid locus string"):
                Locus.from_string(s)

    def test_locus_from_string_invalid_start_end_order(self):
        """Test from_string method when start > end (should fail during Locus creation)."""
        locus_str = "chr1:200-100,+"
        
        # from_string will parse successfully, but Locus validation should catch it
        with pytest.raises(ValueError, match="Start coordinate cannot be greater than end coordinate"):
            Locus.from_string(locus_str)

    def test_locus_from_string_zero_start(self):
        """Test from_string method with zero start coordinate (should fail)."""
        locus_str = "chr1:0-100,+"
        
        # from_string will parse successfully, but Locus validation should catch it
        with pytest.raises(ValueError, match="Start coordinate cannot be less than 1"):
            Locus.from_string(locus_str)

    def test_locus_from_string_whitespace_in_coordinates(self):
        """Test from_string method with strings containing whitespace in coordinates."""
        # Only strings with whitespace in numeric coordinates will fail isdigit() check
        invalid_strings = [
            "chr1: 100-200,+",      # Space after colon, before start
            "chr1:100 -200,+",      # Space after start
            "chr1:100- 200,+",      # Space after dash
            "chr1:100-200 ,+",      # Space after end, before comma
        ]

        for s in invalid_strings:
            with pytest.raises(ValueError, match="Invalid locus string"):
                Locus.from_string(s)

