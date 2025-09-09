"""Tests for the Locus class."""

import pytest
from ...src import Locus


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


