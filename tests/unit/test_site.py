"""Tests for the Site class."""

import pytest
from unittest.mock import Mock

from GenomeUtils import Site, Locus, GenomeElement


class ConcreteSite(Site):
    """Concrete implementation of Site for testing purposes."""
    
    def get_site_type(self) -> str:
        """Return the type of site."""
        return "test_site"


class TestSite:
    """Test cases for the Site class."""

    def test_site_initialization_with_all_params(self):
        """Test Site initialization with all parameters."""
        genome_mock = Mock()
        parent_mock = Mock(spec=GenomeElement)
        
        site = ConcreteSite(
            chr="chr1",
            start=1000,
            end=2000,
            strand="+",
            sequence="ATCGATCGATCG",
            id="test_site",
            parent=parent_mock,
            genome=genome_mock,
            custom_attr="test_value"
        )
        
        assert site.id == "test_site"
        assert site.locus.chr == "chr1"
        assert site.locus.start == 1000
        assert site.locus.end == 2000
        assert site.locus.strand == "+"
        assert site.sequence == "ATCGATCGATCG"
        assert site._parent == parent_mock
        assert site._genome == genome_mock
        assert site.custom_attr == "test_value"

    def test_site_initialization_with_minimal_params(self):
        """Test Site initialization with minimal required parameters."""
        site = ConcreteSite(
            chr="chr2",
            start=500,
            end=600,
            strand="-",
            sequence="GGGGCCCCAAAA"
        )
        
        assert site.id == "chr2:500-600,-"  # Auto-generated from locus
        assert site.locus.chr == "chr2"
        assert site.locus.start == 500
        assert site.locus.end == 600
        assert site.locus.strand == "-"
        assert site.sequence == "GGGGCCCCAAAA"
        assert site._parent is None
        assert site._genome is None

    def test_site_initialization_auto_id_generation(self):
        """Test that Site auto-generates ID from locus when not provided."""
        site = ConcreteSite(
            chr="chrX",
            start=12345,
            end=12445,
            strand="+",
            sequence="ATCG"
        )
        
        expected_id = "chrX:12345-12445,+"
        assert site.id == expected_id

    def test_site_initialization_explicit_id_overrides_auto_generation(self):
        """Test that explicit ID overrides auto-generation."""
        site = ConcreteSite(
            chr="chrY",
            start=1,
            end=100,
            strand="-",
            sequence="TTTT",
            id="custom_site_id"
        )
        
        assert site.id == "custom_site_id"
        # But locus should still be correct
        assert site.locus.chr == "chrY"
        assert site.locus.start == 1
        assert site.locus.end == 100
        assert site.locus.strand == "-"

    def test_site_locus_creation(self):
        """Test that Site correctly creates Locus object."""
        site = ConcreteSite(
            chr="chr3",
            start=2000,
            end=3000,
            strand="+",
            sequence="AAAA"
        )
        
        assert isinstance(site.locus, Locus)
        assert site.locus.chr == "chr3"
        assert site.locus.start == 2000
        assert site.locus.end == 3000
        assert site.locus.strand == "+"

    def test_site_invalid_coordinates_raises_error(self):
        """Test that invalid coordinates raise ValueError."""
        with pytest.raises(ValueError, match="Start coordinate cannot be greater than end coordinate"):
            ConcreteSite(
                chr="chr1",
                start=2000,
                end=1000,  # Invalid: start > end
                strand="+",
                sequence="ATCG"
            )

    def test_site_invalid_start_coordinate_raises_error(self):
        """Test that invalid start coordinate raises ValueError."""
        with pytest.raises(ValueError, match="Start coordinate cannot be less than 1"):
            ConcreteSite(
                chr="chr1",
                start=0,  # Invalid: start < 1
                end=1000,
                strand="+",
                sequence="ATCG"
            )

    def test_site_repr_method(self):
        """Test Site __repr__ method."""
        site = ConcreteSite(
            chr="chr1",
            start=1000,
            end=2000,
            strand="+",
            sequence="ATCGATCGATCG",
            id="test_site"
        )
        
        expected = "ConcreteSite(id='test_site', locus=Locus(chr1:1000-2000, strand=+), sequence='ATCGATCGATCG')"
        assert repr(site) == expected

    def test_site_repr_method_with_auto_generated_id(self):
        """Test Site __repr__ method with auto-generated ID."""
        site = ConcreteSite(
            chr="chr2",
            start=500,
            end=600,
            strand="-",
            sequence="GGGG"
        )
        
        expected = "ConcreteSite(id='chr2:500-600,-', locus=Locus(chr2:500-600, strand=-), sequence='GGGG')"
        assert repr(site) == expected

    def test_site_repr_method_with_long_sequence(self):
        """Test Site __repr__ method with long sequence."""
        long_sequence = "ATCGATCGATCG" * 10  # 120 characters
        site = ConcreteSite(
            chr="chr1",
            start=1000,
            end=2000,
            strand="+",
            sequence=long_sequence,
            id="long_site"
        )
        
        expected = f"ConcreteSite(id='long_site', locus=Locus(chr1:1000-2000, strand=+), sequence='{long_sequence}')"
        assert repr(site) == expected

    def test_site_inheritance_from_genome_element(self):
        """Test that Site properly inherits from GenomeElement."""
        site = ConcreteSite(
            chr="chr1",
            start=1000,
            end=2000,
            strand="+",
            sequence="ATCG"
        )
        
        assert isinstance(site, GenomeElement)
        assert hasattr(site, 'id')
        assert hasattr(site, 'locus')
        assert hasattr(site, '_parent')
        assert hasattr(site, '_genome')

    def test_site_kwargs_passed_to_parent(self):
        """Test that additional kwargs are passed to parent class."""
        site = ConcreteSite(
            chr="chr1",
            start=1000,
            end=2000,
            strand="+",
            sequence="ATCG",
            custom_attribute="test_value",
            another_attr=42
        )
        
        assert site.custom_attribute == "test_value"
        assert site.another_attr == 42

    def test_site_empty_sequence(self):
        """Test Site with empty sequence."""
        site = ConcreteSite(
            chr="chr1",
            start=1000,
            end=1000,  # Single position
            strand="+",
            sequence=""
        )
        
        assert site.sequence == ""
        assert site.locus.start == 1000
        assert site.locus.end == 1000

    def test_site_different_strand_values(self):
        """Test Site with different strand values."""
        # Positive strand
        site_plus = ConcreteSite(
            chr="chr1",
            start=1000,
            end=2000,
            strand="+",
            sequence="ATCG"
        )
        assert site_plus.locus.strand == "+"
        
        # Negative strand
        site_minus = ConcreteSite(
            chr="chr1",
            start=1000,
            end=2000,
            strand="-",
            sequence="CGAT"
        )
        assert site_minus.locus.strand == "-"

    def test_site_with_different_chromosome_formats(self):
        """Test Site with different chromosome naming formats."""
        # Numeric chromosome
        site1 = ConcreteSite(chr="1", start=1000, end=2000, strand="+", sequence="ATCG")
        assert site1.locus.chr == "1"
        
        # Chr prefix
        site2 = ConcreteSite(chr="chr1", start=1000, end=2000, strand="+", sequence="ATCG")
        assert site2.locus.chr == "chr1"
        
        # Sex chromosomes
        site3 = ConcreteSite(chr="X", start=1000, end=2000, strand="+", sequence="ATCG")
        assert site3.locus.chr == "X"
        
        site4 = ConcreteSite(chr="chrY", start=1000, end=2000, strand="+", sequence="ATCG")
        assert site4.locus.chr == "chrY"
        
        # Mitochondrial
        site5 = ConcreteSite(chr="MT", start=1000, end=2000, strand="+", sequence="ATCG")
        assert site5.locus.chr == "MT"


    def test_site_sequence_attribute_access(self):
        """Test that sequence attribute can be accessed and modified."""
        site = ConcreteSite(
            chr="chr1",
            start=1000,
            end=2000,
            strand="+",
            sequence="ATCGATCGATCG"
        )
        
        # Access sequence
        assert site.sequence == "ATCGATCGATCG"
        
        # Modify sequence
        with pytest.raises(AttributeError, match="object has no setter"):
            site.sequence = "GGGGCCCCAAAA"

    def test_site_with_none_values_for_optional_params(self):
        """Test Site initialization with explicit None values for optional parameters."""
        site = ConcreteSite(
            chr="chr1",
            start=1000,
            end=2000,
            strand="+",
            sequence="ATCG",
            id=None,  # Should trigger auto-generation
            parent=None,
            genome=None
        )
        
        assert site.id == "chr1:1000-2000,+"  # Auto-generated
        assert site._parent is None
        assert site._genome is None

    def test_site_abstract_method_implementation(self):
        """Test that concrete implementations must implement abstract methods."""
        site = ConcreteSite(
            chr="chr1",
            start=1000,
            end=2000,
            strand="+",
            sequence="ATCG"
        )
        
        assert site.get_site_type() == "test_site"
