"""Tests for the Exon class."""

import pytest
from unittest.mock import Mock
from Bio.Seq import Seq

from GenomeUtils.Genome import Exon, Locus


class TestExon:
    """Test cases for the Exon class."""

    def create_mock_transcript(self, chr_id="chr1", start=1000, end=5000, strand="+"):
        """Helper to create a mock transcript."""
        transcript = Mock()
        transcript.chr = chr_id
        transcript.start = start
        transcript.end = end
        transcript.strand = strand
        transcript.sequence = Seq("ATCGATCGATCGAAATTTGGGCCCTTTTAAAA")
        return transcript

    def test_exon_creation(self):
        """Test basic exon creation."""
        transcript = self.create_mock_transcript()
        genome_mock = Mock()
        
        exon = Exon(
            id="ENSE00000001",
            chr="chr1",
            start=1100,
            end=1200,
            strand="+",
            transcript=transcript,
            genome=genome_mock
        )
        
        assert exon.id == "ENSE00000001"
        assert exon.start == 1100
        assert exon.end == 1200
        assert exon.strand == "+"
        assert exon.chr == "chr1"
        assert exon._parent == transcript
        assert exon._genome == genome_mock
        assert len(exon) == 101  # 1-based inclusive

    def test_exon_creation_with_kwargs(self):
        """Test exon creation with additional attributes."""
        transcript = self.create_mock_transcript("chr2", 2000, 6000, "-")
        genome_mock = Mock()
        
        exon = Exon(
            id="ENSE00000002",
            chr="chr2",
            start=3000,
            end=3100,
            strand="-",
            transcript=transcript,
            genome=genome_mock,
            exon_number=1,
            phase=0
        )
        
        assert exon.exon_number == 1
        assert exon.phase == 0
        assert exon.strand == "-"

    def test_exon_locus_creation(self):
        """Test that exon creates correct locus."""
        transcript = self.create_mock_transcript("chr3", 1000, 5000, "+")
        genome_mock = Mock()
        
        exon = Exon(
            id="ENSE00000003",
            chr="chr3",
            start=2000,
            end=2500,
            strand="+",
            transcript=transcript,
            genome=genome_mock
        )
        
        expected_locus = Locus("chr3", 2000, 2500, "+")
        assert exon.locus == expected_locus

    def test_exon_get_transcript(self):
        """Test get_transcript method."""
        transcript = self.create_mock_transcript()
        genome_mock = Mock()
        
        exon = Exon(
            id="ENSE00000001",
            chr="chr1",
            start=1100,
            end=1200,
            strand="+",
            transcript=transcript,
            genome=genome_mock
        )
        
        assert exon.get_transcript() == transcript
        assert exon.get_transcript() is exon._parent

    def test_exon_sequence_property(self):
        """Test exon sequence property."""
        transcript = self.create_mock_transcript()
        genome_mock = Mock()
        
        # Create multiple exons to test proper indexing
        exon1 = Exon(
            id="ENSE00000001",
            chr="chr1",
            start=1100,
            end=1109,  # 10 bp
            strand="+",
            transcript=transcript,
            genome=genome_mock
        )
        
        exon2 = Exon(
            id="ENSE00000002",
            chr="chr1",
            start=1200,
            end=1219,  # 20 bp
            strand="+",
            transcript=transcript,
            genome=genome_mock
        )
        
        # Mock the transcript's exons list
        transcript.exons = [exon1, exon2]
        
        # Test first exon - should get first 10 bp of transcript sequence
        sequence1 = exon1.sequence
        expected1 = transcript.sequence[0:10]
        assert str(sequence1) == str(expected1)
        
        # Test second exon - should get next 20 bp of transcript sequence
        sequence2 = exon2.sequence
        expected2 = transcript.sequence[10:30]
        assert str(sequence2) == str(expected2)

    def test_exon_sequence_property_not_in_transcript(self):
        """Test exon sequence when exon is not in transcript's exon list."""
        transcript = self.create_mock_transcript()
        genome_mock = Mock()
        
        exon = Exon(
            id="ENSE00000001",
            chr="chr1",
            start=1100,
            end=1200,
            strand="+",
            transcript=transcript,
            genome=genome_mock
        )
        
        # Mock transcript with empty exons list
        transcript.exons = []
        
        # Should return empty string when exon not found
        sequence = exon.sequence
        assert str(sequence) == ""

    def test_exon_sequence_property_complex_scenario(self):
        """Test exon sequence with multiple exons in realistic scenario."""
        transcript = self.create_mock_transcript()
        genome_mock = Mock()
        
        # Create transcript with longer sequence for testing
        transcript.sequence = Seq("ATCGATCGATCGAAATTTGGGCCCTTTTAAAACCCCGGGGTTTTAAAACCCCGGGG")  # 56 bp
        
        # Create three exons
        exon1 = Exon(
            id="ENSE00000001",
            chr="chr1",
            start=1100,
            end=1119,  # 20 bp
            strand="+",
            transcript=transcript,
            genome=genome_mock
        )
        
        exon2 = Exon(
            id="ENSE00000002", 
            chr="chr1",
            start=1200,
            end=1214,  # 15 bp
            strand="+",
            transcript=transcript,
            genome=genome_mock
        )
        
        exon3 = Exon(
            id="ENSE00000003",
            chr="chr1",
            start=1300,
            end=1320,  # 21 bp
            strand="+",
            transcript=transcript,
            genome=genome_mock
        )
        
        # Mock the transcript's exons list
        transcript.exons = [exon1, exon2, exon3]
        
        # Test middle exon
        sequence2 = exon2.sequence
        expected2 = transcript.sequence[20:35]  # positions 20-34 (15 bp)
        assert str(sequence2) == str(expected2)

    def test_exon_inheritance_from_genome_element(self):
        """Test that Exon properly inherits from GenomeElement."""
        transcript = self.create_mock_transcript("chr5", 10000, 20000, "-")
        genome_mock = Mock()
        
        exon = Exon(
            id="ENSE00000005",
            chr="chr5",
            start=12000,
            end=12100,
            strand="-",
            transcript=transcript,
            genome=genome_mock
        )
        
        # Should have GenomeElement properties
        assert exon.id == "ENSE00000005"
        assert exon.chr == "chr5"
        assert exon.start == 12000
        assert exon.end == 12100
        assert exon.strand == "-"
        assert len(exon) == 101
        
        # Should have proper repr
        expected_repr = "Exon(id='ENSE00000005', locus=Locus(chr5:12000-12100, strand=-))"
        assert repr(exon) == expected_repr

    def test_exon_equality(self):
        """Test exon equality based on ID and locus."""
        transcript1 = self.create_mock_transcript("chr1")
        transcript2 = self.create_mock_transcript("chr2")
        genome_mock = Mock()
        
        exon1 = Exon(
            id="SAME_ID",
            chr="chr1",
            start=1100,
            end=1200,
            strand="+",
            transcript=transcript1,
            genome=genome_mock
        )
        
        exon2 = Exon(
            id="SAME_ID",
            chr="chr1",
            start=1100,      # Same coordinates
            end=1200,
            strand="+",      # Same strand
            transcript=transcript1,  # Same transcript (same chr)
            genome=genome_mock
        )
        
        exon3 = Exon(
            id="SAME_ID",
            chr="chr2",
            start=2000,      # Different coordinates
            end=2100,
            strand="-",      # Different strand
            transcript=transcript2,  # Different transcript (different chr)
            genome=genome_mock
        )
        
        exon4 = Exon(
            id="DIFFERENT_ID",
            chr="chr1",
            start=1100,
            end=1200,
            strand="+",
            transcript=transcript1,
            genome=genome_mock
        )
        
        assert exon1 == exon2  # Same ID and same locus
        assert exon1 != exon3  # Same ID but different locus
        assert exon1 != exon4  # Different ID (even though same locus)

    def test_exon_with_different_strands(self):
        """Test exon creation with different strand orientations."""
        transcript_pos = self.create_mock_transcript("chr1", 1000, 5000, "+")
        transcript_neg = self.create_mock_transcript("chr1", 1000, 5000, "-")
        genome_mock = Mock()
        
        exon_pos = Exon(
            id="ENSE00000001",
            chr="chr1",
            start=2000,
            end=2100,
            strand="+",
            transcript=transcript_pos,
            genome=genome_mock
        )
        
        exon_neg = Exon(
            id="ENSE00000002",
            chr="chr1",
            start=2000,
            end=2100,
            strand="-",
            transcript=transcript_neg,
            genome=genome_mock
        )
        
        assert exon_pos.strand == "+"
        assert exon_neg.strand == "-"
        assert exon_pos.chr == exon_neg.chr
        assert len(exon_pos) == len(exon_neg)

    def test_exon_chromosome_consistency(self):
        """Test that exon's chromosome matches transcript's chromosome."""
        transcript = self.create_mock_transcript("chr22", 10000, 20000, "+")
        genome_mock = Mock()
        
        exon = Exon(
            id="ENSE00000022",
            chr="chr22",
            start=15000,
            end=15100,
            strand="+",
            transcript=transcript,
            genome=genome_mock
        )
        
        # Exon's chromosome should match transcript's chromosome
        assert exon.chr == transcript.chr
        assert exon.chr == "chr22"
        assert exon.locus.chr == "chr22"

    def test_exon_multiple_in_transcript(self):
        """Test multiple exons belonging to the same transcript."""
        transcript = self.create_mock_transcript("chr1", 1000, 5000, "+")
        genome_mock = Mock()
        
        exons = []
        for i in range(3):
            exon = Exon(
                id=f"ENSE0000000{i+1}",
                chr="chr1",
                start=1500 + i * 500,
                end=1600 + i * 500,
                strand="+",
                transcript=transcript,
                genome=genome_mock,
                exon_number=i+1
            )
            exons.append(exon)
        
        # All exons should belong to the same transcript
        for exon in exons:
            assert exon.get_transcript() == transcript
            assert exon.chr == "chr1"
            assert exon.strand == "+"
        
        # Each exon should have different coordinates but same transcript
        assert exons[0].start == 1500
        assert exons[1].start == 2000  
        assert exons[2].start == 2500
        
        assert all(exon.exon_number == i+1 for i, exon in enumerate(exons))

    def test_exon_standalone_creation(self):
        """Test creating an exon without transcript or genome dependencies."""
        exon = Exon(
            id="ENSE00000001",
            chr="chr1",
            start=1100,
            end=1200,
            strand="+"
        )
        
        assert exon.id == "ENSE00000001"
        assert exon.chr == "chr1"
        assert exon.start == 1100
        assert exon.end == 1200
        assert exon.strand == "+"
        assert exon._parent is None
        assert exon._genome is None
        assert len(exon) == 101

    def test_exon_standalone_with_kwargs(self):
        """Test creating standalone exon with additional attributes."""
        exon = Exon(
            id="ENSE00000002",
            chr="chr2",
            start=3000,
            end=3100,
            strand="-",
            exon_number=1,
            phase=0,
            rank=1
        )
        
        assert exon.exon_number == 1
        assert exon.phase == 0
        assert exon.rank == 1
        assert exon._parent is None
        assert exon._genome is None

    def test_exon_standalone_get_transcript_returns_none(self):
        """Test that get_transcript returns None for standalone exon."""
        exon = Exon(
            id="ENSE00000003",
            chr="chr1",
            start=1100,
            end=1200,
            strand="+"
        )
        
        # Should return None since no parent is set
        assert exon._parent is None

    def test_exon_standalone_sequence_property_fails(self):
        """Test that sequence property fails gracefully for standalone exon."""
        exon = Exon(
            id="ENSE00000004",
            chr="chr1",
            start=1100,
            end=1200,
            strand="+"
        )
        
        # Should raise AttributeError since no transcript is set (get_transcript returns None)
        with pytest.raises(AttributeError):
            _ = exon.sequence

    def test_exon_standalone_equality_and_hashing(self):
        """Test equality and hashing for standalone exons."""
        exon1 = Exon(
            id="SAME_ID",
            chr="chr1",
            start=1100,
            end=1200,
            strand="+"
        )
        
        exon2 = Exon(
            id="SAME_ID",
            chr="chr1",
            start=1100,
            end=1200,
            strand="+"
        )
        
        exon3 = Exon(
            id="SAME_ID",
            chr="chr1",
            start=2000,  # Different coordinates
            end=2100,
            strand="+"
        )
        
        assert exon1 == exon2  # Same ID and locus
        assert exon1 != exon3  # Same ID but different locus
        assert hash(exon1) == hash(exon2)
        assert hash(exon1) != hash(exon3)

    def test_exon_standalone_inheritance_from_genome_element(self):
        """Test that standalone Exon properly inherits from GenomeElement."""
        exon = Exon(
            id="ENSE00000005",
            chr="chr5",
            start=12000,
            end=12100,
            strand="-"
        )
        
        # Should have GenomeElement properties
        assert exon.id == "ENSE00000005"
        assert exon.chr == "chr5"
        assert exon.start == 12000
        assert exon.end == 12100
        assert exon.strand == "-"
        assert len(exon) == 101
        
        # Should have proper repr
        expected_repr = "Exon(id='ENSE00000005', locus=Locus(chr5:12000-12100, strand=-))"
        assert repr(exon) == expected_repr




