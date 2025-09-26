"""Tests for the Transcript class."""

import pytest
from unittest.mock import Mock
from Bio.Seq import Seq

from GenomeUtils.Genome import Transcript, Locus


class TestTranscript:
    """Test cases for the Transcript class."""

    def create_mock_gene(self, chr_id="chr1", start=1000, end=5000, strand="+"):
        """Helper to create a mock gene."""
        gene = Mock()
        gene.chr = chr_id
        gene.start = start
        gene.end = end
        gene.strand = strand
        return gene

    def create_mock_exon(self, start, end, strand="+"):
        """Helper to create a mock exon."""
        exon = Mock()
        exon.start = start
        exon.end = end
        exon.strand = strand
        exon.__len__ = lambda self: end - start + 1
        return exon

    def test_transcript_creation(self):
        """Test basic transcript creation."""
        gene = self.create_mock_gene()
        genome_mock = Mock()
        test_seq = Seq("ATCGATCGATCGAAATTTGGGCCC")
        
        transcript = Transcript(
            id="ENST00000001",
            chr="chr1",
            start=1100,
            end=1500,
            strand="+",
            sequence=test_seq,
            gene=gene,
            genome=genome_mock
        )
        
        assert transcript.id == "ENST00000001"
        assert transcript.start == 1100
        assert transcript.end == 1500
        assert transcript.strand == "+"
        assert transcript.sequence == test_seq
        assert transcript._parent == gene
        assert transcript._genome == genome_mock
        assert len(transcript) == 401  # 1-based inclusive genomic coordinates

    def test_transcript_creation_with_kwargs(self):
        """Test transcript creation with additional attributes."""
        gene = self.create_mock_gene()
        genome_mock = Mock()
        test_seq = Seq("ATCGATCGATCG")
        
        transcript = Transcript(
            id="ENST00000002",
            chr="chr3",
            start=2000,
            end=3000,
            strand="-",
            sequence=test_seq,
            gene=gene,
            genome=genome_mock,
            transcript_type="protein_coding",
            support_level=1
        )
        
        assert transcript.transcript_type == "protein_coding"
        assert transcript.support_level == 1

    def test_transcript_locus_creation(self):
        """Test that transcript creates correct locus."""
        gene = self.create_mock_gene("chr3", 1000, 5000, "-")
        genome_mock = Mock()
        test_seq = Seq("ATCGATCGATCG")
        
        transcript = Transcript(
            id="ENST00000003",
            chr="chr3",
            start=2000,
            end=3000,
            strand="-",
            sequence=test_seq,
            gene=gene,
            genome=genome_mock
        )
        
        expected_locus = Locus("chr3", 2000, 3000, "-")
        assert transcript.locus == expected_locus

    def test_transcript_exons_property(self):
        """Test that exons property returns children list."""
        gene = self.create_mock_gene()
        genome_mock = Mock()
        test_seq = Seq("ATCGATCGATCG")
        
        transcript = Transcript(
            id="ENST00000001",
            chr="chr1",
            start=1100,
            end=1500,
            strand="+",
            sequence=test_seq,
            gene=gene,
            genome=genome_mock
        )
        
        # Initially empty
        assert transcript.exons == []
        assert transcript.exons is transcript._children

    def test_transcript_add_exon_positive_strand(self):
        """Test adding exons to transcript on positive strand (sorted ascending)."""
        gene = self.create_mock_gene()
        genome_mock = Mock()
        genome_mock.is_indexed = True
        test_seq = Seq("ATCGATCGATCG")
        
        transcript = Transcript(
            id="ENST00000001",
            chr="chr1",
            start=1000,
            end=2000,
            strand="+",
            sequence=test_seq,
            gene=gene,
            genome=genome_mock
        )
        
        exon1 = self.create_mock_exon(1500, 1600, "+")
        exon2 = self.create_mock_exon(1100, 1200, "+")  # Earlier exon added later
        exon3 = self.create_mock_exon(1700, 1800, "+")
        
        transcript.add_exon(exon1)
        transcript.add_exon(exon2)  # Should be inserted at beginning
        transcript.add_exon(exon3)
        
        # Should be sorted by start position for positive strand
        assert len(transcript.exons) == 3
        assert transcript.exons[0] == exon2  # start=1100
        assert transcript.exons[1] == exon1  # start=1500
        assert transcript.exons[2] == exon3  # start=1700
        assert genome_mock.is_indexed == False

    def test_transcript_add_exon_negative_strand(self):
        """Test adding exons to transcript on negative strand (sorted descending)."""
        gene = self.create_mock_gene("chr1", 1000, 2000, "-")
        genome_mock = Mock()
        genome_mock.is_indexed = True
        test_seq = Seq("ATCGATCGATCG")
        
        transcript = Transcript(
            id="ENST00000001",
            chr="chr1",
            start=1000,
            end=2000,
            strand="-",
            sequence=test_seq,
            gene=gene,
            genome=genome_mock
        )
        
        exon1 = self.create_mock_exon(1500, 1600, "-")
        exon2 = self.create_mock_exon(1100, 1200, "-")  # Earlier exon added later
        exon3 = self.create_mock_exon(1700, 1800, "-")
        
        transcript.add_exon(exon1)
        transcript.add_exon(exon2)  # Should be inserted at end for negative strand
        transcript.add_exon(exon3)
        
        # Should be sorted descending by start position for negative strand
        assert len(transcript.exons) == 3
        assert transcript.exons[0] == exon3  # start=1700
        assert transcript.exons[1] == exon1  # start=1500
        assert transcript.exons[2] == exon2  # start=1100

    def test_transcript_get_gene(self):
        """Test get_gene method."""
        gene = self.create_mock_gene()
        genome_mock = Mock()
        test_seq = Seq("ATCGATCGATCG")
        
        transcript = Transcript(
            id="ENST00000001",
            chr="chr1",
            start=1100,
            end=1500,
            strand="+",
            sequence=test_seq,
            gene=gene,
            genome=genome_mock
        )
        
        assert transcript.get_gene() == gene
        assert transcript.get_gene() is transcript._parent

    def test_transcript_exon_intervals(self):
        """Test exon_intervals method."""
        gene = self.create_mock_gene()
        genome_mock = Mock()
        test_seq = Seq("ATCGATCGATCG")
        
        transcript = Transcript(
            id="ENST00000001",
            chr="chr1",
            start=1000,
            end=2000,
            strand="+",
            sequence=test_seq,
            gene=gene,
            genome=genome_mock
        )
        
        exon1 = self.create_mock_exon(1100, 1200, "+")
        exon2 = self.create_mock_exon(1500, 1600, "+")
        exon3 = self.create_mock_exon(1700, 1800, "+")
        
        transcript.add_exon(exon1)
        transcript.add_exon(exon2)
        transcript.add_exon(exon3)
        
        intervals = transcript.exon_intervals()
        expected = [(1100, 1200), (1500, 1600), (1700, 1800)]
        assert intervals == expected

    def test_transcript_to_genomic_pos_single_point_positive_strand(self):
        """Test transcript to genomic coordinate conversion for single point on positive strand."""
        gene = self.create_mock_gene("chr1", 1000, 2000, "+")
        genome_mock = Mock()
        test_seq = Seq("A" * 250)  # 250 bp transcript
        
        transcript = Transcript(
            id="ENST00000001",
            chr="chr1",
            start=1000,
            end=2000,
            strand="+",
            sequence=test_seq,
            gene=gene,
            genome=genome_mock
        )
        
        # Add exons: [1100-1199] (100bp), [1300-1349] (50bp), [1500-1599] (100bp)
        exon1 = self.create_mock_exon(1100, 1199, "+")  # 100 bp
        exon2 = self.create_mock_exon(1300, 1349, "+")  # 50 bp  
        exon3 = self.create_mock_exon(1500, 1599, "+")  # 100 bp
        
        transcript.add_exon(exon1)
        transcript.add_exon(exon2) 
        transcript.add_exon(exon3)
        
        # Test position in first exon (transcript pos 50 -> genomic pos 1150)
        locus = transcript.transcript_to_genomic_pos(50)
        assert locus.chr == "chr1"
        assert locus.start == 1150
        assert locus.end == 1150
        assert locus.strand == "+"
        
        # Test position in second exon (transcript pos 120 -> genomic pos 1320)
        locus = transcript.transcript_to_genomic_pos(120)
        assert locus.start == 1320
        assert locus.end == 1320
        
        # Test position in third exon (transcript pos 170 -> genomic pos 1520)
        locus = transcript.transcript_to_genomic_pos(170)
        assert locus.start == 1520
        assert locus.end == 1520

    def test_transcript_to_genomic_pos_single_point_negative_strand(self):
        """Test transcript to genomic coordinate conversion for single point on negative strand."""
        gene = self.create_mock_gene("chr1", 1000, 2000, "-")
        genome_mock = Mock()
        test_seq = Seq("A" * 250)  # 250 bp transcript
        
        transcript = Transcript(
            id="ENST00000001",
            chr="chr1",
            start=1000,
            end=2000,
            strand="-",
            sequence=test_seq,
            gene=gene,
            genome=genome_mock
        )
        
        # Add exons for negative strand: [1500-1599] (100bp), [1300-1349] (50bp), [1100-1199] (100bp)
        # These will be sorted in transcriptional order (descending for negative strand)
        exon1 = self.create_mock_exon(1500, 1599, "-")  # 100 bp (first in transcript)
        exon2 = self.create_mock_exon(1300, 1349, "-")  # 50 bp (second in transcript)
        exon3 = self.create_mock_exon(1100, 1199, "-")  # 100 bp (third in transcript)
        
        transcript.add_exon(exon1)
        transcript.add_exon(exon2)
        transcript.add_exon(exon3)
        
        # Test position in first exon (transcript pos 50 -> genomic pos 1549)
        locus = transcript.transcript_to_genomic_pos(50)
        assert locus.chr == "chr1"
        assert locus.start == 1549
        assert locus.end == 1549
        assert locus.strand == "-"

    def test_transcript_to_genomic_pos_range_within_exon(self):
        """Test transcript to genomic coordinate conversion for range within single exon."""
        gene = self.create_mock_gene("chr1", 1000, 2000, "+")
        genome_mock = Mock()
        test_seq = Seq("A" * 250)
        
        transcript = Transcript(
            id="ENST00000001",
            chr="chr1",
            start=1000,
            end=2000,
            strand="+",
            sequence=test_seq,
            gene=gene,
            genome=genome_mock
        )
        
        exon1 = self.create_mock_exon(1100, 1199, "+")  # 100 bp
        transcript.add_exon(exon1)
        
        # Test range within first exon (transcript pos 10-19 -> genomic pos 1110-1119)
        locus = transcript.transcript_to_genomic_pos(10, 20)
        assert locus.chr == "chr1"
        assert locus.start == 1110
        assert locus.end == 1119
        assert locus.strand == "+"

    def test_transcript_to_genomic_pos_range_spanning_exons(self):
        """Test transcript to genomic coordinate conversion for range spanning multiple exons."""
        gene = self.create_mock_gene("chr1", 1000, 2000, "+")
        genome_mock = Mock()
        test_seq = Seq("A" * 250)
        
        transcript = Transcript(
            id="ENST00000001",
            chr="chr1",
            start=1000,
            end=2000,
            strand="+",
            sequence=test_seq,
            gene=gene,
            genome=genome_mock
        )
        
        exon1 = self.create_mock_exon(1100, 1199, "+")  # 100 bp
        exon2 = self.create_mock_exon(1300, 1349, "+")  # 50 bp
        transcript.add_exon(exon1)
        transcript.add_exon(exon2)
        
        # Test range spanning exons (transcript pos 90-120 spans both exons)
        loci = transcript.transcript_to_genomic_pos(90, 120)
        assert isinstance(loci, list)
        assert len(loci) == 2
        
        # First locus: end of exon1 (transcript pos 90-99 -> genomic pos 1190-1199)
        assert loci[0].chr == "chr1"
        assert loci[0].start == 1190
        assert loci[0].end == 1199
        
        # Second locus: start of exon2 (transcript pos 100-119 -> genomic pos 1300-1319)
        assert loci[1].chr == "chr1"
        assert loci[1].start == 1300
        assert loci[1].end == 1319

    def test_transcript_to_genomic_pos_out_of_bounds(self):
        """Test error handling for out of bounds coordinates."""
        gene = self.create_mock_gene()
        genome_mock = Mock()
        test_seq = Seq("A" * 100)  # 100 bp transcript
        
        transcript = Transcript(
            id="ENST00000001",
            chr="chr1",
            start=1100,
            end=1500,
            strand="+",
            sequence=test_seq,
            gene=gene,
            genome=genome_mock
        )
        
        # Test out of bounds single position
        with pytest.raises(ValueError, match="Transcript position .* is out of bounds"):
            transcript.transcript_to_genomic_pos(500)  # Beyond transcript genomic length
        
        # Test out of bounds range
        with pytest.raises(ValueError, match="Transcript positions .* are out of bounds"):
            transcript.transcript_to_genomic_pos(50, 500)  # End beyond transcript length

    def test_transcript_length_calculation(self):
        """Test that transcript length matches sequence length."""
        gene = self.create_mock_gene()
        genome_mock = Mock()
        test_seq = Seq("ATCGATCGATCGAAATTTGGGCCCTTTTAAAA")  # 32 bp
        
        transcript = Transcript(
            id="ENST00000001",
            chr="chr1",
            start=1100,
            end=1500,
            strand="+",
            sequence=test_seq,
            gene=gene,
            genome=genome_mock
        )
        
        assert len(transcript.sequence) == 32
        # Note: len(transcript) returns genomic coordinate length, not sequence length

    def test_transcript_inheritance_from_genome_element(self):
        """Test that Transcript properly inherits from GenomeElement."""
        gene = self.create_mock_gene("chr5", 10000, 20000, "-")
        genome_mock = Mock()
        test_seq = Seq("ATCGATCGATCG")
        
        transcript = Transcript(
            id="ENST00000005",
            chr="chr5",
            start=12000,
            end=15000,
            strand="-",
            sequence=test_seq,
            gene=gene,
            genome=genome_mock
        )
        
        # Should have GenomeElement properties
        assert transcript.id == "ENST00000005"
        assert transcript.chr == "chr5"
        assert transcript.start == 12000
        assert transcript.end == 15000
        assert transcript.strand == "-"
        
        # Should have proper repr
        expected_repr = "Transcript(id='ENST00000005', locus=Locus(chr5:12000-15000, strand=-))"
        assert repr(transcript) == expected_repr

    def test_transcript_coordinate_conversion_edge_cases(self):
        """Test edge cases in coordinate conversion."""
        gene = self.create_mock_gene("chr1", 1000, 2000, "+")
        genome_mock = Mock()
        test_seq = Seq("A" * 100)  # 100 bp transcript
        
        transcript = Transcript(
            id="ENST00000001",
            chr="chr1",
            start=1000,
            end=2000,
            strand="+",
            sequence=test_seq,
            gene=gene,
            genome=genome_mock
        )
        
        exon1 = self.create_mock_exon(1100, 1199, "+")  # 100 bp
        transcript.add_exon(exon1)
        
        # Test position 0 (start of transcript)
        locus = transcript.transcript_to_genomic_pos(0)
        assert locus.start == 1100
        assert locus.end == 1100
        
        # Test last position (end of transcript)
        locus = transcript.transcript_to_genomic_pos(99)
        assert locus.start == 1199
        assert locus.end == 1199
        
        # Test range covering entire transcript
        locus = transcript.transcript_to_genomic_pos(0, 100)
        assert locus.start == 1100
        assert locus.end == 1199

    def test_transcript_standalone_creation(self):
        """Test creating a transcript without gene or genome dependencies."""
        test_seq = Seq("ATCGATCGATCGAAATTTGGGCCC")
        
        transcript = Transcript(
            id="ENST00000001",
            chr="chr1",
            start=1100,
            end=1500,
            strand="+",
            sequence=test_seq
        )
        
        assert transcript.id == "ENST00000001"
        assert transcript.chr == "chr1"
        assert transcript.start == 1100
        assert transcript.end == 1500
        assert transcript.strand == "+"
        assert transcript.sequence == test_seq
        assert transcript._parent is None
        assert transcript._genome is None
        assert len(transcript) == 401

    def test_transcript_standalone_with_kwargs(self):
        """Test creating standalone transcript with additional attributes."""
        test_seq = Seq("ATCGATCGATCG")
        
        transcript = Transcript(
            id="ENST00000002",
            chr="chr2",
            start=2000,
            end=3000,
            strand="-",
            sequence=test_seq,
            transcript_type="protein_coding",
            support_level=1,
            biotype="mRNA"
        )
        
        assert transcript.transcript_type == "protein_coding"
        assert transcript.support_level == 1
        assert transcript.biotype == "mRNA"
        assert transcript._parent is None
        assert transcript._genome is None

    def test_transcript_standalone_exons_functionality(self):
        """Test that exon management works for standalone transcripts."""
        test_seq = Seq("ATCGATCGATCG")
        
        transcript = Transcript(
            id="ENST00000003",
            chr="chr1",
            start=1000,
            end=2000,
            strand="+",
            sequence=test_seq
        )
        
        # Initially empty
        assert transcript.exons == []
        
        # Can't add exons without genome (would cause error)
        # But the exons property should work
        assert len(transcript.exons) == 0

    def test_transcript_standalone_get_gene_returns_none(self):
        """Test that get_gene returns None for standalone transcript."""
        test_seq = Seq("ATCGATCGATCG")
        
        transcript = Transcript(
            id="ENST00000004",
            chr="chr1",
            start=1100,
            end=1500,
            strand="+",
            sequence=test_seq
        )
        
        # Should return None since no parent is set
        assert transcript._parent is None

    def test_transcript_standalone_coordinate_conversion_with_mock_exons(self):
        """Test coordinate conversion for standalone transcript with manually added exons."""
        test_seq = Seq("A" * 250)  # 250 bp transcript
        
        transcript = Transcript(
            id="ENST00000005",
            chr="chr1",
            start=1000,
            end=2000,
            strand="+",
            sequence=test_seq
        )
        
        # Manually create mock exons for testing
        exon1 = Mock()
        exon1.start = 1100
        exon1.end = 1199
        exon1.__len__ = Mock(return_value=100)
        
        exon2 = Mock()
        exon2.start = 1300
        exon2.end = 1349
        exon2.__len__ = Mock(return_value=50)
        
        # Manually add to children (bypassing add_exon to avoid genome dependency)
        transcript._children = [exon1, exon2]
        
        # Test coordinate conversion
        locus = transcript.transcript_to_genomic_pos(50)
        assert locus.chr == "chr1"
        assert locus.start == 1150
        assert locus.end == 1150
        assert locus.strand == "+"

    def test_transcript_standalone_equality_and_hashing(self):
        """Test equality and hashing for standalone transcripts."""
        test_seq = Seq("ATCGATCGATCG")
        
        transcript1 = Transcript(
            id="SAME_ID",
            chr="chr1",
            start=1100,
            end=1500,
            strand="+",
            sequence=test_seq
        )
        
        transcript2 = Transcript(
            id="SAME_ID",
            chr="chr1",
            start=1100,
            end=1500,
            strand="+",
            sequence=test_seq
        )
        
        transcript3 = Transcript(
            id="SAME_ID",
            chr="chr1",
            start=2000,  # Different coordinates
            end=2500,
            strand="+",
            sequence=test_seq
        )
        
        assert transcript1 == transcript2  # Same ID and locus
        assert transcript1 != transcript3  # Same ID but different locus
        assert hash(transcript1) == hash(transcript2)
        assert hash(transcript1) != hash(transcript3)
