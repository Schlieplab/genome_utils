#!/usr/bin/env python
"""
Filename: test_transcript.py
Author: Arash Ayat
Copyright: 2026, Alexander Schliep
Version: 0.1.3
Description: Unit tests for the Transcript class.
License: LGPL-3.0-or-later
"""

from unittest.mock import Mock

import pytest
from Bio.Seq import Seq

from GenomeUtils.Genome import Locus, Transcript, TranscriptLocus


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

    def create_mock_exon(self, start, end, strand="+", chr_id="chr1"):
        """Helper to create a mock exon."""
        exon = Mock()
        exon.chr = chr_id
        exon.start = start
        exon.end = end
        exon.strand = strand
        exon.__len__ = lambda self: end - start + 1
        return exon

    def create_segment_transcript(self, strand="+", sequence_length=250):
        """Create a three-exon transcript for coordinate conversion tests."""
        transcript = Transcript(
            id="ENST00000001",
            chr="chr1",
            start=1000,
            end=2000,
            strand=strand,
            sequence=Seq("A" * sequence_length),
            gene=self.create_mock_gene("chr1", 1000, 2000, strand),
            genome=Mock(),
        )
        for start, end in [(1100, 1199), (1300, 1349), (1500, 1599)]:
            transcript.add_exon(self.create_mock_exon(start, end, strand))
        return transcript

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
        assert transcript.exons is transcript._exons

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

    def test_to_genomic_loci_inside_one_exon_positive_strand(self):
        transcript = self.create_segment_transcript()

        loci = transcript.to_genomic_loci(
            TranscriptLocus(transcript.id, 10, 20)
        )

        assert isinstance(loci, list)
        assert loci == [Locus("chr1", 1110, 1119, "+")]

    def test_to_genomic_loci_crosses_two_exons_positive_strand(self):
        transcript = self.create_segment_transcript()

        loci = transcript.to_genomic_loci(
            TranscriptLocus(transcript.id, 90, 120)
        )

        assert loci == [
            Locus("chr1", 1190, 1199, "+"),
            Locus("chr1", 1300, 1319, "+"),
        ]

    def test_to_genomic_loci_crosses_more_than_two_exons(self):
        transcript = self.create_segment_transcript()

        loci = transcript.to_genomic_loci(
            TranscriptLocus(transcript.id, 90, 170)
        )

        assert loci == [
            Locus("chr1", 1190, 1199, "+"),
            Locus("chr1", 1300, 1349, "+"),
            Locus("chr1", 1500, 1519, "+"),
        ]

    def test_to_genomic_loci_negative_strand_uses_transcript_order(self):
        transcript = self.create_segment_transcript(strand="-")

        loci = transcript.to_genomic_loci(
            TranscriptLocus(transcript.id, 90, 170)
        )

        assert loci == [
            Locus("chr1", 1500, 1509, "-"),
            Locus("chr1", 1300, 1349, "-"),
            Locus("chr1", 1180, 1199, "-"),
        ]

    def test_to_genomic_loci_entire_spliced_sequence(self):
        transcript = self.create_segment_transcript()

        loci = transcript.to_genomic_loci(
            TranscriptLocus(transcript.id, 0, len(transcript.sequence))
        )

        assert loci == [
            Locus("chr1", 1100, 1199, "+"),
            Locus("chr1", 1300, 1349, "+"),
            Locus("chr1", 1500, 1599, "+"),
        ]

    @pytest.mark.parametrize(
        ("strand", "first", "last"),
        [
            ("+", Locus("chr1", 1100, 1100, "+"), Locus("chr1", 1599, 1599, "+")),
            ("-", Locus("chr1", 1599, 1599, "-"), Locus("chr1", 1100, 1100, "-")),
        ],
    )
    def test_to_genomic_loci_first_and_last_nucleotide(self, strand, first, last):
        transcript = self.create_segment_transcript(strand=strand)

        assert transcript.to_genomic_loci(
            TranscriptLocus(transcript.id, 0, 1)
        ) == [first]
        assert transcript.to_genomic_loci(
            TranscriptLocus(transcript.id, 249, 250)
        ) == [last]

    def test_to_genomic_loci_rejects_end_beyond_sequence(self):
        transcript = self.create_segment_transcript()

        with pytest.raises(ValueError, match="out of bounds"):
            transcript.to_genomic_loci(TranscriptLocus(transcript.id, 0, 251))

    def test_to_genomic_loci_rejects_mismatched_transcript_id(self):
        transcript = self.create_segment_transcript()

        with pytest.raises(ValueError, match="must match"):
            transcript.to_genomic_loci(TranscriptLocus("different", 0, 1))

    def test_to_genomic_loci_rejects_transcript_without_exons(self):
        transcript = Transcript(
            id="ENST00000001",
            chr="chr1",
            start=1000,
            end=2000,
            strand="+",
            sequence=Seq("A"),
        )

        with pytest.raises(ValueError, match="at least one exon"):
            transcript.to_genomic_loci(TranscriptLocus(transcript.id, 0, 1))

    def test_to_genomic_loci_rejects_exon_sequence_length_mismatch(self):
        transcript = self.create_segment_transcript(sequence_length=249)

        with pytest.raises(ValueError, match="Sum of exon lengths"):
            transcript.to_genomic_loci(TranscriptLocus(transcript.id, 0, 1))

    @pytest.mark.parametrize(
        ("exon_chr", "exon_strand", "message"),
        [
            ("chr2", "+", "chromosome"),
            ("chr1", "-", "strand"),
        ],
    )
    def test_to_genomic_loci_rejects_exon_on_wrong_locus(
        self,
        exon_chr,
        exon_strand,
        message,
    ):
        transcript = Transcript(
            id="ENST00000001",
            chr="chr1",
            start=1000,
            end=2000,
            strand="+",
            sequence=Seq("A" * 100),
            genome=Mock(),
        )
        transcript.add_exon(
            self.create_mock_exon(1100, 1199, exon_strand, exon_chr)
        )

        with pytest.raises(ValueError, match=message):
            transcript.to_genomic_loci(TranscriptLocus(transcript.id, 0, 1))

    def test_to_genomic_loci_rejects_overlapping_exons(self):
        transcript = Transcript(
            id="ENST00000001",
            chr="chr1",
            start=1000,
            end=2000,
            strand="+",
            sequence=Seq("A" * 202),
            genome=Mock(),
        )
        transcript.add_exon(self.create_mock_exon(1100, 1200))
        transcript.add_exon(self.create_mock_exon(1200, 1300))

        with pytest.raises(ValueError, match="must not overlap"):
            transcript.to_genomic_loci(TranscriptLocus(transcript.id, 0, 1))

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

    def test_to_genomic_loci_single_exon_edge_cases(self):
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
        loci = transcript.to_genomic_loci(TranscriptLocus(transcript.id, 0, 1))
        assert loci == [Locus("chr1", 1100, 1100, "+")]
        
        # Test last position (end of transcript)
        loci = transcript.to_genomic_loci(TranscriptLocus(transcript.id, 99, 100))
        assert loci == [Locus("chr1", 1199, 1199, "+")]
        
        # Test range covering entire transcript
        loci = transcript.to_genomic_loci(TranscriptLocus(transcript.id, 0, 100))
        assert loci == [Locus("chr1", 1100, 1199, "+")]

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
        test_seq = Seq("A" * 150)
        
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
        exon1.chr = "chr1"
        exon1.start = 1100
        exon1.end = 1199
        exon1.strand = "+"
        exon1.__len__ = Mock(return_value=100)
        
        exon2 = Mock()
        exon2.chr = "chr1"
        exon2.start = 1300
        exon2.end = 1349
        exon2.strand = "+"
        exon2.__len__ = Mock(return_value=50)
        
        # Manually add to exons (bypassing add_exon to avoid genome dependency)
        transcript._exons = [exon1, exon2]
        
        # Test coordinate conversion
        loci = transcript.to_genomic_loci(TranscriptLocus(transcript.id, 50, 51))
        assert loci == [Locus("chr1", 1150, 1150, "+")]

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
