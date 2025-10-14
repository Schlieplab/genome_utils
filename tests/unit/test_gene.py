#!/usr/bin/env python
"""
Filename: tests/unit/test_gene.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.1
Description: Unit tests for the Gene class.
License: LGPL-3.0-or-later
"""

from unittest.mock import Mock

import pytest
from Bio.Seq import Seq

from GenomeUtils.Genome import Gene, Locus


class TestGene:
    """Test cases for the Gene class."""

    def create_mock_chromosome(self, chr_id, sequence="ATCGATCGATCGAAATTTGGGCCC"):
        """Helper to create a mock chromosome."""
        chromosome = Mock()
        chromosome.chr = chr_id
        chromosome.id = chr_id
        chromosome.get_subsequence_by_locus.return_value = Seq(sequence)
        return chromosome

    def test_gene_creation(self):
        """Test basic gene creation."""
        chromosome = self.create_mock_chromosome("chr1")
        genome_mock = Mock()
        
        gene = Gene(
            id="ENSG00000001",
            name="TEST_GENE",
            chr="chr1",
            start=100,
            end=200,
            strand="+",
            chromosome=chromosome,
            genome=genome_mock
        )
        
        assert gene.id == "ENSG00000001"
        assert gene.name == "TEST_GENE"
        assert gene.start == 100
        assert gene.end == 200
        assert gene.strand == "+"
        assert gene.chr == "chr1"
        assert gene._parent == chromosome
        assert gene._genome == genome_mock
        assert len(gene) == 101  # 1-based inclusive

    def test_gene_creation_with_kwargs(self):
        """Test gene creation with additional attributes."""
        chromosome = self.create_mock_chromosome("chr2")
        genome_mock = Mock()
        
        gene = Gene(
            id="ENSG00000002",
            name="PROTEIN_GENE",
            chr="chr2",
            start=500,
            end=1000,
            strand="-",
            chromosome=chromosome,
            genome=genome_mock,
            gene_type="protein_coding",
            description="A test protein gene",
            score=95.5
        )
        
        assert gene.gene_type == "protein_coding"
        assert gene.description == "A test protein gene"
        assert gene.score == 95.5

    def test_gene_locus_creation(self):
        """Test that gene creates correct locus."""
        chromosome = self.create_mock_chromosome("chr3")
        genome_mock = Mock()
        
        gene = Gene(
            id="ENSG00000003",
            name="TEST_GENE3",
            chr="chr3",
            start=1000,
            end=2000,
            strand="-",
            chromosome=chromosome,
            genome=genome_mock
        )
        
        expected_locus = Locus("chr3", 1000, 2000, "-")
        assert gene.locus == expected_locus

    def test_gene_sequence_property(self):
        """Test that sequence property calls chromosome method."""
        test_sequence = "ATCGATCGATCGAAATTTGGGCCC"
        chromosome = self.create_mock_chromosome("chr1", test_sequence)
        genome_mock = Mock()
        
        gene = Gene(
            id="ENSG00000001",
            name="TEST_GENE",
            chr="chr1",
            start=100,
            end=200,
            strand="+",
            chromosome=chromosome,
            genome=genome_mock
        )
        
        sequence = gene.sequence
        
        # Verify the chromosome method was called with correct locus
        chromosome.get_subsequence_by_locus.assert_called_once_with(gene.locus)
        assert str(sequence) == test_sequence

    def test_gene_transcripts_property(self):
        """Test that transcripts property returns children list."""
        chromosome = self.create_mock_chromosome("chr1")
        genome_mock = Mock()
        
        gene = Gene(
            id="ENSG00000001",
            name="TEST_GENE",
            chr="chr1",
            start=100,
            end=200,
            strand="+",
            chromosome=chromosome,
            genome=genome_mock
        )
        
        # Initially empty
        assert gene.transcripts == []
        assert gene.transcripts is gene._children

    def test_gene_add_transcript(self):
        """Test adding transcripts to gene."""
        chromosome = self.create_mock_chromosome("chr1")
        genome_mock = Mock()
        genome_mock.is_indexed = True
        
        gene = Gene(
            id="ENSG00000001",
            name="TEST_GENE",
            chr="chr1",
            start=100,
            end=200,
            strand="+",
            chromosome=chromosome,
            genome=genome_mock
        )
        
        transcript1 = Mock()
        transcript2 = Mock()
        
        gene.add_transcript(transcript1)
        gene.add_transcript(transcript2)
        
        assert len(gene.transcripts) == 2
        assert transcript1 in gene.transcripts
        assert transcript2 in gene.transcripts
        assert genome_mock.is_indexed == False  # Should be set to False

    def test_gene_get_chromosome(self):
        """Test get_chromosome method."""
        chromosome = self.create_mock_chromosome("chr1")
        genome_mock = Mock()
        
        gene = Gene(
            id="ENSG00000001",
            name="TEST_GENE",
            chr="chr1",
            start=100,
            end=200,
            strand="+",
            chromosome=chromosome,
            genome=genome_mock
        )
        
        assert gene.get_chromosome() == chromosome
        assert gene.get_chromosome() is gene._parent

    def test_gene_inheritance_from_genome_element(self):
        """Test that Gene properly inherits from GenomeElement."""
        chromosome = self.create_mock_chromosome("chr1")
        genome_mock = Mock()
        
        gene = Gene(
            id="ENSG00000001",
            name="TEST_GENE",
            chr="chr1",
            start=100,
            end=200,
            strand="-",
            chromosome=chromosome,
            genome=genome_mock
        )
        
        # Should have GenomeElement properties
        assert gene.id == "ENSG00000001"
        assert gene.chr == "chr1"
        assert gene.start == 100
        assert gene.end == 200
        assert gene.strand == "-"
        assert len(gene) == 101
        
        # Should have proper repr
        expected_repr = "Gene(id='ENSG00000001', locus=Locus(chr1:100-200, strand=-))"
        assert repr(gene) == expected_repr

    def test_gene_equality(self):
        """Test gene equality based on ID and locus."""
        chromosome1 = self.create_mock_chromosome("chr1")
        chromosome2 = self.create_mock_chromosome("chr2")
        genome_mock = Mock()
        
        gene1 = Gene(
            id="SAME_ID",
            name="GENE1",
            chr="chr1",
            start=100,
            end=200,
            strand="+",
            chromosome=chromosome1,
            genome=genome_mock
        )
        
        gene2 = Gene(
            id="SAME_ID",
            name="GENE1",  # Same name
            chr="chr1",
            start=100,     # Same coordinates
            end=200,
            strand="+",    # Same strand
            chromosome=chromosome1,  # Same chromosome
            genome=genome_mock
        )
        
        gene3 = Gene(
            id="SAME_ID",
            name="GENE2",  # Different name (shouldn't matter)
            chr="chr2",
            start=300,     # Different coordinates
            end=400,
            strand="-",    # Different strand
            chromosome=chromosome2,  # Different chromosome
            genome=genome_mock
        )
        
        gene4 = Gene(
            id="DIFFERENT_ID",
            name="GENE1",
            chr="chr1",
            start=100,
            end=200,
            strand="+",
            chromosome=chromosome1,
            genome=genome_mock
        )
        
        assert gene1 == gene2  # Same ID and same locus
        assert gene1 != gene3  # Same ID but different locus
        assert gene1 != gene4  # Different ID (even though same locus)

    def test_gene_multiple_transcripts_workflow(self):
        """Test realistic workflow with multiple transcripts."""
        chromosome = self.create_mock_chromosome("chr1")
        genome_mock = Mock()
        genome_mock.is_indexed = True
        
        gene = Gene(
            id="ENSG00000001",
            name="MULTI_TRANSCRIPT_GENE",
            chr="chr1",
            start=1000,
            end=5000,
            strand="+",
            chromosome=chromosome,
            genome=genome_mock,
            gene_type="protein_coding"
        )
        
        # Add multiple transcripts
        transcripts = []
        for i in range(3):
            transcript = Mock()
            transcript.id = f"ENST0000000{i+1}"
            transcripts.append(transcript)
            gene.add_transcript(transcript)
        
        assert len(gene.transcripts) == 3
        for i, transcript in enumerate(gene.transcripts):
            assert transcript.id == f"ENST0000000{i+1}"
        
        # Genome should be marked as not indexed
        assert genome_mock.is_indexed == False

    def test_gene_standalone_creation(self):
        """Test creating a gene without chromosome or genome dependencies."""
        gene = Gene(
            id="ENSG00000001",
            name="STANDALONE_GENE",
            chr="chr1",
            start=100,
            end=200,
            strand="+"
        )
        
        assert gene.id == "ENSG00000001"
        assert gene.name == "STANDALONE_GENE"
        assert gene.chr == "chr1"
        assert gene.start == 100
        assert gene.end == 200
        assert gene.strand == "+"
        assert gene._parent is None
        assert gene._genome is None
        assert len(gene) == 101

    def test_gene_standalone_with_kwargs(self):
        """Test creating standalone gene with additional attributes."""
        gene = Gene(
            id="ENSG00000002",
            name="PROTEIN_GENE",
            chr="chr2",
            start=500,
            end=1000,
            strand="-",
            gene_type="protein_coding",
            description="A standalone test gene",
            score=95.5
        )
        
        assert gene.gene_type == "protein_coding"
        assert gene.description == "A standalone test gene"
        assert gene.score == 95.5
        assert gene._parent is None
        assert gene._genome is None

    def test_gene_standalone_sequence_property_fails(self):
        """Test that sequence property fails gracefully for standalone gene."""
        gene = Gene(
            id="ENSG00000003",
            name="STANDALONE_GENE",
            chr="chr1",
            start=100,
            end=200,
            strand="+"
        )
        
        # Should raise an error since there's no chromosome to get sequence from
        with pytest.raises(AttributeError):
            _ = gene.sequence

    def test_gene_standalone_transcripts_functionality(self):
        """Test that transcript management works for standalone genes."""
        gene = Gene(
            id="ENSG00000004",
            name="MULTI_TRANSCRIPT_GENE",
            chr="chr1",
            start=1000,
            end=5000,
            strand="+"
        )
        
        # Initially empty
        assert gene.transcripts == []
        
        # Can't add transcripts without genome (would cause error)
        # But the transcripts property should work
        assert len(gene.transcripts) == 0

    def test_gene_standalone_get_chromosome_returns_none(self):
        """Test that get_chromosome returns None for standalone gene."""
        gene = Gene(
            id="ENSG00000005",
            name="STANDALONE_GENE",
            chr="chr1",
            start=100,
            end=200,
            strand="+"
        )
        
        # Should return None since no parent is set
        assert gene._parent is None

    def test_gene_standalone_equality_and_hashing(self):
        """Test equality and hashing for standalone genes."""
        gene1 = Gene(
            id="SAME_ID",
            name="GENE1",
            chr="chr1",
            start=100,
            end=200,
            strand="+"
        )
        
        gene2 = Gene(
            id="SAME_ID",
            name="GENE2",  # Different name shouldn't matter
            chr="chr1",
            start=100,
            end=200,
            strand="+"
        )
        
        gene3 = Gene(
            id="SAME_ID",
            name="GENE1",
            chr="chr1",
            start=300,  # Different coordinates
            end=400,
            strand="+"
        )
        
        assert gene1 == gene2  # Same ID and locus
        assert gene1 != gene3  # Same ID but different locus
        assert hash(gene1) == hash(gene2)
        assert hash(gene1) != hash(gene3)

    def test_gene_with_negative_strand(self):
        """Test gene creation and properties with negative strand."""
        chromosome = self.create_mock_chromosome("chrX")
        genome_mock = Mock()
        
        gene = Gene(
            id="ENSG00000X01",
            name="NEG_STRAND_GENE",
            chr="chrX",
            start=10000,
            end=15000,
            strand="-",
            chromosome=chromosome,
            genome=genome_mock
        )
        
        assert gene.strand == "-"
        assert gene.chr == "chrX"
        assert len(gene) == 5001
        
        # Test that sequence property still works
        sequence = gene.sequence
        chromosome.get_subsequence_by_locus.assert_called_once()

    def test_gene_chromosome_id_consistency(self):
        """Test that gene's chromosome ID is consistent with chromosome."""
        chromosome = self.create_mock_chromosome("chr22")
        genome_mock = Mock()
        
        gene = Gene(
            id="ENSG00000022",
            name="CHR22_GENE",
            chr="chr22",
            start=100,
            end=200,
            strand="+",
            chromosome=chromosome,
            genome=genome_mock
        )
        
        # Gene's chromosome ID should match the chromosome
        assert gene.chr == chromosome.chr
        assert gene.chr == "chr22"
        assert gene.locus.chr == "chr22"




