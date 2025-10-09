#!/usr/bin/env python
"""
Filename: test_chromosome.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.0
Description: Unit tests for the Chromosome class.
License: LGPL-3.0-or-later
"""

from unittest.mock import Mock

import pytest
from Bio.Seq import Seq

from GenomeUtils.Genome import Chromosome, Locus


class TestChromosome:
    """Test cases for the Chromosome class."""

    def create_mock_seq_index(self, seq_dict):
        """Helper to create a mock SeqIO.index with given sequences."""
        mock_index = {}
        
        for key, seq_str in seq_dict.items():
            mock_record = Mock()
            mock_record.seq = Seq(seq_str)
            mock_index[key] = mock_record
        
        return mock_index

    def test_chromosome_creation(self):
        """Test basic chromosome creation."""
        test_seq = "ATCGATCGATCG"
        seq_index = self.create_mock_seq_index({"chr1": test_seq})
        genome_mock = Mock()
        
        chromosome = Chromosome(
            id="chr1",
            seq_index=seq_index,
            genome=genome_mock
        )
        
        assert chromosome.id == "chr1"
        assert chromosome._seq_index == seq_index
        assert chromosome._genome == genome_mock
        assert len(chromosome) == len(test_seq)
        assert chromosome.chr == "chr1"
        assert chromosome.start == 1
        assert chromosome.end == len(test_seq)
        assert chromosome.strand == "+"

    def test_chromosome_creation_with_explicit_length(self):
        """Test chromosome creation with explicitly provided length."""
        test_seq = "ATCGATCGATCG"
        seq_index = self.create_mock_seq_index({"chr2": test_seq})
        genome_mock = Mock()
        
        chromosome = Chromosome(
            id="chr2",
            seq_index=seq_index,
            genome=genome_mock,
            length=100  # Different from actual sequence length
        )
        
        assert len(chromosome) == 100
        assert chromosome.end == 100

    def test_chromosome_creation_with_kwargs(self):
        """Test chromosome creation with additional attributes."""
        test_seq = "ATCGATCGATCG"
        seq_index = self.create_mock_seq_index({"chr3": test_seq})
        genome_mock = Mock()
        
        chromosome = Chromosome(
            id="chr3",
            seq_index=seq_index,
            genome=genome_mock,
            assembly="GRCh38",
            species="Homo sapiens"
        )
        
        assert chromosome.assembly == "GRCh38"
        assert chromosome.species == "Homo sapiens"

    def test_chromosome_sequence_property(self):
        """Test that sequence property returns the correct sequence."""
        test_seq = "ATCGATCGATCGAAATTTGGGCCC"
        seq_index = self.create_mock_seq_index({"chr1": test_seq})
        genome_mock = Mock()
        
        chromosome = Chromosome(
            id="chr1",
            seq_index=seq_index,
            genome=genome_mock
        )
        
        assert str(chromosome.sequence) == test_seq

    def test_chromosome_add_gene(self):
        """Test adding genes to chromosome."""
        test_seq = "ATCGATCGATCG"
        seq_index = self.create_mock_seq_index({"chr1": test_seq})
        genome_mock = Mock()
        genome_mock.is_indexed = True
        
        chromosome = Chromosome(
            id="chr1",
            seq_index=seq_index,
            genome=genome_mock
        )
        
        gene_mock1 = Mock()
        gene_mock2 = Mock()
        
        chromosome.add_gene(gene_mock1)
        chromosome.add_gene(gene_mock2)
        
        assert len(chromosome.genes) == 2
        assert gene_mock1 in chromosome.genes
        assert gene_mock2 in chromosome.genes
        assert chromosome._genome.is_indexed == False  # Should be set to False

    def test_chromosome_genes_property(self):
        """Test that genes property returns children list."""
        test_seq = "ATCGATCGATCG"
        seq_index = self.create_mock_seq_index({"chr1": test_seq})
        genome_mock = Mock()
        
        chromosome = Chromosome(
            id="chr1",
            seq_index=seq_index,
            genome=genome_mock
        )
        
        # Initially empty
        assert chromosome.genes == []
        
        # Add some genes
        gene_mock = Mock()
        chromosome.add_gene(gene_mock)
        
        assert chromosome.genes == [gene_mock]
        assert chromosome.genes is chromosome._children

    def test_chromosome_get_subsequence_by_locus_positive_strand(self):
        """Test getting subsequence for positive strand."""
        test_seq = "ATCGATCGATCGAAATTTGGGCCC"
        seq_index = self.create_mock_seq_index({"chr1": test_seq})
        genome_mock = Mock()
        
        chromosome = Chromosome(
            id="chr1",
            seq_index=seq_index,
            genome=genome_mock
        )
        
        # Test positive strand subsequence
        locus = Locus("chr1", 3, 8, "+")  # Should get "CGATCG"
        subseq = chromosome.get_subsequence_by_locus(locus)
        
        expected = test_seq[2:8]  # 0-based indexing: positions 2-7
        assert str(subseq) == expected

    def test_chromosome_get_subsequence_by_locus_negative_strand(self):
        """Test getting subsequence for negative strand."""
        test_seq = "ATCGATCGATCGAAATTTGGGCCC"
        seq_index = self.create_mock_seq_index({"chr1": test_seq})
        genome_mock = Mock()
        
        chromosome = Chromosome(
            id="chr1",
            seq_index=seq_index,
            genome=genome_mock
        )
        
        # Test negative strand subsequence
        locus = Locus("chr1", 3, 8, "-")
        subseq = chromosome.get_subsequence_by_locus(locus)
        
        # Should get reverse complement of "CGATCG" -> "CGATCG"
        forward_seq = test_seq[2:8]
        expected = str(Seq(forward_seq).reverse_complement())
        assert str(subseq) == expected

    def test_chromosome_get_subsequence_wrong_chromosome(self):
        """Test error when locus chromosome doesn't match."""
        test_seq = "ATCGATCGATCG"
        seq_index = self.create_mock_seq_index({"chr1": test_seq})
        genome_mock = Mock()
        
        chromosome = Chromosome(
            id="chr1",
            seq_index=seq_index,
            genome=genome_mock
        )
        
        # Locus for different chromosome
        locus = Locus("chr2", 3, 8, "+")
        
        with pytest.raises(ValueError, match="The Locus does not belong to this chromosome"):
            chromosome.get_subsequence_by_locus(locus)

    def test_chromosome_get_subsequence_out_of_bounds(self):
        """Test error when locus end is out of bounds."""
        test_seq = "ATCGATCGATCG"  # Length 12
        seq_index = self.create_mock_seq_index({"chr1": test_seq})
        genome_mock = Mock()
        
        chromosome = Chromosome(
            id="chr1",
            seq_index=seq_index,
            genome=genome_mock
        )
        
        # Locus that extends beyond sequence
        locus = Locus("chr1", 10, 20, "+")  # End at 20, but sequence is only 12 bp
        
        with pytest.raises(ValueError, match="End coordinate .* is out of bounds"):
            chromosome.get_subsequence_by_locus(locus)

    def test_chromosome_get_subsequence_invalid_strand(self):
        """Test error with invalid strand."""
        test_seq = "ATCGATCGATCG"
        seq_index = self.create_mock_seq_index({"chr1": test_seq})
        genome_mock = Mock()
        
        chromosome = Chromosome(
            id="chr1",
            seq_index=seq_index,
            genome=genome_mock
        )
        
        # Create locus with invalid strand (bypass Locus validation for testing)
        locus = Locus("chr1", 3, 8, "+")
        # Manually change strand to invalid value
        object.__setattr__(locus, 'strand', 'invalid')
        
        with pytest.raises(ValueError, match="Invalid strand"):
            chromosome.get_subsequence_by_locus(locus)

    def test_chromosome_get_subsequence_full_sequence(self):
        """Test getting the full sequence."""
        test_seq = "ATCGATCGATCG"
        seq_index = self.create_mock_seq_index({"chr1": test_seq})
        genome_mock = Mock()
        
        chromosome = Chromosome(
            id="chr1",
            seq_index=seq_index,
            genome=genome_mock
        )
        
        # Locus spanning the entire chromosome
        locus = Locus("chr1", 1, len(test_seq), "+")
        subseq = chromosome.get_subsequence_by_locus(locus)
        
        assert str(subseq) == test_seq

    def test_chromosome_get_subsequence_single_base(self):
        """Test getting a single base."""
        test_seq = "ATCGATCGATCG"
        seq_index = self.create_mock_seq_index({"chr1": test_seq})
        genome_mock = Mock()
        
        chromosome = Chromosome(
            id="chr1",
            seq_index=seq_index,
            genome=genome_mock
        )
        
        # Single base locus
        locus = Locus("chr1", 5, 5, "+")  # Should get the 5th base (G)
        subseq = chromosome.get_subsequence_by_locus(locus)
        
        expected = test_seq[4]  # 0-based index 4
        assert str(subseq) == expected

    def test_chromosome_inheritance_from_genome_element(self):
        """Test that Chromosome properly inherits from GenomeElement."""
        test_seq = "ATCGATCGATCG"
        seq_index = self.create_mock_seq_index({"chr1": test_seq})
        genome_mock = Mock()
        
        chromosome = Chromosome(
            id="chr1",
            seq_index=seq_index,
            genome=genome_mock
        )
        
        # Should have GenomeElement properties
        assert chromosome.id == "chr1"
        assert chromosome.chr == "chr1"
        assert chromosome.start == 1
        assert chromosome.end == len(test_seq)
        assert chromosome.strand == "+"
        assert len(chromosome) == len(test_seq)
        
        # Should have proper repr
        expected_repr = f"Chromosome(id='chr1', locus=Locus(chr1:1-{len(test_seq)}, strand=+))"
        assert repr(chromosome) == expected_repr

    def test_chromosome_standalone_creation(self):
        """Test creating a chromosome without genome dependencies."""
        test_seq = "ATCGATCGATCG"
        seq_index = self.create_mock_seq_index({"chr1": test_seq})
        
        chromosome = Chromosome(
            id="chr1",
            seq_index=seq_index
        )
        
        assert chromosome.id == "chr1"
        assert chromosome._seq_index == seq_index
        assert chromosome._genome is None
        assert len(chromosome) == len(test_seq)
        assert chromosome.chr == "chr1"
        assert chromosome.start == 1
        assert chromosome.end == len(test_seq)
        assert chromosome.strand == "+"

    def test_chromosome_standalone_with_kwargs(self):
        """Test creating standalone chromosome with additional attributes."""
        test_seq = "ATCGATCGATCG"
        seq_index = self.create_mock_seq_index({"chr3": test_seq})
        
        chromosome = Chromosome(
            id="chr3",
            seq_index=seq_index,
            assembly="GRCh38",
            species="Homo sapiens",
            build="hg38"
        )
        
        assert chromosome.assembly == "GRCh38"
        assert chromosome.species == "Homo sapiens"
        assert chromosome.build == "hg38"
        assert chromosome._genome is None

    def test_chromosome_standalone_sequence_functionality(self):
        """Test that sequence functionality works for standalone chromosome."""
        test_seq = "ATCGATCGATCGAAATTTGGGCCC"
        seq_index = self.create_mock_seq_index({"chr1": test_seq})
        
        chromosome = Chromosome(
            id="chr1",
            seq_index=seq_index
        )
        
        assert str(chromosome.sequence) == test_seq
        
        # Test subsequence functionality
        locus = Locus("chr1", 3, 8, "+")
        subseq = chromosome.get_subsequence_by_locus(locus)
        expected = test_seq[2:8]  # 0-based indexing
        assert str(subseq) == expected

    def test_chromosome_standalone_add_gene_fails(self):
        """Test that add_gene fails gracefully for standalone chromosome."""
        test_seq = "ATCGATCGATCG"
        seq_index = self.create_mock_seq_index({"chr1": test_seq})
        
        chromosome = Chromosome(
            id="chr1",
            seq_index=seq_index
        )
        
        gene_mock = Mock()
        
        # This should fail because there's no genome to set is_indexed on
        with pytest.raises(AttributeError):
            chromosome.add_gene(gene_mock)

    def test_chromosome_standalone_genes_property(self):
        """Test that genes property works for standalone chromosome."""
        test_seq = "ATCGATCGATCG"
        seq_index = self.create_mock_seq_index({"chr1": test_seq})
        
        chromosome = Chromosome(
            id="chr1",
            seq_index=seq_index
        )
        
        # Initially empty
        assert chromosome.genes == []
        assert chromosome.genes is chromosome._children

    def test_chromosome_standalone_with_explicit_length(self):
        """Test standalone chromosome creation with explicitly provided length."""
        test_seq = "ATCGATCGATCG"
        seq_index = self.create_mock_seq_index({"chr2": test_seq})
        
        chromosome = Chromosome(
            id="chr2",
            seq_index=seq_index,
            length=100  # Different from actual sequence length
        )
        
        assert len(chromosome) == 100
        assert chromosome.end == 100
        assert chromosome._genome is None

    def test_chromosome_standalone_equality_and_hashing(self):
        """Test equality and hashing for standalone chromosomes."""
        test_seq1 = "ATCGATCGATCG"
        test_seq2 = "ATCGATCGATCG"
        test_seq3 = "GGGGCCCCAAAA"
        seq_index1 = self.create_mock_seq_index({"chr1": test_seq1})
        seq_index2 = self.create_mock_seq_index({"chr1": test_seq2})
        seq_index3 = self.create_mock_seq_index({"chr2": test_seq3})
        
        chromosome1 = Chromosome(id="chr1", seq_index=seq_index1)
        chromosome2 = Chromosome(id="chr1", seq_index=seq_index2)
        chromosome3 = Chromosome(id="chr2", seq_index=seq_index3)  # Different ID
        
        assert chromosome1 == chromosome2  # Same ID and locus
        assert chromosome1 != chromosome3  # Different ID
        assert hash(chromosome1) == hash(chromosome2)
        assert hash(chromosome1) != hash(chromosome3)
