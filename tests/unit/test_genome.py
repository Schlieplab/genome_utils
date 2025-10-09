#!/usr/bin/env python
"""
Filename: test_genome.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.0
Description: Unit tests for the Genome class.
License: LGPL-3.0-or-later
"""

from unittest.mock import Mock

import pytest
from Bio.Seq import Seq

from GenomeUtils.Genome import Chromosome, Exon, Gene, Genome, Locus, Transcript


class TestGenome:
    """Test cases for the Genome class."""

    def test_genome_creation(self):
        """Test basic genome creation."""
        genome = Genome(
            id="test_genome",
            species="Homo sapiens",
            name="Test Genome"
        )
        
        assert genome.id == "test_genome"
        assert genome.species == "Homo sapiens"
        assert genome.name == "Test Genome"
        assert genome._chromosomes == {}
        assert genome._genes == {}
        assert genome._transcripts == {}
        assert genome._exons == {}
        assert genome.is_indexed == False

    def test_genome_creation_with_kwargs(self):
        """Test genome creation with additional attributes."""
        genome = Genome(
            id="human_genome",
            species="Homo sapiens",
            name="Human Reference Genome",
            assembly="GRCh38",
            version="p14",
            release=110
        )
        
        assert genome.assembly == "GRCh38"
        assert genome.version == "p14"
        assert genome.release == 110

    def test_genome_repr(self):
        """Test genome string representation."""
        genome = Genome(
            id="test_genome",
            species="Homo sapiens",
            name="Test Genome"
        )
        
        expected = "Genome(id='test_genome', species='Homo sapiens', name='Test Genome')"
        assert repr(genome) == expected

    def test_genome_add_chromosome(self):
        """Test adding chromosomes to genome."""
        genome = Genome("test", "Homo sapiens", "name")
        
        # Create mock chromosomes
        chr1 = Mock(spec=Chromosome)
        chr1.id = "chr1"
        chr2 = Mock(spec=Chromosome)
        chr2.id = "chr2"
        
        genome.add_chromosome(chr1)
        genome.add_chromosome(chr2)
        
        assert len(genome._chromosomes) == 2
        assert genome._chromosomes["chr1"] == chr1
        assert genome._chromosomes["chr2"] == chr2
        assert genome.is_indexed == False

    def test_genome_add_chromosome_duplicate_id(self):
        """Test error when adding chromosome with duplicate ID."""
        genome = Genome("test", "Homo sapiens", "name")
        
        chr1 = Mock(spec=Chromosome)
        chr1.id = "chr1"
        chr_duplicate = Mock(spec=Chromosome)
        chr_duplicate.id = "chr1"
        
        genome.add_chromosome(chr1)
        
        with pytest.raises(ValueError, match="Chromosome with ID 'chr1' already exists"):
            genome.add_chromosome(chr_duplicate)

    def test_genome_chromosomes_property(self):
        """Test chromosomes property returns list of chromosomes."""
        genome = Genome("test", "Homo sapiens", "name")
        
        # Initially empty
        assert genome.chromosomes == []
        
        # Add chromosomes
        chr1 = Mock(spec=Chromosome)
        chr1.id = "chr1"
        chr2 = Mock(spec=Chromosome)
        chr2.id = "chr2"
        
        genome.add_chromosome(chr1)
        genome.add_chromosome(chr2)
        
        chromosomes = genome.chromosomes
        assert len(chromosomes) == 2
        assert chr1 in chromosomes
        assert chr2 in chromosomes

    def test_genome_index_simple(self):
        """Test genome indexing with simple structure."""
        genome = Genome("test", "Homo sapiens", "name")
        
        # Create chromosome
        chromosome = Mock(spec=Chromosome)
        chromosome.id = "chr1"
        
        # Create gene
        gene = Mock(spec=Gene)
        gene.id = "GENE001"
        
        # Create transcript
        transcript = Mock(spec=Transcript)
        transcript.id = "TRANS001"
        
        # Create exon
        exon = Mock(spec=Exon)
        exon.id = "EXON001"
        
        # Set up relationships
        chromosome.genes = [gene]
        gene.transcripts = [transcript]
        transcript.exons = [exon]
        
        genome.add_chromosome(chromosome)
        genome.index()
        
        assert genome.is_indexed == True
        assert genome._genes["GENE001"] == gene
        assert genome._transcripts["TRANS001"] == transcript
        assert genome._exons["EXON001"] == exon

    def test_genome_index_complex(self):
        """Test genome indexing with complex structure."""
        genome = Genome("test", "Homo sapiens", "name")
        
        # Create multiple chromosomes
        chr1 = Mock(spec=Chromosome)
        chr1.id = "chr1"
        chr2 = Mock(spec=Chromosome)
        chr2.id = "chr2"
        
        # Create genes
        gene1 = Mock(spec=Gene)
        gene1.id = "GENE001"
        gene2 = Mock(spec=Gene)
        gene2.id = "GENE002"
        gene3 = Mock(spec=Gene)
        gene3.id = "GENE003"
        
        # Create transcripts
        transcript1 = Mock(spec=Transcript)
        transcript1.id = "TRANS001"
        transcript2 = Mock(spec=Transcript)
        transcript2.id = "TRANS002"
        transcript3 = Mock(spec=Transcript)
        transcript3.id = "TRANS003"
        transcript4 = Mock(spec=Transcript)
        transcript4.id = "TRANS004"
        
        # Create exons
        exon1 = Mock(spec=Exon)
        exon1.id = "EXON001"
        exon2 = Mock(spec=Exon)
        exon2.id = "EXON002"
        exon3 = Mock(spec=Exon)
        exon3.id = "EXON003"
        
        # Set up relationships
        chr1.genes = [gene1, gene2]
        chr2.genes = [gene3]
        
        gene1.transcripts = [transcript1, transcript2]
        gene2.transcripts = [transcript3]
        gene3.transcripts = [transcript4]
        
        transcript1.exons = [exon1]
        transcript2.exons = [exon2]
        transcript3.exons = [exon3]
        transcript4.exons = []  # No exons
        
        genome.add_chromosome(chr1)
        genome.add_chromosome(chr2)
        genome.index()
        
        # Check indexing
        assert len(genome._genes) == 3
        assert len(genome._transcripts) == 4
        assert len(genome._exons) == 3
        
        assert genome._genes["GENE001"] == gene1
        assert genome._transcripts["TRANS003"] == transcript3
        assert genome._exons["EXON002"] == exon2

    def test_genome_chromosome_by_id(self):
        """Test getting chromosome by ID."""
        genome = Genome("test", "Homo sapiens", "name")
        
        chromosome = Mock(spec=Chromosome)
        chromosome.id = "chr1"
        genome.add_chromosome(chromosome)
        
        # Should work before indexing
        result = genome.chromosome_by_id("chr1")
        assert result == chromosome
        
        # Should raise error for non-existent chromosome
        with pytest.raises(ValueError, match="Chromosome with ID 'chr_nonexistent' not found"):
            genome.chromosome_by_id("chr_nonexistent")

    def test_genome_gene_by_id_before_indexing(self):
        """Test error when getting gene by ID before indexing."""
        genome = Genome("test", "Homo sapiens", "name")
        
        with pytest.raises(RuntimeError, match="The genome is not indexed"):
            genome.gene_by_id("GENE001")

    def test_genome_gene_by_id_after_indexing(self):
        """Test getting gene by ID after indexing."""
        genome = Genome("test", "Homo sapiens", "name")
        
        # Set up genome structure
        chromosome = Mock(spec=Chromosome)
        chromosome.id = "chr1"
        
        gene = Mock(spec=Gene)
        gene.id = "GENE001"
        
        chromosome.genes = [gene]
        gene.transcripts = []
        
        genome.add_chromosome(chromosome)
        genome.index()
        
        # Should work after indexing
        result = genome.gene_by_id("GENE001")
        assert result == gene
        
        # Should raise error for non-existent gene
        with pytest.raises(ValueError, match="Gene with ID 'GENE_NONEXISTENT' not found"):
            genome.gene_by_id("GENE_NONEXISTENT")

    def test_genome_transcript_by_id(self):
        """Test getting transcript by ID."""
        genome = Genome("test", "Homo sapiens", "name")
        
        # Set up genome structure
        chromosome = Mock(spec=Chromosome)
        chromosome.id = "chr1"
        
        gene = Mock(spec=Gene)
        gene.id = "GENE001"
        
        transcript = Mock(spec=Transcript)
        transcript.id = "TRANS001"
        
        chromosome.genes = [gene]
        gene.transcripts = [transcript]
        transcript.exons = []
        
        genome.add_chromosome(chromosome)
        
        # Should fail before indexing
        with pytest.raises(RuntimeError, match="The genome is not indexed"):
            genome.transcript_by_id("TRANS001")
        
        genome.index()
        
        # Should work after indexing
        result = genome.transcript_by_id("TRANS001")
        assert result == transcript
        
        # Should raise error for non-existent transcript
        with pytest.raises(ValueError, match="Transcript with ID 'TRANS_NONEXISTENT' not found"):
            genome.transcript_by_id("TRANS_NONEXISTENT")

    def test_genome_exon_by_id(self):
        """Test getting exon by ID."""
        genome = Genome("test", "Homo sapiens", "name")
        
        # Set up genome structure
        chromosome = Mock(spec=Chromosome)
        chromosome.id = "chr1"
        
        gene = Mock(spec=Gene)
        gene.id = "GENE001"
        
        transcript = Mock(spec=Transcript)
        transcript.id = "TRANS001"
        
        exon = Mock(spec=Exon)
        exon.id = "EXON001"
        
        chromosome.genes = [gene]
        gene.transcripts = [transcript]
        transcript.exons = [exon]
        
        genome.add_chromosome(chromosome)
        
        # Should fail before indexing
        with pytest.raises(RuntimeError, match="The genome is not indexed"):
            genome.exon_by_id("EXON001")
        
        genome.index()
        
        # Should work after indexing
        result = genome.exon_by_id("EXON001")
        assert result == exon
        
        # Should raise error for non-existent exon
        with pytest.raises(ValueError, match="Exon with ID 'EXON_NONEXISTENT' not found"):
            genome.exon_by_id("EXON_NONEXISTENT")

    def test_genome_properties_empty(self):
        """Test genome properties when empty."""
        genome = Genome("test", "Homo sapiens", "name")
        
        assert genome.genes == []
        assert genome.transcripts == []
        assert genome.exons == []

    def test_genome_properties_after_indexing(self):
        """Test genome properties after indexing."""
        genome = Genome("test", "Homo sapiens", "name")
        
        # Set up structure
        chromosome = Mock(spec=Chromosome)
        chromosome.id = "chr1"
        
        gene1 = Mock(spec=Gene)
        gene1.id = "GENE001"
        gene2 = Mock(spec=Gene)
        gene2.id = "GENE002"
        
        transcript1 = Mock(spec=Transcript)
        transcript1.id = "TRANS001"
        transcript2 = Mock(spec=Transcript)
        transcript2.id = "TRANS002"
        
        exon1 = Mock(spec=Exon)
        exon1.id = "EXON001"
        
        chromosome.genes = [gene1, gene2]
        gene1.transcripts = [transcript1]
        gene2.transcripts = [transcript2]
        transcript1.exons = [exon1]
        transcript2.exons = []
        
        genome.add_chromosome(chromosome)
        genome.index()
        
        genes = genome.genes
        transcripts = genome.transcripts
        exons = genome.exons
        
        assert len(genes) == 2
        assert len(transcripts) == 2
        assert len(exons) == 1
        
        assert gene1 in genes
        assert gene2 in genes
        assert transcript1 in transcripts
        assert transcript2 in transcripts
        assert exon1 in exons

    def test_genome_get_sequence_by_locus_before_indexing(self):
        """Test error when getting sequence by locus before indexing."""
        genome = Genome("test", "Homo sapiens", "name")
        locus = Locus("chr1", 100, 200)
        
        with pytest.raises(RuntimeError, match="The genome is not indexed"):
            genome.get_sequence_by_locus(locus)

    def test_genome_get_sequence_by_locus_after_indexing(self):
        """Test getting sequence by locus after indexing."""
        genome = Genome("test", "Homo sapiens", "name")
        
        # Create chromosome with mock sequence method
        chromosome = Mock(spec=Chromosome)
        chromosome.id = "chr1"
        test_sequence = Seq("ATCGATCGATCG")
        chromosome.get_subsequence_by_locus.return_value = test_sequence
        chromosome.genes = []
        
        genome.add_chromosome(chromosome)
        genome.index()
        
        locus = Locus("chr1", 100, 200)
        result = genome.get_sequence_by_locus(locus)
        
        chromosome.get_subsequence_by_locus.assert_called_once_with(locus)
        assert result == test_sequence

    def test_genome_indexing_resets_flag(self):
        """Test that adding elements resets the indexing flag."""
        genome = Genome("test", "Homo sapiens", "name")
        
        chromosome = Mock(spec=Chromosome)
        chromosome.id = "chr1"
        chromosome.genes = []
        
        genome.add_chromosome(chromosome)
        genome.index()
        assert genome.is_indexed == True
        
        # Adding another chromosome should reset index flag
        chromosome2 = Mock(spec=Chromosome)
        chromosome2.id = "chr2"
        genome.add_chromosome(chromosome2)
        assert genome.is_indexed == False

    def test_genome_multiple_indexing(self):
        """Test that genome can be indexed multiple times."""
        genome = Genome("test", "Homo sapiens", "name")
        
        chromosome = Mock(spec=Chromosome)
        chromosome.id = "chr1"
        
        gene = Mock(spec=Gene)
        gene.id = "GENE001"
        
        chromosome.genes = [gene]
        gene.transcripts = []
        
        genome.add_chromosome(chromosome)
        
        # Index first time
        genome.index()
        assert genome.is_indexed == True
        assert len(genome._genes) == 1
        
        # Index second time (should work)
        genome.index()
        assert genome.is_indexed == True
        assert len(genome._genes) == 1

    def test_genome_empty_indexing(self):
        """Test indexing genome with no chromosomes."""
        genome = Genome("test", "Homo sapiens", "name")
        
        genome.index()
        
        assert genome.is_indexed == True
        assert len(genome._genes) == 0
        assert len(genome._transcripts) == 0
        assert len(genome._exons) == 0




