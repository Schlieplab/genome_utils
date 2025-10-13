#!/usr/bin/env python
"""
Filename: tests/unit/test_genome_builder.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.0
Description: Unit tests for the GenomeBuilder class.
License: LGPL-3.0-or-later
"""

import gzip
from unittest.mock import Mock, patch

import pytest

from GenomeUtils.genome.builder import BuilderStateError, _strip_version
from GenomeUtils.Genome import Genome, GenomeBuilder


class TestGenomeBuilder:
    """Test cases for the GenomeBuilder class."""

    def test_genome_builder_initialization(self):
        """Test basic GenomeBuilder initialization."""
        builder = GenomeBuilder(
            id="test_genome",
            species="Homo sapiens",
            name="Test Genome",
            assembly="Test_v1"
        )
        
        assert builder._genome.id == "test_genome"
        assert builder._genome.species == "Homo sapiens"
        assert builder._genome.name == "Test Genome"
        assert builder._genome.assembly == "Test_v1"
        assert builder._separate_scaffolds == True  # Default
        assert builder._scaffold_genome is not None
        assert builder._scaffold_genome.id == "test_genome_scaffolds"
        assert builder._scaffold_genome.name == "Test Genome (Scaffolds)"

    def test_genome_builder_initialization_no_scaffolds(self):
        """Test GenomeBuilder initialization without scaffold separation."""
        builder = GenomeBuilder(
            id="test_genome",
            species="Homo sapiens", 
            name="Test Genome",
            separate_scaffolds=False
        )
        
        assert builder._separate_scaffolds == False
        assert builder._scaffold_genome is None

    def test_genome_builder_initialization_custom_chromosomes(self):
        """Test GenomeBuilder initialization with custom main chromosomes."""
        custom_chroms = ["chr1", "chr2", "chrX"]
        builder = GenomeBuilder(
            id="test_genome",
            species="Homo sapiens",
            name="Test Genome",
            main_chromosomes=custom_chroms
        )
        
        assert builder._main_chromosomes == set(custom_chroms)

    def test_genome_builder_default_chromosomes(self):
        """Test that default main chromosomes include standard human chromosomes."""
        builder = GenomeBuilder(
            id="test_genome",
            species="Homo sapiens",
            name="Test Genome"
        )
        
        # Should include 1-22, X, Y, M, MT and chr-prefixed versions
        expected_chroms = {str(i) for i in range(1, 23)} | {'X', 'Y', 'M', 'MT'}
        expected_chroms.update({f'chr{c}' for c in expected_chroms})
        
        assert builder._main_chromosomes == expected_chroms

    def test_set_chromosome_filter_success(self):
        """Test setting chromosome filter before DNA loading."""
        builder = GenomeBuilder("test", "Homo sapiens", "name")
        result = builder.set_chromosome_filter(["chr1", "chr2"])
        
        assert result is builder  # Test method chaining
        assert builder._chromosome_filter == {"chr1", "chr2"}

    def test_set_chromosome_filter_after_dna_fails(self, temp_dir, sample_fasta_content):
        """Test that setting chromosome filter after DNA loading raises error."""
        # Create test FASTA file
        fasta_file = temp_dir / "test.fa"
        fasta_file.write_text(sample_fasta_content)
        
        builder = GenomeBuilder("test", "Homo sapiens", "name")
        builder.with_dna_fasta(fasta_file)
        
        with pytest.raises(BuilderStateError, match="Cannot set chromosome filter after with_dna_fasta"):
            builder.set_chromosome_filter(["chr1"])

    def test_with_dna_fasta_success(self, temp_dir, sample_fasta_content):
        """Test successful DNA FASTA loading."""
        fasta_file = temp_dir / "test.fa"
        fasta_file.write_text(sample_fasta_content)
        
        builder = GenomeBuilder("test", "Homo sapiens", "name", separate_scaffolds=False)
        result = builder.with_dna_fasta(fasta_file)
        
        assert result is builder  # Test method chaining
        assert len(builder._genome.chromosomes) == 2
        assert builder._genome.chromosome_by_id("chr1") is not None
        assert builder._genome.chromosome_by_id("chr2") is not None

    def test_with_dna_fasta_gzipped(self, temp_dir, sample_fasta_content):
        """Test DNA FASTA loading from gzipped file."""
        fasta_file = temp_dir / "test.fa.gz"
        with gzip.open(fasta_file, 'wt') as f:
            f.write(sample_fasta_content)
        
        builder = GenomeBuilder("test", "Homo sapiens", "name", separate_scaffolds=False)
        builder.with_dna_fasta(fasta_file)
        
        assert len(builder._genome.chromosomes) == 2
        # Check that extracted file was created
        extracted_file = temp_dir / "test.fa"
        assert extracted_file.exists()

    def test_with_dna_fasta_chromosome_filter(self, temp_dir, sample_fasta_content):
        """Test DNA FASTA loading with chromosome filter."""
        fasta_file = temp_dir / "test.fa"
        fasta_file.write_text(sample_fasta_content)
        
        builder = GenomeBuilder("test", "Homo sapiens", "name", separate_scaffolds=False)
        builder.set_chromosome_filter(["chr1"])
        builder.with_dna_fasta(fasta_file)
        
        assert len(builder._genome.chromosomes) == 1
        assert builder._genome.chromosome_by_id("chr1") is not None
        # chr2 should not exist, so this should raise ValueError
        with pytest.raises(ValueError):
            builder._genome.chromosome_by_id("chr2")

    def test_with_dna_fasta_scaffold_separation(self, temp_dir):
        """Test DNA FASTA loading with scaffold separation."""
        # Create FASTA with main chromosomes and scaffolds
        fasta_content = """>chr1
ATCGATCGATCGAAATTTGGGCCC
>scaffold_123
GGGGCCCCAAAATTTTGGGGCCCC
>chrM
TTTTAAAACCCCGGGGTTTTAAAA
"""
        fasta_file = temp_dir / "test.fa"
        fasta_file.write_text(fasta_content)
        
        builder = GenomeBuilder("test", "Homo sapiens", "name", separate_scaffolds=True)
        builder.with_dna_fasta(fasta_file)
        
        # chr1 and chrM should be in main genome
        assert len(builder._genome.chromosomes) == 2
        assert builder._genome.chromosome_by_id("chr1") is not None
        assert builder._genome.chromosome_by_id("chrM") is not None
        
        # scaffold_123 should be in scaffold genome
        assert len(builder._scaffold_genome.chromosomes) == 1
        assert builder._scaffold_genome.chromosome_by_id("scaffold_123") is not None

    def test_with_dna_fasta_already_called(self, temp_dir, sample_fasta_content):
        """Test that calling with_dna_fasta twice raises error."""
        fasta_file = temp_dir / "test.fa"
        fasta_file.write_text(sample_fasta_content)
        
        builder = GenomeBuilder("test", "Homo sapiens", "name")
        builder.with_dna_fasta(fasta_file)
        
        with pytest.raises(BuilderStateError, match="with_dna_fasta\\(\\) has already been called"):
            builder.with_dna_fasta(fasta_file)

    def test_with_cdna_fasta_success(self, temp_dir, sample_cdna_content):
        """Test successful cDNA FASTA loading."""
        cdna_file = temp_dir / "cdna.fa"
        cdna_file.write_text(sample_cdna_content)
        
        builder = GenomeBuilder("test", "Homo sapiens", "name")
        result = builder.with_cdna_fasta(cdna_file)
        
        assert result is builder  # Test method chaining
        assert len(builder._cdna_records) == 2
        assert "TRANS001" in builder._cdna_records
        assert "TRANS002" in builder._cdna_records

    def test_with_cdna_fasta_gzipped(self, temp_dir, sample_cdna_content):
        """Test cDNA FASTA loading from gzipped file."""
        cdna_file = temp_dir / "cdna.fa.gz"
        with gzip.open(cdna_file, 'wt') as f:
            f.write(sample_cdna_content)
        
        builder = GenomeBuilder("test", "Homo sapiens", "name")
        builder.with_cdna_fasta(cdna_file)
        
        assert len(builder._cdna_records) == 2

    def test_with_cdna_fasta_already_called(self, temp_dir, sample_cdna_content):
        """Test that calling with_cdna_fasta twice raises error."""
        cdna_file = temp_dir / "cdna.fa"
        cdna_file.write_text(sample_cdna_content)
        
        builder = GenomeBuilder("test", "Homo sapiens", "name")
        builder.with_cdna_fasta(cdna_file)
        
        with pytest.raises(BuilderStateError, match="with_cdna_fasta\\(\\) has already been called"):
            builder.with_cdna_fasta(cdna_file)

    def test_with_gtf_file_prerequisites(self, temp_dir, sample_gtf_content):
        """Test GTF file loading prerequisites."""
        gtf_file = temp_dir / "test.gtf"
        gtf_file.write_text(sample_gtf_content)
        
        builder = GenomeBuilder("test", "Homo sapiens", "name")
        
        # Should fail without DNA FASTA
        with pytest.raises(BuilderStateError, match="Must call with_dna_fasta\\(\\) before with_gtf_file"):
            builder.with_gtf_file(gtf_file)

    def test_with_gtf_file_prerequisites_no_cdna(self, temp_dir, sample_fasta_content, sample_gtf_content):
        """Test GTF file loading without cDNA."""
        fasta_file = temp_dir / "test.fa"
        fasta_file.write_text(sample_fasta_content)
        gtf_file = temp_dir / "test.gtf"
        gtf_file.write_text(sample_gtf_content)
        
        builder = GenomeBuilder("test", "Homo sapiens", "name")
        builder.with_dna_fasta(fasta_file)
        
        # Should fail without cDNA FASTA
        with pytest.raises(BuilderStateError, match="Must call with_cdna_fasta\\(\\) before with_gtf_file"):
            builder.with_gtf_file(gtf_file)

    @patch('GenomeUtils.genome.builder.gffutils')
    def test_with_gtf_file_success(self, mock_gffutils, temp_dir, sample_fasta_content, sample_cdna_content, sample_gtf_content):
        """Test successful GTF file processing."""
        # Setup files
        fasta_file = temp_dir / "test.fa"
        fasta_file.write_text(sample_fasta_content)
        cdna_file = temp_dir / "cdna.fa"
        cdna_file.write_text(sample_cdna_content)
        gtf_file = temp_dir / "test.gtf"
        gtf_file.write_text(sample_gtf_content)
        
        # Mock gffutils database
        mock_db = Mock()
        
        # Create proper mock results for database queries
        def mock_execute(query):
            mock_result = Mock()
            if "count(*)" in query:
                mock_result.fetchone.return_value = (1,)
                return mock_result
            elif "featuretype = 'gene'" in query:
                return [("GENE001", "chr1", 100, 400, "+", '{"gene_id": ["GENE001"], "gene_name": ["GENE_1"]}')]
            elif "featuretype = 'transcript'" in query:
                return [("TRANS001", 120, 380, "+", '{"gene_id": ["GENE001"], "transcript_id": ["TRANS001"]}')]
            elif "featuretype = 'exon'" in query:
                return [("EXON001", "chr1", 120, 200, "+", '{"transcript_id": ["TRANS001"], "exon_id": ["EXON001"]}')]
            return []
        
        mock_db.conn.execute.side_effect = mock_execute
        mock_gffutils.FeatureDB.return_value = mock_db
        mock_gffutils.create_db.return_value = mock_db
        
        builder = GenomeBuilder("test", "Homo sapiens", "name", separate_scaffolds=False)
        result = (builder
                 .with_dna_fasta(fasta_file)
                 .with_cdna_fasta(cdna_file)
                 .with_gtf_file(gtf_file))
        
        assert result is builder  # Test method chaining
        assert len(builder._genes_map) == 1
        assert len(builder._transcripts_map) == 1

    def test_with_gtf_file_already_called(self, temp_dir, sample_fasta_content, sample_cdna_content, sample_gtf_content):
        """Test that calling with_gtf_file twice raises error."""
        # Setup files
        fasta_file = temp_dir / "test.fa"
        fasta_file.write_text(sample_fasta_content)
        cdna_file = temp_dir / "cdna.fa"
        cdna_file.write_text(sample_cdna_content)
        gtf_file = temp_dir / "test.gtf"
        gtf_file.write_text(sample_gtf_content)
        
        builder = GenomeBuilder("test", "Homo sapiens", "name")
        builder.with_dna_fasta(fasta_file)
        builder.with_cdna_fasta(cdna_file)
        
        # Mock the GTF processing to set _genes_map
        builder._genes_map = {"GENE001": Mock()}
        
        with pytest.raises(BuilderStateError, match="with_gtf_file\\(\\) has already been called"):
            builder.with_gtf_file(gtf_file)

    def test_build_without_gtf(self):
        """Test that build fails without GTF data."""
        builder = GenomeBuilder("test", "Homo sapiens", "name")
        
        with pytest.raises(BuilderStateError, match="Cannot build Genome.*GTF data is missing"):
            builder.build()

    @patch('GenomeUtils.genome.builder.gffutils')
    def test_build_success_single_genome(self, mock_gffutils, temp_dir, sample_fasta_content, sample_cdna_content, sample_gtf_content):
        """Test successful build returning single genome."""
        # Setup files and mocks
        fasta_file = temp_dir / "test.fa"
        fasta_file.write_text(sample_fasta_content)
        cdna_file = temp_dir / "cdna.fa" 
        cdna_file.write_text(sample_cdna_content)
        gtf_file = temp_dir / "test.gtf"
        gtf_file.write_text(sample_gtf_content)
        
        # Mock gffutils database
        mock_db = Mock()
        
        # Create proper mock results for database queries
        def mock_execute(query):
            mock_result = Mock()
            if "count(*)" in query:
                mock_result.fetchone.return_value = (0,)
                return mock_result
            return []
        
        mock_db.conn.execute.side_effect = mock_execute
        mock_gffutils.FeatureDB.return_value = mock_db
        mock_gffutils.create_db.return_value = mock_db
        
        builder = GenomeBuilder("test", "Homo sapiens", "name", separate_scaffolds=False)
        builder.with_dna_fasta(fasta_file)
        builder.with_cdna_fasta(cdna_file)
        builder.with_gtf_file(gtf_file)
        
        # Mock at least one gene to satisfy build requirements
        builder._genes_map = {"GENE001": Mock()}
        
        genome = builder.build()
        
        assert isinstance(genome, Genome)
        assert genome.id == "test"
        assert genome.is_indexed == True

    @patch('GenomeUtils.genome.builder.gffutils')
    def test_build_success_with_scaffolds(self, mock_gffutils, temp_dir, sample_fasta_content, sample_cdna_content, sample_gtf_content):
        """Test successful build returning tuple with scaffolds."""
        # Setup files and mocks
        fasta_file = temp_dir / "test.fa"
        fasta_file.write_text(sample_fasta_content)
        cdna_file = temp_dir / "cdna.fa"
        cdna_file.write_text(sample_cdna_content)
        gtf_file = temp_dir / "test.gtf"
        gtf_file.write_text(sample_gtf_content)
        
        # Mock gffutils database
        mock_db = Mock()
        
        # Create proper mock results for database queries
        def mock_execute(query):
            mock_result = Mock()
            if "count(*)" in query:
                mock_result.fetchone.return_value = (0,)
                return mock_result
            return []
        
        mock_db.conn.execute.side_effect = mock_execute
        mock_gffutils.FeatureDB.return_value = mock_db
        mock_gffutils.create_db.return_value = mock_db
        
        builder = GenomeBuilder("test", "Homo sapiens", "name", separate_scaffolds=True)
        builder.with_dna_fasta(fasta_file)
        builder.with_cdna_fasta(cdna_file)
        builder.with_gtf_file(gtf_file)
        
        # Mock at least one gene to satisfy build requirements
        builder._genes_map = {"GENE001": Mock()}
        
        result = builder.build()
        
        assert isinstance(result, tuple)
        main_genome, scaffold_genome = result
        assert isinstance(main_genome, Genome)
        assert isinstance(scaffold_genome, Genome)
        assert main_genome.id == "test"
        assert scaffold_genome.id == "test_scaffolds"

    def test_method_chaining(self, temp_dir, sample_fasta_content, sample_cdna_content):
        """Test that all methods support method chaining."""
        fasta_file = temp_dir / "test.fa"
        fasta_file.write_text(sample_fasta_content)
        cdna_file = temp_dir / "cdna.fa"
        cdna_file.write_text(sample_cdna_content)
        
        builder = GenomeBuilder("test", "Homo sapiens", "name")
        
        # Test chaining
        chained_builder = (builder
                          .set_chromosome_filter(["chr1", "chr2"])
                          .with_dna_fasta(fasta_file)
                          .with_cdna_fasta(cdna_file))
        
        assert chained_builder is builder

    def test_memory_offload(self):
        """Test that memory is properly offloaded after build."""
        builder = GenomeBuilder("test", "Homo sapiens", "name")
        
        # Setup some data
        builder._cdna_records = {"test": Mock()}
        builder._genes_map = {"test": Mock()}
        builder._transcripts_map = {"test": Mock()}
        
        builder._offload_memory()
        
        assert len(builder._cdna_records) == 0
        assert len(builder._genes_map) == 0
        assert len(builder._transcripts_map) == 0

    def test_strip_version_function(self):
        """Test the _strip_version utility function."""
        assert _strip_version("NC_000001.11") == "NC_000001"
        assert _strip_version("ENST00000123456.2") == "ENST00000123456"
        assert _strip_version("no_version") == "no_version"
        assert _strip_version("multiple.dots.here") == "multiple"

    def test_builder_state_error(self):
        """Test BuilderStateError custom exception."""
        error = BuilderStateError("Test error message")
        assert str(error) == "Test error message"
        assert isinstance(error, Exception)


class TestGenomeBuilderIntegration:
    """Integration tests for GenomeBuilder with real file processing."""

    @pytest.mark.integration
    def test_full_builder_workflow_minimal(self, temp_dir):
        """Test complete builder workflow with minimal data."""
        # Create minimal test files
        fasta_content = """>chr1
ATCGATCGATCGAAATTTGGGCCCTTTTAAAACCCCGGGGTTTTAAAACCCCGGGG
"""
        cdna_content = """>TRANS001
ATCGATCGATCGAAATTTGGGCCCTTTTAAAACCCCGGGG
"""
        gtf_content = """chr1	test	gene	10	50	.	+	.	gene_id "GENE001"; gene_name "TEST_GENE";
chr1	test	transcript	15	45	.	+	.	gene_id "GENE001"; transcript_id "TRANS001";
chr1	test	exon	15	25	.	+	.	gene_id "GENE001"; transcript_id "TRANS001"; exon_id "EXON001";
chr1	test	exon	35	45	.	+	.	gene_id "GENE001"; transcript_id "TRANS001"; exon_id "EXON002";
"""
        
        fasta_file = temp_dir / "test.fa"
        fasta_file.write_text(fasta_content)
        cdna_file = temp_dir / "cdna.fa"
        cdna_file.write_text(cdna_content)
        gtf_file = temp_dir / "test.gtf"
        gtf_file.write_text(gtf_content)
        
        # Build genome
        builder = GenomeBuilder(
            id="test_genome",
            species="Homo sapiens",
            name="Test Genome",
            separate_scaffolds=False
        )
        
        genome = (builder
                 .set_chromosome_filter(["chr1"])
                 .with_dna_fasta(fasta_file)
                 .with_cdna_fasta(cdna_file)
                 .with_gtf_file(gtf_file)
                 .build())
        
        # Verify genome structure
        assert isinstance(genome, Genome)
        assert len(genome.chromosomes) == 1
        assert genome.chromosome_by_id("chr1") is not None
        
        chromosome = genome.chromosome_by_id("chr1")
        assert len(chromosome.genes) == 1
        
        gene = chromosome.genes[0]
        assert gene.id == "GENE001"
        assert gene.name == "TEST_GENE"
        
        assert len(gene.transcripts) == 1
        transcript = gene.transcripts[0]
        assert transcript.id == "TRANS001"
        
        assert len(transcript.exons) == 2
        assert transcript.exons[0].id == "EXON001"
        assert transcript.exons[1].id == "EXON002"
        
        # Verify exon-transcript relationship
        exon1 = transcript.exons[0]
        exon2 = transcript.exons[1]
        assert transcript in exon1.get_transcripts()
        assert transcript in exon2.get_transcripts()
        
        # Verify exon-gene relationship
        assert exon1.get_gene() == gene
        assert exon2.get_gene() == gene

    @pytest.mark.integration
    def test_exon_shared_across_transcripts(self, temp_dir):
        """Test that exons can be shared across multiple transcripts (alternative splicing)."""
        # Create test files with alternative transcripts sharing exons
        fasta_content = """>chr1
ATCGATCGATCGAAATTTGGGCCCTTTTAAAACCCCGGGGTTTTAAAACCCCGGGG
"""
        cdna_content = """>TRANS001
ATCGATCGATCGAAATTTGGGCCCTTTTAAAACCCCGGGG
>TRANS002
ATCGATCGATCGAAATTTAAAACCCCGGGG
"""
        # Two transcripts of the same gene sharing EXON001
        gtf_content = """chr1	test	gene	10	50	.	+	.	gene_id "GENE001"; gene_name "TEST_GENE";
chr1	test	transcript	15	45	.	+	.	gene_id "GENE001"; transcript_id "TRANS001";
chr1	test	exon	15	25	.	+	.	gene_id "GENE001"; transcript_id "TRANS001"; exon_id "EXON001";
chr1	test	exon	35	45	.	+	.	gene_id "GENE001"; transcript_id "TRANS001"; exon_id "EXON002";
chr1	test	transcript	15	40	.	+	.	gene_id "GENE001"; transcript_id "TRANS002";
chr1	test	exon	15	25	.	+	.	gene_id "GENE001"; transcript_id "TRANS002"; exon_id "EXON001";
chr1	test	exon	30	40	.	+	.	gene_id "GENE001"; transcript_id "TRANS002"; exon_id "EXON003";
"""
        
        fasta_file = temp_dir / "test.fa"
        fasta_file.write_text(fasta_content)
        cdna_file = temp_dir / "cdna.fa"
        cdna_file.write_text(cdna_content)
        gtf_file = temp_dir / "test.gtf"
        gtf_file.write_text(gtf_content)
        
        # Build genome
        builder = GenomeBuilder(
            id="test_genome",
            species="Homo sapiens",
            name="Test Genome",
            separate_scaffolds=False
        )
        
        genome = (builder
                 .set_chromosome_filter(["chr1"])
                 .with_dna_fasta(fasta_file)
                 .with_cdna_fasta(cdna_file)
                 .with_gtf_file(gtf_file)
                 .build())
        
        # Get the gene and transcripts
        gene = genome.gene_by_id("GENE001")
        assert len(gene.transcripts) == 2
        
        transcript1 = genome.transcript_by_id("TRANS001")
        transcript2 = genome.transcript_by_id("TRANS002")
        
        # Get EXON001 which should be shared
        exon001 = genome.exon_by_id("EXON001")
        
        # EXON001 should be in both transcripts
        assert exon001 in transcript1.exons
        assert exon001 in transcript2.exons
        
        # EXON001 should have both transcripts in its transcripts list
        assert len(exon001.get_transcripts()) == 2
        assert transcript1 in exon001.get_transcripts()
        assert transcript2 in exon001.get_transcripts()
        
        # Both transcripts should have the same gene
        assert exon001.get_gene() == gene
        
        # Verify other exons are specific to their transcripts
        exon002 = genome.exon_by_id("EXON002")
        exon003 = genome.exon_by_id("EXON003")
        
        assert len(exon002.get_transcripts()) == 1
        assert transcript1 in exon002.get_transcripts()
        
        assert len(exon003.get_transcripts()) == 1
        assert transcript2 in exon003.get_transcripts()