#!/usr/bin/env python
"""
Filename: tests/conftest.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.3
Description: Configuration and fixtures for genome_utils tests.
License: LGPL-3.0-or-later
"""

import pytest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from GenomeUtils.Genome import Locus, Chromosome, Transcript, Exon
from GenomeUtils.Genome import Gene, Genome




@pytest.fixture
def temp_dir():
    """Provide a temporary directory that gets cleaned up after the test."""
    temp_dir = Path(tempfile.mkdtemp())
    yield temp_dir
    shutil.rmtree(temp_dir)


@pytest.fixture
def sample_locus():
    """Provide a sample Locus object for testing."""
    return Locus("chr1", 1000, 2000, "+")


@pytest.fixture
def sample_sequence():
    """Provide a sample DNA sequence for testing."""
    return Seq("ATCGATCGATCGAAATTTGGGCCCTTTTAAAACCCCGGGGTTTTAAAACCCCGGGG")


@pytest.fixture
def mock_seq_index():
    """Provide a mock SeqIO index for testing."""
    def create_index(seq_dict):
        mock_index = {}
        for key, seq_str in seq_dict.items():
            mock_record = Mock()
            mock_record.seq = Seq(seq_str)
            mock_index[key] = mock_record
        return mock_index
    return create_index


@pytest.fixture
def sample_chromosome(mock_seq_index, sample_sequence):
    """Provide a sample Chromosome object for testing."""
    genome_mock = Mock()
    seq_index = mock_seq_index({"chr1": str(sample_sequence)})
    
    chromosome = Chromosome(
        id="chr1",
        seq_index=seq_index,
        genome=genome_mock
    )
    return chromosome


@pytest.fixture
def sample_gene(sample_chromosome):
    """Provide a sample Gene object for testing."""
    genome_mock = Mock()
    
    gene = Gene(
        id="ENSG00000001",
        name="TEST_GENE",
        start=1100,
        end=1500,
        strand="+",
        chromosome=sample_chromosome,
        genome=genome_mock,
        gene_type="protein_coding"
    )
    return gene


@pytest.fixture
def sample_transcript(sample_gene, sample_sequence):
    """Provide a sample Transcript object for testing."""
    genome_mock = Mock()
    
    transcript = Transcript(
        id="ENST00000001",
        start=1100,
        end=1400,
        strand="+",
        sequence=sample_sequence[:300],  # 300 bp transcript
        gene=sample_gene,
        genome=genome_mock,
        transcript_type="protein_coding"
    )
    return transcript


@pytest.fixture
def sample_exon(sample_transcript):
    """Provide a sample Exon object for testing."""
    genome_mock = Mock()
    
    exon = Exon(
        id="ENSE00000001",
        chr="chr1",
        start=1150,
        end=1250,
        strand="+",
        transcripts=[sample_transcript],
        genome=genome_mock,
        exon_number=1
    )
    return exon


@pytest.fixture
def sample_genome():
    """Provide a sample Genome object for testing."""
    return Genome(
        id="test_genome",
        species="Homo sapiens",
        name="Test Genome",
        assembly="Test_v1"
    )


@pytest.fixture
def complex_genome_structure(mock_seq_index):
    """Provide a complex genome structure for integration testing."""
    # Create genome
    genome = Genome(
        id="test_genome",
        species="Homo sapiens",
        name="Test Genome"
    )
    
    # Create sequences
    sequences = {
        "chr1": "ATCGATCGATCGAAATTTGGGCCCTTTTAAAACCCCGGGGTTTTAAAACCCCGGGG" * 10,  # 560 bp
        "chr2": "GGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCC" * 8,   # 448 bp
    }
    
    # Create chromosomes
    seq_index1 = mock_seq_index({"chr1": sequences["chr1"]})
    seq_index2 = mock_seq_index({"chr2": sequences["chr2"]})
    
    chr1 = Chromosome(id="chr1", seq_index=seq_index1, genome=genome)
    chr2 = Chromosome(id="chr2", seq_index=seq_index2, genome=genome)
    
    genome.add_chromosome(chr1)
    genome.add_chromosome(chr2)
    
    # Create genes
    gene1 = Gene(
        id="GENE001", name="GENE_1", chr="chr1", start=100, end=400, strand="+",
        chromosome=chr1, genome=genome, gene_type="protein_coding"
    )
    gene2 = Gene(
        id="GENE002", name="GENE_2", chr="chr2", start=200, end=350, strand="-",
        chromosome=chr2, genome=genome, gene_type="lncRNA"
    )
    
    chr1.add_gene(gene1)
    chr2.add_gene(gene2)
    
    # Create transcripts
    transcript1 = Transcript(
        id="TRANS001", chr="chr1", start=120, end=380, strand="+",
        sequence=Seq(sequences["chr1"][119:380]),  # 261 bp
        gene=gene1, genome=genome
    )
    transcript2 = Transcript(
        id="TRANS002", chr="chr2", start=220, end=330, strand="-",
        sequence=Seq(sequences["chr2"][219:330]),  # 111 bp
        gene=gene2, genome=genome
    )
    
    gene1.add_transcript(transcript1)
    gene2.add_transcript(transcript2)
    
    # Create exons
    exon1 = Exon(
        id="EXON001", chr="chr1", start=120, end=200, strand="+",
        transcripts=[transcript1], genome=genome
    )
    exon2 = Exon(
        id="EXON002", chr="chr1", start=300, end=380, strand="+",
        transcripts=[transcript1], genome=genome
    )
    exon3 = Exon(
        id="EXON003", chr="chr2", start=220, end=330, strand="-",
        transcripts=[transcript2], genome=genome
    )
    
    transcript1.add_exon(exon1)
    transcript1.add_exon(exon2)
    transcript2.add_exon(exon3)
    
    # Index the genome
    genome.index()
    
    return genome


@pytest.fixture
def sample_fasta_content():
    """Provide sample FASTA content for testing."""
    return """>chr1
ATCGATCGATCGAAATTTGGGCCCTTTTAAAACCCCGGGGTTTTAAAACCCCGGGG
GGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCC
>chr2
GGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCC
ATCGATCGATCGAAATTTGGGCCCTTTTAAAACCCCGGGGTTTTAAAACCCCGGGG
"""


@pytest.fixture
def sample_gtf_content():
    """Provide sample GTF content for testing."""
    return """chr1	test	gene	100	400	.	+	.	gene_id "GENE001"; gene_name "GENE_1"; gene_biotype "protein_coding";
chr1	test	transcript	120	380	.	+	.	gene_id "GENE001"; transcript_id "TRANS001"; gene_name "GENE_1";
chr1	test	exon	120	200	.	+	.	gene_id "GENE001"; transcript_id "TRANS001"; exon_id "EXON001";
chr1	test	exon	300	380	.	+	.	gene_id "GENE001"; transcript_id "TRANS001"; exon_id "EXON002";
chr2	test	gene	200	350	.	-	.	gene_id "GENE002"; gene_name "GENE_2"; gene_biotype "lncRNA";
chr2	test	transcript	220	330	.	-	.	gene_id "GENE002"; transcript_id "TRANS002"; gene_name "GENE_2";
chr2	test	exon	220	330	.	-	.	gene_id "GENE002"; transcript_id "TRANS002"; exon_id "EXON003";
"""


@pytest.fixture
def sample_cdna_content():
    """Provide sample cDNA FASTA content for testing."""
    return """>TRANS001
ATCGATCGATCGAAATTTGGGCCCTTTTAAAACCCCGGGGTTTTAAAACCCCGGGG
GGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCC
ATCGATCGATCGAAATTTGGGCCCTTTTAAAACCCCGGGGTTTTAAAACCCCGGGG
GGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCC
ATCGATCGATCGAAATTT
>TRANS002
GGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCC
ATCGATCGATCGAAATTTGGGCCCTTTTAAAA
"""


# Markers for different test categories
pytest_plugins = []

def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line("markers", "unit: Unit tests")
    config.addinivalue_line("markers", "integration: Integration tests")
    config.addinivalue_line("markers", "slow: Slow tests that might download data")
    config.addinivalue_line("markers", "network: Tests that require network access")



