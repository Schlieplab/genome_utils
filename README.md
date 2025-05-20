# Genome Utils

A Python package for working with genomic data, providing classes and utilities for handling genes, transcripts, exons, and genome sequences.

## Features

- Load and parse genome annotations from GTF files
- Handle genomic sequences from primary assembly FASTA files
- Manage transcript sequences from FASTA files
- Work with genes, transcripts, and exons with their genomic coordinates
- Support for transcript support levels and biotypes
- Efficient sequence retrieval and coordinate mapping
- Support for both forward and reverse strand sequences


## Usage

### Basic Usage

```python
from genome_utils import Genome

# Initialize a genome with reference data
genome = Genome(
    reference_name="GRCm38",
    annotation_version="113",
    gtf_path="path/to/annotation.gtf",
    transcript_fasta_paths=["path/to/transcripts.fa"],
    primary_assembly_path="path/to/primary_assembly.fa"
)

# Index the genome (parse annotations and sequences)
genome.index()

# Get a gene by ID
gene = genome.gene_by_id("ENSMUSG00000000001")

# Get gene sequence
sequence = gene.get_sequence()

# Get a transcript by ID
transcript = genome.transcript_by_id("ENSMUST00000000001")

# Get transcript sequence
transcript_sequence = transcript.sequence

# Get chromosomal positions for transcript coordinates
positions = transcript.get_chromosomal_positions(
    positions=[100, 200, 300],
    window_length=20
)
```

### Working with Genes

```python
# Get all genes
all_genes = genome.genes

# Get transcripts for a gene
transcripts = gene.transcripts

# Get transcripts filtered by support level (1-5)
high_quality_transcripts = gene.get_transcripts_by_support_level(max_level=2)
```

### Working with Transcripts

```python
# Get exons for a transcript
exons = transcript.exons

# Get exon intervals
intervals = transcript.exon_intervals

# Get the exon containing a specific position
exon = transcript.get_exon_by_position(position=100)

# Get a subsequence from the transcript
subsequence = transcript.get_subsequence(start_pos=100, length=50)
```

### Working with Sites

```python
from genome_utils import Site, TargetSite

# Create a basic site
site = Site(sequence="ATCG", chromosomal_position="1:100-104:+")

# Create a target site with additional information
target_site = TargetSite(
    sequence="ATCG",
    chromosomal_position="1:100-104:+",
    gene_id="ENSMUSG00000000001",
    dG=-10.5,
    oligo_dG=-8.2,
    pedersen_steady_state=0.75
)
```

## Dependencies

- Biopython
- Python 3.7+
