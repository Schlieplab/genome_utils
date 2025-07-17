# genome-utils

A modern, object-oriented Python library for working with genomic data.

`genome-utils` provides a clean, intuitive, and efficient object model for representing and manipulating genomic features like genes, transcripts, and exons. It is designed to be both powerful for complex bioinformatics tasks and easy to use for everyday scripting.

## Key Features

- **Object-Oriented Design:** Represents genomic features as a clear hierarchy of Python objects (`Genome` > `Chromosome` > `Gene` > `Transcript` > `Exon`).
- **Performance-Aware:** Utilizes lazy loading for sequence data from FASTA files, ensuring low memory usage even with large genomes.
- **Modern & Pythonic:** Built with modern Python features, including dataclasses for robust value objects (`Locus`) and comprehensive type hinting for clarity and editor support.
- **Intuitive API:** Access data in a natural way, from iterating over all genes in a genome to mapping coordinates between transcripts and chromosomes.
- **Flexible & Extensible:** Easily attach arbitrary metadata to any genomic feature.

## Core Concepts

The library is built around a few central classes:

-   `Genome`: The top-level container for an entire genome assembly. It manages chromosomes and provides fast, indexed lookups for all genomic features.
-   `Chromosome`: Represents a chromosome, providing access to its full sequence and the genes it contains.
-   `Gene`: Represents a gene, which is a collection of transcripts.
-   `Transcript`: Represents a specific transcript of a gene, containing its exons and its spliced mRNA sequence.
-   `Exon`: Represents an exon, the fundamental unit of a transcript.
-   `Locus`: An immutable dataclass that precisely defines a genomic coordinate range (`chromosome_id`, `start`, `end`, `strand`).

## Installation

```bash
pip install .
```
*(Assuming a `setup.py` or `pyproject.toml` is present for local installation).*

## Quickstart: Example Usage

Here's how to model a simple gene and explore its features.

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from genome_utils import Genome, Chromosome, Gene, Transcript, Exon

# 1. Set up a Chromosome with its sequence
# (In a real scenario, this comes from a FASTA file)
chr1_seq = SeqRecord(Seq("AGCATGATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"), id="chr1")
chromosome = Chromosome(seq_record=chr1_seq)

# 2. Create a Genome object and add the chromosome
genome = Genome(id="hg38_toy", species="Homo sapiens", name="Toy Genome")
genome.add_chromosome(chromosome)

# 3. Build the feature hierarchy (Gene -> Transcript -> Exon)
# In a real application, you would parse this from a GFF/GTF file.

# Create a Gene
gene = Gene(id="GENE001", 
            name="MYGENE", 
            start=5, 
            end=35, 
            strand='+', 
            chromosome=chromosome)
chromosome.add_gene(gene)

# Create a Transcript for the Gene
# The sequence is the final spliced mRNA sequence
transcript = Transcript(
    id="TRANSCRIPT001",
    start=5,
    end=35,
    strand='+',
    sequence=Seq("CATGATGCATGCATGCATGCATGCATGC"), # Spliced sequence
    gene=gene
)
gene.add_transcript(transcript)

# Add Exons to the Transcript
exon1 = Exon(id="EXON001", start=5, end=15, strand='+', transcript=transcript)
exon2 = Exon(id="EXON002", start=25, end=35, strand='+', transcript=transcript)
transcript.add_exon(exon1)
transcript.add_exon(exon2)

# 4. Index the genome for fast lookups
genome.index()

# 5. Now, you can easily access your data

# Get a gene by its ID
my_gene = genome.gene_by_id("GENE001")
print(f"Found Gene: {my_gene.name}")
# > Found Gene: MYGENE

# Get the gene's pre-mRNA sequence from the chromosome
print(f"Gene Sequence: {my_gene.sequence}")
# > Gene Sequence: CATGATGCATGCATGCATGCATGCATGCATG

# Get the transcript's spliced mRNA sequence
my_transcript = genome.transcript_by_id("TRANSCRIPT001")
print(f"Transcript Sequence: {my_transcript.sequence}")
# > Transcript Sequence: CATGATGCATGCATGCATGCATGCATGC

# Map a position from transcript coordinates to genomic coordinates
# Where is the 15th base of the spliced mRNA located on the chromosome?
locus = my_transcript.get_locus_from_transcript_position(15)
print(f"Position 15 in transcript maps to: {locus}")
# > Position 15 in transcript maps to: Locus(chr1:29-29 +)

```

## Contributing

Contributions are welcome! Please feel free to submit a pull request or open an issue for bugs, feature requests, or improvements. 