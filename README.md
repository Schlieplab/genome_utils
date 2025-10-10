# GenomeUtils

A Python library for working with annotated genomes. We developed GenomeUtils as an alternative/replacement for the no longer maintained pyensembl.

Object-oriented model for representing genomic features: genomes, chromosomes, genes, transcripts, and exons.

## Features

- Object model: `Genome` > `Chromosome` > `Gene` > `Transcript` > `Exon` (+ `Locus`)
- Builder workflow: `GenomeBuilder` assembles a `Genome` from FASTA (DNA, cDNA) and GTF
- Indexed lookups, optional scaffold separation, streaming/gzip handling
- Downloader utilities: Fetch Ensembl DNA, cDNA, and GTF assets with `EnsemblGenomeDownloader`

## Installation
You can install GenomeUtils via pip with the following command:
```bash
pip install GenomeUtils
```

Requires Python >= 3.10. Dependencies that will be installed automatically by pip are: `biopython`, `gffutils`, `requests`, `gget`.

## Quickstart

### 1) Download and build a genome (complete workflow)

```python
from pathlib import Path
from GenomeUtils.Downloaders import EnsemblGenomeDownloader
from GenomeUtils.Genome import GenomeBuilder

# Download Ensembl assets
downloader = EnsemblGenomeDownloader(
    assembly_id="GRCh38",
    ensembl_release=109,
    species="homo_sapiens",
    genomes_root_dir=Path("./data/genomes"),
)

files = downloader.download()
print(files)  # { 'dna': Path(...), 'cdna': Path(...), 'annotation': Path(...) }


# Build genome from downloaded files
# The builder automatically uses species-appropriate chromosomes:
# Human: 1-22,X,Y,M,MT | Mouse: 1-19,X,Y,M,MT
genome, scaffold_genome = (
    GenomeBuilder(id="GRCh38", species="Homo sapiens", name="Human")
      .with_dna_fasta(files['dna'])
      .with_cdna_fasta(files['cdna'])
      .with_gtf_file(files['annotation'])
      .build()
)

# For other species:
# mouse_genome = GenomeBuilder(id="GRCm39", species="Mus musculus", name="Mouse")...

# Access features
chromosome = genome.chromosome_by_id("1")
first_gene = chromosome.genes[0]
print(first_gene.id, first_gene.name)

# Fast lookups (after build() the genome is indexed)
print(genome.gene_by_id(first_gene.id))
```

### 2) Build a genome from existing files

```python
from pathlib import Path
from GenomeUtils.Genome import GenomeBuilder

# Prepare input files (can be .gz):
dna_fasta = Path("/path/to/genome.dna.fa.gz")
cdna_fasta = Path("/path/to/genome.cdna.fa.gz")
gtf_file  = Path("/path/to/annotations.gtf.gz")

builder = GenomeBuilder(
    id="hg38",
    species="Homo sapiens",
    name="Human Reference Genome",
    separate_scaffolds=False,  # set True to split non-main scaffolds
)

# Optional: limit to specific chromosomes (must be called before with_dna_fasta)
builder.set_chromosome_filter(["chr1", "chr2", "chrX"])  # or ["1","2","X"]

genome, _ = (
    builder
      .with_dna_fasta(dna_fasta)
      .with_cdna_fasta(cdna_fasta)
      .with_gtf_file(gtf_file)
      .build()
)

# Access features
chromosome = genome.chromosome_by_id("chr1")
first_gene = chromosome.genes[0]
print(first_gene.id, first_gene.name)

# Fast lookups (after build() the genome is indexed)
print(genome.gene_by_id(first_gene.id))
```

### 3) Minimal toy example (no files)

```python
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from GenomeUtils.Genome import Genome, Chromosome, Gene, Transcript, Exon

# Create a tiny in-memory genome
genome = Genome(id="toy", species="Test species", name="Toy Genome")
chr1_seq = SeqRecord(Seq("AGCATGATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"), id="chr1")
chromosome = Chromosome("chr1", seq_index={"chr1": chr1_seq}, genome=genome, length=len(chr1_seq.seq))

genome.add_chromosome(chromosome)

gene = Gene(id="GENE001", chr=chromosome, name="MYGENE", start=5, end=35, strand='+', genome=genome)
chromosome.add_gene(gene)

transcript = Transcript(
    id="TRANSCRIPT001",
    chr=chromosome,
    start=5,
    end=35,
    strand='+',
    sequence=Seq("CATGATGCATGCATGCATGCATGCATGC"),
    gene=gene,
    genome=genome,
)

gene.add_transcript(transcript)

transcript.add_exon(Exon(id="EXON001", chr=chromosome, start=5, end=15, strand='+', transcript=transcript, genome=genome))
transcript.add_exon(Exon(id="EXON002", chr=chromosome, start=25, end=35, strand='+', transcript=transcript, genome=genome))

genome.index()
print(genome.gene_by_id("GENE001").name)
```

### 4) Species-specific examples

```python
from pathlib import Path
from GenomeUtils.Downloaders import EnsemblGenomeDownloader
from GenomeUtils.Genome import GenomeBuilder

# Human genome (uses chromosomes 1-22, X, Y, M, MT)
human_genome, _ = GenomeBuilder(
    id="GRCh38", 
    species="Homo sapiens", 
    name="Human Reference Genome"
).with_dna_fasta(human_dna).with_gtf_file(human_gtf).build()

# Mouse genome (uses chromosomes 1-19, X, Y, M, MT)  
mouse_genome, _ = GenomeBuilder(
    id="GRCm39", 
    species="Mus musculus", 
    name="Mouse Reference Genome"
).with_dna_fasta(mouse_dna).with_gtf_file(mouse_gtf).build()


# Override default chromosomes if needed
custom_genome, _ = GenomeBuilder(
    id="custom", 
    species="Custom species", 
    name="Custom Genome",
    main_chromosomes=["chr1", "chr2", "chrX"]  # Only these chromosomes
).with_dna_fasta(custom_dna).with_gtf_file(custom_gtf).build()
```


## Project status

Early-stage library. APIs may evolve.

## Contributing

Issues and PRs are welcome.

Copyright 2025, Alexander Schliep