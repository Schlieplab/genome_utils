Getting started
===============

The :mod:`GenomeUtils` package provides an object-oriented toolkit for downloading,
building, and exploring genomes. This guide walks you through installation and
common workflows so you can start integrating GenomeUtils into your data pipelines.

Installation
------------

The package is published on PyPI. Install it with pip (Python 3.10+):

.. code-block:: bash

   pip install GenomeUtils

Core concepts
-------------

GenomeUtils models genomic information with a hierarchy of Python classes:

- :class:`GenomeUtils.Genome.Genome` contains chromosomes and maintains fast lookup indexes.
- :class:`GenomeUtils.Genome.Chromosome` stores sequence references and gene collections.
- :class:`GenomeUtils.Genome.Gene`, :class:`GenomeUtils.Genome.Transcript`, and :class:`GenomeUtils.Genome.Exon`
    represent individual genomic features.
- :class:`GenomeUtils.Genome.GenomeBuilder` orchestrates parsing FASTA and GTF files to build genomes.
- :class:`GenomeUtils.Downloaders.EnsemblGenomeDownloader` fetches genome assets from Ensembl.

Complete workflow example
-------------------------

The snippet below downloads Ensembl resources and builds an indexed genome.

.. code-block:: python

    from pathlib import Path
    from GenomeUtils.Downloaders import EnsemblGenomeDownloader
    from GenomeUtils.Genome import GenomeBuilder

   downloader = EnsemblGenomeDownloader(
       assembly_id="GRCh38",
       ensembl_release=109,
       species="homo_sapiens",
       genomes_root_dir=Path("./data/genomes"),
   )

   files = downloader.download()

   human_genome, scaffold_genome = (
       GenomeBuilder(id="GRCh38", species="Homo sapiens", name="Human")
       .with_dna_fasta(files["dna"])
       .with_cdna_fasta(files["cdna"])
       .with_gtf_file(files["annotation"])
       .build()
   )

   chromosome = human_genome.chromosome_by_id("1")
   first_gene = chromosome.genes[0]
   print(first_gene.id, first_gene.name)
   print(human_genome.gene_by_id(first_gene.id))

Building from existing files
----------------------------

If you already have FASTA and GTF files on disk, pass them directly to the builder.

.. code-block:: python

    from pathlib import Path
    from GenomeUtils.Genome import GenomeBuilder

   dna_fasta = Path("/path/to/genome.dna.fa.gz")
   cdna_fasta = Path("/path/to/genome.cdna.fa.gz")
   gtf_file = Path("/path/to/annotations.gtf.gz")

   builder = GenomeBuilder(
       id="hg38",
       species="Homo sapiens",
       name="Human Reference Genome",
       separate_scaffolds=False,
   )

   builder.set_chromosome_filter(["chr1", "chr2", "chrX"])

   genome, _ = (
       builder
       .with_dna_fasta(dna_fasta)
       .with_cdna_fasta(cdna_fasta)
       .with_gtf_file(gtf_file)
       .build()
   )

   chromosome = genome.chromosome_by_id("chr1")
   first_gene = chromosome.genes[0]
   print(first_gene.id, first_gene.name)
   print(genome.gene_by_id(first_gene.id))

Minimal in-memory example
-------------------------

For unit tests or demonstrations, you can construct entire genomes in memory.

.. code-block:: python

    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from GenomeUtils.Genome import Genome, Chromosome, Gene, Transcript, Exon

    # Create a tiny in-memory genome
    genome = Genome(id="toy", species="Test species", name="Toy Genome")
    chr1_seq = SeqRecord(Seq("AGCATGATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"), id="chr1")
    chromosome = Chromosome("chr1", seq_index={"chr1": chr1_seq}, genome=genome, length=len(chr1_seq.seq))

    genome.add_chromosome(chromosome)

    gene = Gene(id="GENE001", chr=chromosome.id, name="MYGENE", start=5, end=35, strand='+', genome=genome, chromosome=chromosome)
    chromosome.add_gene(gene)

    transcript = Transcript(
        id="TRANSCRIPT001",
        chr=chromosome.id,
        start=5,
        end=35,
        strand='+',
        sequence=Seq("CATGATGCATGCATGCATGCATGCATGC"),
        gene=gene,
        genome=genome,
    )

    gene.add_transcript(transcript)

    Exon(id="EXON001", chr=chromosome.id, start=5, end=15, strand='+', gene=gene, genome=genome).add_to_transcript(transcript)
    Exon(id="EXON002", chr=chromosome.id, start=25, end=35, strand='+', gene=gene, genome=genome).add_to_transcript(transcript)


    genome.index()
    print(genome.gene_by_id("GENE001").name)

Next steps
----------

- Browse the :doc:`API reference <api/index>` for module-level documentation.
- Review the repository README for release notes and contribution guidelines.
