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

To work from a local clone instead, install the project in editable mode:

.. code-block:: bash

   pip install -e .

Core concepts
-------------

GenomeUtils models genomic information with a hierarchy of Python classes:

- :class:`GenomeUtils.genome.Genome` contains chromosomes and maintains fast lookup indexes.
- :class:`GenomeUtils.genome.Chromosome` stores sequence references and gene collections.
- :class:`GenomeUtils.genome.Gene`, :class:`GenomeUtils.genome.Transcript`, and :class:`GenomeUtils.genome.Exon`
    represent individual genomic features.
- :class:`GenomeUtils.genome.GenomeBuilder` orchestrates parsing FASTA and GTF files to build genomes.
- :class:`GenomeUtils.Downloaders.EnsemblGenomeDownloader` fetches genome assets from Ensembl.

Complete workflow example
-------------------------

The snippet below downloads Ensembl resources and builds an indexed genome.

.. code-block:: python

    from pathlib import Path
    from GenomeUtils.Downloaders import EnsemblGenomeDownloader
    from GenomeUtils.genome import GenomeBuilder

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
    from GenomeUtils.genome import GenomeBuilder

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
    from GenomeUtils.genome import Genome
    from GenomeUtils.genome import Chromosome
    from GenomeUtils.genome import Gene
    from GenomeUtils.genome import Transcript
    from GenomeUtils.genome import Exon

   genome = Genome(id="toy", species="Test species", name="Toy Genome")
   chr1_seq = SeqRecord(Seq("AGCATGATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"), id="chr1")
   chromosome = Chromosome("chr1", seq_index={"chr1": chr1_seq}, genome=genome, length=len(chr1_seq.seq))
   genome.add_chromosome(chromosome)

   gene = Gene(id="GENE001", chr=chromosome, name="MYGENE", start=5, end=35, strand="+", genome=genome)
   chromosome.add_gene(gene)

   transcript = Transcript(
       id="TRANSCRIPT001",
       chr=chromosome,
       start=5,
       end=35,
       strand="+",
       sequence=Seq("CATGATGCATGCATGCATGCATGCATGC"),
       gene=gene,
       genome=genome,
   )

   gene.add_transcript(transcript)
   transcript.add_exon(Exon(id="EXON001", chr=chromosome, start=5, end=15, strand="+", transcript=transcript, genome=genome))
   transcript.add_exon(Exon(id="EXON002", chr=chromosome, start=25, end=35, strand="+", transcript=transcript, genome=genome))

   genome.index()
   assert genome.gene_by_id("GENE001").name == "MYGENE"

Next steps
----------

- Browse the :doc:`API reference <api/index>` for module-level documentation.
- Review the repository README for release notes and contribution guidelines.
