GenomeUtils.Genome
==================

.. currentmodule:: GenomeUtils.Genome


.. list-table:: Class Descriptions
   :header-rows: 0
   :widths: 25 75

   * - :py:class:`Genome`
     - Top-level container that manages chromosomes, genes, transcripts, and exons with indexing utilities.
   * - :py:class:`Gene`
     - Locus-bound gene that owns transcripts and provides access to the genomic sequence.
   * - :py:class:`Transcript`
     - Transcript with ordered exons, canonical sequence, and helpers for coordinate conversion.
   * - :py:class:`Exon`
     - Exon segment attached to a transcript and able to derive its nucleotide sequence.
   * - :py:class:`Chromosome`
     - Lazy-loaded chromosome wrapper that exposes sequence slices via loci and tracks genes.
   * - :py:class:`Locus`
     - Immutable utility for 1-based inclusive genomic coordinates with overlap/containment helpers.
   * - :py:class:`GenomeElement`
     - Abstract base class that unifies shared behavior for loci-based genome entities.
   * - :py:class:`GenomeBuilder`
     - Fluent builder that assembles a :py:class:`Genome` from FASTA and GTF inputs.


.. rubric:: Classes


.. automodule:: GenomeUtils.Genome
   :members:
   :undoc-members:
   :show-inheritance: