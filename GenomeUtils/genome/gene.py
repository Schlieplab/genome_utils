#!/usr/bin/env python
"""
Filename: GenomeUtils/genome/gene.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.3
Description: This file defines the Gene class, representing a biological gene.
License: LGPL-3.0-or-later
"""

from __future__ import annotations

from typing import cast, TYPE_CHECKING

from Bio.Seq import Seq

from .genome_element import GenomeElement
from .locus import Locus, Strand


if TYPE_CHECKING:
    from .chromosome import Chromosome
    from .transcript import Transcript
    from .genome import Genome

class Gene(GenomeElement):
    """Represents a gene."""

    name: str

    def __init__(self, 
                 id: str, 
                 name: str,
                 chr: str,
                 start: int, 
                 end: int, 
                 strand: Strand,
                 chromosome: Chromosome | None = None,
                 genome: Genome | None = None,
                 **kwargs):
        """
        Initializes a Gene object.

        Args:
            id: The ID of the gene.
            name: The name of the gene.
            chr: The chromosome identifier (e.g., 'chr1', '1', 'X').
            start: The genomic start position of the gene in chromosome.
            end: The genomic end position of the gene in chromosome.
            strand: The strand in which the gene is oriented.
            chromosome: The `Chromosome` object that the gene is on. Optional, defaults to None.
            genome: The `Genome` object in which the gene is located. Optional, defaults to None.
            kwargs: Additional keyword arguments.
        """
        self.name = name
        locus = Locus(chr, start, end, strand)
        super().__init__(id, locus, chromosome, genome, **kwargs)
    
    @property
    def sequence(self) -> Seq:
        """
        Returns the pre-mRNA sequence of the gene.
        
        The pre-mRNA (precursor mRNA) is the complete genomic sequence spanning
        from the gene's start to end position, including all introns and exons.
        
        Returns:
            The pre-mRNA sequence as a Bio.Seq object.
        """
        return self.get_chromosome().get_subsequence_by_locus(self.locus)
    
    @property
    def transcripts(self) -> list[Transcript]:
        return cast("list[Transcript]", self._children)

    def add_transcript(self, transcript: "Transcript"):
        """Add a transcript to the gene."""
        self._children.append(transcript)
        cast("Genome", self._genome).is_indexed = False
    
    def get_chromosome(self) -> "Chromosome":
        """Returns the `Chromosome` object that this gene is on."""
        return cast("Chromosome", self.parent)
