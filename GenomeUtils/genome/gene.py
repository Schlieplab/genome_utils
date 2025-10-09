#!/usr/bin/env python
"""
Filename: GenomeUtils/genome/gene.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.0
Description: This file defines the Gene class, representing a biological gene.
License: LGPL-3.0-or-later
"""

from __future__ import annotations

from typing import List, Literal, TYPE_CHECKING

from Bio.Seq import Seq

from .genome_element import GenomeElement
from .locus import Locus


if TYPE_CHECKING:
    from .chromosome import Chromosome
    from .transcript import Transcript
    from .genome import Genome

class Gene(GenomeElement):
    """Represents a gene."""

    def __init__(self, 
                 id: str, 
                 name: str,
                 chr: str,
                 start: int, 
                 end: int, 
                 strand: Literal["+", "-"], 
                 chromosome: "Chromosome" = None, 
                 genome: "Genome" = None,
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
        return self.get_chromosome().get_subsequence_by_locus(self.locus)
    
    @property
    def transcripts(self) -> List["Transcript"]:
        return self._children

    def add_transcript(self, transcript: "Transcript"):
        """Add a transcript to the gene."""
        self._children.append(transcript)
        self._genome.is_indexed = False
    
    def get_chromosome(self) -> "Chromosome":
        """Returns the `Chromosome` object that this gene is on."""
        return self.parent
