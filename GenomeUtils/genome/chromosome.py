#!/usr/bin/env python
"""
Filename: GenomeUtils/genome/chromosome.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.0
Description: This file defines the Chromosome class, representing a biological chromosome.
License: LGPL-3.0-or-later
"""

from __future__ import annotations

from typing import List, TYPE_CHECKING

from Bio import SeqIO
from Bio.Seq import Seq

from .genome_element import GenomeElement
from .locus import Locus


if TYPE_CHECKING:
    from .gene import Gene
    from .genome import Genome

class Chromosome(GenomeElement):
    """Represents a chromosome, with sequence data loaded from file on demand."""

    def __init__(self, 
                 id: str, 
                 seq_index: SeqIO.index, 
                 genome: "Genome" = None,
                 length: int = None,
                 **kwargs):
        """
        Initializes a Chromosome object.

        Args:
            id: The ID of the chromosome.
            seq_index: The `Bio.SeqIO.index` including the sequence of the chromosome.
            genome: The `Genome` object in which the chromosome is located. Optional, defaults to None.
            length: The length of the chromosome. If not provided, it will be inferred from the sequence index.
            **kwargs: Additional keyword arguments.
        """
        self._seq_index = seq_index
        if length is None:
            length = len(self._seq_index[id].seq)
            
        super().__init__(id, Locus(id, 1, length, "+"), genome=genome, **kwargs)
        
    def add_gene(self, gene: "Gene"):
        self._children.append(gene) 
        self._genome.is_indexed = False
        
    @property
    def genes(self) -> List["Gene"]:
        return self._children
        
    @property
    def sequence(self) -> Seq:
        return self._seq_index[self.id].seq
    
    def get_subsequence_by_locus(self, locus: Locus) -> Seq:
        """
        Returns a subsequence of the chromosome for a given Locus.
        """        
        if locus.chr != self.chr:
            raise ValueError(f"The Locus does not belong to this chromosome: locus.chr={locus.chr} != self.chr={self.chr}")
        
        sequence = self.sequence
        if locus.end > len(sequence):
            raise ValueError(f"End coordinate ({locus.end}) is out of bounds for chromosome '{self.id}' (length: {len(sequence)}).")
        
        sequence_slice = sequence[locus.start - 1:locus.end]
        if locus.strand == '+':
            return sequence_slice
        elif locus.strand == '-':
            return sequence_slice.reverse_complement()
        else:
            raise ValueError(f"Invalid strand: {locus.strand}")
        
        
        
        
        
