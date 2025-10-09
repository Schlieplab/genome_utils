#!/usr/bin/env python
"""
Filename: GenomeUtils/genome/transcript.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.0
Description: This file defines the Transcript class, representing a biological transcript.
License: LGPL-3.0-or-later
"""

from __future__ import annotations

from typing import List, Optional, Tuple, TYPE_CHECKING, Union

from Bio.Seq import Seq

from .genome_element import GenomeElement
from .locus import Locus


if TYPE_CHECKING:
    from .exon import Exon
    from .gene import Gene
    from .genome import Genome

class Transcript(GenomeElement):
    """Represents a transcript."""

    def __init__(self, 
                 id: str, 
                 chr: str,
                 start: int, 
                 end: int, 
                 strand: str, 
                 sequence: Seq,
                 gene: "Gene" = None, 
                 genome: "Genome" = None,
                 **kwargs):
        """
        Initializes a Transcript object.

        Args:
            id: The ID of the transcript.
            chr: The chromosome identifier (e.g., 'chr1', '1', 'X').
            start: The genomic start position of the transcript in chromosome.
            end: The genomic end position of the transcript in chromosome.
            strand: The strand in which the transcript is oriented.
            sequence: The sequence of the transcript.
            gene: The `Gene` object that the transcript is associated with. Optional, defaults to None.
            genome: The `Genome` object in which the transcript is located. Optional, defaults to None.
            kwargs: Additional keyword arguments.
        """
        self._sequence = sequence
        locus = Locus(chr, start, end, strand)
        super().__init__(id, locus, gene, genome, **kwargs)
        
    @property
    def sequence(self) -> Seq:
        return self._sequence
    
    
    @property
    def exons(self) -> List["Exon"]:
        return self._children

    def add_exon(self, exon: "Exon"):
        """Add an `Exon` to the transcript in a sorted manner."""
        pos = 0
        # For '+' strand, sort ascending by start coordinate.
        # For '-' strand, sort descending by start coordinate (transcriptional order).
        if self.strand == "+":
            while pos < len(self._children) and self._children[pos].start < exon.start:
                pos += 1
        else:  # self.strand == "-"
            while pos < len(self._children) and self._children[pos].start > exon.start:
                pos += 1
        self._children.insert(pos, exon)
        
        self._genome.is_indexed = False

    
    def get_gene(self) -> "Gene":
        """Returns the `Gene` object that this transcript is associated with."""
        return self.parent
    
    def exon_intervals(self) -> List[Tuple[int, int]]:
        """Get the exon intervals for this transcript."""
        return [(exon.start, exon.end) for exon in self.exons]
    
    def transcript_to_genomic_pos(self, start: int, end: Optional[int] = None) -> Union[Locus, List[Locus], None]:
        """
        Converts a 0-based, half-open transcript coordinate (or range) to a
        1-based, inclusive genomic coordinate (or list of Locus objects).

        Args:
            start: The 0-based start position on the transcript.
            end: The optional 0-based end position on the transcript. If None, a single
                 point is converted. If provided, the range is [start, end).

        Returns:
            - A Locus object for a single point or for a range within a single exon.
            - A list of Locus objects if the range spans multiple exons.
            - None if a single point maps to no location; an empty list for a range.
        """
        is_single_point = end is None
        if is_single_point:
            end = start + 1

        if not (0 <= start < end <= len(self)):
            if is_single_point:
                raise ValueError(f"Transcript position {start} is out of bounds.")
            else:
                raise ValueError(f"Transcript positions [{start}, {end}) are out of bounds.")

        genomic_loci = []
        transcript_pos = 0
        
        exons_in_order = self.exons

        for exon in exons_in_order:
            exon_len = len(exon)

            overlap_start = max(start, transcript_pos)
            overlap_end = min(end, transcript_pos + exon_len)

            if overlap_start < overlap_end:
                start_in_exon = overlap_start - transcript_pos
                end_in_exon = overlap_end - transcript_pos
                
                if self.strand == '+':
                    genomic_start = exon.start + start_in_exon
                    genomic_end = exon.start + end_in_exon - 1
                else:  # self.strand == "-"
                    genomic_end = exon.end - start_in_exon
                    genomic_start = exon.end - (end_in_exon - 1)
                
                genomic_loci.append(Locus(self.chr, genomic_start, genomic_end, self.strand))

            transcript_pos += exon_len
            
            if transcript_pos >= end:
                break

        if len(genomic_loci) == 1:
            return genomic_loci[0]
        
        return genomic_loci

