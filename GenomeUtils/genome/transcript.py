#!/usr/bin/env python
"""
Filename: GenomeUtils/genome/transcript.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.3
Description: This file defines the Transcript class, representing a biological transcript.
License: LGPL-3.0-or-later
"""

from __future__ import annotations

from typing import cast, TYPE_CHECKING

from Bio.Seq import Seq

from .genome_element import GenomeElement
from .locus import Locus, Strand


if TYPE_CHECKING:
    from .exon import Exon
    from .gene import Gene
    from .genome import Genome

class Transcript(GenomeElement):
    """Represents a transcript."""

    _sequence: Seq
    _exons: list[Exon]

    def __init__(self, 
                 id: str, 
                 chr: str,
                 start: int, 
                 end: int, 
                 strand: Strand,
                 sequence: Seq,
                 gene: Gene | None = None,
                 genome: Genome | None = None,
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
        self._exons = []
        locus = Locus(chr, start, end, strand)
        super().__init__(id, locus, gene, genome, **kwargs)
        
    @property
    def sequence(self) -> Seq:
        return self._sequence
    
    
    @property
    def exons(self) -> list[Exon]:
        """Returns the list of exons associated with this transcript."""
        return self._exons

    def add_exon(self, exon: "Exon"):
        """Add an `Exon` to the transcript in a sorted manner."""
        # Only add if not already present
        if exon in self._exons:
            return
            
        pos = 0
        # For '+' strand, sort ascending by start coordinate.
        # For '-' strand, sort descending by start coordinate (transcriptional order).
        if self.strand == "+":
            while pos < len(self._exons) and self._exons[pos].start < exon.start:
                pos += 1
        else:  # self.strand == "-"
            while pos < len(self._exons) and self._exons[pos].start > exon.start:
                pos += 1
        self._exons.insert(pos, exon)
        
        cast("Genome", self._genome).is_indexed = False

    
    def get_gene(self) -> "Gene":
        """Returns the `Gene` object that this transcript is associated with."""
        return cast("Gene", self.parent)
    
    def exon_intervals(self) -> list[tuple[int, int]]:
        """Get the exon intervals for this transcript."""
        return [(exon.start, exon.end) for exon in self.exons]
    
    def segment_to_loci(
        self,
        start: int,
        end: int,
    ) -> list[Locus]:
        """Convert a transcript segment to genomic loci.

        Args:
            start: The 0-based, inclusive segment start on the spliced transcript.
            end: The 0-based, exclusive segment end on the spliced transcript.

        Returns:
            One 1-based, inclusive locus per contributing exon, in 5'-to-3'
            transcript order.
        """
        if not self.exons:
            raise ValueError("Transcript must have at least one exon.")

        for exon in self.exons:
            if exon.chr != self.chr:
                raise ValueError("Exon chromosome must match transcript chromosome.")
            if exon.strand != self.strand:
                raise ValueError("Exon strand must match transcript strand.")

        exons_by_position = sorted(self.exons, key=lambda exon: exon.start)
        for previous, current in zip(exons_by_position, exons_by_position[1:]):
            if current.start <= previous.end:
                raise ValueError("Transcript exons must not overlap.")

        if sum(len(exon) for exon in self.exons) != len(self.sequence):
            raise ValueError("Sum of exon lengths must equal transcript sequence length.")

        if not 0 <= start < end <= len(self.sequence):
            raise ValueError(f"Transcript segment [{start}, {end}) is out of bounds.")

        genomic_loci: list[Locus] = []
        transcript_pos = 0

        exons_in_order = sorted(
            self.exons,
            key=lambda exon: exon.start,
            reverse=self.strand == "-",
        )

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

        return genomic_loci
