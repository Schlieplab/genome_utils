#!/usr/bin/env python
"""
Filename: GenomeUtils/genome/exon.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.0
Description: This file defines the Exon class, representing a biological exon.
License: LGPL-3.0-or-later
"""

from __future__ import annotations

from typing import Literal, TYPE_CHECKING

from Bio.Seq import Seq

from .genome_element import GenomeElement
from .locus import Locus


if TYPE_CHECKING:
    from .transcript import Transcript
    from .genome import Genome
    
class Exon(GenomeElement):
    """Represents an exon."""
    def __init__(self, 
                 id: str, 
                 chr: str,
                 start: int, 
                 end: int, 
                 strand: Literal["+", "-"], 
                 transcript: "Transcript" = None, 
                 genome: "Genome" = None,
                 **kwargs):
        """
        Initializes an Exon object.

        Args:
            id: The ID of the exon.
            chr: The chromosome identifier (e.g., 'chr1', '1', 'X').
            start: The genomic start position of the exon in transcript.
            end: The genomic end position of the exon in transcript.
            strand: The strand in which the exon is oriented.
            transcript: The `Transcript` object that the exon belongs to. Optional, defaults to None.
            genome: The `Genome` object in which the exon is located. Optional, defaults to None.
            kwargs: Additional keyword arguments.
        """
        locus = Locus(chr, start, end, strand)
        super().__init__(id, locus, transcript, genome, **kwargs)

    
    def get_transcript(self) -> "Transcript":
        """Returns the `Transcript` object that the exon belongs to."""
        return self.parent
    
    @property
    def sequence(self) -> Seq:
        transcript = self.get_transcript()
        
        try:
            exon_index = transcript.exons.index(self)
        except ValueError:
            return "" 
            
        start_in_transcript = sum(len(exon) for exon in transcript.exons[:exon_index])
        end_in_transcript = start_in_transcript + len(self)
        
        return transcript.sequence[start_in_transcript:end_in_transcript]
    
    