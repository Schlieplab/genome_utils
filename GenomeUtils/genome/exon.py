#!/usr/bin/env python
"""
Filename: GenomeUtils/genome/exon.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.1
Description: This file defines the Exon class, representing a biological exon.
License: LGPL-3.0-or-later
"""

from __future__ import annotations

from typing import Literal, TYPE_CHECKING, List

from Bio.Seq import Seq

from .genome_element import GenomeElement
from .locus import Locus


if TYPE_CHECKING:
    from .gene import Gene
    from .genome import Genome
    from .transcript import Transcript
    
class Exon(GenomeElement):
    """Represents an exon."""
    def __init__(self, 
                 id: str, 
                 chr: str,
                 start: int, 
                 end: int, 
                 strand: Literal["+", "-"], 
                 gene: "Gene" = None, 
                 transcripts: List["Transcript"] = None,
                 genome: "Genome" = None,
                 sequence: Seq = None,
                 **kwargs):
        """
        Initializes an Exon object.

        Args:
            id: The ID of the exon.
            chr: The chromosome identifier (e.g., 'chr1', '1', 'X').
            start: The genomic start position of the exon in transcript.
            end: The genomic end position of the exon in transcript.
            strand: The strand in which the exon is oriented.
            gene: The `Gene` object that the exon belongs to. Optional, defaults to None.
            transcripts: The `Transcript` objects that the exon belongs to. Optional, defaults to None.
            genome: The `Genome` object in which the exon is located. Optional, defaults to None.
            sequence: The sequence of the exon. Optional, takes sequence from transcript if not provided.
            kwargs: Additional keyword arguments.
        """
        self._transcripts = transcripts if transcripts is not None else []
        self._sequence = sequence
        locus = Locus(chr, start, end, strand)
        super().__init__(id, locus, gene, genome, **kwargs)

    
    def get_transcripts(self) -> List["Transcript"]:
        """Returns the `Transcript` object that the exon belongs to."""
        if len(self._transcripts) == 0:
            raise AttributeError("Exon is not associated with any transcripts.")
        return self._transcripts
    
    def add_to_transcript(self, transcript: "Transcript"):
        """Add the exon to the transcript."""
        if transcript not in self._transcripts:
            self._transcripts.append(transcript)
        if self not in transcript.exons:
            transcript.add_exon(self)
    
    def get_gene(self) -> "Gene":
        """Returns the `Gene` object that the exon belongs to."""
        return self.parent
    
    @property
    def sequence(self) -> Seq:
        if self._sequence:
            return self._sequence

        transcripts = self.get_transcripts()

        if len(transcripts) == 0:
            return None

        transcript = transcripts[0]
        exon_index = transcript.exons.index(self)
 
        start_in_transcript = sum(len(exon) for exon in transcript.exons[:exon_index])
        end_in_transcript = start_in_transcript + len(self)
        
        return transcript.sequence[start_in_transcript:end_in_transcript]
    
    