from __future__ import annotations
from typing import List, Tuple
from Bio.Seq import Seq
from .genome_element import GenomeElement
from .exon import Exon
from .locus import Locus
from .gene import Gene
from .chromosome import Chromosome

class Transcript(GenomeElement):
    """Represents a transcript."""

    def __init__(self, 
                 id: str, 
                 start: int, 
                 end: int, 
                 strand: str, 
                 sequence: Seq,
                 gene: Gene, 
                 **kwargs):
        locus = Locus(gene.chromosome_id, start, end, strand)
        super().__init__(id, locus, gene, **kwargs)
        self.sequence = sequence
    @property
    def exons(self) -> List[Exon]:
        """Returns the list of exons (children) for this transcript."""
        return self._children

    def add_exon(self, exon: Exon):
        """Add an exon to the transcript."""
        exon._parent = self
        self._children.append(exon) 
    
    def __len__(self) -> int:
        return len(self.sequence)
    
    def get_gene(self) -> Gene:
        """Returns the gene that this transcript is on."""
        return self._parent
    
    @property
    def exon_intervals(self) -> List[Tuple[int, int]]:
        """Get the exon intervals for this transcript."""
        return [(exon.start, exon.end) for exon in self.exons]
    
    def get_locus_from_transcript_position(self, transcript_pos: int) -> Locus:
        """
        Converts a 1-based position within the transcript's spliced sequence to a genomic Locus.

        Args:
            transcript_pos: The 1-based position within the transcript's sequence.

        Returns:
            A Locus object representing the specific genomic coordinate.
            
        Raises:
            ValueError: If the transcript has no exons or the position is out of bounds.
        """
        if not self.exons:
            raise ValueError("Transcript has no exons to map coordinates from.")

        if not (1 <= transcript_pos <= len(self)):
            raise ValueError(f"Position {transcript_pos} is out of bounds for transcript of length {len(self)}")

        if self.strand == '+':
            sorted_exons = sorted(self.exons, key=lambda e: e.start)
        else:
            sorted_exons = sorted(self.exons, key=lambda e: e.start, reverse=True)

        transcript_base_count = 0
        for exon in sorted_exons:
            exon_len = len(exon)
            if transcript_base_count + exon_len >= transcript_pos:
                offset_in_exon = transcript_pos - transcript_base_count
                
                if self.strand == '+':
                    genomic_coord = exon.start + offset_in_exon - 1
                else:
                    genomic_coord = exon.end - offset_in_exon + 1
                
                return Locus(self.chromosome_id, genomic_coord, genomic_coord, self.strand)
            
            transcript_base_count += exon_len 
    
    