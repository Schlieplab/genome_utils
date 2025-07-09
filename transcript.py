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
        """
        Initializes a Transcript object.

        Args:
            id: The ID of the transcript.
            start: The genomic start position of the transcript in chromosome.
            end: The genomic end position of the transcript in chromosome.
            strand: The strand in which the transcript is oriented.
            sequence: The sequence of the transcript.
            gene: The gene that the transcript is associated with.
            kwargs: Additional keyword arguments.
        """
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
    
    def get_locus_from_transcript_position(self, transcript_start_pos: int, transcript_end_pos: int = None) -> Locus:
        """
        Converts a 1-based position or window within the transcript's spliced sequence 
        to a genomic Locus.

        Note: If the window spans an intron, this method returns a single Locus
        that covers the entire genomic region from the start to the end, including
        the intron.

        Args:
            transcript_start_pos: The 1-based start position within the transcript's sequence.
            transcript_end_pos: Optional 1-based end position for a window. If not provided,
                                a single-base Locus is returned.

        Returns:
            A Locus object representing the specific genomic coordinate or window.
            
        Raises:
            ValueError: If the transcript has no exons or positions are out of bounds.
        """
        if transcript_end_pos is None:
            transcript_end_pos = transcript_start_pos
        
        if not self.exons:
            raise ValueError("Transcript has no exons to map coordinates from.")

        # Validate positions
        if not (1 <= transcript_start_pos <= len(self)):
            raise ValueError(f"Start position {transcript_start_pos} is out of bounds for transcript of length {len(self)}")
        if not (1 <= transcript_end_pos <= len(self)):
            raise ValueError(f"End position {transcript_end_pos} is out of bounds for transcript of length {len(self)}")
        if transcript_start_pos > transcript_end_pos:
            raise ValueError(f"Start position {transcript_start_pos} cannot be greater than end position {transcript_end_pos}.")

        # Sort exons by genomic start position, reversing for negative strand
        if self.strand == '+':
            sorted_exons = sorted(self.exons, key=lambda e: e.start)
        else:
            sorted_exons = sorted(self.exons, key=lambda e: e.start, reverse=True)

        # Find the genomic coordinates for the start and end of the window
        genomic_start_coord = None
        genomic_end_coord = None
        
        transcript_base_count = 0
        for exon in sorted_exons:
            exon_len = len(exon)
            
            # Check if the start of the window is in this exon
            if genomic_start_coord is None and transcript_base_count + exon_len >= transcript_start_pos:
                offset_in_exon = transcript_start_pos - transcript_base_count
                if self.strand == '+':
                    genomic_start_coord = exon.start + offset_in_exon - 1
                else:
                    genomic_start_coord = exon.end - offset_in_exon + 1

            # Check if the end of the window is in this exon
            if genomic_end_coord is None and transcript_base_count + exon_len >= transcript_end_pos:
                offset_in_exon = transcript_end_pos - transcript_base_count
                if self.strand == '+':
                    genomic_end_coord = exon.start + offset_in_exon - 1
                else:
                    genomic_end_coord = exon.end - offset_in_exon + 1
            
            # If both are found, we can stop iterating
            if genomic_start_coord is not None and genomic_end_coord is not None:
                break
                
            transcript_base_count += exon_len
        
        # Determine the final genomic start and end for the Locus
        final_start = min(genomic_start_coord, genomic_end_coord)
        final_end = max(genomic_start_coord, genomic_end_coord)
        
        return Locus(self.chromosome_id, final_start, final_end, self.strand)
    
    