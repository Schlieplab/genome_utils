from __future__ import annotations
from typing import TYPE_CHECKING, Literal
from .genome_element import GenomeElement
from .locus import Locus
from Bio.Seq import Seq

if TYPE_CHECKING:
    from .transcript import Transcript
    from .genome import Genome
    
class Exon(GenomeElement):
    """Represents an exon."""
    def __init__(self, 
                 id: str, 
                 start: int, 
                 end: int, 
                 strand: Literal["+", "-"], 
                 transcript: "Transcript", 
                 genome: "Genome",
                 **kwargs):
        """
        Initializes an Exon object.

        Args:
            id: The ID of the exon.
            start: The genomic start position of the exon in transcript.
            end: The genomic end position of the exon in transcript.
            strand: The strand in which the exon is oriented.
            transcript: The `Transcript` object that the exon belongs to.
            genome: The `Genome` object in which the exon is located.
            kwargs: Additional keyword arguments.
        """
        locus = Locus(transcript.chr, start, end, strand)
        super().__init__(id, locus, transcript, genome, **kwargs)

    
    def get_transcript(self) -> "Transcript":
        """Returns the `Transcript` object that the exon belongs to."""
        return self._parent
    
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
    
    