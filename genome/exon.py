from __future__ import annotations
from typing import TYPE_CHECKING
from .genome_element import GenomeElement
from .locus import Locus

if TYPE_CHECKING:
    from .transcript import Transcript
    from .genome import Genome
    
class Exon(GenomeElement):
    """Represents an exon."""
    def __init__(self, 
                 id: str, 
                 start: int, 
                 end: int, 
                 strand: str, 
                 transcript: "Transcript", 
                 genome: "Genome",
                 **kwargs):
        locus = Locus(transcript.chromosome_id, start, end, strand)
        super().__init__(id, locus, transcript, genome, **kwargs)

    
    def get_transcript(self) -> "Transcript":
        """Returns the transcript that the exon belongs to."""
        return self._parent
    
    @property
    def sequence(self) -> str:
        """Returns the sequence of the exon."""
        transcript = self.get_transcript()
        
        # Find the start position of this exon within the transcript's spliced sequence
        # by summing the lengths of all preceding exons.
        try:
            exon_index = transcript.exons.index(self)
        except ValueError:
            # This should not happen if the exon is properly associated with its transcript.
            return "" 
            
        start_in_transcript = sum(len(exon) for exon in transcript.exons[:exon_index])
        end_in_transcript = start_in_transcript + len(self)
        
        return transcript.sequence[start_in_transcript:end_in_transcript]
    
    