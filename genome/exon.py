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

    @property
    def exon_id(self) -> str:
        """Returns the ID of the exon."""
        return self.id
    
    def get_transcript(self) -> "Transcript":
        """Returns the transcript that the exon belongs to."""
        return self._parent
    
    @property
    def sequence(self) -> str:
        """Returns the sequence of the exon."""
        transcript_start_pos = self.get_transcript().genomic_to_transcript_pos(self.start)
        transcript_end_pos = self.get_transcript().genomic_to_transcript_pos(self.end)
        return self.get_transcript().sequence[transcript_start_pos:transcript_end_pos]
    