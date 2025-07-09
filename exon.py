from __future__ import annotations
from .genome_element import GenomeElement
from .locus import Locus
from .transcript import Transcript

class Exon(GenomeElement):
    """Represents an exon."""
    def __init__(self, 
                 id: str, 
                 start: int, 
                 end: int, 
                 strand: str, 
                 transcript: Transcript, 
                 **kwargs):
        locus = Locus(transcript.chromosome_id, start, end, strand)
        super().__init__(id, locus, transcript, **kwargs)

    @property
    def exon_id(self) -> str:
        """Returns the ID of the exon."""
        return self.id