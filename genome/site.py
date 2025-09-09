from .locus import Locus
from .genome_element import GenomeElement
from abc import abstractmethod
from typing import TYPE_CHECKING, Literal
from Bio.Seq import Seq

if TYPE_CHECKING:
    from .genome import Genome

class Site(GenomeElement):
    """Abstract base class for genomic sites."""
    def __init__(self, 
                 chr: str, 
                 start: int, 
                 end: int, 
                 strand: Literal["+", "-"], 
                 sequence: Seq,
                 id: str = None,
                 parent: "GenomeElement" = None,
                 genome: "Genome" = None,
                 **kwargs):
        """
        Initializes a Site object.
        
        Args:
            chr: The chromosome of the site.
            start: The start position of the site.
            end: The end position of the site.
            strand: The strand of the site.
            sequence: The sequence of the site.
            id: The ID of the site.
            parent: The parent of the site.
            genome: The genome of the site.
            kwargs: Additional keyword arguments.
        """
        locus = Locus(chr, start, end, strand)
        
        if id is None:
            id = str(locus)
            
        self._sequence = sequence
        super().__init__(id, locus, parent, genome, **kwargs)
        

    @property
    def sequence(self) -> Seq:
        return self._sequence
    

    def __repr__(self):
        return f"{self.__class__.__name__}(id='{self.id}', locus={self.locus!r}, sequence='{self.sequence}')"
    