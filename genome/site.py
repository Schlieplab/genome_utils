from .locus import Locus
from .genome_element import GenomeElement
from abc import abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .genome import Genome

class Site(GenomeElement):
    """Abstract base class for genomic sites."""
    def __init__(self, 
                 chr: str, 
                 start: int, 
                 end: int, 
                 strand: str, 
                 sequence: str,
                 id: str = None,
                 parent: "GenomeElement" = None,
                 genome: "Genome" = None,
                 **kwargs):
        
        locus = Locus(chr, start, end, strand)
        
        if id is None:
            id = str(locus)
            
        self.sequence = sequence
        super().__init__(id, locus, parent, genome, **kwargs)
        

    def __repr__(self):
        return f"{self.__class__.__name__}(id='{self.id}', locus={self.locus!r}, sequence='{self.sequence}')"
    
    @abstractmethod
    def get_site_type(self) -> str:
        """Return the type of site. Must be implemented by subclasses."""
        pass
    