from __future__ import annotations
from abc import ABC
from typing import Optional, List, Any, Dict, TYPE_CHECKING

from .locus import Locus

if TYPE_CHECKING:
    from .genome import Genome

class GenomeElement(ABC):
    """Abstract base class for genomic elements."""

    def __init__(self, id: str, locus: Locus,
                 parent: Optional[GenomeElement] = None,
                 genome: "Genome" = None,
                 **kwargs):
        
        self.id = id
        self.locus = locus
        self._parent = parent
        self._children: List[GenomeElement] = []
        self._genome: "Genome" = genome

        for key, value in kwargs.items():
            # Unpack single-item lists to save memory
            if isinstance(value, list) and len(value) == 1:
                setattr(self, key, value[0])
            else:
                setattr(self, key, value)
    
    
    @property
    def chromosome_id(self) -> str:
        return self.locus.chromosome_id
    
    @property
    def start(self) -> int:
        return self.locus.start

    @property
    def end(self) -> int:
        return self.locus.end

    @property
    def strand(self) -> str:
        return self.locus.strand

    def __len__(self) -> int:
        return len(self.locus)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(id='{self.id}', locus={self.locus!r})"
    
    def __eq__(self, other: GenomeElement) -> bool:
        return self.id == other.id

