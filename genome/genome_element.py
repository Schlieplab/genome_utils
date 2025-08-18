from __future__ import annotations
from abc import ABC
from typing import Optional, List, Any, Dict

from .locus import Locus


class GenomeElement(ABC):
    """Abstract base class for genomic elements."""

    def __init__(self, id: str, locus: Locus,
                 parent: Optional[GenomeElement] = None, **kwargs):
        self._attributes: Dict[str, Any] = {
            key: value[0] if isinstance(value, list) and len(value) == 1 else value
            for key, value in kwargs.items()
        }
        self.id = id
        self.locus = locus
        self._parent = parent
        self._children: List[GenomeElement] = []
    
    
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

    def __getattr__(self, name: str) -> Any:
        """Allow direct access to attributes in the attributes dictionary."""
        try:
            return self._attributes[name]
        except KeyError:
            raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")
    
    def __setattr__(self, name: str, value: Any):
        """Allow setting attributes. Explicitly defined attributes are set normally. New, dynamic attributes are stored in the 'attributes' dictionary."""
        if name in self.__dict__ or name in self.__class__.__dict__ or name == '_attributes' or not hasattr(self, '_attributes'):
            super().__setattr__(name, value)
        else:
            # Unpack single-item lists to save memory.
            if isinstance(value, list) and len(value) == 1:
                self._attributes[name] = value[0]
            else:
                self._attributes[name] = value

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(id='{self.id}', locus={self.locus!r})"
    
    def __eq__(self, other: GenomeElement) -> bool:
        return self.id == other.id
    
    def __hash__(self) -> int:
        return hash(self.id)
