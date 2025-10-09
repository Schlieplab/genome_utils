#!/usr/bin/env python
"""
Filename: GenomeUtils/genome/genome_element.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.0
Description: This file defines the abstract base class for all genome elements.
License: LGPL-3.0-or-later
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import List, Optional, TYPE_CHECKING

from Bio.Seq import Seq

from .locus import Locus


if TYPE_CHECKING:
    from .genome import Genome

class GenomeElement(ABC):
    """Abstract base class for genomic elements (e.g. chromosomes, genes, transcripts, exons, etc.)."""

    def __init__(self, 
                 id: str, 
                 locus: Locus,
                 parent: Optional[GenomeElement] = None,
                 genome: "Genome" = None,
                 **kwargs):
        """
        Initializes the GenomeElement.
        
        Args:
            id: The identifier for the genome element.
            locus: The locus of the genome element.
            parent: The parent of the genome element.
            genome: The genome of the genome element.
            **kwargs: Additional attributes for the genome element.
        """
        self.id = id
        self.locus = locus
        self._parent = parent
        self._children: List[GenomeElement] = []
        self._genome: "Genome" = genome

        for key, value in kwargs.items():
            setattr(self, key, value)
    
    
    @property
    def chr(self) -> str:
        return self.locus.chr
    
    @property
    def start(self) -> int:
        return self.locus.start

    @property
    def end(self) -> int:
        return self.locus.end
    
    @property
    def parent(self) -> GenomeElement:
        """Returns the parent of the genome element."""
        if self._parent is None:
            raise AttributeError("Parent `GenomeElement` is not set.")
        return self._parent

    @property
    def strand(self) -> str:
        return self.locus.strand

    def __len__(self) -> int:
        return len(self.locus)

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(id='{self.id}', locus={self.locus!r})"
    
    def __eq__(self, other: GenomeElement) -> bool:
        if not isinstance(other, GenomeElement):
            return False
        return self.id == other.id and self.locus == other.locus
    
    def __hash__(self) -> int:
        return hash((self.id, self.locus))
    
    @property
    @abstractmethod
    def sequence(self) -> Seq:
        pass

