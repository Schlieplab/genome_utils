#!/usr/bin/env python
"""
Filename: GenomeUtils/genome/locus.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.0
Description: This file defines the Locus class, representing a genomic location.
License: LGPL-3.0-or-later
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal


@dataclass(frozen=True, order=True)
class Locus:
    """Represents a 1-based inclusive genomic coordinates on a chromosome."""
    chr: str
    start: int
    end: int
    strand: Literal["+", "-"] = "+"

    def __post_init__(self):
        """Validate coordinates after initialization."""
        if self.start > self.end:
            raise ValueError("Start coordinate cannot be greater than end coordinate.")
        if self.start < 1:
            raise ValueError("Start coordinate cannot be less than 1.")

    def __len__(self) -> int:
        """Return the length of the locus."""
        return self.end - self.start + 1

    def __repr__(self):
        return f"{self.__class__.__name__}({self.chr}:{self.start}-{self.end}, strand={self.strand})"
    
    def __str__(self):
        return f"{self.chr}:{self.start}-{self.end},{self.strand}"

    def overlaps(self, other: Locus) -> bool:
        """Check if this locus overlaps with another."""
        if self.chr != other.chr:
            return False
        return self.end >= other.start and self.start <= other.end

    def contains(self, other: Locus) -> bool:
        """Check if this locus completely contains another."""
        if self.chr != other.chr:
            return False
        return self.start <= other.start and self.end >= other.end 
    

    