from __future__ import annotations
from dataclasses import dataclass
from typing import Literal

@dataclass(frozen=True, order=True)
class Locus:
    """Represents a locus on a chromosome."""
    chromosome_id: str
    start: int
    end: int
    strand: Literal["+", "-"] = "+"

    def __post_init__(self):
        """Validate coordinates after initialization."""
        if self.start > self.end:
            raise ValueError("Start coordinate cannot be greater than end coordinate.")

    def __len__(self) -> int:
        """Return the length of the locus."""
        return self.end - self.start + 1

    def __repr__(self):
        return f"{self.__class__.__name__}({self.chromosome_id}:{self.start}-{self.end} {self.strand})"

    def overlaps(self, other: Locus) -> bool:
        """Check if this locus overlaps with another."""
        if self.chromosome_id != other.chromosome_id:
            return False
        return self.end >= other.start and self.start <= other.end

    def contains(self, other: Locus) -> bool:
        """Check if this locus completely contains another."""
        if self.chromosome_id != other.chromosome_id:
            return False
        return self.start <= other.start and self.end >= other.end 