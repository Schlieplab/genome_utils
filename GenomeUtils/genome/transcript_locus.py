#!/usr/bin/env python
"""
Filename: GenomeUtils/genome/transcript_locus.py
Author: Arash Ayat
Copyright: 2026, Alexander Schliep
Version: 0.1.3
Description: This file defines a locus on a spliced transcript sequence.
License: LGPL-3.0-or-later
"""

from dataclasses import dataclass


@dataclass(frozen=True, order=True)
class TranscriptLocus:
    """Represent a 0-based, half-open interval on a transcript sequence."""

    transcript_id: str
    start: int
    end: int

    def __post_init__(self):
        if not self.transcript_id:
            raise ValueError("Transcript ID cannot be empty.")
        if self.start < 0:
            raise ValueError("Start coordinate cannot be less than zero.")
        if self.end <= self.start:
            raise ValueError("End coordinate must be greater than start coordinate.")

    def __len__(self) -> int:
        return self.end - self.start
