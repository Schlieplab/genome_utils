from __future__ import annotations
from typing import List, Optional, Tuple, TYPE_CHECKING
from Bio.Seq import Seq
from .genome_element import GenomeElement
from .locus import Locus
import os
import gzip
import logging

if TYPE_CHECKING:
    from .chromosome import Chromosome
    from .transcript import Transcript
    from .genome import Genome

class Gene(GenomeElement):
    """Represents a gene."""

    def __init__(self, 
                 id: str, 
                 name: str,
                 start: int, 
                 end: int, 
                 strand: str, 
                 chromosome: "Chromosome", 
                 genome: "Genome",
                 **kwargs):
        """
        Initializes a Gene object.

        Args:
            id: The ID of the gene.
            name: The name of the gene.
            start: The genomic start position of the gene in chromosome.
            end: The genomic end position of the gene in chromosome.
            strand: The strand in which the gene is oriented.
            chromosome: The `Chromosome` object that the gene is on.
            genome: The `Genome` object in which the gene is located.
            kwargs: Additional keyword arguments.
        """
        self.name = name
        locus = Locus(chromosome.chr, start, end, strand)
        
        super().__init__(id, locus, chromosome, genome, **kwargs)
    
    @property
    def sequence(self) -> Seq:
        return self.get_chromosome().get_subsequence_by_locus(self.locus)
    
    @property
    def transcripts(self) -> List["Transcript"]:
        return self._children

    def add_transcript(self, transcript: "Transcript"):
        """Add a transcript to the gene."""
        self._children.append(transcript)
        self._genome.is_indexed = False
    
    def get_chromosome(self) -> "Chromosome":
        """Returns the `Chromosome` object that this gene is on."""
        return self._parent

    

    
