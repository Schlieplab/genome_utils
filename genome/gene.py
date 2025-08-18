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
            chromosome: The chromosome that the gene is on.
            kwargs: Additional keyword arguments.
        """
        locus = Locus(chromosome.chromosome_id, start, end, strand)
        super().__init__(id, locus, chromosome, genome, **kwargs)
        self.name = name


    @property
    def transcripts(self) -> List["Transcript"]:
        """Returns the list of transcripts (children) for this gene."""
        return self._children

    def add_transcript(self, transcript: "Transcript"):
        """Add a transcript to the gene."""
        transcript._parent = self
        self._children.append(transcript)
        self._genome.is_indexed = False
    
    def get_chromosome(self) -> "Chromosome":
        """Returns the chromosome that this gene is on."""
        return self._parent
    
    @property
    def sequence(self) -> Optional[Seq]:
        """
        Get the pre-mRNA sequence for this gene if available.
        
        Returns:
            Optional[Seq]: The pre-mRNA sequence as a Bio.Seq object if available, None otherwise.
        """
        if self._parent is None:
            return None
        return self.get_chromosome().get_subsequence(self.start, self.end, self.strand)

    

    
