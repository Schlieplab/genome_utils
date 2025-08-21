from __future__ import annotations
from typing import List, TYPE_CHECKING
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path

from .genome_element import GenomeElement
from .locus import Locus

if TYPE_CHECKING:
    from .gene import Gene
    from .genome import Genome

class Chromosome(GenomeElement):
    """Represents a chromosome, with sequence data loaded on demand."""

    def __init__(self, 
                 id: str, 
                 start: int, 
                 end: int, 
                 strand: str, 
                 seq_index: SeqIO.index, 
                 genome: "Genome",
                 **kwargs):
        super().__init__(id, Locus(id, start, end, strand), genome=genome, **kwargs)
        self._seq_index = seq_index

    @property
    def genes(self) -> List["Gene"]:
        """Returns the list of genes (children) for this chromosome."""
        return self._children

    def add_gene(self, gene: "Gene"):
        gene._parent = self
        self._children.append(gene) 
        self._genome.is_indexed = False
        
    @property
    def sequence(self) -> Seq:
        return str(self._seq_index[self.id].seq)
    
    
    def get_subsequence(self, locus: Locus) -> str:
        """
        Returns a subsequence of the chromosome for a given Locus.
        """
        if locus.chromosome_id != self.id:
            raise ValueError(f"Locus is for chromosome '{locus.chromosome_id}', but this is chromosome '{self.id}'.")

        if locus.start < 1:
            raise ValueError(f"Start coordinate cannot be less than 1: {locus.start}")
        if locus.end > len(self.sequence):
            raise ValueError(f"End coordinate ({locus.end}) is out of bounds for chromosome '{self.id}' (length: {len(self.sequence)}).")
        
        # The sequence is stored as a string, convert it to a Seq object.
        sequence_slice = Seq(self.sequence[locus.start - 1:locus.end])
        
        if locus.strand == '+':
            return str(sequence_slice)
        elif locus.strand == '-':
            return str(sequence_slice.reverse_complement())
        else:
            raise ValueError(f"Invalid strand: {locus.strand}")
        
        
        
        
        
