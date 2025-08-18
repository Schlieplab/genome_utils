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

    def __init__(self, id: str, start: int, end: int, strand: str, seq_index: SeqIO.index, fasta_path: Path, genome: "Genome"):
        super().__init__(id, Locus(id, start, end, strand), genome=genome)
        self._seq_index = seq_index
        self._fasta_path = fasta_path

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
        return self._seq_index[self.id].seq
    
    def __getitem__(self, key: slice) -> Seq:
        return self.sequence[key]
    
    def get_subsequence(self, start: int, end: int, strand: str = "+") -> Seq:

        if start > end:
            raise ValueError(f"Start coordinate cannot be greater than end coordinate: {start} > {end}")
        if start < 1:
            raise ValueError(f"Start coordinate cannot be less than 1: {start}")
        if end > len(self.sequence):
            raise ValueError(f"End coordinate cannot be greater than the length of the chromosome: {end} > {len(self.sequence)}")
        
        if strand == '+':
            return self.sequence[start - 1:end]
        elif strand == '-':
            return self.sequence[start - 1:end].reverse_complement()
        else:
            raise ValueError(f"Invalid strand: {strand}")
        
    def __getstate__(self):
        """Prepare the object for pickling."""
        state = self.__dict__.copy()
        # Remove the unpickleable entry, which is the file handle
        del state['_seq_index']
        return state

    def __setstate__(self, state):
        """Restore the object after unpickling."""
        self.__dict__.update(state)
        # Re-create the seq_index from the fasta path
        self._seq_index = SeqIO.index(str(self._fasta_path), "fasta")
        
        
        
