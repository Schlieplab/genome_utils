from __future__ import annotations
from typing import List, TYPE_CHECKING
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .genome_element import GenomeElement
from .locus import Locus

if TYPE_CHECKING:
    from .gene import Gene


class Chromosome(GenomeElement):
    """Represents a chromosome, with sequence data loaded on demand."""

    def __init__(self, seq_record: SeqRecord, **kwargs):
        """
        Initializes a Chromosome from a SeqRecord proxy.
        The ID, start, end, and strand are derived from the record.
        Args:
            seq_record: The SeqRecord proxy for the chromosome.
            kwargs: Additional keyword arguments.
        """
        id = seq_record.id
        start = 1
        end = len(seq_record)
        strand = '+'
        locus = Locus(id, start, end, strand)
        super().__init__(id, locus, **kwargs)
        self._seq_record = seq_record


    @property
    def genes(self) -> List["Gene"]:
        """Returns the list of genes (children) for this chromosome."""
        return self._children

    def get_subsequence(self, start: int, end: int, strand: str) -> Seq:
        """
        Efficiently retrieves a subsequence from disk without loading
        the entire chromosome into memory.
        Coordinates are 1-based and inclusive.
        """
        if start > end:
            raise ValueError("Start must be less than end")
        if start < 1:
            raise ValueError("Start must be greater than 0")
        if end > len(self._seq_record):
            raise ValueError("End must be less than the length of the chromosome")
        if strand not in ['+', '-']:
            raise ValueError("Strand must be + or -")
        if strand == '+':
            return self._seq_record[start - 1:end].seq
        elif strand == '-':
            return self._seq_record[start - 1:end].reverse_complement().seq
        else:
            return None

    @property
    def sequence(self) -> Seq:
        """
        Returns the full sequence for this chromosome.
        Warning: This loads the entire sequence into memory and can be
        resource-intensive for large chromosomes.
        """
        return self._seq_record.seq

    def add_gene(self, gene: "Gene"):
        gene._parent = self
        self._children.append(gene) 