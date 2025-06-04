from typing import List, Optional, TYPE_CHECKING
from Bio.Seq import Seq
from .transcript import Transcript

class Gene:
    """Class representing a gene with transcripts."""
    def __init__(self, 
                 gene_id: str, 
                 gene_name: str, 
                 chromosome: str, 
                 start: int, 
                 end: int, 
                 strand: str, 
                 biotype: Optional[str] = None,
                 ) -> None:
        """
        Initialize the Gene object.

        Parameters:
            - gene_id: The ID of the gene.
            - gene_name: The name of the gene.
            - chromosome: The chromosome of the gene.
            - start: The start position of the gene in 1-based genomic coordinates.
            - end: The end position of the gene in 1-based genomic coordinates.
            - strand: The strand of the gene (+ or -).
            - biotype: The biotype of the gene.
        """ 
        
        self.gene_id: str = gene_id
        self.gene_name: str = gene_name
        self.chromosome: str = chromosome
        self.start: int = start
        self.end: int = end
        self.strand: str = strand
        self.biotype: Optional[str] = biotype
        self.transcripts: List[Transcript] = []
        self._pre_mrna_sequence: Optional[str] = None
    
    def __str__(self) -> str:
        return f"Gene(gene_id={self.gene_id}, gene_name={self.gene_name}, chromosome={self.chromosome}, start={self.start}, end={self.end}, strand={self.strand}, biotype={self.biotype})"
    
    @property
    def length(self) -> int:
        """Get the length of the gene based on its genomic coordinates."""
        return self.end - self.start + 1
    
    def __len__(self) -> int:
        """Return the length of the gene."""
        return self.length
    
    def add_transcript(self, transcript: Transcript) -> None:
        """Add a transcript to this gene."""
        self.transcripts.append(transcript)
    
    
    @property
    def pre_mrna_sequence(self) -> Optional[str]:
        """
        Get the pre-mRNA sequence for this gene if available.
        
        Returns:
            Optional[str]: The pre-mRNA sequence if available, None otherwise.
        """
        return self._pre_mrna_sequence
    
    @pre_mrna_sequence.setter
    def pre_mrna_sequence(self, sequence: str) -> None:
        """
        Set the pre-mRNA sequence for this gene.
        
        Parameters:
            sequence (str): The pre-mRNA sequence to set.
        """
        self._pre_mrna_sequence = sequence
    
