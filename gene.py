from typing import List, Optional, TYPE_CHECKING
from Bio.Seq import Seq
from .transcript import Transcript

class Gene:
    """Class representing a gene with transcripts."""
    def __init__(self, gene_id: str, gene_name: str, chromosome: str, 
                 start: int, end: int, strand: str, biotype: Optional[str] = None,) -> None:
        self.gene_id: str = gene_id
        self.gene_name: str = gene_name
        self.chromosome: str = chromosome
        self.start: int = start  # 1-based genomic coordinates
        self.end: int = end      # 1-based genomic coordinates
        self.strand: str = strand
        self.biotype: Optional[str] = biotype
        self.transcripts: List[Transcript] = []
    
    def add_transcript(self, transcript: Transcript) -> None:
        """Add a transcript to this gene."""
        self.transcripts.append(transcript)
    
    def get_transcripts_by_support_level(self, max_level: Optional[int] = None) -> List['Transcript']:
        """
        Get transcripts filtered by support level.
        
        Args:
            max_level (Optional[int]): Maximum support level to include (1-5), 
                                      or None to include all transcripts
            
        Returns:
            List[Transcript]: List of transcripts with support level less than or equal to max_level,
                             or all transcripts if max_level is None
        """
        if max_level is None:
            return self.transcripts
        
        return [t for t in self.transcripts 
                if t.support_level is not None and t.support_level <= max_level]
    
