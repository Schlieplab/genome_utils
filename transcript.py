from typing import Dict, List, Optional, Tuple
from Bio.Seq import Seq


from .exon import Exon

class Transcript:
    """Class representing a transcript with exons."""
    def __init__(self, 
                 transcript_id: str, 
                 gene_id: str, 
                 chromosome: str, 
                 start: int, 
                 end: int, 
                 strand: str, 
                 biotype: Optional[str] = None,
                 support_level: Optional[int] = None
                 ) -> None:
        """
        Initialize the Transcript object.

        Parameters:
            - transcript_id: The ID of the transcript.
            - gene_id: The ID of the gene.
            - chromosome: The chromosome of the transcript.
            - start: The start position of the transcript in 1-based genomic coordinates.
            - end: The end position of the transcript in 1-based genomic coordinates.
            - strand: The strand of the transcript (+ or -).
            - biotype: The biotype of the transcript.
            - support_level: The support level of the transcript.
        """
        
        self.transcript_id: str = transcript_id
        self.gene_id: str = gene_id
        self.chromosome: str = chromosome
        self.start: int = start  
        self.end: int = end     
        self.strand: str = strand
        self.biotype: Optional[str] = biotype
        self.support_level: Optional[int] = support_level  
        self.exons: List['Exon'] = []
        self._sequence: Optional[str] = None
        self._genomic_coordinate_map: Optional[Dict[int, int]] = None
    
    def __str__(self) -> str:
        return f"Transcript(transcript_id={self.transcript_id}, gene_id={self.gene_id}, chromosome={self.chromosome}, start={self.start}, end={self.end}, strand={self.strand}, biotype={self.biotype}, support_level={self.support_level})"
    
    @property
    def length(self) -> int:
        """Get the total length of the transcript (sum of exon lengths)."""
        return len(self._sequence)
    
    def __len__(self) -> int:
        """Return the total length of the transcript."""
        return self.length
    
    @property
    def sequence(self) -> Optional[str]:
        """Get the transcript sequence."""
        return self._sequence
    
    @sequence.setter
    def sequence(self, seq: str) -> None:
        """Set the transcript sequence."""
        self._sequence = seq
        
    def add_exon(self, exon: 'Exon') -> None:
        """Add an exon to this transcript."""
        self.exons.append(exon)
        # Keep exons sorted by position
        self.exons.sort(key=lambda e: e.start)
        # Reset the coordinate map since exon structure changed
        self._genomic_coordinate_map = None
        
    @property
    def exon_intervals(self) -> List[Tuple[int, int]]:
        """Get the exon intervals for this transcript."""
        return [(exon.start, exon.end) for exon in self.exons]
    
    @property
    def genomic_coordinate_map(self) -> Dict[int, int]:
        """
        Get the cached genomic coordinate map, building it if necessary.
        
        Returns:
            Dict[int, int]: Dictionary mapping transcript positions to genomic positions
        """
        if self._genomic_coordinate_map is None:
            self._genomic_coordinate_map = self._build_genomic_coordinate_map()
        return self._genomic_coordinate_map
    
    def _build_genomic_coordinate_map(self) -> Dict[int, int]:
        """
        Internal method to build a mapping from transcript positions to genomic positions.
        
        Returns:
            Dict[int, int]: Dictionary mapping transcript positions to genomic positions
        """
        if not self.exons:
            return {}
            
        mapping: Dict[int, int] = {}
        transcript_pos: int = 1  # 1-based position in transcript
        
        if self.strand == '+':
            # Forward strand: process exons in genomic order (5' to 3')
            sorted_exons = sorted(self.exons, key=lambda e: e.start)
            
            for exon in sorted_exons:
                for genomic_coord in range(exon.start, exon.end + 1):
                    mapping[transcript_pos] = genomic_coord
                    transcript_pos += 1
        else:
            # Reverse strand: process exons in reverse genomic order (5' to 3' for transcript)
            sorted_exons = sorted(self.exons, key=lambda e: e.start, reverse=True)
            
            for exon in sorted_exons:
                for genomic_coord in range(exon.end, exon.start - 1, -1):
                    mapping[transcript_pos] = genomic_coord
                    transcript_pos += 1
                    
        return mapping
    

    
    def get_exon_by_position(self, position: int) -> Optional['Exon']:
        """
        Get the exon containing a specific position in the transcript.
        
        Parameters:
            position (int): 1-based position in the transcript
            
        Returns:
            Optional[Exon]: The exon containing the position, or None if not found
        """
        if position < 1 or not self.exons:
            return None
            
        current_pos: int = 1
        
        # Sort exons based on strand direction
        if self.strand == '+':
            exon_order = sorted(self.exons, key=lambda e: e.start)
        else:  # reverse strand
            exon_order = sorted(self.exons, key=lambda e: e.start, reverse=True)
            
        for exon in exon_order:
            exon_length = exon.end - exon.start + 1
            if current_pos <= position < current_pos + exon_length:
                return exon
            current_pos += exon_length
            
        return None
        
    def get_subsequence(self, start_pos: int, length: int) -> Optional[str]:
        """
        Get a subsequence from this transcript starting at the specified position.
        
        Parameters:
            start_pos (int): 1-based start position in the transcript
            length (int): Length of the subsequence to return
            
        Returns:
            Optional[str]: The subsequence of specified length, or None if invalid position/length
        """
        if not self._sequence:
            return None
            
        # Convert to 0-based indexing for Python string operations
        idx = start_pos - 1
        
        # Check bounds
        if idx < 0 or idx + length > len(self._sequence):
            return None
            
        return self._sequence[idx:idx + length]

    def get_chromosomal_position_in_chrmosome(self, position: int) -> Optional[int]:
        """
        Get chromosomal position for a single position within a chromosome.
        
        Parameters:
            position (int): 1-based position within the transcript
        
        Returns:
            Optional[str]: Chromosomal coordinates in format "chrom:start-end:strand"
        """
        
        mapping = self.genomic_coordinate_map
        if not mapping: 
            return None
            
        return mapping.get(position)
    
    def get_chromosomal_window(
        self, position: int, window_length: int = 1
        ) -> Optional[str]:
            """
            Get chromosomal position for a window within a transcript.
            
            Parameters:
                position (int): 1-based position within the transcript
                window_length (int): Length of the window starting at the position (default: 1)
                
            Returns:
                Optional[str]: Chromosomal coordinates in format "chrom:start-end:strand"
            """
            if window_length < 1:
                return None

            start_transcript_pos = position
            end_transcript_pos = position + window_length - 1

            start_genomic = self.get_chromosomal_position_in_chrmosome(start_transcript_pos)
            end_genomic = self.get_chromosomal_position_in_chrmosome(end_transcript_pos)

            if start_genomic is None or end_genomic is None:
                return None

            genomic_interval_start = min(start_genomic, end_genomic)
            genomic_interval_end = max(start_genomic, end_genomic)

            return f"{self.chromosome}:{genomic_interval_start}-{genomic_interval_end}:{self.strand}" 

    def to_gtf_entry(self, source: str = "custom") -> str:
        """
        Returns a GTF formatted string for this transcript.
        'source' is customizable. Attributes follow GTF specification.
        """
        attributes = [
            f'gene_id "{self.gene_id}"',
            f'transcript_id "{self.transcript_id}"'
        ]
        if self.biotype:
            attributes.append(f'transcript_biotype "{self.biotype}"')
        if self.support_level is not None:
            attributes.append(f'transcript_support_level "{self.support_level}"') 
        
        attributes_str = "; ".join(attributes) + ";"
        
        return f"{self.chromosome}\t{source}\ttranscript\t{self.start}\t{self.end}\t.\t{self.strand}\t.\t{attributes_str}" 