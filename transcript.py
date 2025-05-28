from typing import Dict, List, Optional, Tuple
from Bio.Seq import Seq


from .exon import Exon

class Transcript:
    """Class representing a transcript with exons."""
    def __init__(self, transcript_id: str, gene_id: str, chromosome: str, 
                 start: int, end: int, strand: str, biotype: Optional[str] = None,
                 support_level: Optional[int] = None) -> None:
        self.transcript_id: str = transcript_id
        self.gene_id: str = gene_id
        self.chromosome: str = chromosome
        self.start: int = start  # 1-based genomic coordinates
        self.end: int = end      # 1-based genomic coordinates
        self.strand: str = strand
        self.biotype: Optional[str] = biotype
        self.support_level: Optional[int] = support_level  # 1-5 or None
        self.exons: List['Exon'] = []
        self._sequence: Optional[str] = None
        self._genomic_coordinate_map: Optional[Dict[int, int]] = None
    
    @property
    def length(self) -> int:
        """Get the total length of the transcript (sum of exon lengths)."""
        return sum(exon.end - exon.start + 1 for exon in self.exons)
    
    def add_exon(self, exon: 'Exon') -> None:
        """Add an exon to this transcript."""
        self.exons.append(exon)
        # Keep exons sorted by position
        self.exons.sort(key=lambda e: e.start)
        # Reset the coordinate map since exon structure changed
        self._genomic_coordinate_map = None
    
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
                exon_length = exon.end - exon.start + 1
                for offset in range(exon_length):
                    mapping[transcript_pos] = exon.start + offset
                    transcript_pos += 1
        else:
            # Reverse strand: process exons in reverse genomic order (5' to 3' for transcript)
            sorted_exons = sorted(self.exons, key=lambda e: e.start, reverse=True)
            
            for exon in sorted_exons:
                exon_length = exon.end - exon.start + 1
                for offset in range(exon_length):
                    mapping[transcript_pos] = exon.end - offset
                    transcript_pos += 1
                    
        return mapping
    

    
    def get_exon_by_position(self, position: int) -> Optional['Exon']:
        """
        Get the exon containing a specific position in the transcript.
        
        Args:
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
        
        Args:
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
        
    def get_chromosomal_positions(self, positions: List[int], window_length: int) -> List[Optional[str]]:
        """
        Get chromosomal coordinates for specified sequences within this transcript.
        
        Args:
            positions (List[int]): List of 1-based positions within the transcript
            window_length (int): Length of the window/sequence at each position
            
        Returns:
            List[Optional[str]]: Chromosomal coordinates in format "chrom:start-end:strand"
        """
        # Build mapping from transcript to genomic coordinates
        mapping = self.genomic_coordinate_map
        if not mapping:
            return [None] * len(positions)
            
        results: List[Optional[str]] = []
        
        for pos in positions:
            try:
                # Check if both start and end positions are mapped
                start_genomic = mapping.get(pos)
                end_genomic = mapping.get(pos + window_length - 1)
                
                if start_genomic is None or end_genomic is None:
                    results.append(None)
                    continue
                    
                # Sort the coordinates in ascending order regardless of strand
                start, end = min(start_genomic, end_genomic), max(start_genomic, end_genomic)
                
                # Format result string
                result = f"{self.chromosome}:{start}-{end}:{self.strand}"
                results.append(result)
                
            except Exception as e:
                results.append(None)
                
        return results
    
    def get_chromosomal_position(
        self, position: int, window_length: int
        ) -> Optional[str]:
            """
            Get chromosomal position for a single position within a transcript.
            
            Args:
                position (int): 1-based position within the transcript
                window_length (int): Length of the window starting at the position
                
            Returns:
                Optional[str]: Chromosomal coordinates in format "chrom:start-end:strand"
            """
            results = self.get_chromosomal_positions(
                positions=[position],
                window_length=window_length
            )
            return results[0] if results else None 