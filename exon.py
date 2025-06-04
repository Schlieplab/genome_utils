class Exon:
    """Class representing an exon within a transcript."""
    def __init__(self, exon_id: str, start: int, end: int, transcript_id: str) -> None:
        self.exon_id: str = exon_id
        self.start: int = start 
        self.end: int = end
        self.transcript_id: str = transcript_id 
        
    def __str__(self) -> str:
        return f"Exon(exon_id={self.exon_id}, start={self.start}, end={self.end}, transcript_id={self.transcript_id})"