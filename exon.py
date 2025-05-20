class Exon:
    """Class representing an exon within a transcript."""
    def __init__(self, exon_id: str, start: int, end: int, transcript_id: str) -> None:
        self.exon_id: str = exon_id
        self.start: int = start  # 1-based genomic coordinates
        self.end: int = end      # 1-based genomic coordinates
        self.transcript_id: str = transcript_id 