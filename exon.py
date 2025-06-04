class Exon:
    """Class representing an exon within a transcript."""
    def __init__(self, exon_id: str, start: int, end: int, transcript_id: str) -> None:
        self.exon_id: str = exon_id
        self.start: int = start 
        self.end: int = end
        self.transcript_id: str = transcript_id 
        
    def __str__(self) -> str:
        return f"Exon(exon_id={self.exon_id}, start={self.start}, end={self.end}, transcript_id={self.transcript_id})"

    def to_gtf_entry(self, chromosome: str, strand: str, gene_id: str, source: str = "custom") -> str:
        """
        Returns a GTF formatted string for this exon.
        Requires chromosome, strand, and gene_id from the parent transcript.
        'source' is customizable. Attributes follow GTF specification.
        """
        attributes = [
            f'gene_id "{gene_id}"',
            f'transcript_id "{self.transcript_id}"',
            f'exon_id "{self.exon_id}"'
        ]
        
        # Standard GTF attributes for exons can also include exon_number.
        # This would require exon numbering logic within the Transcript class when adding exons.
        # For now, keeping it to gene_id, transcript_id, and exon_id.
        # if hasattr(self, 'exon_number'):
        #     attributes.append(f'exon_number "{self.exon_number}"')

        attributes_str = "; ".join(attributes) + ";"
        
        return f"{chromosome}\t{source}\texon\t{self.start}\t{self.end}\t.\t{strand}\t.\t{attributes_str}"