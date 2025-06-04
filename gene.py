from typing import List, Optional, Tuple
from Bio.Seq import Seq
from .transcript import Transcript
import os
import gzip
import logging


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
    
    def to_gtf_entry(self, source: str = "custom") -> str:
        """
        Returns a GTF formatted string for this gene.
        The 'source' field in GTF is customizable.
        Attributes are formatted according to GTF specification.
        Example: gene_id "ENSG00000223972"; gene_version "5"; ...
        """
        attributes = [
            f'gene_id "{self.gene_id}"',
            f'gene_name "{self.gene_name}"'
        ]
        if self.biotype:
            attributes.append(f'gene_biotype "{self.biotype}"')
        
        attributes_str = "; ".join(attributes) + ";" 
        
        return f"{self.chromosome}\t{source}\tgene\t{self.start}\t{self.end}\t.\t{self.strand}\t.\t{attributes_str}"
    
    def export_gene_data(self, output_dir: str, source_tag: str = "gene_export") -> Tuple[str, str]:
        """
        Exports this gene's data to FASTA (cDNA) and GTF files.

        FASTA Naming: {gene_id}.cdna.fa.gz
        GTF Naming:   {gene_id}.gtf.gz

        Args:
            output_dir: The directory to save the exported files.
            source_tag: The source tag to use in the GTF file (column 2).

        Returns:
            A tuple containing the paths to the exported FASTA and GTF files.
        """
        os.makedirs(output_dir, exist_ok=True)

        # --- Prepare FASTA file --- 
        fasta_filename = f"{self.gene_id}.cdna.fa.gz"
        output_fasta_path = os.path.join(output_dir, fasta_filename)

        logging.info(f"Exporting cDNA FASTA for gene {self.gene_id} to: {output_fasta_path}")
        transcripts_written_count = 0
        with gzip.open(output_fasta_path, 'wt') as f_fasta_out:
            for transcript in self.transcripts:
                if transcript.sequence:
                    f_fasta_out.write(f">{transcript.transcript_id} gene_id={self.gene_id}\\n{transcript.sequence}\\n")
                    transcripts_written_count += 1
                else:
                    logging.debug(f"Transcript {transcript.transcript_id} for gene {self.gene_id} has no sequence. Not written to FASTA.")
        logging.info(f"Wrote {transcripts_written_count} transcripts for gene {self.gene_id} to {output_fasta_path}")

        # --- Prepare GTF file --- 
        gtf_filename = f"{self.gene_id}.gtf.gz"
        output_gtf_path = os.path.join(output_dir, gtf_filename)
        
        logging.info(f"Exporting GTF for gene {self.gene_id} to: {output_gtf_path}")
        gtf_entries_written = 0
        with gzip.open(output_gtf_path, 'wt') as f_gtf_out:
            f_gtf_out.write(self.to_gtf_entry(source=source_tag) + "\\n")
            gtf_entries_written += 1

            for transcript in self.transcripts:
                f_gtf_out.write(transcript.to_gtf_entry(source=source_tag) + "\\n")
                gtf_entries_written += 1
                
                sorted_exons = sorted(transcript.exons, key=lambda ex: ex.start)
                for exon in sorted_exons:
                    f_gtf_out.write(exon.to_gtf_entry(chromosome=transcript.chromosome, 
                                                          strand=transcript.strand, 
                                                          gene_id=self.gene_id, 
                                                          source=source_tag) + "\\n")
                    gtf_entries_written += 1
        logging.info(f"Wrote {gtf_entries_written} GTF entries for gene {self.gene_id} to {output_gtf_path}")
        
        return output_fasta_path, output_gtf_path
    
