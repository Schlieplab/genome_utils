from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import gzip
import re
from typing import Dict, List, Optional, Union, Any
import logging

from .gene import Gene
from .transcript import Transcript
from .exon import Exon

class Genome:
    def __init__(self, reference_name: str, annotation_version: Optional[str] = None,
                 gtf_path: Optional[str] = None, 
                 transcript_fasta_paths: Optional[Union[str, List[str]]] = None,
                 primary_assembly_path: Optional[str] = None) -> None:
        """
        A Genome class using Biopython to load and manipulate genome sequences.

        Parameters:
        - reference_name: E.g. 'GRCm38'
        - annotation_version: Annotation version/release (e.g. 113)
        - gtf_path: Path to the GTF file
        - transcript_fasta_paths: Path(s) to transcript FASTA file(s)
        - primary_assembly_path: Path to the primary assembly FASTA file
        """
        
        self.reference_name: str = reference_name
        self.annotation_version: Optional[str] = annotation_version
        self.gtf_path: Optional[str] = gtf_path
        self.primary_assembly_path: Optional[str] = primary_assembly_path
        
        # Handle either a single path or a list of paths
        if transcript_fasta_paths is not None:
            if isinstance(transcript_fasta_paths, list):
                self.transcript_fasta_paths: List[str] = transcript_fasta_paths
            else:
                self.transcript_fasta_paths: List[str] = [transcript_fasta_paths]
        else:
            self.transcript_fasta_paths: List[str] = []
        
        # Initialize data structures
        self._genes: Dict[str, Gene] = {}  # gene_id -> Gene
        self._transcripts: Dict[str, Transcript] = {}  # transcript_id -> Transcript
        self._transcript_sequences: Dict[str, str] = {}  # transcript_id -> sequence
        self._exons: Dict[str, Exon] = {}  # exon_id -> Exon
        self._indexed: bool = False
    
    def index(self, overwrite: bool = False) -> None:
        """
        Parse annotation files and build indices for genes, transcripts, exons, and sequences.
        
        Args:
            overwrite: If True, rebuild indices even if they already exist
        """
        if self._indexed and not overwrite:
            return
        
        # Parse GTF file to build gene/transcript/exon data structures
        self._parse_gtf()
        
        # Parse transcript FASTA files to get sequences
        self._parse_transcript_fasta()
        
        self._indexed = True
        
    def _parse_gtf(self) -> None:
        """Parse GTF file to extract gene, transcript, and exon information."""
        if not self.gtf_path or not os.path.exists(self.gtf_path):
            raise FileNotFoundError(f"GTF file not found: {self.gtf_path}")
        
        # Determine if file is gzipped
        is_gzipped: bool = self.gtf_path.endswith('.gz')
        open_func: Any = gzip.open if is_gzipped else open
        
        # Track current gene and transcript when parsing
        current_gene: Optional[Gene] = None
        current_transcript: Optional[Transcript] = None
        
        with open_func(self.gtf_path, 'rt') as gtf:
            for line in gtf:
                # Skip comments/headers
                if line.startswith('#'):
                    continue
                
                fields: List[str] = line.strip().split('\t')
                if len(fields) < 9:  # GTF has at least 9 fields
                    continue
                
                # Extract feature fields
                seqname, source, feature_type, start, end, score, strand, frame, attributes = fields
                
                # Skip if not gene, transcript, or exon
                if feature_type not in ['gene', 'transcript', 'exon']:
                    continue
                
                # Parse attributes
                attr_dict: Dict[str, str] = {}
                for attr in attributes.split(';'):
                    attr = attr.strip()
                    if not attr:
                        continue
                    try:
                        key, value = attr.split(' ', 1)
                        attr_dict[key] = value.strip('"')
                    except ValueError:
                        pass
                
                # Process different feature types
                if feature_type == 'gene':
                    gene_id: Optional[str] = attr_dict.get('gene_id')
                    gene_name: str = attr_dict.get('gene_name', gene_id)
                    biotype: Optional[str] = attr_dict.get('gene_biotype', attr_dict.get('biotype', None))
                    
                    if gene_id:
                        gene = Gene(
                            gene_id=gene_id,
                            gene_name=gene_name,
                            chromosome=seqname,
                            start=int(start),
                            end=int(end),
                            strand=strand,
                            biotype=biotype,
                        )
                        self._genes[gene_id] = gene
                        current_gene = gene
                
                elif feature_type == 'transcript' and current_gene:
                    transcript_id: Optional[str] = attr_dict.get('transcript_id')
                    gene_id: Optional[str] = attr_dict.get('gene_id')
                    biotype: Optional[str] = attr_dict.get('transcript_biotype', attr_dict.get('biotype', None))
                    
                    # Extract support level - could be 'transcript_support_level' or 'tsl'
                    support_level_str: Optional[str] = attr_dict.get('transcript_support_level', 
                                                                    attr_dict.get('tsl', None))
                    support_level: Optional[int] = None
                    
                    # Parse support level if it exists
                    if support_level_str:
                        # Sometimes formatted as "1 (assigned)", so extract just the number
                        match = re.match(r'^(\d+)', support_level_str)
                        if match:
                            try:
                                support_level = int(match.group(1))
                            except ValueError:
                                pass  # Keep as None if conversion fails
                    
                    if transcript_id and gene_id:
                        transcript = Transcript(
                            transcript_id=transcript_id,
                            gene_id=gene_id,
                            chromosome=seqname,
                            start=int(start),
                            end=int(end),
                            strand=strand,
                            biotype=biotype,
                            support_level=support_level
                        )
                        self._transcripts[transcript_id] = transcript
                        
                        # Add to gene if it exists
                        if gene_id in self._genes:
                            self._genes[gene_id].add_transcript(transcript)
                        
                        current_transcript = transcript
                
                elif feature_type == 'exon' and current_transcript:
                    exon_id: str = attr_dict.get('exon_id', f"{seqname}:{start}-{end}")
                    transcript_id: Optional[str] = attr_dict.get('transcript_id')
                    
                    if exon_id and transcript_id and transcript_id in self._transcripts:
                        exon = Exon(
                            exon_id=exon_id,
                            start=int(start),
                            end=int(end),
                            transcript_id=transcript_id
                        )
                        self._exons[exon_id] = exon
                        self._transcripts[transcript_id].add_exon(exon)
    
    def _parse_transcript_fasta(self) -> None:
        """Parse transcript FASTA files to extract sequences."""
        for fasta_path in self.transcript_fasta_paths:
            if not os.path.exists(fasta_path):
                logging.warning(f"Warning: FASTA file not found: {fasta_path}")
                continue
            
            # Determine if file is gzipped
            is_gzipped: bool = fasta_path.endswith('.gz')
            open_func: Any = gzip.open if is_gzipped else open
            
            # Parse FASTA file
            with open_func(fasta_path, 'rt') as fasta_file:
                for record in SeqIO.parse(fasta_file, 'fasta'):
                    # Extract transcript ID from header (may need adjustment based on format)
                    header_parts: List[str] = record.id.split('|')
                    if len(header_parts) > 1:
                        # Try to extract a clean transcript ID
                        transcript_id: str = header_parts[0].split('.')[0]
                    else:
                        transcript_id: str = record.id.split('.')[0]
                    
                    # Store sequence if we have this transcript
                    if transcript_id in self._transcripts:
                        sequence: str = str(record.seq)
                        self._transcript_sequences[transcript_id] = sequence
                        self._transcripts[transcript_id].sequence = sequence
    
    def gene_by_id(self, gene_id: str) -> Gene:
        """Get a gene by its ID."""
        if not self._indexed:
            self.index()
        
        if gene_id in self._genes:
            return self._genes[gene_id]
        
        raise ValueError(f"Gene not found with ID: {gene_id}")
    
    def transcript_by_id(self, transcript_id: str) -> Transcript:
        """Get a transcript by its ID."""
        if not self._indexed:
            self.index()
        
        if transcript_id in self._transcripts:
            return self._transcripts[transcript_id]
        
        raise ValueError(f"Transcript not found with ID: {transcript_id}")
    
    def transcripts(self) -> List[Transcript]:
        """Get all transcripts."""
        if not self._indexed:
            self.index()
        
        return list(self._transcripts.values())
    
    @property
    def genes(self) -> List[Gene]:
        """Get all genes."""
        if not self._indexed:
            self.index()
        
        return list(self._genes.values())
    
    def get_sequence_for_transcript_id(self, transcript_id: str) -> str:
        """Get the sequence for a transcript by ID."""
        if not self._indexed:
            self.index()
        
        if transcript_id in self._transcript_sequences:
            return self._transcript_sequences[transcript_id]
        
        raise ValueError(f"No sequence found for transcript: {transcript_id}")
    
    def get_transcript_subsequence(self, transcript_id: str, position: int, length: int) -> Optional[str]:
        """
        Get a subsequence from a transcript by ID.
        
        Args:
            transcript_id (str): The ID of the transcript
            position (int): 1-based position within the transcript
            length (int): Length of the subsequence to extract
            
        Returns:
            Optional[str]: The requested subsequence, or None if not available
        """
        if not self._indexed:
            self.index()
            
        # Get the base ID without version if present
        base_id = transcript_id.split('.')[0]
        
        try:
            transcript = self.transcript_by_id(base_id)
            return transcript.get_subsequence(position, length)
        except ValueError:
            return None
    
    def get_chromosomal_positions(
        self, transcript_id: str, positions: List[int], window_length: int
    ) -> List[Optional[str]]:
        """
        Get chromosomal positions for multiple positions within a transcript.
        
        Args:
            transcript_id (str): The ID of the transcript
            positions (List[int]): List of 1-based positions within the transcript
            window_length (int): Length of the window at each position
            
        Returns:
            List[Optional[str]]: List of chromosomal coordinates in format "chrom:start-end:strand"
        """
        if not self._indexed:
            self.index()
            
        # Get the base ID without version if present
        base_id = transcript_id.split('.')[0]
        
        try:
            transcript = self.transcript_by_id(base_id)
            return transcript.get_chromosomal_positions(positions, window_length)
        except ValueError:
            # Return None for each position if transcript not found
            return [None] * len(positions)
    
    def get_chromosomal_position(
        self, transcript_id: str, position: int, window_length: int
    ) -> Optional[str]:
        """
        Get chromosomal position for a single position within a transcript.
        
        Args:
            transcript_id (str): The ID of the transcript
            position (int): 1-based position within the transcript
            window_length (int): Length of the window starting at the position
            
        Returns:
            Optional[str]: Chromosomal coordinates in format "chrom:start-end:strand"
        """
        results = self.get_chromosomal_positions(
            transcript_id=transcript_id,
            positions=[position],
            window_length=window_length
        )
        return results[0] if results else None
    
    def get_exon_at_transcript_position(self, transcript_id: str, position: int) -> Optional[Exon]:
        """
        Get the exon containing a specific position in a transcript.
        
        Args:
            transcript_id (str): The ID of the transcript
            position (int): 1-based position within the transcript
            
        Returns:
            Optional[Exon]: The exon containing the position, or None if not found
        """
        if not self._indexed:
            self.index()
            
        # Get the base ID without version if present
        base_id = transcript_id.split('.')[0]
        
        try:
            transcript = self.transcript_by_id(base_id)
            return transcript.get_exon_by_position(position)
        except ValueError:
            return None
    
    def get_sequence_from_primary_assembly(self, chromosome: str, start: int, end: int) -> Optional[str]:
        """
        Get a sequence from the primary assembly by loading only the needed chromosome.
        
        Args:
            chromosome (str): Chromosome name
            start (int): 1-based start position
            end (int): 1-based end position
            
        Returns:
            Optional[str]: The requested sequence, or None if not available
        """
        if not self.primary_assembly_path or not os.path.exists(self.primary_assembly_path):
            return None
            
        # Determine if file is gzipped
        is_gzipped: bool = self.primary_assembly_path.endswith('.gz')
        open_func: Any = gzip.open if is_gzipped else open
        
        # Convert to 0-based indexing for Python string operations
        start_idx = start - 1
        end_idx = end
        
        with open_func(self.primary_assembly_path, 'rt') as fasta_file:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                if record.id == chromosome:
                    # Check bounds
                    if start_idx < 0 or end_idx > len(record.seq):
                        return None
                    return str(record.seq[start_idx:end_idx])
        
        return None  # Chromosome not found 

    def extract_premrna_sequences(self, output_path: str) -> None:
        """
        Extract pre-mRNA sequences for each gene from the primary assembly and save them in FASTA format.
        
        Args:
            output_path (str): Path to save the FASTA file containing pre-mRNA sequences
        """
        if not self.primary_assembly_path or not os.path.exists(self.primary_assembly_path):
            raise FileNotFoundError(f"Primary assembly file not found: {self.primary_assembly_path}")
            
        if not self._indexed:
            self.index()
            
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        
        # Create FASTA records for each gene
        records = []
        for gene in self.genes:
            # Get the sequence from the primary assembly
            sequence = self.get_sequence_from_primary_assembly(
                chromosome=gene.chromosome,
                start=gene.start,
                end=gene.end
            )
            
            if sequence is not None:
                # If gene is on reverse strand, reverse complement the sequence
                if gene.strand == '-':
                    sequence = str(Seq(sequence).reverse_complement())
                
                # Create FASTA record
                record = SeqRecord(
                    seq=Seq(sequence),
                    id=f"{gene.gene_id}|{gene.gene_name}",
                    description=f"pre-mRNA sequence for gene {gene.gene_name} ({gene.gene_id}) on {gene.chromosome}:{gene.start}-{gene.end}:{gene.strand}"
                )
                records.append(record)
        
        # Write records to FASTA file
        with open(output_path, 'w') as output_handle:
            SeqIO.write(records, output_handle, "fasta")
            
        logging.info(f"Extracted {len(records)} pre-mRNA sequences to {output_path}") 