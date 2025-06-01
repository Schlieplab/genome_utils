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
                        pass # Ignore malformed attributes
                
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
                
                elif feature_type == 'transcript':
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
                            support_level=support_level,
                        )
                        self._transcripts[transcript_id] = transcript
                        
                        # Add to gene if it exists
                        if gene_id in self._genes:
                            self._genes[gene_id].add_transcript(transcript)
                        
                
                elif feature_type == 'exon':
                    
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

    def extract_genome_premrna_sequences(self, output_path: str, force: bool = False, 
                                 exclude_genes: Optional[Union[str, List[str]]] = None) -> str:
        """
        Extract pre-mRNA sequences for each gene from the primary assembly and save them in FASTA format.
        Uses chunked reading to process one chromosome at a time for memory efficiency.
        
        Args:
            output_path (str): Path to save the FASTA file containing pre-mRNA sequences.
                If exclude_genes is provided, the filename will be modified to reflect the exclusions.
            force (bool): If True, overwrite existing file. If False, skip if file exists.
            exclude_genes (Optional[Union[str, List[str]]]): Gene ID(s) to exclude from extraction.
                Can be a single gene ID string or a list of gene IDs.
        """
        if not self.primary_assembly_path or not os.path.exists(self.primary_assembly_path):
            raise FileNotFoundError(f"Primary assembly file not found: {self.primary_assembly_path}")
            
        if not self._indexed:
            self.index()
            
        # Process exclude_genes parameter and modify output path accordingly
        exclude_set = set()
        modified_output_path = output_path
        
        if exclude_genes is not None:
            if isinstance(exclude_genes, str):
                exclude_set.add(exclude_genes)
                # For single gene exclusion, include the gene ID in the filename
                base_path = output_path.replace('.fa.gz', '')
                modified_output_path = f"{base_path}.exclude_{exclude_genes}.fa"
            else:
                exclude_set.update(exclude_genes)
                if len(exclude_genes) > 0:
                    # For multiple exclusions, include the count in the filename
                    base_path = output_path.replace('.fa.gz', '')
                    modified_output_path = f"{base_path}.exclude_{len(exclude_genes)}_genes.fa"
        
        # If no exclusions, ensure .all. is in the filename
        if not exclude_set:
            if '.all.' not in modified_output_path:
                modified_output_path = modified_output_path.replace('.fa.gz', '.all.fa')
        
        # Check if file exists and handle accordingly
        if os.path.exists(modified_output_path):
            if not force and not (isinstance(exclude_genes, list) and len(exclude_genes) > 1):
                logging.info(f"Pre-mRNA sequences file already exists at {modified_output_path}. Skipping extraction.")
                return modified_output_path
            else:
                if isinstance(exclude_genes, list) and len(exclude_genes) > 1:
                    logging.info(f"Multiple genes excluded - overwriting existing file at {modified_output_path}")
                else:
                    logging.info(f"Overwriting existing pre-mRNA sequences file at {modified_output_path}")
            
        # Create output directory if it doesn't exist
        os.makedirs(os.path.dirname(modified_output_path), exist_ok=True)
        
        # Group genes by chromosome for efficient processing, excluding specified genes
        genes_by_chromosome: Dict[str, List[Gene]] = {}
        excluded_count = 0
        for gene in self.genes:
            if gene.gene_id in exclude_set:
                excluded_count += 1
                continue
            if gene.chromosome not in genes_by_chromosome:
                genes_by_chromosome[gene.chromosome] = []
            genes_by_chromosome[gene.chromosome].append(gene)
        
        if excluded_count > 0:
            logging.info(f"Excluding {excluded_count} genes from pre-mRNA extraction")
        
        # Determine if file is gzipped
        is_gzipped: bool = self.primary_assembly_path.endswith('.gz')
        open_func: Any = gzip.open if is_gzipped else open
        
        # Process one chromosome at a time
        records = []
        with open_func(self.primary_assembly_path, 'rt') as fasta_file:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                chromosome = record.id
                if chromosome not in genes_by_chromosome:
                    continue
                
                # Get all genes for this chromosome
                chromosome_genes = genes_by_chromosome[chromosome]
                chromosome_seq = str(record.seq)
                
                # Process all genes on this chromosome
                for gene in chromosome_genes:
                    # Convert to 0-based indexing for Python string operations
                    start_idx = gene.start - 1
                    end_idx = gene.end
                    
                    # Check bounds
                    if start_idx < 0 or end_idx > len(chromosome_seq):
                        logging.warning(f"Gene {gene.gene_id} coordinates out of bounds for chromosome {chromosome}")
                        continue
                    
                    # Extract sequence
                    sequence = chromosome_seq[start_idx:end_idx]
                    
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
                
                # Write records in batches to avoid memory buildup
                if len(records) >= 1000:
                    with open(modified_output_path, 'a') as output_handle:
                        SeqIO.write(records, output_handle, "fasta")
                    records = []
        
        # Write any remaining records
        if records:
            with open(modified_output_path, 'a') as output_handle:
                SeqIO.write(records, output_handle, "fasta")
            
        logging.info(f"Extracted pre-mRNA sequences to {modified_output_path}")
        return modified_output_path

    def extract_premrna_sequences_per_gene(self, gene_ids: Union[str, List[str]], 
                                      output_path: Optional[str] = None) -> Dict[str, str]:
        """
        Extract pre-mRNA sequences for specific genes and optionally save them to a FASTA file.
        
        Args:
            gene_ids (Union[str, List[str]]): Gene ID(s) to extract sequences for.
                Can be a single gene ID string or a list of gene IDs.
            output_path (Optional[str]): If provided, save the sequences to this FASTA file.
                If not provided, only return the sequences without saving to file.
                
        Returns:
            Dict[str, str]: Dictionary mapping gene IDs to their sequences.
        """
        if not self.primary_assembly_path or not os.path.exists(self.primary_assembly_path):
            raise FileNotFoundError(f"Primary assembly file not found: {self.primary_assembly_path}")
            
        if not self._indexed:
            self.index()
            
        # Convert single gene ID to list
        if isinstance(gene_ids, str):
            gene_ids = [gene_ids]
            
        # Get genes and group by chromosome
        genes_by_chromosome: Dict[str, List[Gene]] = {}
        gene_map: Dict[str, Gene] = {}
        
        for gene_id in gene_ids:
            try:
                gene = self.gene_by_id(gene_id)
                if gene.chromosome not in genes_by_chromosome:
                    genes_by_chromosome[gene.chromosome] = []
                genes_by_chromosome[gene.chromosome].append(gene)
                gene_map[gene_id] = gene
            except ValueError:
                logging.warning(f"Gene not found: {gene_id}")
                continue
        
        if not genes_by_chromosome:
            logging.warning("No valid genes found to extract")
            return {}
            
        # Determine if file is gzipped
        is_gzipped: bool = self.primary_assembly_path.endswith('.gz')
        open_func: Any = gzip.open if is_gzipped else open
        
        # Process one chromosome at a time
        sequences: Dict[str, str] = {}
        records = []
        
        with open_func(self.primary_assembly_path, 'rt') as fasta_file:
            for record in SeqIO.parse(fasta_file, 'fasta'):
                chromosome = record.id
                if chromosome not in genes_by_chromosome:
                    continue
                
                # Get all genes for this chromosome
                chromosome_genes = genes_by_chromosome[chromosome]
                chromosome_seq = str(record.seq)
                
                # Process all genes on this chromosome
                for gene in chromosome_genes:
                    # Convert to 0-based indexing for Python string operations
                    start_idx = gene.start - 1
                    end_idx = gene.end
                    
                    # Check bounds
                    if start_idx < 0 or end_idx > len(chromosome_seq):
                        logging.warning(f"Gene {gene.gene_id} coordinates out of bounds for chromosome {chromosome}")
                        continue
                    
                    # Extract sequence
                    sequence = chromosome_seq[start_idx:end_idx]
                    
                    # If gene is on reverse strand, reverse complement the sequence
                    if gene.strand == '-':
                        sequence = str(Seq(sequence).reverse_complement())
                    
                    # Store sequence in gene object and return dictionary
                    gene.pre_mrna_sequence = sequence
                    sequences[gene.gene_id] = sequence
                    
                    # Create FASTA record if output path is provided
                    if output_path:
                        output_path = output_path.replace('.fa.gz', '.fa')
                        record = SeqRecord(
                            seq=Seq(sequence),
                            id=f"{gene.gene_id}|{gene.gene_name}",
                            description=f"pre-mRNA sequence for gene {gene.gene_name} ({gene.gene_id}) on {gene.chromosome}:{gene.start}-{gene.end}:{gene.strand}"
                        )
                        records.append(record)
        
        # Write records if output path is provided
        if output_path and records:
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
            with open(output_path, 'w') as output_handle:
                SeqIO.write(records, output_handle, "fasta")
            logging.info(f"Extracted sequences for {len(sequences)} genes to {output_path}")
        
        return sequences 