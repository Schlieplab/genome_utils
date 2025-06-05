from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import gzip
import re
from typing import Dict, List, Optional, Union, Any, Iterator, Tuple
import logging
from datetime import datetime
from .gene import Gene
from .transcript import Transcript
from .exon import Exon

class Genome:
    """
    A Genome class using Biopython to load and manipulate genome sequences from Ensembl.
    """
    def __init__(self, 
                 reference_name: str, 
                 e_release: Optional[str] = None,
                 gtf_path: Optional[str] = None, 
                 transcript_fasta_paths: Optional[Union[str, List[str]]] = None,
                 primary_assembly_path: Optional[str] = None,
                 tsl_to_keep: Optional[List[Optional[int]]] = None,
                 biotype_to_keep: Optional[List[str]] = None,
                 verbose: bool = False
                 ) -> None:
        """
        Initialize the Genome object.

        Parameters:
        - reference_name: E.g. 'GRCm38' or 'GRCh38'
        - e_release: The Ensembl release version. 
        - gtf_path: Path to the GTF file
        - transcript_fasta_paths: Path(s) to transcript FASTA file(s)
        - primary_assembly_path: Path to the primary assembly FASTA file
        - tsl_to_keep: Optional list of TSL values to keep (e.g., [1, 2, None]). 
                       Transcripts not matching these TSLs will be excluded.
                       If None, all transcripts are loaded.
        - biotype_to_keep: Optional list of gene biotypes to keep (e.g., ['protein_coding', 'lncRNA']).
                       If None, all genes are loaded.
        - verbose: If True, print verbose output.
        """
        
        self.reference_name: str = reference_name
        self.e_release: Optional[str] = e_release
        self.gtf_path: Optional[str] = gtf_path
        self.primary_assembly_path: Optional[str] = primary_assembly_path
        self.biotype_to_keep: Optional[List[str]] = biotype_to_keep
        self.tsl_to_keep: Optional[List[Optional[int]]] = tsl_to_keep
        self.verbose: bool = verbose
        
        self._genes: Dict[str, Gene] = {}  # gene_id -> Gene
        self._transcripts: Dict[str, Transcript] = {}  # transcript_id -> Transcript
        self._exons: Dict[str, Exon] = {}  # exon_id -> Exon
        self._indexed: bool = False
        
        # Handle either a single path or a list of paths
        if transcript_fasta_paths is not None:
            if isinstance(transcript_fasta_paths, list):
                self.transcript_fasta_paths: List[str] = transcript_fasta_paths
            else:
                self.transcript_fasta_paths: List[str] = [transcript_fasta_paths]
        else:
            if self.verbose:
                logging.warning("No transcript FASTA paths provided. No transcripts will be loaded.")
            self.transcript_fasta_paths: List[str] = []

    def __str__(self) -> str:
        return f"Genome(reference_name={self.reference_name}, e_release={self.e_release}, gtf_path={self.gtf_path}, primary_assembly_path={self.primary_assembly_path}, tsl_to_keep={self.tsl_to_keep}, biotype_to_keep={self.biotype_to_keep}, verbose={self.verbose})"
    
    def index(self, overwrite: bool = False) -> None:
        """
        Parse annotation files and build indices for genes, transcripts, exons, and sequences.
        
        Parameters:
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
                (seqname, source, 
                 feature_type, start, 
                 end, score, 
                 strand, frame, 
                 attributes) = fields
                
                if feature_type not in ['gene', 'transcript', 'exon']:
                    continue
                
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
                    
                    if self.biotype_to_keep and biotype not in self.biotype_to_keep:
                        if self.verbose:
                            logging.debug(f"Skipping non-protein-coding gene {gene_id} (biotype: {biotype}) during GTF parsing because biotype_to_keep is True.")
                        continue
                    
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
                    
                    if self.tsl_to_keep is not None and support_level not in self.tsl_to_keep:
                        logging.debug(f"Skipping transcript {transcript_id} in gene {gene_id} due to TSL filtering (TSL: {support_level}, Allowed: {self.tsl_to_keep})")
                        continue
                    
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
                        
                        if gene_id in self._genes:
                            self._genes[gene_id].add_transcript(transcript)
                        
                
                elif feature_type == 'exon':
                    
                    exon_id: str = attr_dict.get('exon_id', f"{seqname}:{start}-{end}:{strand}")
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
    
    @property
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
    
    def get_sequence_from_primary_assembly(self, chromosome: str, start: int, end: int) -> Optional[str]:
        """
        Get a sequence from the primary assembly by loading only the needed chromosome.
        
        Parameters:
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

    def export_pre_mrna_sequences(self, output_dir: str, exclude_genes: Optional[Union[str, List[str]]] = None, force_overwrite: bool = False) -> Optional[str]:
        """
        Exports pre-mRNA sequences for genes to a FASTA file, processed chromosome by chromosome.

        Parameters:
            output_dir: The directory to save the exported FASTA file.
            exclude_genes: Optional gene ID(s) to exclude.

        Returns:
            Optional[str]: The path to the exported FASTA file, or None if an error occurs or no sequences are written.
        
        Raises:
            FileNotFoundError: If the primary assembly file is not found.
        """
        if not self.primary_assembly_path or not os.path.exists(self.primary_assembly_path):
            raise FileNotFoundError(f"Primary assembly file not found: {self.primary_assembly_path}")

        if not self._indexed:
            self.index()

        os.makedirs(output_dir, exist_ok=True)
        e_release_str = self.e_release if self.e_release else "no_release"

        exclude_set = set()
        exclusion_suffix = ""
        log_message_exclusion_details = " (all genes)"

        if exclude_genes:
            if isinstance(exclude_genes, str):
                exclude_set.add(exclude_genes)
                exclusion_suffix = f".excluding_{exclude_genes}"
                log_message_exclusion_details = f" (excluding gene '{exclude_genes}')"
            else:
                exclude_set.update(exclude_genes)
                if len(exclude_genes) == 1:
                    exclusion_suffix = f".excluding_{exclude_genes[0]}"
                    log_message_exclusion_details = f" (excluding gene '{exclude_genes[0]}')"
                else:
                    exclusion_suffix = f".excluding_{len(exclude_genes)}_genes"
                    log_message_exclusion_details = f" (excluding {len(exclude_genes)} genes)"
        
        if self.biotype_to_keep:
            if len(self.biotype_to_keep) == 1:
                biotype_suffix = f".{self.biotype_to_keep[0]}"
            else:
                biotype_suffix = f".{'_'.join(self.biotype_to_keep)}"
        else:
            biotype_suffix = ""
        
        if biotype_suffix == "" and exclusion_suffix == "":
            fasta_filename = f"{self.reference_name}.{e_release_str}.premrna.all.fa"
        else:
            fasta_filename = f"{self.reference_name}.{e_release_str}.premrna{biotype_suffix}{exclusion_suffix}.fa"
        output_fasta_path = os.path.join(output_dir, fasta_filename)

        if os.path.exists(output_fasta_path) and not force_overwrite:
            if self.verbose:
                logging.info(f"Using existing FASTA file: {output_fasta_path}. Export will be skipped.")
            else:
                logging.debug(f"Using existing FASTA file: {output_fasta_path}. Export will be skipped.")
            return output_fasta_path
        elif os.path.exists(output_fasta_path) and force_overwrite:
            if self.verbose:
                logging.info(f"Overwriting existing file: {output_fasta_path}")
            else:
                logging.debug(f"Overwriting existing file: {output_fasta_path}")
            os.remove(output_fasta_path)
        
        if self.verbose:    
            logging.info(f"Exporting pre-mRNA sequences to: {output_fasta_path}{log_message_exclusion_details}")
        else:
            logging.debug(f"Exporting pre-mRNA sequences to: {output_fasta_path}{log_message_exclusion_details}")

        genes_by_chromosome: Dict[str, List[Gene]] = {}
        processed_gene_ids = set()
        all_genes_list = self.genes
        
        # Sort genes first by chromosome, then by start position for ordered processing
        sorted_genes = sorted(all_genes_list, key=lambda g: (g.chromosome, g.start))

        excluded_count = 0
        genes_to_process_count = 0
        for gene in sorted_genes:
            if gene.gene_id in exclude_set:
                excluded_count += 1
                continue
            
            if gene.gene_id in processed_gene_ids:
                continue
            processed_gene_ids.add(gene.gene_id)

            if gene.chromosome not in genes_by_chromosome:
                genes_by_chromosome[gene.chromosome] = []
            genes_by_chromosome[gene.chromosome].append(gene)
            genes_to_process_count +=1

        if excluded_count > 0:
            if self.verbose:
                logging.info(f"Excluded {excluded_count} genes from pre-mRNA sequence export.")
            else:
                logging.debug(f"Excluded {excluded_count} genes from pre-mRNA sequence export.")

        if genes_to_process_count == 0:
            logging.warning(f"No genes found to process for pre-mRNA export after exclusions. Output file '{output_fasta_path}' will be empty.")
            

        is_gzipped_assembly: bool = self.primary_assembly_path.endswith('.gz')
        open_func_assembly: Any = gzip.open if is_gzipped_assembly else open
        
        sequences_written_count = 0
        records_batch = []

        try:
            with open(output_fasta_path, 'wt') as f_fasta_out:
                with open_func_assembly(self.primary_assembly_path, 'rt') as fasta_file_handle:
                    for assembly_record in SeqIO.parse(fasta_file_handle, 'fasta'):
                        chromosome_id = assembly_record.id
                        if chromosome_id not in genes_by_chromosome:
                            continue

                        chromosome_genes_to_process = genes_by_chromosome[chromosome_id]
                        chromosome_seq_str = str(assembly_record.seq)

                        for gene in chromosome_genes_to_process: # These are already filtered from exclude_set
                            start_idx = gene.start - 1 # Convert to 0-based
                            end_idx = gene.end

                            if not (0 <= start_idx < len(chromosome_seq_str) and 0 <= end_idx <= len(chromosome_seq_str) and start_idx < end_idx):
                                logging.warning(
                                    f"Gene {gene.gene_id} coordinates ({gene.start}-{gene.end}) "
                                    f"out of bounds for chromosome {chromosome_id} (length: {len(chromosome_seq_str)}). Skipping."
                                )
                                continue
                            
                            sequence_str = chromosome_seq_str[start_idx:end_idx]

                            if gene.strand == '-':
                                sequence_str = str(Seq(sequence_str).reverse_complement())
                            
                            seq_record = SeqRecord(
                                seq=Seq(sequence_str),
                                id=f"{gene.gene_id}|{gene.gene_name}", 
                                description=f"pre-mRNA sequence for gene {gene.gene_name} ({gene.gene_id}) on {gene.chromosome}:{gene.start}-{gene.end}:{gene.strand}"
                            )
                            records_batch.append(seq_record)
                            sequences_written_count += 1
                        
                        if len(records_batch) > 0:
                            SeqIO.write(records_batch, f_fasta_out, "fasta")
                            records_batch = [] 
                

            if sequences_written_count > 0:
                if self.verbose:
                    logging.info(f"Successfully exported {sequences_written_count} pre-mRNA sequences to {output_fasta_path}")
                else:
                    logging.debug(f"Successfully exported {sequences_written_count} pre-mRNA sequences to {output_fasta_path}")
                return output_fasta_path
            else:
                logging.warning(f"No pre-mRNA sequences were written to {output_fasta_path}. This could be due to all genes being excluded or other issues.")
                
                if os.path.exists(output_fasta_path) and os.path.getsize(output_fasta_path) == 0:
                    try:
                        os.remove(output_fasta_path)
                        logging.info(f"Removed empty pre-mRNA FASTA file: {output_fasta_path}")
                    except OSError as oe:
                        logging.error(f"Error removing empty pre-mRNA FASTA file {output_fasta_path}: {oe}")
                return None

        except Exception as e:
            logging.error(f"An error occurred during pre-mRNA export to {output_fasta_path}: {e}")
            
            if os.path.exists(output_fasta_path):
                try:
                    os.remove(output_fasta_path)
                    logging.info(f"Removed potentially incomplete pre-mRNA FASTA file due to error: {output_fasta_path}")
                except OSError as oe:
                    logging.error(f"Error removing incomplete pre-mRNA FASTA file {output_fasta_path}: {oe}")
            return None

    def extract_premrna_sequences_per_gene(self, gene_ids: Union[str, List[str]], 
                                      output_path: Optional[str] = None) -> Dict[str, str]:
        """
        Extract pre-mRNA sequences for specific genes and optionally save them to a FASTA file.
        
        Parameters:
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
    
    def export_genome_data(self, output_dir: str, exclude_ids: Optional[List[str]] = None, force_overwrite: bool = False) -> Tuple[str, str]:
        """
        Exports genome data to FASTA (cDNA) and GTF files.

        FASTA Naming: {reference_name}.{e_release}.cdna.{tsl_string}.{exclude_str}.fa.gz
        GTF Naming:   {reference_name}.{e_release}.gtf.{exclude_str}.gz

        Args:
            output_dir: The directory to save the exported files.
            exclude_ids: An optional list of gene, transcript, or exon IDs to exclude.

        Returns:
            A tuple containing the paths to the exported FASTA and GTF files.
        """
        os.makedirs(output_dir, exist_ok=True)
        exclude_ids_set = set(exclude_ids) if exclude_ids else set()
        # --- Prepare FASTA file --- 
        
        e_release_str = self.e_release if self.e_release else "no_release"

        # Biotype suffix
        biotype_suffix = ""
        if self.biotype_to_keep:
            if len(self.biotype_to_keep) == 1:
                biotype_suffix = f".{self.biotype_to_keep[0]}"
            else:
                biotype_suffix = f".{'_'.join(self.biotype_to_keep)}"
        
        # Exclusion suffix
        exclusion_suffix = ""
        log_message_exclusion_details = " (all IDs)"
        if exclude_ids_set:
            if len(exclude_ids_set) == 1:
                single_excluded_id = list(exclude_ids_set)[0]
                exclusion_suffix = f".excluding_{single_excluded_id}"
                log_message_exclusion_details = f" (excluding ID '{single_excluded_id}')"
            else:
                exclusion_suffix = f".excluding_{len(exclude_ids_set)}_ids"
                log_message_exclusion_details = f" (excluding {len(exclude_ids_set)} IDs)"

        if self.tsl_to_keep:
            tsl_str = f".tsl{'_'.join(str(tsl) if tsl is not None else 'NA' for tsl in self.tsl_to_keep)}"
        else:
            tsl_str = ""
        
        if biotype_suffix == "" and tsl_str == "" and exclusion_suffix == "":
            fasta_filename = f"{self.reference_name}.{e_release_str}.cdna.all.fa.gz"
        else:
            fasta_filename = f"{self.reference_name}.{e_release_str}.cdna{biotype_suffix}{tsl_str}{exclusion_suffix}.fa.gz"
        output_fasta_path = os.path.join(output_dir, fasta_filename)
        fasta_file_exists = os.path.exists(output_fasta_path)

        if fasta_file_exists and not force_overwrite:
            if self.verbose:
                logging.info(f"Using existing FASTA file: {output_fasta_path}. Export will be skipped.")
            else:
                logging.debug(f"Using existing FASTA file: {output_fasta_path}. Export will be skipped.")
        else:
            if fasta_file_exists and force_overwrite: # Implies force_overwrite is True
                if self.verbose:
                    logging.info(f"Overwriting existing FASTA file: {output_fasta_path}")
                else:
                    logging.debug(f"Overwriting existing FASTA file: {output_fasta_path}")
                os.remove(output_fasta_path)
            # else: # File does not exist, or did exist and was removed by force_overwrite

            logging.debug(f"Exporting cDNA FASTA to: {output_fasta_path}{log_message_exclusion_details}")
            
            transcripts_written_count = 0
            with gzip.open(output_fasta_path, 'wt') as f_fasta_out:
                for transcript in self.transcripts: 
                    if transcript.transcript_id in exclude_ids_set or transcript.gene_id in exclude_ids_set:
                        continue
                    if transcript.sequence:
                        f_fasta_out.write(f">{transcript.transcript_id} gene_id={transcript.gene_id}\n{transcript.sequence}\n")
                        transcripts_written_count += 1
                    else:
                        logging.debug(f"Transcript {transcript.transcript_id} (gene {transcript.gene_id}) has no sequence. Not written to FASTA.")
            
            if self.verbose:
                logging.info(f"Wrote {transcripts_written_count} transcripts to {output_fasta_path}")
            else:
                logging.debug(f"Wrote {transcripts_written_count} transcripts to {output_fasta_path}")

        # --- Prepare GTF file --- 
        gtf_filename = f"{self.reference_name}.{e_release_str}{biotype_suffix}{exclusion_suffix}.gtf.gz"
        output_gtf_path = os.path.join(output_dir, gtf_filename)
        gtf_file_exists = os.path.exists(output_gtf_path)
        
        if gtf_file_exists and not force_overwrite:
            if self.verbose:
                logging.info(f"Using existing GTF file: {output_gtf_path}. Export will be skipped.")
            else:
                logging.debug(f"Using existing GTF file: {output_gtf_path}. Export will be skipped.")
        else:
            if gtf_file_exists and force_overwrite: # Implies force_overwrite is True
                if self.verbose:
                    logging.info(f"Overwriting existing GTF file: {output_gtf_path}")
                else:
                    logging.debug(f"Overwriting existing GTF file: {output_gtf_path}")
                os.remove(output_gtf_path)
            # else: # File does not exist or was removed

            logging.debug(f"Exporting GTF to: {output_gtf_path}{log_message_exclusion_details}")
            
            gtf_entries_written = 0
            with gzip.open(output_gtf_path, 'wt') as f_gtf_out:
                
                f_gtf_out.write(f"#!genome-build {self.reference_name}\n")
                f_gtf_out.write(f"#!genome-version {e_release_str}\n")
                f_gtf_out.write(f"#!genebuild-last-updated {datetime.now().strftime('%Y-%m-%d')}\n")

                source_tag = "custom_export" # Or use self.gtf_path if that was the original source

                for gene in self.genes: # self.genes are already biotype filtered
                    if gene.gene_id in exclude_ids_set:
                        continue
                    f_gtf_out.write(gene.to_gtf_entry(source=source_tag) + "\n")
                    gtf_entries_written += 1
                    for transcript in gene.transcripts: # These are TSL/biotype filtered if gene.transcripts come from self._transcripts
                        if transcript.transcript_id in exclude_ids_set or transcript.gene_id in exclude_ids_set: # Redundant check for gene_id if already checked for gene
                            continue
                        f_gtf_out.write(transcript.to_gtf_entry(source=source_tag) + "\n")
                        gtf_entries_written += 1
                        # Exons should be sorted by start position for GTF
                        sorted_exons = sorted(transcript.exons, key=lambda ex: ex.start)
                        for exon in sorted_exons:
                            if exon.exon_id in exclude_ids_set or exon.transcript_id in exclude_ids_set or transcript.gene_id in exclude_ids_set:
                                continue
                            # Exon needs chromosome, strand, gene_id from its transcript for GTF
                            f_gtf_out.write(exon.to_gtf_entry(chromosome=transcript.chromosome, 
                                                          strand=transcript.strand, 
                                                          gene_id=transcript.gene_id, 
                                                          source=source_tag) + "\n")
                            gtf_entries_written += 1
            
            if self.verbose:
                logging.info(f"Wrote {gtf_entries_written} entries to {output_gtf_path}")
            else:
                logging.debug(f"Wrote {gtf_entries_written} entries to {output_gtf_path}")
        
        return output_fasta_path, output_gtf_path