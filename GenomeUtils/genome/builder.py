#!/usr/bin/env python
"""
Filename: GenomeUtils/genome/builder.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.0
Description: This file contains the GenomeBuilder class for constructing genome objects.
License: LGPL-3.0-or-later
"""

from __future__ import annotations

import gzip
import json
import logging
import shutil
from pathlib import Path
from typing import Dict, Optional

import gffutils
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from tqdm import tqdm

from .chromosome import Chromosome
from .exon import Exon
from .gene import Gene
from .genome import Genome
from .transcript import Transcript


def _get_default_chromosomes_for_species(species: str) -> set[str]:
    """
    Returns the default main chromosomes for common species.
    
    Args:
        species: Species name (case-insensitive)
        
    Returns:
        Set of chromosome identifiers including both with and without 'chr' prefix
    """
    species_lower = species.lower()
    
    if any(term in species_lower for term in ['human', 'homo sapiens', 'homo_sapiens']):
        # Human: 1-22, X, Y, M, MT
        standard_set = {str(i) for i in range(1, 23)} | {'X', 'Y', 'M', 'MT'}
    elif any(term in species_lower for term in ['mouse', 'mice', 'mus musculus', 'mus_musculus']):
        # Mouse: 1-19, X, Y, M, MT
        standard_set = {str(i) for i in range(1, 20)} | {'X', 'Y', 'M', 'MT'}
    else:
        # Default to human if species not recognized
        raise ValueError(f"Species '{species}' not recognized. Please use a supported species.")
    
    # Return both with and without 'chr' prefix
    return set(standard_set).union({f'chr{c}' for c in standard_set})


class BuilderStateError(Exception):
    """Custom exception for GenomeBuilder state errors."""
    pass


def _strip_version(seq_id: str) -> str:
    """Removes version numbers from a sequence ID (e.g., 'NC_000001.11' -> 'NC_000001')."""
    seq_id_parts = seq_id.split('.')
    return seq_id_parts[0] if len(seq_id_parts) > 1 else seq_id


class GenomeBuilder:
    """Constructs a Genome object from various file formats.

    This builder simplifies the process of assembling a complete Genome object
    by handling the parsing and integration of DNA sequences, cDNA sequences,
    and gene annotations from standard bioinformatics files.

    The correct order of operations is:

    1. with_dna_fasta()
    2. with_cdna_fasta()
    3. with_gtf_file()
    4. build()

    Example::

        builder = GenomeBuilder(id="hg38", species="homo_sapiens", name="Human Reference Genome")
        genome = (
            builder.with_dna_fasta(Path("path/to/dna.fa"))
            .with_cdna_fasta(Path("path/to/cdna.fa"))
            .with_gtf_file(Path("path/to/annotations.gtf"))
            .build()
        )
    """

    def __init__(self, 
                 id: str, 
                 species: str, 
                 name: str, 
                 main_chromosomes: Optional[list[str]] = None, 
                 separate_scaffolds: bool = True, 
                 **kwargs):
        """
        Initializes the GenomeBuilder.

        Args:
            id: The ID of the genome.
            species: The species of the genome.
            name: The name of the genome.
            main_chromosomes: A list of chromosome IDs to be considered as the main set.
                              If None, defaults to species-appropriate chromosomes 
                              (Human: 1-22,X,Y,M,MT; Mouse: 1-19,X,Y,M,MT).
            separate_scaffolds: If True, separates scaffold chromosomes into a second Genome object.
                                The `build()` method will then return a tuple: (main_genome, scaffold_genome).
            kwargs: Additional attributes for the Genome object.
        """
        self._genome = Genome(id, species, name, **kwargs)
        self._cdna_records: Dict[str, SeqRecord] = {}
        self._genes_map: Dict[str, Gene] = {}
        self._transcripts_map: Dict[str, Transcript] = {}
        self._chromosome_filter = None
        self._separate_scaffolds = separate_scaffolds
        self._scaffold_genome: Optional[Genome] = None

        if main_chromosomes is None:
            # Use species-dependent default chromosomes
            self._main_chromosomes = _get_default_chromosomes_for_species(species)
        else:
            self._main_chromosomes = set(main_chromosomes)

        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(self.__class__.__name__)

        if self._separate_scaffolds:
            self.logger.info("Scaffold separation enabled. `build()` will return (main_genome, scaffold_genome).")
            self._scaffold_genome = Genome(
                id=f"{id}_scaffolds",
                species=species,
                name=f"{name} (Scaffolds)",
                **kwargs
            )

    def set_chromosome_filter(self, chromosomes: list[str]) -> "GenomeBuilder":
        """
        Set a filter to only include specified chromosomes.
        """
        if self._genome.chromosomes:
            raise BuilderStateError("Cannot set chromosome filter after with_dna_fasta() has been called.")
        
        self._chromosome_filter = set(chromosomes)
        self.logger.info(f"Chromosome filter set to: {self._chromosome_filter}")
        return self

    def with_dna_fasta(self, dna_fasta_path: Path) -> "GenomeBuilder":
        """
        Loads chromosome sequences from a genomic DNA FASTA file.
        This must be the first step in the build process.
        """
        if self._genome.chromosomes:
            raise BuilderStateError("with_dna_fasta() has already been called.")

        dna_file_to_use = dna_fasta_path

        if str(dna_fasta_path).endswith('.gz'):
            extracted_path = dna_fasta_path.with_suffix('')
            if extracted_path.exists():
                dna_file_to_use = extracted_path
            else:
                self.logger.info(f"Extracting gzipped DNA FASTA to: {extracted_path}")
                with gzip.open(dna_fasta_path, 'rt') as gz_in:
                    with open(extracted_path, 'w') as f_out:
                        shutil.copyfileobj(gz_in, f_out)
                dna_file_to_use = extracted_path

        self.logger.info(f"Loading DNA sequences from {dna_file_to_use}...")

        dna_records = SeqIO.index(str(dna_file_to_use), "fasta")
        
        for record in SeqIO.parse(dna_file_to_use, "fasta"):
            if self._chromosome_filter and record.id not in self._chromosome_filter:
                continue
            
            chromosome = Chromosome(record.id, dna_records, genome=self._genome, length=len(record.seq))

            if self._separate_scaffolds and record.id not in self._main_chromosomes:
                if self._scaffold_genome:
                    self._scaffold_genome.add_chromosome(chromosome)
            else:
                self._genome.add_chromosome(chromosome)

        self.logger.info(f"Loaded {len(self._genome.chromosomes)} main chromosomes.")
        if self._scaffold_genome:
            self.logger.info(f"Loaded {len(self._scaffold_genome.chromosomes)} scaffold chromosomes.")
        return self

    def with_cdna_fasta(self, cdna_fasta_path: Path) -> "GenomeBuilder":
        """
        Loads transcript sequences from a cDNA FASTA file.
        """
        if self._cdna_records:
            raise BuilderStateError("with_cdna_fasta() has already been called.")
        
        self.logger.info(f"Loading cDNA sequences from {cdna_fasta_path}...")
        
        open_func = gzip.open if str(cdna_fasta_path).endswith('.gz') else open
        with open_func(cdna_fasta_path, "rt") as handle:
            self._cdna_records = SeqIO.to_dict(SeqIO.parse(handle, "fasta"), key_function=lambda x: _strip_version(x.id))

        self.logger.info(f"Loaded {len(self._cdna_records)} cDNA sequences.")
        return self

    def with_gtf_file(self, gtf_path: Path) -> "GenomeBuilder":
        """
        Parses a GTF file to build the gene-transcript-exon hierarchy.
        `with_dna_fasta()` and `with_cdna_fasta()` must be called before this method.
        """

        if not self._genome.chromosomes:
            raise BuilderStateError("Must call with_dna_fasta() before with_gtf_file().")
        if not self._cdna_records:
            raise BuilderStateError("Must call with_cdna_fasta() before with_gtf_file().")
        if self._genes_map:
            raise BuilderStateError("with_gtf_file() has already been called.")

        self.logger.info(f"Processing annotations from {gtf_path}...")

        gtf_db_path = gtf_path.with_suffix('.db')

        if gtf_db_path.exists():
            self.logger.info(f"Loading existing gffutils database: {gtf_db_path}")
            try:
                db = gffutils.FeatureDB(str(gtf_db_path))
            except Exception as e:
                self.logger.warning(f"Error loading existing gffutils database: {e}. Creating new database.")
                gtf_db_path.unlink()
                db = gffutils.create_db(str(gtf_path), 
                                        dbfn=str(gtf_db_path), 
                                        keep_order=False, 
                                        merge_strategy='error', 
                                        id_spec={'gene': 'gene_id', 'transcript': 'transcript_id'}, 
                                        disable_infer_genes=True, 
                                        disable_infer_transcripts=True)
        else:
            self.logger.info(f"Database not found. Creating new database at: {gtf_db_path}")
            gtf_file_to_use = gtf_path
            
            if str(gtf_path).endswith('.gz'):
                extracted_path = gtf_path.with_suffix('')
                
                if extracted_path.exists():
                    self.logger.info(f"Using existing extracted GTF file: {extracted_path}")
                    gtf_file_to_use = extracted_path
                else:
                    self.logger.info(f"Extracting gzipped GTF file to: {extracted_path}")
                    with gzip.open(gtf_path, 'rt') as gz_file:
                        with open(extracted_path, 'w') as out_file:
                            out_file.write(gz_file.read())
                    gtf_file_to_use = extracted_path
            
            db = gffutils.create_db(
                    str(gtf_file_to_use),
                    dbfn=str(gtf_db_path),
                    keep_order=False,
                    merge_strategy='error',
                    id_spec={'gene': 'gene_id', 'transcript': 'transcript_id'},
                    disable_infer_genes=True,
                    disable_infer_transcripts=True
            )

        logging.info(f"GTF database created at: {gtf_db_path}")
        
        self._create_genes(db)

        self._create_transcripts(db)

        self._create_exons(db)

        self.logger.info(f"Successfully parsed and linked {len(self._genes_map)} genes "
                         f"and {len(self._transcripts_map)} transcripts.")
        return self

    def _create_genes(self, db: gffutils.FeatureDB):
        """Creates Gene objects from the GTF database."""
        query = "SELECT id, seqid, start, end, strand, attributes FROM features WHERE featuretype = 'gene'"
        
        count_query = "SELECT count(*) FROM features WHERE featuretype = 'gene'"
        total_genes = db.conn.execute(count_query).fetchone()[0]


        for g_id, seqid, start, end, strand, attributes_json in tqdm(db.conn.execute(query), total=total_genes, desc="Creating genes"):
            if self._chromosome_filter and seqid not in self._chromosome_filter:
                continue
            
            try:
                chromosome = self._genome.chromosome_by_id(seqid)
            except ValueError:
                if self._scaffold_genome:
                    chromosome = self._scaffold_genome.chromosome_by_id(seqid)
                else:
                    self.logger.warning(f"Chromosome '{seqid}' for gene '{g_id}' not found. Skipping gene.")
                    continue
            
            if not chromosome:
                self.logger.warning(f"Chromosome '{seqid}' for gene '{g_id}' not found. Skipping gene.")
                continue

            try:
                attributes = json.loads(attributes_json)

                gene_names = attributes.pop('gene_name', attributes.pop('gene', [g_id]))
                gene_name = gene_names[0]
                attributes['gene_synonyms'] = gene_names[1:]

                attributes = {k: v for k, v in attributes.items() 
                            if not (k.startswith('exon') or k.startswith('transcript'))}
                
                gene_id = attributes.pop('gene_id', [g_id])[0]
                attributes = {k.replace('gene_', ''): v for k, v in attributes.items()}
                attributes = {k: (v[0] if isinstance(v, list) and len(v) == 1 else v) for k, v in attributes.items()}

                gene = Gene(id=gene_id, name=gene_name, chr=chromosome.id, start=start,
                            end=end, strand=strand, chromosome=chromosome,
                            genome=self._genome,
                            **attributes)
                chromosome.add_gene(gene)
                self._genes_map[g_id] = gene

            except Exception as e:
                self.logger.warning(f"Error processing gene '{g_id}': {e}. Skipping.")

    def _create_transcripts(self, db: gffutils.FeatureDB):
        """Creates Transcript objects and links them to genes."""
        query = "SELECT id, start, end, strand, attributes FROM features WHERE featuretype = 'transcript'"
        
        count_query = "SELECT count(*) FROM features WHERE featuretype = 'transcript'"
        total_transcripts = db.conn.execute(count_query).fetchone()[0]

        for t_id, start, end, strand, attributes_json in tqdm(db.conn.execute(query), total=total_transcripts, desc="Creating transcripts"):
            attributes = json.loads(attributes_json)

            gene_id = attributes.pop('gene_id', attributes.pop('gene', [None]))[0]
            transcript_id = attributes.pop('transcript_id', [t_id])[0]  
            
            attributes = {k: v for k, v in attributes.items() 
                            if not (k.startswith('exon') or k.startswith('gene'))}
            attributes = {k.replace('transcript_', ''): v for k, v in attributes.items()}
            attributes = {k: (v[0] if isinstance(v, list) and len(v) == 1 else v) for k, v in attributes.items()}
            if gene_id and gene_id in self._genes_map:
                gene = self._genes_map[gene_id]
                sequence = self._cdna_records.pop(transcript_id, SeqRecord(Seq(""))).seq
                transcript = Transcript(id=transcript_id, chr=gene.chr, start=start, end=end, strand=strand,
                                        sequence=sequence, gene=gene, genome=self._genome, **attributes)
                gene.add_transcript(transcript)
                self._transcripts_map[t_id] = transcript
            else:
                self.logger.warning(f"Gene '{gene_id}' for transcript '{t_id}' not found. Skipping transcript.")

    def _create_exons(self, db: gffutils.FeatureDB):
        """Creates Exon objects and links them to transcripts."""
        query = "SELECT id, seqid, start, end, strand, attributes FROM features WHERE featuretype = 'exon'"
        
        count_query = "SELECT count(*) FROM features WHERE featuretype = 'exon'"
        total_exons = db.conn.execute(count_query).fetchone()[0]

        for e_id, seqid, start, end, strand, attributes_json in tqdm(db.conn.execute(query), total=total_exons, desc="Creating exons"):
            if self._chromosome_filter and seqid not in self._chromosome_filter:
                continue

            attributes = json.loads(attributes_json)

            transcript_id = attributes.pop('transcript_id', [None])[0]
            exon_id = attributes.pop('exon_id', [e_id])[0]
            
            attributes = {k: v for k, v in attributes.items() 
                            if not (k.startswith('transcript') or k.startswith('gene'))}
            attributes = {k.replace('exon_', ''): v for k, v in attributes.items()}
            attributes = {k: (v[0] if isinstance(v, list) and len(v) == 1 else v) for k, v in attributes.items()}
            if transcript_id and transcript_id in self._transcripts_map:
                transcript = self._transcripts_map[transcript_id]
                exon = Exon(id=exon_id, chr=transcript.chr, start=start, end=end, strand=strand, transcript=transcript,
                            genome=self._genome,
                            **attributes)
                transcript.add_exon(exon)
            else:
                self.logger.warning(f"Transcript '{transcript_id}' for exon '{e_id}' not found. Skipping exon.")

    def build(self) -> Genome | tuple[Genome, Genome]:
        """
        Finalizes the Genome object by creating an index for fast lookups.
        """
        if not self._genes_map:
            raise BuilderStateError("Cannot build Genome. GTF data is missing. "
                                    "Please call with_gtf_file() before build().")
        
        self._genome.index()
        if self._scaffold_genome:
            self.logger.info("Indexing scaffold genome for fast lookups...")
            self._scaffold_genome.index()

        self.logger.info("Genome construction complete.")
        
        self._offload_memory()
        

        if self._scaffold_genome:
            return self._genome, self._scaffold_genome
        
        return self._genome

    def _offload_memory(self):
        """Clears large data structures from memory after the build is complete."""
        self._cdna_records.clear()
        self._genes_map.clear()
        self._transcripts_map.clear()
