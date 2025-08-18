from __future__ import annotations
from pathlib import Path
from typing import Dict, Optional
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gffutils
import logging
import gzip
import shutil
from io import StringIO
import time
import pickle
from .genome import Genome
from .chromosome import Chromosome
from .gene import Gene
from .transcript import Transcript
from .exon import Exon
from .locus import Locus

class BuilderStateError(Exception):
    """Custom exception for GenomeBuilder state errors."""
    pass


def _strip_version(seq_id: str) -> str:
    """Removes version numbers from a sequence ID (e.g., 'NC_000001.11' -> 'NC_000001')."""
    seq_id_parts = seq_id.split('.')
    return seq_id_parts[0], seq_id_parts[1] if len(seq_id_parts) > 1 else None


class GenomeBuilder:
    """
    Constructs a Genome object from various file formats.

    This builder simplifies the process of assembling a complete Genome object
    by handling the parsing and integration of DNA sequences, cDNA sequences,
    and gene annotations from standard bioinformatics files.

    The correct order of operations is:
    1. with_dna_fasta()
    2. with_cdna_fasta() (optional)
    3. with_gtf_file()
    4. build()

    Example:
        builder = GenomeBuilder(id="hg38", species="Homo sapiens", name="Human Reference Genome")
        genome = (
            builder.with_dna_fasta(Path("path/to/dna.fa"))
            .with_cdna_fasta(Path("path/to/cdna.fa"))
            .with_gtf_file(Path("path/to/annotations.gtf"))
            .build()
        )
    """

    def __init__(self, id: str, species: str, name: str, 
                 main_chromosomes: Optional[list[str]] = None, 
                 separate_scaffolds: bool = False, 
                 output_dir: Path = Path('./data'), **kwargs):
        """
        Initializes the GenomeBuilder.

        Args:
            id: The ID of the genome.
            species: The species of the genome.
            name: The name of the genome.
            main_chromosomes: A list of chromosome IDs to be considered as the main set.
                              If None, defaults to human standard chromosomes (1-22, X, Y, M, MT).
            separate_scaffolds: If True, separates scaffold chromosomes into a second Genome object.
                                The `build()` method will then return a tuple: (main_genome, scaffold_genome).
            output_dir: The directory to save the pickled genome file. 
                        If None, defaults to ./data
            kwargs: Additional attributes for the Genome object.
        """
        self._genome = Genome(id, species, name, **kwargs)
        self._cdna_records: Dict[str, SeqRecord] = {}
        self._genes_map: Dict[str, Gene] = {}
        self._transcripts_map: Dict[str, Transcript] = {}
        self._chromosome_filter = None
        self._separate_scaffolds = separate_scaffolds
        self._scaffold_genome: Optional[Genome] = None
        self._output_dir = output_dir

        if main_chromosomes is None:
            # Default to standard human chromosomes
            standard_set = {str(i) for i in range(1, 23)} | {'X', 'Y', 'M', 'MT'}
            self._main_chromosomes = standard_set.union({f'chr{c}' for c in standard_set})
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

        self.logger.info(f"Loading DNA sequences from {dna_fasta_path}...")

        dna_file_to_use = dna_fasta_path

        if str(dna_fasta_path).endswith('.gz'):
            extracted_path = dna_fasta_path.with_suffix('')
            if extracted_path.exists():
                self.logger.info(f"Using existing extracted DNA FASTA file: {extracted_path}")
                dna_file_to_use = extracted_path
            else:
                self.logger.info(f"Extracting gzipped DNA FASTA to: {extracted_path}")
                with gzip.open(dna_fasta_path, 'rt') as gz_in:
                    with open(extracted_path, 'w') as f_out:
                        shutil.copyfileobj(gz_in, f_out)
                dna_file_to_use = extracted_path

        dna_records = SeqIO.index(str(dna_file_to_use), "fasta")
        
        for seq_id in dna_records:
            if self._chromosome_filter and seq_id not in self._chromosome_filter:
                continue
            
            chromosome = Chromosome(seq_id, 1, len(dna_records[seq_id]), '+', dna_records, fasta_path=dna_file_to_use, genome=self._genome)

            if self._separate_scaffolds and seq_id not in self._main_chromosomes:
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
        
        def parse_cdna():
            for record in SeqIO.parse(handle, "fasta"):
                record.id = _strip_version(record.id)
                yield record

        if str(cdna_fasta_path).endswith('.gz'):
            self.logger.info(f"Reading gzipped cDNA FASTA file: {cdna_fasta_path}")
            with gzip.open(cdna_fasta_path, "rt") as handle:
                self._cdna_records = SeqIO.to_dict(parse_cdna())
        else:
            with open(cdna_fasta_path, "rt") as handle:
                self._cdna_records = SeqIO.to_dict(parse_cdna())

        self.logger.info(f"Loaded {len(self._cdna_records)} cDNA sequences.")
        return self

    def with_gtf_file(self, gtf_path: Path) -> "GenomeBuilder":
        """
        Parses a GTF file to build the gene-transcript-exon hierarchy.
        `with_dna_fasta()` must be called before this method.
        """
        if not self._genome.chromosomes:
            raise BuilderStateError("Must call with_dna_fasta() before with_gtf_file().")
        if self._genes_map:
            raise BuilderStateError("with_gtf_file() has already been called.")

        self.logger.info(f"Processing annotations from {gtf_path}...")

        gtf_db_path = gtf_path.with_suffix('.db')

        if gtf_db_path.exists():
            self.logger.info(f"Loading existing gffutils database: {gtf_db_path}")
            db = gffutils.FeatureDB(str(gtf_db_path))
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
        
        start_time = time.time()
        self._create_genes(db)
        self.logger.info(f"Created genes in {time.time() - start_time:.2f} seconds")

        start_time = time.time()
        self._create_transcripts(db)
        self.logger.info(f"Created transcripts in {time.time() - start_time:.2f} seconds")

        start_time = time.time()
        self._create_exons(db)
        self.logger.info(f"Created exons in {time.time() - start_time:.2f} seconds")

        self.logger.info(f"Successfully parsed and linked {len(self._genes_map)} genes, "
                         f"{len(self._transcripts_map)} transcripts.")
        return self

    def _create_genes(self, db: gffutils.FeatureDB):
        """Creates Gene objects from the GTF database."""
        for g in db.features_of_type('gene'):
            if self._chromosome_filter and g.chrom not in self._chromosome_filter:
                continue

            chromosome = self._genome.chromosomes.get(g.chrom)
            if not chromosome and self._scaffold_genome:
                chromosome = self._scaffold_genome.chromosomes.get(g.chrom)

            if not chromosome:
                self.logger.warning(f"Chromosome '{g.chrom}' for gene '{g.id}' not found. Skipping gene.")
                continue

            try:
                attributes = dict(g.attributes)

                gene_names = attributes.pop('gene_name', attributes.pop('gene', [g.id]))
                gene_name = gene_names.pop(0) if isinstance(gene_names, list) else gene_names
                attributes['gene_synonyms'] = gene_names

                # Remove exon and transcript related attributes
                attributes = {k: v for k, v in attributes.items() 
                            if not (k.startswith('exon') or k.startswith('transcript'))}
                gene_id = attributes.pop('gene_id', [g.id])[0]
                
                gene = Gene(id=gene_id, name=gene_name, start=g.start,
                            end=g.end, strand=g.strand, chromosome=chromosome,
                            genome=self._genome,
                            **attributes)
                chromosome.add_gene(gene)
                self._genes_map[g.id] = gene

            except KeyError:
                self.logger.warning(f"Chromosome '{g.chrom}' for gene '{g.id}' not found in FASTA. Skipping gene.")
                raise KeyError(f"Chromosome '{g.chrom}' for gene '{g.id}' not found in FASTA. Skipping gene.")

    def _create_transcripts(self, db: gffutils.FeatureDB):
        """Creates Transcript objects and links them to genes."""
        for t in db.features_of_type('transcript'):
            attributes = dict(t.attributes)

            gene_id = attributes.pop('gene_id', attributes.pop('gene', [None]))[0]
            transcript_id = attributes.pop('transcript_id', [t.id])[0]  
            # Remove exon and gene related attributes
            attributes = {k: v for k, v in attributes.items() 
                            if not (k.startswith('exon') or k.startswith('gene'))}
            
            if gene_id and gene_id in self._genes_map:
                gene = self._genes_map[gene_id]
                sequence = self._cdna_records.pop(transcript_id, SeqRecord(Seq(""))).seq
                transcript = Transcript(id=transcript_id, start=t.start, end=t.end, strand=t.strand,
                                        sequence=sequence, gene=gene, genome=self._genome, **attributes)
                gene.add_transcript(transcript)
                self._transcripts_map[t.id] = transcript
            else:
                self.logger.warning(f"Gene '{gene_id}' for transcript '{t.id}' not found. Skipping transcript.")

    def _create_exons(self, db: gffutils.FeatureDB):
        """Creates Exon objects and links them to transcripts."""
        for e in db.features_of_type('exon'):
            if self._chromosome_filter and e.chrom not in self._chromosome_filter:
                continue

            attributes = dict(e.attributes)

            transcript_id = attributes.pop('transcript_id', [None])[0]
            exon_id = attributes.pop('exon_id', [e.id])[0]
            # Remove transcript and gene related attributes
            attributes = {k: v for k, v in attributes.items() 
                            if not (k.startswith('transcript') or k.startswith('gene'))}
            
            if transcript_id and transcript_id in self._transcripts_map:
                transcript = self._transcripts_map[transcript_id]
                exon = Exon(id=exon_id, start=e.start, end=e.end, strand=e.strand, transcript=transcript,
                            genome=self._genome,
                            **attributes)
                transcript.add_exon(exon)
            else:
                self.logger.warning(f"Transcript '{transcript_id}' for exon '{e.id}' not found. Skipping exon.")

    def build(self, pickle_genome: bool = False) -> Genome | tuple[Genome, Genome]:
        """
        Finalizes the Genome object by creating an index for fast lookups.
        
        Args:
            pickle_genome: If True, saves the final genome object(s) to a pickle file.
        """
        if not self._genes_map:
            raise BuilderStateError("Cannot build Genome. GTF data is missing. "
                                    "Please call with_gtf_file() before build().")
        
        self.logger.info("Indexing genome for fast lookups...")
        self._genome.index()
        if self._scaffold_genome:
            self.logger.info("Indexing scaffold genome for fast lookups...")
            self._scaffold_genome.index()

        self.logger.info("Genome construction complete.")
        
        self._offload_memory()

        if pickle_genome:
            output_path = self._output_dir / f"{self._genome.species}.{self._genome.id}.pkl"
            output_path.parent.mkdir(parents=True, exist_ok=True)
            self.logger.info(f"Saving genome to {output_path}...")
            if self._scaffold_genome:
                with open(output_path, "wb") as f:
                    pickle.dump((self._genome, self._scaffold_genome), f)
            else:
                with open(output_path, "wb") as f:
                    pickle.dump(self._genome, f)
            self.logger.info("Genome saved successfully.")

        if self._scaffold_genome:
            return self._genome, self._scaffold_genome
        
        return self._genome

    def _offload_memory(self):
        """Clears large data structures from memory after the build is complete."""
        self.logger.info("Offloading builder memory...")
        self._cdna_records.clear()
        self._genes_map.clear()
        self._transcripts_map.clear()
        
        
        self.logger.info("Memory offload complete.")
        
    @staticmethod
    def load_from_file(file_path: Path) -> Genome | tuple[Genome, Genome]:
        """
        Loads a Genome object (or a tuple of Genome objects) from a pickle file.

        Args:
            file_path: The path to the pickle file.

        Returns:
            The loaded Genome object or a tuple of (main_genome, scaffold_genome).
        """
        with open(file_path, "rb") as f:
            return pickle.load(f) 