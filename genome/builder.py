from __future__ import annotations
from pathlib import Path
from typing import Dict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gffutils
import logging
import gzip
import shutil
from io import StringIO
import time
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
    return seq_id.split('.')[0]


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

    def __init__(self, id: str, species: str, name: str, **kwargs):
        """
        Initializes the GenomeBuilder.

        Args:
            id: The ID of the genome.
            species: The species of the genome.
            name: The name of the genome.
            kwargs: Additional attributes for the Genome object.
        """
        self._genome = Genome(id, species, name, **kwargs)
        self._cdna_records: Dict[str, SeqRecord] = {}
        self._genes_map: Dict[str, Gene] = {}
        self._transcripts_map: Dict[str, Transcript] = {}
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(self.__class__.__name__)

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
            chromosome = Chromosome(seq_id, 1, len(dna_records[seq_id]), '+', dna_records)
            self._genome.add_chromosome(chromosome)

        self.logger.info(f"Loaded {len(self._genome.chromosomes)} chromosomes.")
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
            try:
                chromosome = self._genome.chromosomes[g.chrom]
                attributes = dict(g.attributes)
                gene_name = attributes.pop('gene_name', [g.id])[0]
                attributes.pop('gene_id', None)  # Already used for id
                gene = Gene(id=g.id, name=gene_name, start=g.start,
                            end=g.end, strand=g.strand, chromosome=chromosome,
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
            gene_id = attributes.pop('gene_id', [None])[0]
            attributes.pop('transcript_id', None)  # Already used for id
            if gene_id and gene_id in self._genes_map:
                gene = self._genes_map[gene_id]
                sequence = str(self._cdna_records.pop(t.id, SeqRecord(Seq(""))).seq)
                transcript = Transcript(id=t.id, start=t.start, end=t.end, strand=t.strand,
                                        sequence=sequence, gene=gene, **attributes)
                gene.add_transcript(transcript)
                self._transcripts_map[t.id] = transcript
            else:
                self.logger.warning(f"Gene '{gene_id}' for transcript '{t.id}' not found. Skipping transcript.")

    def _create_exons(self, db: gffutils.FeatureDB):
        """Creates Exon objects and links them to transcripts."""
        for e in db.features_of_type('exon'):
            attributes = dict(e.attributes)
            transcript_id = attributes.pop('transcript_id', [None])[0]
            attributes.pop('gene_id', None)  # Redundant in exon context
            if transcript_id and transcript_id in self._transcripts_map:
                transcript = self._transcripts_map[transcript_id]
                exon = Exon(id=e.id, start=e.start, end=e.end, strand=e.strand, transcript=transcript,
                            **attributes)
                transcript.add_exon(exon)
            else:
                self.logger.warning(f"Transcript '{transcript_id}' for exon '{e.id}' not found. Skipping exon.")

    def build(self) -> Genome:
        """
        Finalizes the Genome object by creating an index for fast lookups.
        """
        if not self._genes_map:
            raise BuilderStateError("Cannot build Genome. GTF data is missing. "
                                    "Please call with_gtf_file() before build().")
        
        self.logger.info("Indexing genome for fast lookups...")
        self._genome.index()
        self.logger.info("Genome construction complete.")
        
        self._offload_memory()
        
        return self._genome

    def _offload_memory(self):
        """Clears large data structures from memory after the build is complete."""
        self.logger.info("Offloading builder memory...")
        self._cdna_records.clear()
        self._genes_map.clear()
        self._transcripts_map.clear()
        
        
        self.logger.info("Memory offload complete.") 