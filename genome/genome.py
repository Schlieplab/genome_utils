from __future__ import annotations
from typing import Dict, Iterator, Any, List

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .chromosome import Chromosome
from .gene import Gene
from .transcript import Transcript
from .exon import Exon
from .locus import Locus


class Genome:
    """Represents a collection of chromosomes, managing sequence data from a FASTA file."""

    def __init__(self, id: str, species: str, name: str, **kwargs):
        """
        Initializes the Genome by indexing a FASTA file.
        Chromosomes are automatically created based on the sequences in the file.
        Args:
            id: The ID of the genome.
            species: The species of the genome.
            name: The name of the genome.
            kwargs: Additional keyword arguments.
        """
        self._attributes: Dict[str, Any] = {
            key: value[0] if isinstance(value, list) and len(value) == 1 else value
            for key, value in kwargs.items()
        }
        self.id = id
        self.species = species
        self.name = name
        self._chromosomes: Dict[str, Chromosome] = {}
        self._genes: Dict[str, Gene] = {}
        self._transcripts: Dict[str, Transcript] = {}
        self._exons: Dict[str, Exon] = {}
        self.is_indexed: bool = False


    def __repr__(self) -> str:
        """Return a developer-friendly representation of the Genome."""
        return (f"{self.__class__.__name__}("
                f"id='{self.id}', "
                f"species='{self.species}', "
                f"name='{self.name}')")

    def add_chromosome(self, chromosome: Chromosome):
        """Add a chromosome to the genome."""
        if chromosome.id in self._chromosomes:
            raise ValueError(f"Chromosome with ID '{chromosome.id}' already exists.")
        self._chromosomes[chromosome.id] = chromosome
        self.is_indexed = False

    def index(self):
        """
        Creates an index of all genes, transcripts, and exons for fast lookup.
        This method MUST be called after all genomic features have been added.
        """
        for chrom in self._chromosomes.values():
            for gene in chrom.genes:
                self._genes[gene.id] = gene
                for transcript in gene.transcripts:
                    self._transcripts[transcript.id] = transcript
                    for exon in transcript.exons:
                        self._exons[exon.id] = exon
        self.is_indexed = True
    
    def sequence_by_locus(self, locus: Locus) -> Seq:
        """Get a sequence by its locus."""
        if not self.is_indexed:
            raise RuntimeError("The genome is not indexed. Call .index() after adding features.")
        
        return self._chromosomes[locus.chromosome_id].get_subsequence(locus.start, 
                                                                     locus.end, 
                                                                     locus.strand)

    def chromosomes_iter(self) -> Iterator[Chromosome]:
        """Iterate over all chromosomes in the genome."""
        if not self.is_indexed:
            raise RuntimeError("The genome is not indexed. Call .index() after adding features.")
        for chrom in self._chromosomes.values():
            yield chrom

    @property
    def chromosomes(self) -> List[Chromosome]:
        """Get all chromosomes in the genome."""
        return list(self.chromosomes_iter())

    def genes_iter(self) -> Iterator[Gene]:
        """Iterate over all genes in the genome."""
        if not self.is_indexed:
            raise RuntimeError("The genome is not indexed. Call .index() after adding features.")
        for chrom in self._chromosomes.values():
            for gene in chrom.genes:
                yield gene

    @property
    def genes(self) -> List[Gene]:
        """Get all genes in the genome."""
        return list(self.genes_iter())


    def transcripts_iter(self) -> Iterator[Transcript]:
        """Iterate over all transcripts in the genome."""
        if not self.is_indexed:
            raise RuntimeError("The genome is not indexed. Call .index() after adding features.")
        for chrom in self._chromosomes.values():
            for gene in chrom.genes: # Access genes from chromosome, then transcripts from gene
                for transcript in gene.transcripts:
                    yield transcript
    
    @property
    def transcripts(self) -> List[Transcript]:
        """Get all transcripts in the genome."""
        return list(self.transcripts_iter())

    def exons_iter(self) -> Iterator[Exon]:
        """Iterate over all exons in the genome."""
        if not self.is_indexed:
            raise RuntimeError("The genome is not indexed. Call .index() after adding features.")
        for chrom in self._chromosomes.values():
            for gene in chrom.genes:
                for transcript in gene.transcripts:
                    for exon in transcript.exons:
                        yield exon

    @property
    def exons(self) -> List[Exon]:
        """Get all exons in the genome."""
        return list(self.exons_iter())


    def chromosome_by_id(self, chromosome_id: str) -> Chromosome:
        """Get a chromosome by its ID using the index. Raises ValueError if not found."""
        try:
            return self._chromosomes[chromosome_id]
        except KeyError:
            raise ValueError(f"Chromosome with ID '{chromosome_id}' not found.")

    def gene_by_id(self, gene_id: str) -> Gene:
        """Get a gene by its ID using the index. Raises ValueError if not found."""
        if not self.is_indexed:
            raise RuntimeError("The genome is not indexed. Call .index() after adding features.")
        try:
            return self._genes[gene_id]
        except KeyError:
            raise ValueError(f"Gene with ID '{gene_id}' not found.")

    def transcript_by_id(self, transcript_id: str) -> Transcript:
        """Get a transcript by its ID using the index. Raises ValueError if not found."""
        if not self.is_indexed:
            raise RuntimeError("The genome is not indexed. Call .index() after adding features.")
        try:
            return self._transcripts[transcript_id]
        except KeyError:
            raise ValueError(f"Transcript with ID '{transcript_id}' not found.")

    def exon_by_id(self, exon_id: str) -> Exon:
        """Get an exon by its ID using the index. Raises ValueError if not found."""
        if not self.is_indexed:
            raise RuntimeError("The genome is not indexed. Call .index() after adding features.")
        try:
            return self._exons[exon_id]
        except KeyError:
            raise ValueError(f"Exon with ID '{exon_id}' not found.")


    def __getattr__(self, name: str) -> Any:
        """Allow direct access to attributes in the attributes dictionary."""
        try:
            return self._attributes[name]
        except KeyError:
            raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")
        
    def __setattr__(self, name: str, value: Any):
        """Allow setting attributes. Explicitly defined attributes are set normally. New, dynamic attributes are stored in the 'attributes' dictionary."""
        if name in self.__dict__ or name in self.__class__.__dict__ or name == '_attributes' or not hasattr(self, '_attributes'):
            super().__setattr__(name, value)
        else:
            self._attributes[name] = value
            
    