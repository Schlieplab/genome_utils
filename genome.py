from __future__ import annotations
from typing import Dict, Iterator, Any

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from .chromosome import Chromosome
from .gene import Gene
from .transcript import Transcript
from .exon import Exon


class Genome:
    """Represents a collection of chromosomes, managing sequence data from a FASTA file."""

    def __init__(self, id: str, species: str, name: str, **kwargs):
        """
        Initializes the Genome by indexing a FASTA file.
        Chromosomes are automatically created based on the sequences in the file.
        """
        self.id = id
        self.species = species
        self.name = name
        self.attributes: Dict[str, Any] = kwargs
        self.chromosomes: Dict[str, Chromosome] = {}

        self._genes_by_id: Dict[str, Gene] = {}
        self._transcripts_by_id: Dict[str, Transcript] = {}
        self._exons_by_id: Dict[str, Exon] = {}
        self._is_indexed: bool = False


    def __repr__(self) -> str:
        """Return a developer-friendly representation of the Genome."""
        return (f"{self.__class__.__name__}("
                f"id='{self.id}', "
                f"species='{self.species}', "
                f"name='{self.name}')")

    def add_chromosome(self, chromosome: Chromosome):
        """Add a chromosome to the genome."""
        if chromosome.id in self.chromosomes:
            raise ValueError(f"Chromosome with ID '{chromosome.id}' already exists.")
        self.chromosomes[chromosome.id] = chromosome

    def index(self):
        """
        Creates an index of all genes, transcripts, and exons for fast lookup.
        This method should be called after all genomic features have been added.
        """
        for chrom in self.chromosomes.values():
            for gene in chrom.genes:
                self._genes_by_id[gene.id] = gene
                for transcript in gene.transcripts:
                    self._transcripts_by_id[transcript.id] = transcript
                    for exon in transcript.exons:
                        self._exons_by_id[exon.id] = exon
        self._is_indexed = True

    def __getitem__(self, key: str) -> Chromosome:
        """Get a chromosome by its ID."""
        if not self._is_indexed:
            raise RuntimeError("The genome is not indexed. Call .index() after adding features.")
        return self.chromosomes[key]

    @property
    def genes(self) -> Iterator[Gene]:
        """Iterate over all genes in the genome."""
        if not self._is_indexed:
            raise RuntimeError("The genome is not indexed. Call .index() after adding features.")
        for chrom in self.chromosomes.values():
            for gene in chrom.genes:
                yield gene


    @property
    def transcripts(self) -> Iterator[Transcript]:
        """Iterate over all transcripts in the genome."""
        if not self._is_indexed:
            raise RuntimeError("The genome is not indexed. Call .index() after adding features.")
        for chrom in self.chromosomes.values():
            for gene in chrom.genes: # Access genes from chromosome, then transcripts from gene
                for transcript in gene.transcripts:
                    yield transcript

    @property
    def exons(self) -> Iterator[Exon]:
        """Iterate over all exons in the genome."""
        if not self._is_indexed:
            raise RuntimeError("The genome is not indexed. Call .index() after adding features.")
        for chrom in self.chromosomes.values():
            for gene in chrom.genes:
                for transcript in gene.transcripts:
                    for exon in transcript.exons:
                        yield exon

    def chromosome_by_id(self, chromosome_id: str) -> Chromosome:
        """Get a chromosome by its ID using the index. Raises ValueError if not found."""
        if not self._is_indexed:
            raise RuntimeError("The genome is not indexed. Call .index() after adding features.")
        try:
            return self.chromosomes[chromosome_id]
        except KeyError:
            raise ValueError(f"Chromosome with ID '{chromosome_id}' not found.")

    def gene_by_id(self, gene_id: str) -> Gene:
        """Get a gene by its ID using the index. Raises ValueError if not found."""
        if not self._is_indexed:
            raise RuntimeError("The genome is not indexed. Call .index() after adding features.")
        try:
            return self._genes_by_id[gene_id]
        except KeyError:
            raise ValueError(f"Gene with ID '{gene_id}' not found.")

    def transcript_by_id(self, transcript_id: str) -> Transcript:
        """Get a transcript by its ID using the index. Raises ValueError if not found."""
        if not self._is_indexed:
            raise RuntimeError("The genome is not indexed. Call .index() after adding features.")
        try:
            return self._transcripts_by_id[transcript_id]
        except KeyError:
            raise ValueError(f"Transcript with ID '{transcript_id}' not found.")

    def exon_by_id(self, exon_id: str) -> Exon:
        """Get an exon by its ID using the index. Raises ValueError if not found."""
        if not self._is_indexed:
            raise RuntimeError("The genome is not indexed. Call .index() after adding features.")
        try:
            return self._exons_by_id[exon_id]
        except KeyError:
            raise ValueError(f"Exon with ID '{exon_id}' not found.")


    def __iter__(self) -> Iterator[Chromosome]:
        """Iterate over all chromosomes in the genome."""
        if not self._is_indexed:
            raise RuntimeError("The genome is not indexed. Call .index() after adding features.")
        return iter(self.chromosomes.values())
    
    def __getattr__(self, name: str) -> Any:
        """Allow direct access to attributes in the attributes dictionary."""
        try:
            return self.attributes[name]
        except KeyError:
            raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")
        
    def __setattr__(self, name: str, value: Any):
        """Allow setting attributes. Explicitly defined attributes are set normally. New, dynamic attributes are stored in the 'attributes' dictionary."""
        if name in self.__dict__ or name in self.__class__.__dict__ or name == 'attributes' or not hasattr(self, 'attributes'):
            super().__setattr__(name, value)
        else:
            self.attributes[name] = value
            
    