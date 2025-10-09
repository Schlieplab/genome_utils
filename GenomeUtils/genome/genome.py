#!/usr/bin/env python
"""
Filename: GenomeUtils/genome/genome.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.0
Description: This file defines the main Genome class, encapsulating chromosomes, genes, and transcripts.
License: LGPL-3.0-or-later
"""

from __future__ import annotations

from typing import Dict, List

from Bio.Seq import Seq

from .chromosome import Chromosome
from .exon import Exon
from .gene import Gene
from .locus import Locus
from .transcript import Transcript


class Genome:
    """Represents a Genome object, includes a collection of chromosomes, genes, transcripts, and exons."""

    def __init__(self, id: str, species: str, name: str, **kwargs):
        """
        Initializes the Genome object.
        
        Args:
            id: The ID of the genome.
            species: The species of the genome.
            name: The name of the genome.
            kwargs: Additional keyword arguments.
        """
        self.id = id
        self.species = species
        self.name = name

        for key, value in kwargs.items():
            setattr(self, key, value)
        
        self._chromosomes: Dict[str, Chromosome] = {}
        self._genes: Dict[str, Gene] = {}
        self._transcripts: Dict[str, Transcript] = {}
        self._exons: Dict[str, Exon] = {}
        self.is_indexed: bool = False


    def __repr__(self) -> str:
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
    
    def get_sequence_by_locus(self, locus: Locus) -> Seq:
        """Get a sequence by its locus."""
        if not self.is_indexed:
            raise RuntimeError("The genome is not indexed. Call .index() after adding features.")
        
        return self._chromosomes[locus.chr].get_subsequence_by_locus(locus)


    @property
    def chromosomes(self) -> List[Chromosome]:
        """Get all chromosomes in the genome."""
        return list(self._chromosomes.values())


    @property
    def genes(self) -> List[Gene]:
        """Get all genes in the genome."""
        return list(self._genes.values())

    
    @property
    def transcripts(self) -> List[Transcript]:
        """Get all transcripts in the genome."""
        return list(self._transcripts.values())


    @property
    def exons(self) -> List[Exon]:
        """Get all exons in the genome."""
        return list(self._exons.values())


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
