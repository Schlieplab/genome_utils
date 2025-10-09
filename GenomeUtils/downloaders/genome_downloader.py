#!/usr/bin/env python
"""
Filename: GenomeUtils/downloaders/genome_downloader.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.0
Description: This file defines the abstract base class for genome downloaders.
License: LGPL-3.0-or-later
"""

from __future__ import annotations

from pathlib import Path

import gget

from .downloader import Downloader


class EnsemblGenomeDownloader(Downloader):
    """
    Downloads genome data from Ensembl.

    This downloader fetches the download URLs
    for genomic data using `gget`, downloads the files, and stores them in `genomes_root_dir/ensembl/{assembly_id}/{ensembl_release}`.
    """

    def __init__(self, 
                 assembly_id: str, 
                 ensembl_release: int, 
                 species: str, 
                 genomes_root_dir: Path | str = Path('./data/genomes')
                 ):
        """
        Initializes the EnsemblGenomeDownloader.
        
        Args:
            assembly_id: The identifier for the genome assembly (e.g., 'GRCh38').
            ensembl_release: The release number of the Ensembl database.
            species: The scientific name for the species (e.g., 'homo_sapiens').
            genomes_root_dir: The parent directory to store all downloaded genomes. Defaults to './data/genomes'.
        """
        self.ensembl_release = ensembl_release
        self.assembly_id = assembly_id
        self.species = species
        self.genomes_root_dir = Path(genomes_root_dir)
        genome_dir = self.genomes_root_dir / 'ensembl' / assembly_id / str(ensembl_release)
        super().__init__(genome_dir)
    
    def __repr__(self) -> str:
        return (f"{self.__class__.__name__}("
                f"assembly_id={self.assembly_id}, "
                f"ensembl_release={self.ensembl_release}, "
                f"species={self.species}, "
                f"genomes_root_dir={self.genomes_root_dir})")

    def download(self) -> dict[str, Path]:
        """
        Downloads all necessary genome files using gget to retrieve the URLs.

        Returns:
            A dictionary mapping a file type to the local Path.
            Keys are `dna`, `cdna`, and `annotation`.
        """
        gtf_url, cdna_url, dna_url = tuple(
            gget.ref(self.species, 
                     which=["gtf", "cdna", "dna"], 
                     release=self.ensembl_release, 
                     ftp=True, 
                     verbose=False)
        )

        dna_path = self.download_file(dna_url, Path(dna_url).name)
        cdna_path = self.download_file(cdna_url, Path(cdna_url).name)
        annotation_path = self.download_file(gtf_url, Path(gtf_url).name)

        return {
            'dna': dna_path,
            'cdna': cdna_path,
            'annotation': annotation_path,
        } 