#!/usr/bin/env python
"""
Filename: GenomeUtils/downloaders/genome_downloader.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.2
Description: This file defines the abstract base class for genome downloaders.
License: LGPL-3.0-or-later
"""

from __future__ import annotations

from pathlib import Path

from .downloader import Downloader


class EnsemblGenomeDownloader(Downloader):
    """
    Downloads genome data from Ensembl FTP.

    This downloader constructs URLs directly from the Ensembl FTP layout and
    downloads the files, storing them in
    `genomes_root_dir/ensembl/{assembly_id}/{ensembl_release}`.
    """

    FTP_BASE = "https://ftp.ensembl.org/pub"

    def _build_urls(self) -> tuple[str, str, str]:
        """
        Build Ensembl FTP URLs for DNA, cDNA, and Annotation files.

        Uses the standard Ensembl FTP layout:
        - release-N/fasta/{species}/dna/{Species}.{Assembly}.dna.primary_assembly.fa.gz
        - release-N/fasta/{species}/cdna/{Species}.{Assembly}.cdna.all.fa.gz
        - release-N/gtf/{species}/{Species}.{Assembly}.{release}.gtf.gz

        Returns:
            Tuple of (dna_url, cdna_url, gtf_url).
        """
        parts = self.species.split("_")
        species_cap = parts[0].capitalize() + "_" + "_".join(p.lower() for p in parts[1:]) if len(parts) > 1 else parts[0].capitalize()
        release_path = f"release-{self.ensembl_release}"

        dna_filename = f"{species_cap}.{self.assembly_id}.dna.primary_assembly.fa.gz"
        cdna_filename = f"{species_cap}.{self.assembly_id}.cdna.all.fa.gz"
        gtf_filename = f"{species_cap}.{self.assembly_id}.{self.ensembl_release}.gtf.gz"

        dna_url = f"{self.FTP_BASE}/{release_path}/fasta/{self.species}/dna/{dna_filename}"
        cdna_url = f"{self.FTP_BASE}/{release_path}/fasta/{self.species}/cdna/{cdna_filename}"
        gtf_url = f"{self.FTP_BASE}/{release_path}/gtf/{self.species}/{gtf_filename}"

        return dna_url, cdna_url, gtf_url

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
            genomes_root_dir: The parent directory to store all downloaded genomes.
                Defaults to './data/genomes'.
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

    def get_urls(self) -> dict[str, str]:
        """
        Build the Ensembl FTP URLs for DNA, cDNA, and Annotation files.

        Returns:
            A dictionary with keys 'dna', 'cdna', 'annotation' mapping to URLs.
        """
        dna_url, cdna_url, gtf_url = self._build_urls()
        return {
            'dna': dna_url,
            'cdna': cdna_url,
            'annotation': gtf_url,
        }

    def download(self, force: bool = False) -> dict[str, Path]:
        """
        Download DNA, cDNA, and Annotation files from Ensembl FTP.

        Returns:
            A dictionary mapping a file type to the local Path.
            Keys are `dna`, `cdna`, and `annotation`.

        Args:
            force: If True, redownload the files even if they already exist. Defaults to False.
        """
        urls = self.get_urls()

        dna_path = self.download_file(urls['dna'], Path(urls['dna']).name, force=force)
        cdna_path = self.download_file(urls['cdna'], Path(urls['cdna']).name, force=force)
        annotation_path = self.download_file(urls['annotation'], Path(urls['annotation']).name, force=force)

        return {
            'dna': dna_path,
            'cdna': cdna_path,
            'annotation': annotation_path,
        }