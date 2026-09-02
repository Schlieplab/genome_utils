#!/usr/bin/env python
"""
Filename: GenomeUtils/downloaders/genome_downloader.py
Author: Arash Ayat
Copyright: 2026, Alexander Schliep
Version: 0.2.0
Description: This file defines the abstract base class for genome downloaders.
License: LGPL-3.0-or-later
"""

from __future__ import annotations

from pathlib import Path
from typing import Mapping

from ..genome.builder import create_gtf_database
from .downloader import Downloader


class EnsemblGenomeDownloader(Downloader):
    """
    Downloads genome data from Ensembl FTP.

    This downloader constructs URLs directly from the Ensembl FTP layout and
    downloads the files. By default they are stored in
    ``genomes_root_dir/ensembl/{assembly_id}/{ensembl_release}``. Passing
    ``destinations`` instead writes directly to those exact paths without
    creating the default directory hierarchy.
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

    def __init__(
        self,
        assembly_id: str,
        ensembl_release: int,
        species: str,
        genomes_root_dir: Path | str = Path("./data/genomes"),
        destinations: Mapping[str, Path | str] | None = None,
    ):
        """
        Initializes the EnsemblGenomeDownloader.

        Args:
            assembly_id: The identifier for the genome assembly (e.g., 'GRCh38').
            ensembl_release: The release number of the Ensembl database.
            species: The scientific name for the species (e.g., 'homo_sapiens').
            genomes_root_dir: The parent directory to store all downloaded genomes.
                Defaults to './data/genomes'.
            destinations: Exact paths for the downloaded artifacts. ``dna``,
                ``cdna``, and ``annotation`` are required. ``db`` is required
                only when :meth:`download` is called with ``output_db=True``.
        """
        self.ensembl_release = ensembl_release
        self.assembly_id = assembly_id
        self.species = species
        self.genomes_root_dir = Path(genomes_root_dir)
        self.destinations = self._normalize_destinations(destinations)
        genome_dir = (
            self.genomes_root_dir
            / "ensembl"
            / assembly_id
            / str(ensembl_release)
        )
        super().__init__(
            genome_dir,
            create_download_dir=self.destinations is None,
        )

    @staticmethod
    def _normalize_destinations(
        destinations: Mapping[str, Path | str] | None,
    ) -> dict[str, Path] | None:
        """Validate and normalize an explicit destination mapping."""
        if destinations is None:
            return None

        required = {"dna", "cdna", "annotation"}
        allowed = required | {"db"}
        provided = set(destinations)
        missing = required - provided
        unexpected = provided - allowed
        if missing:
            raise ValueError(
                "Missing required destinations: " + ", ".join(sorted(missing))
            )
        if unexpected:
            raise ValueError(
                "Unexpected destinations: " + ", ".join(sorted(unexpected))
            )

        normalized = {key: Path(value) for key, value in destinations.items()}
        resolved_paths = [path.resolve(strict=False) for path in normalized.values()]
        if len(resolved_paths) != len(set(resolved_paths)):
            raise ValueError("Destination paths must be unique after resolution")
        return normalized

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"assembly_id={self.assembly_id}, "
            f"ensembl_release={self.ensembl_release}, "
            f"species={self.species}, "
            f"genomes_root_dir={self.genomes_root_dir})"
        )

    def get_urls(self) -> dict[str, str]:
        """
        Build the Ensembl FTP URLs for DNA, cDNA, and Annotation files.

        Returns:
            A dictionary with keys 'dna', 'cdna', 'annotation' mapping to URLs.
        """
        dna_url, cdna_url, gtf_url = self._build_urls()
        return {
            "dna": dna_url,
            "cdna": cdna_url,
            "annotation": gtf_url,
        }

    def download(
        self,
        force: bool = False,
        output_db: bool = False,
    ) -> dict[str, Path]:
        """
        Download DNA, cDNA, and Annotation files from Ensembl FTP.

        Returns:
            A dictionary mapping a file type to the local Path.
            Keys are `dna`, `cdna`, and `annotation`. When ``output_db`` is
            true, the returned mapping also contains the annotation database
            path under the `db` key.

        Args:
            force: If True, redownload the files even if they already exist. Defaults to False.
            output_db: If True, create a reusable gffutils annotation database
                and include its path in the returned mapping.

                In explicit-destination mode, a ``db`` destination must be
                present when this is true and must be absent when this is false.
        """
        if self.destinations is not None:
            if output_db and "db" not in self.destinations:
                raise ValueError("The 'db' destination is required when output_db=True")
            if not output_db and "db" in self.destinations:
                raise ValueError("The 'db' destination requires output_db=True")

        urls = self.get_urls()

        if self.destinations is not None:
            for destination in self.destinations.values():
                destination.parent.mkdir(parents=True, exist_ok=True)

            paths = {
                key: self.download_file_to(
                    urls[key],
                    self.destinations[key],
                    force=force,
                )
                for key in ("dna", "cdna", "annotation")
            }
            if output_db:
                db_destination = self.destinations["db"]
                create_gtf_database(
                    paths["annotation"],
                    db_path=db_destination,
                    force=force,
                    keep_extracted_gtf=False,
                )
                paths["db"] = db_destination
            return paths

        paths = {
            key: self.download_file(url, Path(url).name, force=force)
            for key, url in urls.items()
        }
        if output_db:
            _, db_path = create_gtf_database(paths["annotation"], force=force)
            paths["db"] = db_path
        return paths
