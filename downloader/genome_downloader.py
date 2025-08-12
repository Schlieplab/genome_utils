from __future__ import annotations
from pathlib import Path
import logging
import gget

from .downloader import Downloader


class GgetEnsemblGenomeDownloader(Downloader):
    """
    Downloads genome data using the `gget` library.

    This downloader leverages `gget` to fetch the download URLs
    for genomic data, and then uses the base Downloader to manage the
    file transfer and caching.
    """

    def __init__(self, assembly_id: str, ensembl_release: int, species: str, genomes_root_dir: Path = Path('./data')):
        """
        Initializes the GgetEnsemblGenomeDownloader.
        
        Args:
            assembly_id: The identifier for the genome assembly (e.g., 'GRCh38').
            species: The scientific name for the species (e.g., 'homo_sapiens').
            genomes_root_dir: The parent directory to store all downloaded genomes.
        """
        self.ensembl_release = ensembl_release
        self.species = species
        
        genome_dir = genomes_root_dir / assembly_id / str(ensembl_release)
        super().__init__(genome_dir)
        self.logger = logging.getLogger(self.__class__.__name__)
    
    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(assembly_id={self.assembly_id}, ensembl_release={self.ensembl_release}, species={self.species}, genomes_root_dir={self.genomes_root_dir})"

    def download(self) -> dict[str, Path]:
        """
        Downloads all necessary genome files using gget.



        Returns:
            A dictionary mapping a file type to the local Path.
            Keys are 'dna', 'cdna', and 'annotation'.
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