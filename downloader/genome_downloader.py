from abc import ABC, abstractmethod
from typing import Dict, Optional
from pathlib import Path

from .downloader import Downloader


class GenomeDownloader(Downloader, ABC):
    """
    Abstract base class for source-specific genome downloaders.

    This class uses the "Template Method" design pattern. It defines the
    skeleton of the download algorithm in the `download` method, and
    defers the implementation of specific steps to its subclasses.
    """
    def __init__(
        self, 
        assembly_id: str, 
        species: str, 
        genomes_root_dir: Path = Path.home() / "genomes"
    ):
        """
        Initializes the source-specific genome downloader.
        
        Args:
            assembly_id: The identifier for the genome assembly (e.g., 'GRCh38').
            species: The scientific name for the species (e.g., 'homo_sapiens').
            genomes_root_dir: The parent directory to store all downloaded genomes.
        """
        self.assembly_id = assembly_id
        self.species = species
        
        genome_dir = genomes_root_dir / assembly_id
        super().__init__(genome_dir)
    
    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(assembly_id={self.assembly_id}, species={self.species}, genome_dir={self.genome_dir})"
        
    def download(self) -> Dict[str, Path]:
        """
        Template method that downloads all necessary genome files.

        This method orchestrates the download process by calling the specific
        download methods for DNA, cDNA, and annotations.

        Returns:
            A dictionary mapping a file type to the local Path.
            Keys are 'dna', 'cdna', and 'annotation'.
        """
        dna_path = self._download_dna()
        cdna_path = self._download_cdna()
        annotation_path = self._download_annotations()

        return {
            'dna': dna_path,
            'cdna': cdna_path,
            'annotation': annotation_path,
        }

    @property
    def genome_dir(self) -> Path:
        return self.cache_dir
     
    @abstractmethod
    def _download_dna(self) -> Path:
        """
        Downloads the primary assembly DNA FASTA file.

        Returns:
            The path to the downloaded DNA FASTA file.
        """
        pass

    @abstractmethod
    def _download_cdna(self) -> Path:
        """
        Downloads the transcript/cDNA FASTA file.

        Returns:
            The path to the downloaded cDNA FASTA file.
        """
        pass

    @abstractmethod
    def _download_annotations(self) -> Path:
        """
        Downloads the genome annotation file (GFF/GTF).

        Returns:
            The path to the downloaded annotation file.
        """
        pass 