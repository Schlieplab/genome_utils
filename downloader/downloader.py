from typing import Optional, Set
from pathlib import Path
import tempfile
import shutil
import requests
import logging


class Downloader:
    """Abstract base class for all downloaders."""
    
    def __init__(self, download_dir: Optional[Path] = None):
        """      
        Initializes the Downloader.
        
        Args:
            download_dir: Directory for storing downloaded files.
                      If None, uses a temporary directory.
        """
        self._is_temp_cache = download_dir is None
        self.download_dir = download_dir or Path(tempfile.mkdtemp())
        self.download_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logging.getLogger(self.__class__.__name__)
        self._created_files: Set[Path] = set()
    
    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(download_dir={self.download_dir})"
    
    def download_file(self, url: str, filename: str = None, force: bool = False) -> Path:
        """
        Download a single file from a URL and saves it in the cache directory.

        Args:
            url: The URL of the file to download.
            filename: The name of the file to be saved in the cache directory.
            force: If True, redownload the file even if it exists. Defaults to False.

        Returns:
            The path to the downloaded file.
        """
        if filename is None:
            filename = url.split('/')[-1]
        destination_path = self.download_dir / filename
        
        if not force and destination_path.exists():
            self.logger.info(f"File '{filename}' already exists in cache. Skipping download.")
            return destination_path

        self.logger.info(f"Downloading '{filename}'...")
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(destination_path, 'wb') as f:
                shutil.copyfileobj(r.raw, f)
        self._created_files.add(destination_path)
        self.logger.info(f"Successfully downloaded '{filename}'.")
        return destination_path
    
    def cleanup(self):
        """
        Clean up created files.
        """
        if self._is_temp_cache:
            if self.download_dir.exists():
                shutil.rmtree(self.download_dir)
        else:
            for path in self._created_files:
                if path.exists():
                    path.unlink() 