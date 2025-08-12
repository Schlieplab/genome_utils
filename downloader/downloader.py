from typing import Optional, Dict, Any
from pathlib import Path
import tempfile
import shutil
import requests
import logging


class Downloader:
    """Abstract base class for all downloaders."""
    
    def __init__(self, cache_dir: Optional[Path] = None):
        """
        Initialize the downloader.
        
        Args:
            cache_dir: Optional directory for caching downloaded files.
                      If None, uses a temporary directory.
        """
        self.cache_dir = cache_dir or Path(tempfile.mkdtemp())
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logging.getLogger(self.__class__.__name__)
    
    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(cache_dir={self.cache_dir})"
    
    def download_file(self, url: str, filename: str, force: bool = False) -> Path:
        """
        Download a single file from a URL and saves it in the cache directory.

        Args:
            url: The URL of the file to download.
            filename: The name of the file to be saved in the cache directory.
            force: If True, redownload the file even if it exists. Defaults to False.

        Returns:
            The path to the downloaded file.
        """
        destination_path = self.cache_dir / filename
        
        if not force and destination_path.exists():
            self.logger.info(f"File '{filename}' already exists in cache. Skipping download.")
            return destination_path

        self.logger.info(f"Downloading '{filename}'...")
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            with open(destination_path, 'wb') as f:
                shutil.copyfileobj(r.raw, f)
        self.logger.info(f"Successfully downloaded '{filename}'.")
        return destination_path
    
    def cleanup(self):
        """Clean up any temporary files/directories."""
        if self.cache_dir.exists():
            shutil.rmtree(self.cache_dir) 