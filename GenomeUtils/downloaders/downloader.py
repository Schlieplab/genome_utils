#!/usr/bin/env python
"""
Filename: GenomeUtils/downloaders/downloader.py
Author: Arash Ayat
Copyright: 2026, Alexander Schliep
Version: 0.2.0
Description: This file defines the base Downloader class for handling file downloads.
License: LGPL-3.0-or-later
"""

from abc import ABC
import logging
from pathlib import Path
import shutil
import tempfile

import requests


class Downloader(ABC):
    """Abstract base class for all downloaders."""

    def __init__(
        self,
        download_dir: Path | None = None,
        *,
        create_download_dir: bool = True,
    ):
        """Initialize the downloader.

        Args:
            download_dir: Directory for storing downloaded files.
                If None, uses a temporary directory.
            create_download_dir: Whether to create ``download_dir`` now.
        """
        self._is_temp_cache = download_dir is None
        self.download_dir = download_dir or Path(tempfile.mkdtemp())
        if create_download_dir:
            self.download_dir.mkdir(parents=True, exist_ok=True)
        self.logger = logging.getLogger(self.__class__.__name__)
        self._created_files: set[Path] = set()

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(download_dir={self.download_dir})"

    def download_file(
        self,
        url: str,
        filename: str | None = None,
        force: bool = False,
    ) -> Path:
        """Download a file into ``download_dir``.

        Args:
            url: The URL of the file to download.
            filename: The name of the file to be saved in the cache directory.
            force: If True, redownload the file even if it exists. Defaults to False.

        Returns:
            The path to the downloaded file.
        """
        if filename is None:
            filename = url.rsplit("/", 1)[-1].split("?", 1)[0]
        return self.download_file_to(
            url,
            self.download_dir / filename,
            force=force,
        )

    def download_file_to(
        self,
        url: str,
        destination: Path | str,
        *,
        force: bool = False,
    ) -> Path:
        """Download ``url`` to an exact caller-supplied destination.

        The download is written to a temporary file in the destination
        directory. On success, that file replaces the destination. On failure,
        the destination is left unchanged.
        """
        destination_path = Path(destination)
        destination_path.parent.mkdir(parents=True, exist_ok=True)
        existed_before_download = destination_path.exists()

        if not force and existed_before_download:
            self.logger.info(
                "File '%s' already exists. Skipping download.",
                destination_path,
            )
            return destination_path

        self.logger.info("Downloading to '%s'...", destination_path)
        temporary_file = tempfile.NamedTemporaryFile(
            mode="wb",
            dir=destination_path.parent,
            prefix=f".{destination_path.name}.",
            suffix=".part",
            delete=False,
        )
        temporary_path = Path(temporary_file.name)
        try:
            with temporary_file:
                with requests.get(url, stream=True) as response:
                    response.raise_for_status()
                    shutil.copyfileobj(response.raw, temporary_file)
            temporary_path.replace(destination_path)
        except Exception:
            temporary_path.unlink(missing_ok=True)
            raise

        if not existed_before_download:
            self._created_files.add(destination_path)
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

    def __del__(self):
        """
        Clean up temporary directory when the instance is garbage collected.
        Only removes the download directory if it was created as a temp dir.
        """
        is_temp_cache = getattr(self, "_is_temp_cache", False)
        download_dir = getattr(self, "download_dir", None)
        if is_temp_cache and download_dir is not None and download_dir.exists():
            try:
                shutil.rmtree(download_dir)
            except OSError:
                pass  # Ignore errors during cleanup (e.g. dir already removed)
