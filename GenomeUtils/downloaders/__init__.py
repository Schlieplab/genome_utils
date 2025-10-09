#!/usr/bin/env python
"""
Filename: __init__.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 1.0
Description: Initialization file for the downloaders package.
License: LGPL-3.0-or-later
"""


from .downloader import Downloader
from .genome_downloader import EnsemblGenomeDownloader

__all__ = [
    "Downloader",
    "EnsemblGenomeDownloader",
] 