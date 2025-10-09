#!/usr/bin/env python
"""
Filename: GenomeUtils/Downloaders.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.0
Description: This file serves as a convenient entry point for downloader classes.
License: LGPL-3.0-or-later
"""


from .downloaders import Downloader, EnsemblGenomeDownloader

__all__ = ["EnsemblGenomeDownloader", "Downloader"]
