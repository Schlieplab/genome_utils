#!/usr/bin/env python
"""
Filename: GenomeUtils/__init__.py
Author: Arash Ayat
Copyright: 2026, Alexander Schliep
Version: 0.2.0
Description: This file is the initialization file for the GenomeUtils package.
License: LGPL-3.0-or-later
"""


from . import Downloaders
from . import Genome

__version__ = "0.2.0"

__all__ = [
    "Genome",
    "Downloaders",
    "__version__",
]
