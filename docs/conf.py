from __future__ import annotations

import os
import sys
from datetime import datetime
from pathlib import Path

try:  # Python 3.11+
    import tomllib  # type: ignore[attr-defined]
except ModuleNotFoundError:  # pragma: no cover - fallback for <3.11
    import tomli as tomllib  # type: ignore[no-redef]

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

project = "GenomeUtils"
author = "Schlieplab"

release = ""
pyproject_path = PROJECT_ROOT / "pyproject.toml"
if pyproject_path.exists():
    with pyproject_path.open("rb") as f:
        pyproject_data = tomllib.load(f)
    release = pyproject_data.get("project", {}).get("version", "")
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
]

autosummary_generate = True

autodoc_mock_imports = [
    "Bio",
    "gffutils",
    "gget",
    "tqdm",
]

autodoc_typehints = "description"
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}

napoleon_google_docstring = True
napoleon_numpy_docstring = True

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_last_updated_fmt = "%Y-%m-%d"
html_title = "GenomeUtils Documentation"
