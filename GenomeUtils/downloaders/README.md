# Downloader

This directory contains modules for downloading data.

## Modules

- `downloader.py`: Provides a general-purpose `Downloader` class for downloading and caching files from URLs.
- `genome_downloader.py`: Contains the `EnsemblGenomeDownloader` class, which is a specialized downloader for fetching genome assemblies (DNA, cDNA, and annotations) from the Ensembl database using the `gget` library.

## Usage Example

Here's how to use `EnsemblGenomeDownloader` to download a human genome assembly:

```python
from pathlib import Path
from genome_utils.downloader import EnsemblGenomeDownloader

# Define genome parameters
assembly_id = 'GRCh38'
ensembl_release = 109
species = 'homo_sapiens'
genomes_root_dir = Path('./data/genomes')

# Initialize the downloader
downloader = EnsemblGenomeDownloader(
    assembly_id=assembly_id,
    ensembl_release=ensembl_release,
    species=species,
    genomes_root_dir=genomes_root_dir
)

# Download the genome data
genome_files = downloader.download()

# The returned `genome_files` dictionary contains paths to the downloaded files
print(genome_files)
```

The output will look something like this:

```
{
'dna':'data/genomes/ensembl/GRCh38/109/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz',
'cdna':'data/genomes/ensembl/GRCh38/109/Homo_sapiens.GRCh38.cdna.all.fa.gz',
'annotation':'data/genomes/ensembl/GRCh38/109/Homo_sapiens.GRCh38.109.gtf.gz'
}
```
