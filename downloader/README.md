# Genome Downloader Package

This package provides a flexible framework for downloading genomic data from various sources using well-established design patterns.

## Package Structure

```
downloader/
├── __init__.py                 # Package exports
├── downloader.py              # Base Downloader ABC
├── genome_downloader.py       # GenomeDownloader ABC with Template Method
├── ensembl_downloader.py      # Ensembl-specific implementation
├── ncbi_downloader.py         # NCBI-specific implementation
├── example_usage.py           # Usage examples and demos
└── README.md                  # This file
```

## Design Patterns Implemented

### 1. Template Method Pattern
- **Location**: `GenomeDownloader.download_genome()`
- **Purpose**: Defines a consistent 6-step algorithm for all genome downloads
- **Steps**:
  1. Download genomic DNA sequences
  2. Download cDNA sequences  
  3. Download annotations (GTF/GFF)
  4. Create Genome object
  5. Parse and add genomic features
  6. Index the genome

### 2. Factory Method Pattern
- **Location**: Each concrete downloader's `_create_genome()` method
- **Purpose**: Creates `Genome` objects with source-specific metadata
- **Benefits**: Encapsulates genome creation logic per data source

### 3. Strategy Pattern
- **Location**: `EnsemblDownloader` vs `NCBIDownloader`
- **Purpose**: Different algorithms for different data sources
- **Benefits**: Easy to add new sources without modifying existing code

### 4. Abstract Base Class Pattern
- **Location**: `Downloader` and `GenomeDownloader`
- **Purpose**: Enforces interface contracts
- **Benefits**: Guarantees all implementations have required methods

## Usage

### Basic Usage

```python
from genome_utils import EnsemblDownloader, NCBIDownloader

# Download from Ensembl (includes both DNA and cDNA automatically)
ensembl = EnsemblDownloader(ensembl_release="109")
genome = ensembl.download_genome(
    genome_id="GRCh38",
    species="homo_sapiens"
)

# Download from NCBI (includes both DNA and RNA automatically)
ncbi = NCBIDownloader()
genome = ncbi.download_genome(
    genome_id="GCF_000001405.40",
    species="Homo sapiens",
    assembly_name="GRCh38.p14"
)

# Both genomes now have cdna_path attribute
print(f"cDNA available at: {genome.cdna_path}")
```

### Adding New Data Sources

To add a new data source, inherit from `GenomeDownloader` and implement four abstract methods:

```python
class CustomDownloader(GenomeDownloader):
    def _download_dna(self, genome_id, species, **kwargs):
        # Download genomic DNA sequences
        pass
    
    def _download_cdna(self, genome_id, species, **kwargs):
        # Download cDNA sequences
        pass
    
    def _download_annotations(self, genome_id, species, **kwargs):
        # Download annotation data  
        pass
    
    def _create_genome(self, genome_id, species, dna_path, cdna_path, **kwargs):
        # Create Genome object with both DNA and cDNA paths
        pass
    
    def _parse_and_add_features(self, genome, annotation_path, **kwargs):
        # Parse annotations and add features
        pass
```

## Benefits

- **Extensible**: Add new data sources easily
- **Maintainable**: Clear separation of concerns
- **Consistent**: Same interface for all sources
- **Robust**: Built-in error handling and cleanup
- **Cacheable**: Automatic file caching to avoid re-downloads

## Examples

See `example_usage.py` for comprehensive examples demonstrating:
- Downloading from different sources (with automatic DNA + cDNA)
- Comparing genome data across sources
- Creating custom downloaders
- Using the design patterns effectively
- Accessing both genomic DNA and cDNA sequences 