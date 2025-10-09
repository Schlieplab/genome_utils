#!/usr/bin/env python
"""
Filename: test_ensembl_genome_downloader.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.0
Description: Unit tests for the EnsemblGenomeDownloader class.
License: LGPL-3.0-or-later
"""

from pathlib import Path
from unittest.mock import Mock, patch

import pytest

from GenomeUtils.Downloaders import EnsemblGenomeDownloader


class TestEnsemblGenomeDownloader:
    """Test cases for the EnsemblGenomeDownloader class."""

    def test_ensembl_downloader_creation_default_dir(self):
        """Test EnsemblGenomeDownloader creation with default directory."""
        downloader = EnsemblGenomeDownloader(
            assembly_id="GRCh38",
            ensembl_release=110,
            species="homo_sapiens"
        )
        
        assert downloader.assembly_id == "GRCh38"
        assert downloader.ensembl_release == 110
        assert downloader.species == "homo_sapiens"
        assert downloader.genomes_root_dir == Path('./data/genomes')
        
        expected_dir = Path('./data/genomes/ensembl/GRCh38/110')
        assert downloader.download_dir == expected_dir

    def test_ensembl_downloader_creation_custom_dir(self):
        """Test EnsemblGenomeDownloader creation with custom directory."""
        custom_root = Path('/custom/genomes')
        
        with patch.object(Path, 'mkdir'):
            downloader = EnsemblGenomeDownloader(
                assembly_id="GRCm39",
                ensembl_release=105,
                species="mus_musculus",
                genomes_root_dir=custom_root
            )
        
        assert downloader.genomes_root_dir == custom_root
        expected_dir = custom_root / 'ensembl' / 'GRCm39' / '105'
        assert downloader.download_dir == expected_dir

    def test_ensembl_downloader_creation_string_path(self):
        """Test EnsemblGenomeDownloader creation with string path."""
        custom_root = '/string/path/genomes'
        
        with patch.object(Path, 'mkdir'):
            downloader = EnsemblGenomeDownloader(
                assembly_id="GRCh37",
                ensembl_release=75,
                species="homo_sapiens",
                genomes_root_dir=custom_root
            )
        
        assert downloader.genomes_root_dir == Path(custom_root)
        expected_dir = Path(custom_root) / 'ensembl' / 'GRCh37' / '75'
        assert downloader.download_dir == expected_dir

    def test_ensembl_downloader_repr(self):
        """Test EnsemblGenomeDownloader string representation."""
        with patch.object(Path, 'mkdir'):
            downloader = EnsemblGenomeDownloader(
                assembly_id="GRCh38",
                ensembl_release=110,
                species="homo_sapiens",
                genomes_root_dir=Path('/test/genomes')
            )
        
        expected = ("EnsemblGenomeDownloader("
                   "assembly_id=GRCh38, "
                   "ensembl_release=110, "
                   "species=homo_sapiens, "
                   "genomes_root_dir=/test/genomes)")
        assert repr(downloader) == expected

    @patch('gget.ref')
    @patch.object(EnsemblGenomeDownloader, 'download_file')
    def test_download_success(self, mock_download_file, mock_gget_ref):
        """Test successful download of genome files."""
        # Setup
        downloader = EnsemblGenomeDownloader(
            assembly_id="GRCh38",
            ensembl_release=110,
            species="homo_sapiens"
        )
        
        # Mock gget.ref response
        mock_urls = (
            "ftp://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz",
            "ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
            "ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        )
        mock_gget_ref.return_value = mock_urls
        
        # Mock download_file returns
        mock_paths = {
            'dna': Path('/path/to/dna.fa.gz'),
            'cdna': Path('/path/to/cdna.fa.gz'),
            'annotation': Path('/path/to/annotation.gtf.gz')
        }
        
        def mock_download_side_effect(url, filename):
            if 'gtf' in url:
                return mock_paths['annotation']
            elif 'cdna' in url:
                return mock_paths['cdna']
            elif 'dna' in url:
                return mock_paths['dna']
        
        mock_download_file.side_effect = mock_download_side_effect
        
        # Execute
        result = downloader.download()
        
        # Assertions
        mock_gget_ref.assert_called_once_with(
            "homo_sapiens",
            which=["gtf", "cdna", "dna"],
            release=110,
            ftp=True,
            verbose=False
        )
        
        assert mock_download_file.call_count == 3
        
        assert result['dna'] == mock_paths['dna']
        assert result['cdna'] == mock_paths['cdna']
        assert result['annotation'] == mock_paths['annotation']

    @patch('gget.ref')
    @patch.object(EnsemblGenomeDownloader, 'download_file')
    def test_download_different_species(self, mock_download_file, mock_gget_ref):
        """Test download for different species."""
        downloader = EnsemblGenomeDownloader(
            assembly_id="GRCm39",
            ensembl_release=105,
            species="mus_musculus"
        )
        
        mock_urls = (
            "ftp://ftp.ensembl.org/pub/release-105/gtf/mus_musculus/Mus_musculus.GRCm39.105.gtf.gz",
            "ftp://ftp.ensembl.org/pub/release-105/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz",
            "ftp://ftp.ensembl.org/pub/release-105/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
        )
        mock_gget_ref.return_value = mock_urls
        
        mock_download_file.return_value = Path('/mock/path')
        
        result = downloader.download()
        
        mock_gget_ref.assert_called_once_with(
            "mus_musculus",
            which=["gtf", "cdna", "dna"],
            release=105,
            ftp=True,
            verbose=False
        )
        
        assert 'dna' in result
        assert 'cdna' in result
        assert 'annotation' in result

    @patch('gget.ref')
    def test_download_gget_error(self, mock_gget_ref):
        """Test error handling when gget.ref fails."""
        downloader = EnsemblGenomeDownloader(
            assembly_id="Invalid",
            ensembl_release=999,
            species="invalid_species"
        )
        
        # Mock gget.ref to raise an exception
        mock_gget_ref.side_effect = Exception("Species not found")
        
        with pytest.raises(Exception, match="Species not found"):
            downloader.download()

    @patch('gget.ref')
    @patch.object(EnsemblGenomeDownloader, 'download_file')
    def test_download_file_error(self, mock_download_file, mock_gget_ref):
        """Test error handling when file download fails."""
        downloader = EnsemblGenomeDownloader(
            assembly_id="GRCh38",
            ensembl_release=110,
            species="homo_sapiens"
        )
        
        mock_urls = (
            "ftp://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/file.gtf.gz",
            "ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/file.fa.gz",
            "ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/file.fa.gz"
        )
        mock_gget_ref.return_value = mock_urls
        
        # Mock download_file to fail
        mock_download_file.side_effect = Exception("Download failed")
        
        with pytest.raises(Exception, match="Download failed"):
            downloader.download()

    @patch('gget.ref')
    @patch.object(EnsemblGenomeDownloader, 'download_file')
    def test_download_filename_extraction(self, mock_download_file, mock_gget_ref):
        """Test that filenames are correctly extracted from URLs."""
        downloader = EnsemblGenomeDownloader(
            assembly_id="GRCh38",
            ensembl_release=110,
            species="homo_sapiens"
        )
        
        mock_urls = (
            "ftp://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz",
            "ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz",
            "ftp://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
        )
        mock_gget_ref.return_value = mock_urls
        
        mock_download_file.return_value = Path('/mock/path')
        
        downloader.download()
        
        # Check that download_file was called with correct filenames
        calls = mock_download_file.call_args_list
        
        # Check GTF call
        gtf_call = next(call for call in calls if 'gtf' in str(call[0][0]))
        assert gtf_call[0][1] == "Homo_sapiens.GRCh38.110.gtf.gz"
        
        # Check cDNA call
        cdna_call = next(call for call in calls if 'cdna' in str(call[0][0]))
        assert cdna_call[0][1] == "Homo_sapiens.GRCh38.cdna.all.fa.gz"
        
        # Check DNA call
        dna_call = next(call for call in calls if 'dna' in str(call[0][0]))
        assert dna_call[0][1] == "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

    def test_inheritance_from_downloader(self):
        """Test that EnsemblGenomeDownloader properly inherits from Downloader."""
        from GenomeUtils.Downloaders import Downloader
        
        downloader = EnsemblGenomeDownloader(
            assembly_id="GRCh38",
            ensembl_release=110,
            species="homo_sapiens"
        )
        
        # Should be instance of base Downloader
        assert isinstance(downloader, Downloader)
        
        # Should have Downloader attributes
        assert hasattr(downloader, 'download_dir')
        assert hasattr(downloader, '_is_temp_cache')
        assert hasattr(downloader, '_created_files')
        assert hasattr(downloader, 'logger')
        
        # Should have Downloader methods
        assert hasattr(downloader, 'download_file')
        assert hasattr(downloader, 'cleanup')

    def test_directory_structure_creation(self):
        """Test that the correct directory structure is created."""
        # Test with various parameters
        test_cases = [
            {
                'assembly_id': 'GRCh38',
                'release': 110,
                'species': 'homo_sapiens',
                'root': Path('/genomes'),
                'expected': Path('/genomes/ensembl/GRCh38/110')
            },
            {
                'assembly_id': 'GRCm39',
                'release': 105,
                'species': 'mus_musculus',
                'root': Path('/data'),
                'expected': Path('/data/ensembl/GRCm39/105')
            },
            {
                'assembly_id': 'TAIR10',
                'release': 57,
                'species': 'arabidopsis_thaliana',
                'root': Path('./genomes'),
                'expected': Path('./genomes/ensembl/TAIR10/57')
            }
        ]
        
        for case in test_cases:
            with patch.object(Path, 'mkdir'):
                downloader = EnsemblGenomeDownloader(
                    assembly_id=case['assembly_id'],
                    ensembl_release=case['release'],
                    species=case['species'],
                    genomes_root_dir=case['root']
                )
                
                assert downloader.download_dir == case['expected']

    @patch('gget.ref')
    @patch.object(EnsemblGenomeDownloader, 'download_file')
    def test_download_return_format(self, mock_download_file, mock_gget_ref):
        """Test that download returns the expected dictionary format."""
        downloader = EnsemblGenomeDownloader(
            assembly_id="GRCh38",
            ensembl_release=110,
            species="homo_sapiens"
        )
        
        mock_urls = ("gtf_url", "cdna_url", "dna_url")
        mock_gget_ref.return_value = mock_urls
        
        mock_paths = [
            Path('/path/to/dna.fa.gz'),
            Path('/path/to/cdna.fa.gz'),
            Path('/path/to/annotation.gtf.gz')
        ]
        mock_download_file.side_effect = mock_paths
        
        result = downloader.download()
        
        # Check result structure
        assert isinstance(result, dict)
        assert set(result.keys()) == {'dna', 'cdna', 'annotation'}
        
        # Check that all values are Path objects
        for value in result.values():
            assert isinstance(value, Path)




