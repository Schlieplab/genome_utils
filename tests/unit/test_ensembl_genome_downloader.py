#!/usr/bin/env python
"""
Filename: tests/unit/test_ensembl_genome_downloader.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.2
Description: Unit tests for the EnsemblGenomeDownloader class.
License: LGPL-3.0-or-later
"""

from pathlib import Path
from unittest.mock import patch

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

    def test_get_urls_homo_sapiens(self):
        """Test URL construction for human genome."""
        downloader = EnsemblGenomeDownloader(
            assembly_id="GRCh38",
            ensembl_release=110,
            species="homo_sapiens"
        )
        urls = downloader.get_urls()

        assert "release-110" in urls['dna']
        assert "homo_sapiens" in urls['dna']
        assert "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" in urls['dna']

        assert "release-110" in urls['cdna']
        assert "Homo_sapiens.GRCh38.cdna.all.fa.gz" in urls['cdna']

        assert "release-110" in urls['annotation']
        assert "Homo_sapiens.GRCh38.110.gtf.gz" in urls['annotation']

    def test_get_urls_mus_musculus(self):
        """Test URL construction for mouse genome."""
        downloader = EnsemblGenomeDownloader(
            assembly_id="GRCm39",
            ensembl_release=105,
            species="mus_musculus"
        )
        urls = downloader.get_urls()

        assert "Mus_musculus.GRCm39" in urls['dna']
        assert "Mus_musculus.GRCm39" in urls['cdna']
        assert "Mus_musculus.GRCm39.105.gtf.gz" in urls['annotation']

    @patch.object(EnsemblGenomeDownloader, 'download_file')
    def test_download_success(self, mock_download_file):
        """Test successful download of genome files."""
        downloader = EnsemblGenomeDownloader(
            assembly_id="GRCh38",
            ensembl_release=110,
            species="homo_sapiens"
        )

        mock_paths = {
            'dna': Path('/path/to/dna.fa.gz'),
            'cdna': Path('/path/to/cdna.fa.gz'),
            'annotation': Path('/path/to/annotation.gtf.gz')
        }

        def mock_download_side_effect(url, filename, force=False):
            if 'gtf' in url:
                return mock_paths['annotation']
            elif 'cdna' in url:
                return mock_paths['cdna']
            elif 'dna' in url:
                return mock_paths['dna']

        mock_download_file.side_effect = mock_download_side_effect

        result = downloader.download()

        assert mock_download_file.call_count == 3
        assert all(call.kwargs == {'force': False} for call in mock_download_file.call_args_list)
        assert result['dna'] == mock_paths['dna']
        assert result['cdna'] == mock_paths['cdna']
        assert result['annotation'] == mock_paths['annotation']

    @patch.object(EnsemblGenomeDownloader, 'download_file')
    def test_download_different_species(self, mock_download_file):
        """Test download for different species."""
        downloader = EnsemblGenomeDownloader(
            assembly_id="GRCm39",
            ensembl_release=105,
            species="mus_musculus"
        )

        mock_download_file.return_value = Path('/mock/path')

        result = downloader.download()

        assert 'dna' in result
        assert 'cdna' in result
        assert 'annotation' in result
        assert mock_download_file.call_count == 3

    @patch.object(EnsemblGenomeDownloader, 'download_file')
    def test_download_file_error(self, mock_download_file):
        """Test error handling when file download fails."""
        downloader = EnsemblGenomeDownloader(
            assembly_id="GRCh38",
            ensembl_release=110,
            species="homo_sapiens"
        )

        mock_download_file.side_effect = Exception("Download failed")

        with pytest.raises(Exception, match="Download failed"):
            downloader.download()

    @patch.object(EnsemblGenomeDownloader, 'download_file')
    def test_download_filename_extraction(self, mock_download_file):
        """Test that filenames are correctly extracted from URLs."""
        downloader = EnsemblGenomeDownloader(
            assembly_id="GRCh38",
            ensembl_release=110,
            species="homo_sapiens"
        )

        mock_download_file.return_value = Path('/mock/path')

        downloader.download()

        calls = mock_download_file.call_args_list

        gtf_call = next(call for call in calls if 'gtf' in str(call[0][0]))
        assert gtf_call[0][1] == "Homo_sapiens.GRCh38.110.gtf.gz"

        cdna_call = next(call for call in calls if 'cdna' in str(call[0][0]))
        assert cdna_call[0][1] == "Homo_sapiens.GRCh38.cdna.all.fa.gz"

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

    @patch.object(EnsemblGenomeDownloader, 'download_file')
    def test_download_return_format(self, mock_download_file):
        """Test that download returns the expected dictionary format."""
        downloader = EnsemblGenomeDownloader(
            assembly_id="GRCh38",
            ensembl_release=110,
            species="homo_sapiens"
        )

        mock_paths = [
            Path('/path/to/dna.fa.gz'),
            Path('/path/to/cdna.fa.gz'),
            Path('/path/to/annotation.gtf.gz')
        ]
        mock_download_file.side_effect = mock_paths

        result = downloader.download()

        assert isinstance(result, dict)
        assert set(result.keys()) == {'dna', 'cdna', 'annotation'}

        for value in result.values():
            assert isinstance(value, Path)

    def test_build_urls(self):
        """Test URL building."""
        downloader = EnsemblGenomeDownloader(
            assembly_id="GRCh38",
            ensembl_release=115,
            species="homo_sapiens"
        )
        dna_url, cdna_url, gtf_url = downloader._build_urls()

        assert "release-115" in dna_url
        assert "fasta/homo_sapiens/dna" in dna_url
        assert "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" in dna_url

        assert "fasta/homo_sapiens/cdna" in cdna_url
        assert "Homo_sapiens.GRCh38.cdna.all.fa.gz" in cdna_url

        assert "gtf/homo_sapiens" in gtf_url
        assert "Homo_sapiens.GRCh38.115.gtf.gz" in gtf_url



