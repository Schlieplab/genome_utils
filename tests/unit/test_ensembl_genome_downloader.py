#!/usr/bin/env python
"""
Filename: tests/unit/test_ensembl_genome_downloader.py
Author: Arash Ayat
Copyright: 2026, Alexander Schliep
Description: Unit tests for the EnsemblGenomeDownloader class.
License: LGPL-3.0-or-later
"""

import gzip
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

    def test_legacy_construction_creates_nested_directory(self, tmp_path):
        root = tmp_path / "genomes"

        EnsemblGenomeDownloader(
            assembly_id="GRCh38",
            ensembl_release=114,
            species="homo_sapiens",
            genomes_root_dir=root,
        )

        assert (root / "ensembl" / "GRCh38" / "114").is_dir()

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

        assert set(urls) == {"dna", "cdna", "annotation"}
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
        urls = downloader.get_urls(outputs=["dna", "cdna", "ncrna", "annotation"])

        assert "Mus_musculus.GRCm39" in urls['dna']
        assert "Mus_musculus.GRCm39" in urls['cdna']
        assert urls["ncrna"] == (
            "https://ftp.ensembl.org/pub/release-105/fasta/mus_musculus/ncrna/"
            "Mus_musculus.GRCm39.ncrna.fa.gz"
        )
        assert "Mus_musculus.GRCm39.105.gtf.gz" in urls['annotation']

    @pytest.mark.parametrize(
        "outputs",
        [
            ["dna", "cdna", "annotation"],
            ["dna", "cdna", "annotation", "ncrna"],
            ["ncrna"],
            ["annotation", "dna"],
        ],
    )
    @patch.object(EnsemblGenomeDownloader, 'download_file')
    def test_download_success(self, mock_download_file, outputs):
        """Test successful download of genome files."""
        downloader = EnsemblGenomeDownloader(
            assembly_id="GRCh38",
            ensembl_release=110,
            species="homo_sapiens"
        )

        mock_paths = {
            'dna': Path('/path/to/dna.fa.gz'),
            'cdna': Path('/path/to/cdna.fa.gz'),
            'ncrna': Path('/path/to/ncrna.fa.gz'),
            'annotation': Path('/path/to/annotation.gtf.gz')
        }

        def mock_download_side_effect(url, filename, force=False):
            if 'gtf' in url:
                return mock_paths['annotation']
            elif 'cdna' in url:
                return mock_paths['cdna']
            elif 'ncrna' in url:
                return mock_paths['ncrna']
            elif 'dna' in url:
                return mock_paths['dna']

        mock_download_file.side_effect = mock_download_side_effect

        result = downloader.download(outputs=outputs)

        assert mock_download_file.call_count == len(outputs)
        assert all(call.kwargs == {'force': False} for call in mock_download_file.call_args_list)
        assert result == {key: mock_paths[key] for key in outputs}
        assert list(result) == outputs

    @patch('GenomeUtils.downloaders.genome_downloader.create_gtf_database')
    @patch.object(EnsemblGenomeDownloader, 'download_file')
    def test_download_can_output_database_path(
        self,
        mock_download_file,
        mock_create_gtf_database,
        tmp_path,
    ):
        """output_db adds the reusable annotation database to returned paths."""
        downloader = EnsemblGenomeDownloader(
            assembly_id="GRCh38",
            ensembl_release=110,
            species="homo_sapiens",
            genomes_root_dir=tmp_path,
        )
        dna_path = tmp_path / "dna.fa.gz"
        cdna_path = tmp_path / "cdna.fa.gz"
        annotation_path = tmp_path / "annotation.gtf.gz"
        db_path = tmp_path / "annotation.gtf.db"
        mock_download_file.side_effect = [dna_path, cdna_path, annotation_path]
        mock_create_gtf_database.return_value = (Mock(), db_path)

        result = downloader.download(force=True, output_db=True)

        mock_create_gtf_database.assert_called_once_with(annotation_path, force=True)
        assert result["db"] == db_path

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
        assert 'ncrna' not in result
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
        urls = downloader.get_urls(outputs=["dna", "cdna", "ncrna", "annotation"])
        dna_url = urls["dna"]
        cdna_url = urls["cdna"]
        ncrna_url = urls["ncrna"]
        gtf_url = urls["annotation"]

        assert "release-115" in dna_url
        assert "fasta/homo_sapiens/dna" in dna_url
        assert "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" in dna_url

        assert "fasta/homo_sapiens/cdna" in cdna_url
        assert "Homo_sapiens.GRCh38.cdna.all.fa.gz" in cdna_url

        assert ncrna_url == (
            "https://ftp.ensembl.org/pub/release-115/fasta/homo_sapiens/ncrna/"
            "Homo_sapiens.GRCh38.ncrna.fa.gz"
        )

        assert "gtf/homo_sapiens" in gtf_url
        assert "Homo_sapiens.GRCh38.115.gtf.gz" in gtf_url

    @pytest.mark.parametrize("outputs", [[], "dna", ["db"], ["unknown"], ["dna", "dna"]])
    def test_invalid_outputs_fail_before_downloading(self, tmp_path, outputs):
        downloader = EnsemblGenomeDownloader(
            "GRCh38", 110, "homo_sapiens", genomes_root_dir=tmp_path,
        )

        with patch.object(downloader, "download_file") as download_file:
            with pytest.raises(ValueError):
                downloader.download(outputs=outputs)

        download_file.assert_not_called()

    def test_database_requires_annotation_output(self, tmp_path):
        downloader = EnsemblGenomeDownloader(
            "GRCh38", 110, "homo_sapiens", genomes_root_dir=tmp_path,
        )

        with patch.object(downloader, "download_file") as download_file:
            with pytest.raises(ValueError, match="'annotation' output is required"):
                downloader.download(outputs=["ncrna"], output_db=True)

        download_file.assert_not_called()

    def test_selected_output_requires_explicit_destination(self, tmp_path):
        output_root = tmp_path / "untouched"
        downloader = EnsemblGenomeDownloader(
            "GRCh38", 110, "homo_sapiens",
            destinations={"dna": output_root / "dna.fa.gz"},
        )

        with patch.object(downloader, "download_file_to") as download_file_to:
            with pytest.raises(ValueError, match="Missing required destinations: ncrna"):
                downloader.download(outputs=["ncrna"])

        download_file_to.assert_not_called()
        assert not output_root.exists()

    @pytest.mark.parametrize("outputs", [None, ["ncrna"]])
    def test_explicit_destinations_download_only_selected_outputs(self, tmp_path, outputs):
        destinations = {
            key: tmp_path / key / "file.gz"
            for key in ("dna", "cdna", "annotation", "ncrna")
        }
        downloader = EnsemblGenomeDownloader(
            "GRCh38", 110, "homo_sapiens", destinations=destinations,
        )
        selected = ["dna", "cdna", "annotation"] if outputs is None else outputs
        expected = {key: destinations[key] for key in selected}

        with patch.object(
            downloader, "download_file_to", side_effect=lambda url, path, **kwargs: path,
        ) as download_file_to:
            result = downloader.download() if outputs is None else downloader.download(outputs=outputs)

        assert result == expected
        assert download_file_to.call_count == len(selected)
        for key, path in destinations.items():
            assert path.parent.exists() == (key in selected)

    @pytest.mark.parametrize("include_ncrna", [False, True])
    def test_explicit_destinations_write_exact_paths_without_legacy_layout(
        self,
        tmp_path,
        sample_gtf_content,
        include_ncrna,
    ):
        """Explicit mode writes only its declared artifacts."""
        destinations = {
            "dna": str(tmp_path / "reference" / "genome.fa.gz"),
            "cdna": tmp_path / "transcripts" / "all.fa.gz",
            "annotation": tmp_path / "annotation" / "genes.gtf.gz",
            "db": tmp_path / "database" / "features.sqlite",
        }
        if include_ncrna:
            destinations["ncrna"] = tmp_path / "transcripts" / "ncrna.fa.gz"
        legacy_root = tmp_path / "legacy-root"
        downloader = EnsemblGenomeDownloader(
            assembly_id="GRCh38",
            ensembl_release=114,
            species="homo_sapiens",
            genomes_root_dir=legacy_root,
            destinations=destinations,
        )
        assert not legacy_root.exists()

        def write_download(_url, destination, *, force=False):
            path = Path(destination)
            if path == Path(destinations["annotation"]):
                with gzip.open(path, "wt") as gtf_file:
                    gtf_file.write(sample_gtf_content)
            else:
                path.write_bytes(b"fasta")
            return path

        with patch.object(
            downloader,
            "download_file_to",
            side_effect=write_download,
        ) as download_file_to:
            outputs = ["dna", "cdna", "annotation"]
            if include_ncrna:
                outputs.append("ncrna")
            paths = downloader.download(force=True, output_db=True, outputs=outputs)

        expected = {key: Path(path) for key, path in destinations.items()}
        assert paths == expected
        assert all(path.exists() for path in expected.values())
        assert download_file_to.call_count == 3 + include_ncrna
        assert not legacy_root.exists()
        assert not expected["annotation"].with_suffix("").exists()

    @pytest.mark.parametrize(
        ("destinations", "message"),
        [
            (
                {
                    "dna": "dna",
                    "cdna": "cdna",
                    "annotation": "annotation",
                    "extra": "extra",
                },
                "Unexpected destinations: extra",
            ),
            (
                {
                    "dna": "same",
                    "cdna": "same",
                    "annotation": "annotation",
                },
                "unique after resolution",
            ),
        ],
    )
    def test_invalid_explicit_destinations_fail_during_construction(
        self,
        tmp_path,
        destinations,
        message,
    ):
        """Invalid mappings fail before the constructor creates directories."""
        root = tmp_path / "untouched"
        destination_paths = {
            key: root / "outputs" / value
            for key, value in destinations.items()
        }
        with pytest.raises(ValueError, match=message):
            EnsemblGenomeDownloader(
                assembly_id="GRCh38",
                ensembl_release=114,
                species="homo_sapiens",
                genomes_root_dir=root,
                destinations=destination_paths,
            )
        assert not root.exists()

    @pytest.mark.parametrize(
        ("include_db", "output_db", "message"),
        [
            (False, True, "required when output_db=True"),
            (True, False, "requires output_db=True"),
        ],
    )
    def test_db_destination_contract_fails_before_downloading(
        self,
        tmp_path,
        include_db,
        output_db,
        message,
    ):
        output_root = tmp_path / "untouched"
        destinations = {
            "dna": output_root / "dna.fa.gz",
            "cdna": output_root / "cdna.fa.gz",
            "annotation": output_root / "annotation.gtf.gz",
        }
        if include_db:
            destinations["db"] = output_root / "annotation.db"
        downloader = EnsemblGenomeDownloader(
            assembly_id="GRCh38",
            ensembl_release=114,
            species="homo_sapiens",
            destinations=destinations,
        )

        with patch.object(downloader, "download_file_to") as download_file_to:
            with pytest.raises(ValueError, match=message):
                downloader.download(output_db=output_db)

        download_file_to.assert_not_called()
        assert not output_root.exists()
