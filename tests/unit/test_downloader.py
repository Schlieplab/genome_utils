"""Focused tests for the generic Downloader exact-destination behavior."""

from io import BytesIO
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import requests

from GenomeUtils.Downloaders import Downloader


def _response(raw: object) -> MagicMock:
    response = MagicMock()
    response.raw = raw
    response.__enter__.return_value = response
    return response


def test_downloader_creates_configured_directory(tmp_path: Path) -> None:
    download_dir = tmp_path / "downloads"

    Downloader(download_dir)

    assert download_dir.is_dir()


def test_download_file_to_reuses_existing_destination(tmp_path: Path) -> None:
    destination = tmp_path / "nested" / "dna.fa.gz"
    destination.parent.mkdir()
    destination.write_bytes(b"existing")
    downloader = Downloader(tmp_path / "cache")

    with patch("GenomeUtils.downloaders.downloader.requests.get") as get:
        result = downloader.download_file_to("https://example.test/dna", destination)

    assert result == destination
    assert destination.read_bytes() == b"existing"
    get.assert_not_called()


def test_download_file_to_is_atomic_on_failure(tmp_path: Path) -> None:
    class FailingStream:
        def __init__(self) -> None:
            self.read_count = 0

        def read(self, _: int = -1) -> bytes:
            self.read_count += 1
            if self.read_count == 1:
                return b"partial"
            raise RuntimeError("stream interrupted")

    destination = tmp_path / "output" / "dna.fa.gz"
    downloader = Downloader(tmp_path / "cache")

    with patch(
        "GenomeUtils.downloaders.downloader.requests.get",
        return_value=_response(FailingStream()),
    ):
        with pytest.raises(RuntimeError, match="stream interrupted"):
            downloader.download_file_to("https://example.test/dna", destination)

    assert not destination.exists()
    assert not list(destination.parent.glob(f".{destination.name}.*.part"))


def test_failed_forced_download_preserves_caller_owned_file(tmp_path: Path) -> None:
    destination = tmp_path / "dna.fa.gz"
    destination.write_bytes(b"caller owned")
    response = _response(BytesIO(b"replacement"))
    response.raise_for_status.side_effect = requests.HTTPError("not found")
    downloader = Downloader(tmp_path / "cache")

    with patch(
        "GenomeUtils.downloaders.downloader.requests.get",
        return_value=response,
    ):
        with pytest.raises(requests.HTTPError, match="not found"):
            downloader.download_file_to(
                "https://example.test/dna",
                destination,
                force=True,
            )

    downloader.cleanup()
    assert destination.read_bytes() == b"caller owned"


def test_cleanup_does_not_remove_overwritten_caller_file(tmp_path: Path) -> None:
    destination = tmp_path / "dna.fa.gz"
    destination.write_bytes(b"old")
    downloader = Downloader(tmp_path / "cache")

    with patch(
        "GenomeUtils.downloaders.downloader.requests.get",
        return_value=_response(BytesIO(b"new")),
    ):
        downloader.download_file_to(
            "https://example.test/dna",
            destination,
            force=True,
        )

    downloader.cleanup()
    assert destination.read_bytes() == b"new"
