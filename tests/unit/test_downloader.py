"""Tests for the Downloader base class."""

import pytest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock, patch, mock_open
import requests

from ...src.downloaders import Downloader


class TestDownloader:
    """Test cases for the Downloader base class."""

    def test_downloader_creation_default_dir(self):
        """Test downloader creation with default temporary directory."""
        with patch('tempfile.mkdtemp') as mock_mkdtemp:
            mock_mkdtemp.return_value = '/tmp/test_dir'
            
            downloader = Downloader()
            
            assert downloader.download_dir == Path('/tmp/test_dir')
            assert downloader._is_temp_cache == True
            mock_mkdtemp.assert_called_once()

    def test_downloader_creation_custom_dir(self):
        """Test downloader creation with custom directory."""
        custom_dir = Path('/custom/download/dir')
        
        with patch.object(Path, 'mkdir') as mock_mkdir:
            downloader = Downloader(download_dir=custom_dir)
            
            assert downloader.download_dir == custom_dir
            assert downloader._is_temp_cache == False
            mock_mkdir.assert_called_once_with(parents=True, exist_ok=True)

    def test_downloader_repr(self):
        """Test downloader string representation."""
        custom_dir = Path('/test/dir')
        
        with patch.object(Path, 'mkdir'):
            downloader = Downloader(download_dir=custom_dir)
            
            expected = "Downloader(download_dir=/test/dir)"
            assert repr(downloader) == expected

    @patch('requests.get')
    @patch('shutil.copyfileobj')
    def test_download_file_success(self, mock_copyfileobj, mock_get):
        """Test successful file download."""
        # Setup
        test_dir = Path('/test/dir')
        
        with patch.object(Path, 'mkdir'):
            downloader = Downloader(download_dir=test_dir)
        
        # Mock response as context manager
        mock_response = Mock()
        mock_response.raise_for_status.return_value = None
        mock_response.raw = Mock()
        mock_response.__enter__ = Mock(return_value=mock_response)
        mock_response.__exit__ = Mock(return_value=None)
        mock_get.return_value = mock_response
        
        # Mock file operations
        mock_file = mock_open()
        
        with patch('builtins.open', mock_file):
            with patch.object(Path, 'exists', return_value=False):
                result = downloader.download_file('http://example.com/file.txt', 'file.txt')
        
        # Assertions
        expected_path = test_dir / 'file.txt'
        assert result == expected_path
        assert expected_path in downloader._created_files
        
        mock_get.assert_called_once_with('http://example.com/file.txt', stream=True)
        mock_response.raise_for_status.assert_called_once()
        mock_copyfileobj.assert_called_once_with(mock_response.raw, mock_file())

    @patch('requests.get')
    @patch('shutil.copyfileobj')
    def test_download_file_auto_filename(self, mock_copyfileobj, mock_get):
        """Test download with automatic filename extraction."""
        test_dir = Path('/test/dir')
        
        with patch.object(Path, 'mkdir'):
            downloader = Downloader(download_dir=test_dir)
        
        # Mock response as context manager
        mock_response = Mock()
        mock_response.raise_for_status.return_value = None
        mock_response.raw = Mock()
        mock_response.__enter__ = Mock(return_value=mock_response)
        mock_response.__exit__ = Mock(return_value=None)
        mock_get.return_value = mock_response
        
        mock_file = mock_open()
        
        with patch('builtins.open', mock_file):
            with patch.object(Path, 'exists', return_value=False):
                result = downloader.download_file('http://example.com/path/auto_file.txt')
        
        expected_path = test_dir / 'auto_file.txt'
        assert result == expected_path

    def test_download_file_already_exists_no_force(self):
        """Test download when file exists and force=False."""
        test_dir = Path('/test/dir')
        
        with patch.object(Path, 'mkdir'):
            downloader = Downloader(download_dir=test_dir)
        
        with patch.object(Path, 'exists', return_value=True):
            with patch('requests.get') as mock_get:
                result = downloader.download_file('http://example.com/file.txt', 'file.txt')
        
        expected_path = test_dir / 'file.txt'
        assert result == expected_path
        mock_get.assert_not_called()  # Should not download if file exists

    @patch('requests.get')
    @patch('shutil.copyfileobj')
    def test_download_file_exists_with_force(self, mock_copyfileobj, mock_get):
        """Test download when file exists and force=True."""
        test_dir = Path('/test/dir')
        
        with patch.object(Path, 'mkdir'):
            downloader = Downloader(download_dir=test_dir)
        
        # Mock response as context manager
        mock_response = Mock()
        mock_response.raise_for_status.return_value = None
        mock_response.raw = Mock()
        mock_response.__enter__ = Mock(return_value=mock_response)
        mock_response.__exit__ = Mock(return_value=None)
        mock_get.return_value = mock_response
        
        mock_file = mock_open()
        
        with patch('builtins.open', mock_file):
            with patch.object(Path, 'exists', return_value=True):
                result = downloader.download_file(
                    'http://example.com/file.txt', 
                    'file.txt', 
                    force=True
                )
        
        expected_path = test_dir / 'file.txt'
        assert result == expected_path
        mock_get.assert_called_once()  # Should download even if file exists

    @patch('requests.get')
    def test_download_file_http_error(self, mock_get):
        """Test download with HTTP error."""
        test_dir = Path('/test/dir')
        
        with patch.object(Path, 'mkdir'):
            downloader = Downloader(download_dir=test_dir)
        
        # Mock response with error as context manager
        mock_response = Mock()
        mock_response.raise_for_status.side_effect = requests.HTTPError("404 Not Found")
        mock_response.__enter__ = Mock(return_value=mock_response)
        mock_response.__exit__ = Mock(return_value=None)
        mock_get.return_value = mock_response
        
        with pytest.raises(requests.HTTPError):
            downloader.download_file('http://example.com/nonexistent.txt', 'file.txt')

    def test_cleanup_temp_cache(self):
        """Test cleanup with temporary cache."""
        with patch('tempfile.mkdtemp') as mock_mkdtemp:
            with patch('shutil.rmtree') as mock_rmtree:
                mock_mkdtemp.return_value = '/tmp/test_dir'
                
                downloader = Downloader()  # Creates temp cache
                
                with patch.object(Path, 'exists', return_value=True):
                    downloader.cleanup()
                
                mock_rmtree.assert_called_once_with(Path('/tmp/test_dir'))

    def test_cleanup_persistent_cache(self):
        """Test cleanup with persistent cache."""
        test_dir = Path('/test/dir')
        
        with patch.object(Path, 'mkdir'):
            downloader = Downloader(download_dir=test_dir)
        
        # Add some files to created_files
        file1 = Path('/test/dir/file1.txt')
        file2 = Path('/test/dir/file2.txt')
        downloader._created_files.add(file1)
        downloader._created_files.add(file2)
        
        with patch.object(Path, 'exists', return_value=True) as mock_exists:
            with patch.object(Path, 'unlink') as mock_unlink:
                downloader.cleanup()
        
        # Should check existence and unlink each file
        assert mock_exists.call_count == 2
        assert mock_unlink.call_count == 2

    def test_cleanup_persistent_cache_nonexistent_files(self):
        """Test cleanup with files that don't exist."""
        test_dir = Path('/test/dir')
        
        with patch.object(Path, 'mkdir'):
            downloader = Downloader(download_dir=test_dir)
        
        # Add some files to created_files
        file1 = Path('/test/dir/file1.txt')
        downloader._created_files.add(file1)
        
        with patch.object(Path, 'exists', return_value=False):
            with patch.object(Path, 'unlink') as mock_unlink:
                downloader.cleanup()
        
        # Should not try to unlink non-existent files
        mock_unlink.assert_not_called()

    def test_created_files_tracking(self):
        """Test that created files are properly tracked."""
        test_dir = Path('/test/dir')
        
        with patch.object(Path, 'mkdir'):
            downloader = Downloader(download_dir=test_dir)
        
        # Initially empty
        assert len(downloader._created_files) == 0
        
        # Mock successful download
        with patch('requests.get') as mock_get:
            with patch('shutil.copyfileobj'):
                with patch('builtins.open', mock_open()):
                    with patch.object(Path, 'exists', return_value=False):
                        mock_response = Mock()
                        mock_response.raise_for_status.return_value = None
                        mock_response.raw = Mock()
                        mock_get.return_value.__enter__.return_value = mock_response
                        
                        result1 = downloader.download_file('http://example.com/file1.txt')
                        result2 = downloader.download_file('http://example.com/file2.txt')
        
        # Should track both files
        assert len(downloader._created_files) == 2
        assert result1 in downloader._created_files
        assert result2 in downloader._created_files

    def test_logger_initialization(self):
        """Test that logger is properly initialized."""
        with patch('tempfile.mkdtemp') as mock_mkdtemp:
            mock_mkdtemp.return_value = '/tmp/test_logger_dir'
            downloader = Downloader()
        
        assert hasattr(downloader, 'logger')
        assert downloader.logger.name == 'Downloader'

    @patch('requests.get')
    @patch('shutil.copyfileobj')
    def test_download_file_complex_url(self, mock_copyfileobj, mock_get):
        """Test download with complex URL."""
        test_dir = Path('/test/dir')
        
        with patch.object(Path, 'mkdir'):
            downloader = Downloader(download_dir=test_dir)
        
        complex_url = 'http://example.com/path/to/file.txt?param=value&other=123'
        
        mock_response = Mock()
        mock_response.raise_for_status.return_value = None
        mock_response.raw = Mock()
        mock_response.__enter__ = Mock(return_value=mock_response)
        mock_response.__exit__ = Mock(return_value=None)
        mock_get.return_value = mock_response
        
        with patch('builtins.open', mock_open()):
            with patch.object(Path, 'exists', return_value=False):
                result = downloader.download_file(complex_url)
        
        # Should extract just the filename
        expected_path = test_dir / 'file.txt'
        assert result == expected_path

    def test_download_dir_creation(self):
        """Test that download directory is created."""
        custom_dir = Path('/custom/new/dir')
        
        with patch.object(Path, 'mkdir') as mock_mkdir:
            Downloader(download_dir=custom_dir)
        
        mock_mkdir.assert_called_once_with(parents=True, exist_ok=True)

    def test_is_temp_cache_flag(self):
        """Test _is_temp_cache flag behavior."""
        # With default directory (should be temp)
        with patch('tempfile.mkdtemp') as mock_mkdtemp:
            mock_mkdtemp.return_value = '/tmp/test_temp_cache_dir'
            downloader1 = Downloader()
            assert downloader1._is_temp_cache == True
        
        # With custom directory (should not be temp)
        with patch.object(Path, 'mkdir'):
            downloader2 = Downloader(download_dir=Path('/custom'))
            assert downloader2._is_temp_cache == False




