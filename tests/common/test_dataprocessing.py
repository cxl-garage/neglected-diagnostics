from unittest.mock import Mock

from app.common.data_processing import get_headers


def test_get_header():
    # empty cases
    assert get_headers([]) == []
    mock_file = Mock()
    mock_file.read.return_value = b""
    assert get_headers([mock_file]) == []

    mock_file1 = Mock()
    mock_file1.read.return_value = b">header3\nACGT\n>header2\nACGT"
    mock_file2 = Mock()
    mock_file2.read.return_value = b">header1\nACGT"
    assert get_headers([mock_file1, mock_file2]) == ["header1", "header2", "header3"]
