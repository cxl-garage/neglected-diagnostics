from unittest.mock import Mock

from app.common.data_processing import get_first_header


def test_first_header():
    assert get_first_header([]) is None

    mock_file = Mock()
    mock_file.read.return_value = b">header\nACGT\nheader2\nACGT"
    assert get_first_header([mock_file]) == "header"
    mock_file.read.assert_called_once()

    mock_file = Mock()
    mock_file.read.return_value = b""
    assert get_first_header([mock_file]) is None
    mock_file.read.assert_called_once()
