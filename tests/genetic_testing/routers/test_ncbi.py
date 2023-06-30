"""Test NCBI Functions.

This module will test NCBI functionalities like search, summary, and fetch.

"""

from unittest.mock import MagicMock, patch

from genetic_testing.routers import ncbi


@patch("Bio.Entrez.esearch")
@patch("Bio.Entrez.Parser.DataHandler.read")
def test_search(mock_parser_read: MagicMock, mock_esearch: MagicMock):
    """Test NCBI search function

    Parameters
    ----------
    mock_parser_read : MagicMock
        A mock object representing the `DataHandler.read` method from the `Bio.Entrez.Parser` module.
    mock_esearch : MagicMock
        A mock object representing the `Bio.Entrez.esearch` function.
    """
    mock_esearch.return_value = {"IdList": ["1", "2", "3"]}
    mock_parser_read.return_value = {"IdList": ["1", "2", "3"]}
    assert ncbi.search("database", "term") == ["1", "2", "3"]
