"""
Sequence Fetching Module

This module implements Phase 1 of the roadmap: Database Integration and Sequence Fetching.
It provides enhanced capabilities for fetching sequences from multiple databases including
NCBI, BOLD, SILVA, and UNITE.

Features:
- Enhanced NCBI integration with taxonomy and marker-specific searches
- BOLD database integration for COI barcoding
- SILVA database integration for rRNA sequences
- UNITE database integration for fungal ITS sequences
- Quality filtering and sequence validation
- Async downloading capabilities
"""

from .async_downloader import AsyncSequenceDownloader, BatchSequenceFetcher
from .bold_client import BOLDClient
from .database_models import (
    DatabaseConfig,
    DatabaseType,
    MarkerType,
    MoleculeType,
    SearchParameters,
    SearchResult,
    SequenceRecord,
)
from .ncbi_client import NCBIClient
from .sequence_filters import SequenceFilter
from .silva_client import SILVAClient
from .unite_client import UNITEClient

__all__ = [
    "NCBIClient",
    "BOLDClient",
    "SILVAClient",
    "UNITEClient",
    "AsyncSequenceDownloader",
    "BatchSequenceFetcher",
    "SequenceRecord",
    "SearchParameters",
    "SearchResult",
    "DatabaseConfig",
    "DatabaseType",
    "MarkerType",
    "MoleculeType",
    "SequenceFilter",
]
