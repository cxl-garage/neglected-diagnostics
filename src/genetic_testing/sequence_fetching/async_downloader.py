"""
Async Batch Downloader

This module provides asynchronous downloading capabilities for batch sequence retrieval
from multiple databases as specified in Phase 1.2 of the roadmap.
"""

import asyncio
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Dict, List, Optional, Union

import aiohttp

from utils.log import _init_logger

from .bold_client import BOLDClient
from .database_models import (
    DatabaseConfig,
    DatabaseType,
    SearchParameters,
    SearchResult,
)
from .ncbi_client import NCBIClient
from .silva_client import SILVAClient
from .unite_client import UNITEClient

logger = _init_logger(__name__)


class AsyncSequenceDownloader:
    """Asynchronous sequence downloader supporting multiple databases."""

    def __init__(self, config: Optional[DatabaseConfig] = None):
        """
        Initialize async downloader.

        Parameters
        ----------
        config : DatabaseConfig, optional
            Database configuration
        """
        self.config = config or DatabaseConfig()

        # Initialize database clients
        self.clients = {
            DatabaseType.NCBI_NUCLEOTIDE: NCBIClient(config),
            DatabaseType.NCBI_GENE: NCBIClient(config),
            DatabaseType.BOLD: BOLDClient(config),
            DatabaseType.SILVA: SILVAClient(config),
            DatabaseType.UNITE: UNITEClient(config),
        }

        # Thread pool for CPU-bound tasks
        self.thread_pool = ThreadPoolExecutor(max_workers=4)

    async def search_multiple_databases(
        self, search_params: List[SearchParameters]
    ) -> Dict[str, SearchResult]:
        """
        Search multiple databases concurrently.

        Parameters
        ----------
        search_params : List[SearchParameters]
            List of search parameters for different databases

        Returns
        -------
        Dict[str, SearchResult]
            Results keyed by database name
        """
        tasks = []

        for params in search_params:
            task = self._search_database_async(params)
            tasks.append(task)

        # Execute all searches concurrently
        results = await asyncio.gather(*tasks, return_exceptions=True)

        # Compile results
        database_results = {}

        for i, result in enumerate(results):
            db_name = search_params[i].database.value

            if isinstance(result, Exception):
                logger.error(f"Error searching {db_name}: {result}")
                database_results[db_name] = SearchResult(
                    query=search_params[i],
                    total_found=0,
                    sequences=[],
                    search_time=0.0,
                    warnings=[str(result)],
                )
            else:
                database_results[db_name] = result

        return database_results

    async def search_taxa_across_databases(
        self, taxa: List[str], databases: List[DatabaseType], max_per_taxon: int = 100
    ) -> Dict[str, Dict[str, SearchResult]]:
        """
        Search for multiple taxa across multiple databases.

        Parameters
        ----------
        taxa : List[str]
            List of taxon names to search
        databases : List[DatabaseType]
            Databases to search
        max_per_taxon : int
            Maximum sequences per taxon per database

        Returns
        -------
        Dict[str, Dict[str, SearchResult]]
            Results nested by taxon then database
        """
        tasks = []

        # Create search tasks for each taxon-database combination
        for taxon in taxa:
            for database in databases:
                params = SearchParameters(
                    search_term=taxon,
                    database=database,
                    max_results=max_per_taxon,
                    taxon=taxon,
                )

                task = self._search_database_async(params)
                tasks.append((taxon, database.value, task))

        # Execute all searches concurrently with rate limiting
        semaphore = asyncio.Semaphore(self.config.max_concurrent_requests)
        limited_tasks = [
            self._rate_limited_task(semaphore, task) for _, _, task in tasks
        ]

        results = await asyncio.gather(*limited_tasks, return_exceptions=True)

        # Organize results by taxon and database
        taxon_results = {taxon: {} for taxon in taxa}

        for i, result in enumerate(results):
            taxon, database, _ = tasks[i]

            if isinstance(result, Exception):
                logger.error(f"Error searching {database} for {taxon}: {result}")
                taxon_results[taxon][database] = SearchResult(
                    query=SearchParameters(
                        search_term=taxon,
                        database=DatabaseType(database),
                        max_results=max_per_taxon,
                        taxon=taxon,
                    ),
                    total_found=0,
                    sequences=[],
                    search_time=0.0,
                    warnings=[str(result)],
                )
            else:
                taxon_results[taxon][database] = result

        return taxon_results

    async def batch_download_sequences(
        self, search_results: List[SearchResult], include_sequence_data: bool = True
    ) -> List[SearchResult]:
        """
        Download sequence data for multiple search results.

        Parameters
        ----------
        search_results : List[SearchResult]
            Search results to download sequences for
        include_sequence_data : bool
            Whether to include actual sequence data

        Returns
        -------
        List[SearchResult]
            Updated search results with sequence data
        """
        if not include_sequence_data:
            return search_results

        tasks = []

        for result in search_results:
            if result.sequences:
                task = self._download_sequences_for_result(result)
                tasks.append(task)
            else:
                # No sequences to download
                tasks.append(asyncio.create_task(self._return_result(result)))

        # Execute downloads concurrently with rate limiting
        semaphore = asyncio.Semaphore(self.config.max_concurrent_requests)
        limited_tasks = [self._rate_limited_task(semaphore, task) for task in tasks]

        updated_results = await asyncio.gather(*limited_tasks, return_exceptions=True)

        # Handle any exceptions
        final_results = []
        for i, result in enumerate(updated_results):
            if isinstance(result, Exception):
                logger.error(f"Error downloading sequences for result {i}: {result}")
                # Return original result with warning
                original_result = search_results[i]
                original_result.warnings.append(
                    f"Sequence download failed: {str(result)}"
                )
                final_results.append(original_result)
            else:
                final_results.append(result)

        return final_results

    async def _search_database_async(self, params: SearchParameters) -> SearchResult:
        """Search a single database asynchronously."""
        client = self.clients.get(params.database)

        if not client:
            raise ValueError(f"Unsupported database: {params.database}")

        # Run the synchronous search in a thread pool
        loop = asyncio.get_event_loop()
        result = await loop.run_in_executor(
            self.thread_pool, client.search_sequences, params
        )

        return result

    async def _download_sequences_for_result(
        self, result: SearchResult
    ) -> SearchResult:
        """Download sequence data for a search result."""
        database = result.query.database
        client = self.clients.get(database)

        if not client or not hasattr(client, "fetch_sequences_with_data"):
            # Client doesn't support sequence download
            return result

        # Run sequence download in thread pool
        loop = asyncio.get_event_loop()
        updated_result = await loop.run_in_executor(
            self.thread_pool,
            client.fetch_sequences_with_data,
            result.query,
            True,  # include_sequences
        )

        return updated_result

    async def _rate_limited_task(self, semaphore: asyncio.Semaphore, task):
        """Execute task with rate limiting."""
        async with semaphore:
            # Add small delay for rate limiting
            await asyncio.sleep(1.0 / self.config.default_rate_limit)
            return await task

    async def _return_result(self, result: SearchResult) -> SearchResult:
        """Simple async wrapper to return a result."""
        return result

    def __del__(self):
        """Cleanup thread pool on deletion."""
        if hasattr(self, "thread_pool"):
            self.thread_pool.shutdown(wait=False)


class BatchSequenceFetcher:
    """High-level interface for batch sequence fetching operations."""

    def __init__(self, config: Optional[DatabaseConfig] = None):
        """
        Initialize batch fetcher.

        Parameters
        ----------
        config : DatabaseConfig, optional
            Database configuration
        """
        self.downloader = AsyncSequenceDownloader(config)

    async def fetch_marker_sequences_for_taxa(
        self,
        taxa: List[str],
        marker: str,
        databases: Optional[List[str]] = None,
        max_per_taxon: int = 100,
    ) -> Dict[str, Dict[str, SearchResult]]:
        """
        Fetch sequences for a specific marker across multiple taxa and databases.

        Parameters
        ----------
        taxa : List[str]
            List of taxon names
        marker : str
            Molecular marker (e.g., 'COI', '16S', 'ITS')
        databases : List[str], optional
            List of database names. If None, uses appropriate defaults for marker.
        max_per_taxon : int
            Maximum sequences per taxon per database

        Returns
        -------
        Dict[str, Dict[str, SearchResult]]
            Results nested by taxon then database
        """
        # Select appropriate databases for marker if not specified
        if not databases:
            databases = self._get_default_databases_for_marker(marker)

        # Convert string names to DatabaseType enums
        db_types = []
        for db_name in databases:
            try:
                db_types.append(DatabaseType(db_name))
            except ValueError:
                logger.warning(f"Unknown database: {db_name}")

        return await self.downloader.search_taxa_across_databases(
            taxa, db_types, max_per_taxon
        )

    async def fetch_comprehensive_dataset(
        self,
        target_taxa: List[str],
        off_target_taxa: List[str],
        marker: str,
        max_per_taxon: int = 500,
    ) -> Dict[str, Dict[str, SearchResult]]:
        """
        Fetch comprehensive dataset including target and off-target sequences.

        Parameters
        ----------
        target_taxa : List[str]
            Target taxonomic groups
        off_target_taxa : List[str]
            Off-target taxonomic groups
        marker : str
            Molecular marker
        max_per_taxon : int
            Maximum sequences per taxon

        Returns
        -------
        Dict[str, Dict[str, SearchResult]]
            Results with 'targets' and 'off_targets' keys
        """
        databases = self._get_default_databases_for_marker(marker)

        # Fetch target and off-target sequences concurrently
        target_task = self.fetch_marker_sequences_for_taxa(
            target_taxa, marker, databases, max_per_taxon
        )

        off_target_task = self.fetch_marker_sequences_for_taxa(
            off_target_taxa, marker, databases, max_per_taxon
        )

        target_results, off_target_results = await asyncio.gather(
            target_task, off_target_task
        )

        return {"targets": target_results, "off_targets": off_target_results}

    def _get_default_databases_for_marker(self, marker: str) -> List[str]:
        """Get default databases for a molecular marker."""
        marker_upper = marker.upper()

        if marker_upper in ["COI", "COX1"]:
            return ["nucleotide", "bold"]
        elif marker_upper in ["16S", "18S", "23S", "28S"]:
            return ["nucleotide", "silva"]
        elif marker_upper in ["ITS", "ITS1", "ITS2"]:
            return ["nucleotide", "unite"]
        else:
            # Default to NCBI nucleotide for unknown markers
            return ["nucleotide"]


# Convenience functions for common use cases


async def fetch_coi_barcodes(
    taxa: List[str], max_per_taxon: int = 100, config: Optional[DatabaseConfig] = None
) -> Dict[str, Dict[str, SearchResult]]:
    """
    Fetch COI barcode sequences from NCBI and BOLD.

    Parameters
    ----------
    taxa : List[str]
        List of taxon names
    max_per_taxon : int
        Maximum sequences per taxon
    config : DatabaseConfig, optional
        Database configuration

    Returns
    -------
    Dict[str, Dict[str, SearchResult]]
        COI sequences by taxon and database
    """
    fetcher = BatchSequenceFetcher(config)
    return await fetcher.fetch_marker_sequences_for_taxa(
        taxa, "COI", ["nucleotide", "bold"], max_per_taxon
    )


async def fetch_rrna_sequences(
    taxa: List[str],
    rrna_type: str = "16S",
    max_per_taxon: int = 100,
    config: Optional[DatabaseConfig] = None,
) -> Dict[str, Dict[str, SearchResult]]:
    """
    Fetch rRNA sequences from NCBI and SILVA.

    Parameters
    ----------
    taxa : List[str]
        List of taxon names
    rrna_type : str
        rRNA type ('16S', '18S', '23S', '28S')
    max_per_taxon : int
        Maximum sequences per taxon
    config : DatabaseConfig, optional
        Database configuration

    Returns
    -------
    Dict[str, Dict[str, SearchResult]]
        rRNA sequences by taxon and database
    """
    fetcher = BatchSequenceFetcher(config)
    return await fetcher.fetch_marker_sequences_for_taxa(
        taxa, rrna_type, ["nucleotide", "silva"], max_per_taxon
    )


async def fetch_its_sequences(
    taxa: List[str],
    its_region: str = "ITS",
    max_per_taxon: int = 100,
    config: Optional[DatabaseConfig] = None,
) -> Dict[str, Dict[str, SearchResult]]:
    """
    Fetch ITS sequences from NCBI and UNITE.

    Parameters
    ----------
    taxa : List[str]
        List of taxon names
    its_region : str
        ITS region ('ITS', 'ITS1', 'ITS2')
    max_per_taxon : int
        Maximum sequences per taxon
    config : DatabaseConfig, optional
        Database configuration

    Returns
    -------
    Dict[str, Dict[str, SearchResult]]
        ITS sequences by taxon and database
    """
    fetcher = BatchSequenceFetcher(config)
    return await fetcher.fetch_marker_sequences_for_taxa(
        taxa, its_region, ["nucleotide", "unite"], max_per_taxon
    )
