"""
BOLD Database Client

This module provides integration with the Barcode of Life Database (BOLD) for COI 
barcode sequences as specified in Phase 1.2 of the roadmap.
"""

import asyncio
import time
from io import StringIO
from typing import Dict, List, Optional
from urllib.parse import urlencode

import aiohttp
import requests
import requests_cache
from Bio import SeqIO

from utils.log import _init_logger

from .database_models import (
    DatabaseConfig,
    DatabaseType,
    MarkerType,
    SearchParameters,
    SearchResult,
    SequenceRecord,
)

logger = _init_logger(__name__)


class BOLDClient:
    """Client for accessing BOLD (Barcode of Life Database) sequences."""

    def __init__(self, config: Optional[DatabaseConfig] = None):
        """
        Initialize BOLD client.

        Parameters
        ----------
        config : DatabaseConfig, optional
            Database configuration
        """
        self.config = config or DatabaseConfig()
        self.base_url = self.config.bold_base_url
        self.rate_limit = self.config.bold_rate_limit

        # Setup caching if enabled
        if self.config.cache_enabled:
            self.session = requests_cache.CachedSession(
                "bold_cache", expire_after=self.config.cache_expire_hours * 3600
            )
        else:
            self.session = requests.Session()

        # Rate limiting
        self.last_request_time = 0

    def search_sequences(self, params: SearchParameters) -> SearchResult:
        """
        Search BOLD database for sequences.

        Parameters
        ----------
        params : SearchParameters
            Search parameters

        Returns
        -------
        SearchResult
            Search results
        """
        start_time = time.time()

        try:
            # Build BOLD API query
            api_params = self._build_bold_query(params)

            # Make API request
            records = self._fetch_bold_records(api_params, params.max_results)

            # Convert to SequenceRecord objects
            sequences = self._parse_bold_records(records)

            search_time = time.time() - start_time

            return SearchResult(
                query=params,
                total_found=len(sequences),
                sequences=sequences,
                search_time=search_time,
                database_version="BOLD Public API",
            )

        except Exception as e:
            logger.error(f"Error searching BOLD database: {e}", exc_info=True)
            return SearchResult(
                query=params,
                total_found=0,
                sequences=[],
                search_time=time.time() - start_time,
                warnings=[f"BOLD search failed: {str(e)}"],
            )

    def fetch_coi_sequences(
        self, taxon: str, max_results: int = 1000, include_sequences: bool = True
    ) -> SearchResult:
        """
        Fetch COI barcode sequences for a specific taxon.

        Parameters
        ----------
        taxon : str
            Target taxon name
        max_results : int
            Maximum number of sequences to retrieve
        include_sequences : bool
            Whether to include actual sequence data

        Returns
        -------
        SearchResult
            COI sequences from BOLD
        """
        params = SearchParameters(
            search_term=taxon,
            database=DatabaseType.BOLD,
            max_results=max_results,
            taxon=taxon,
            marker=MarkerType.COI,
        )

        result = self.search_sequences(params)

        # Fetch actual sequences if requested
        if include_sequences and result.sequences:
            self._add_sequence_data(result.sequences)

        return result

    def _build_bold_query(self, params: SearchParameters) -> Dict[str, str]:
        """Build BOLD API query parameters."""
        query_params = {}

        # Taxonomic search
        if params.taxon:
            query_params["taxon"] = params.taxon

        # Marker filter (BOLD is primarily COI)
        if params.marker and params.marker in [MarkerType.COI, MarkerType.COX1]:
            query_params["marker"] = "COI-5P"

        # Geographic filter
        if params.country:
            query_params["geo"] = params.country

        # Output format
        query_params["format"] = "tsv"

        return query_params

    def _fetch_bold_records(
        self, api_params: Dict[str, str], max_results: int
    ) -> List[Dict]:
        """Fetch records from BOLD API."""
        # Rate limiting
        self._enforce_rate_limit()

        # BOLD specimen data endpoint
        url = f"{self.base_url}specimen"

        try:
            response = self.session.get(url, params=api_params, timeout=30)
            response.raise_for_status()

            # Parse TSV response
            records = self._parse_tsv_response(response.text)

            # Limit results
            if len(records) > max_results:
                records = records[:max_results]

            return records

        except requests.RequestException as e:
            logger.error(f"Error fetching BOLD records: {e}")
            raise

    def _parse_tsv_response(self, tsv_text: str) -> List[Dict]:
        """Parse TSV response from BOLD API."""
        lines = tsv_text.strip().split("\n")
        if len(lines) < 2:
            return []

        # Extract headers
        headers = lines[0].split("\t")
        records = []

        for line in lines[1:]:
            values = line.split("\t")
            if len(values) == len(headers):
                record = dict(zip(headers, values))
                records.append(record)

        return records

    def _parse_bold_records(self, records: List[Dict]) -> List[SequenceRecord]:
        """Convert BOLD records to SequenceRecord objects."""
        sequences = []

        for record in records:
            try:
                # Extract key information
                process_id = record.get("processid", "")
                specimen_id = record.get("sampleid", "")

                sequence_record = SequenceRecord(
                    id=process_id or specimen_id,
                    accession=record.get("genbank_accession", ""),
                    title=self._build_title(record),
                    organism=self._get_organism_name(record),
                    taxonomy_id=None,  # BOLD doesn't provide NCBI taxonomy IDs directly
                    length=self._parse_length(record.get("nucleotides", "0")),
                    database=DatabaseType.BOLD,
                    marker=MarkerType.COI,
                    country=record.get("country", ""),
                    lat_lon=self._format_coordinates(record),
                    specimen_voucher=record.get("voucher", ""),
                    metadata={
                        "bold_id": record.get("bold_id", ""),
                        "institution": record.get("institution_storing", ""),
                        "collection_date": record.get("collection_date", ""),
                        "collectors": record.get("collectors", ""),
                        "identification": record.get("identification", ""),
                        "identification_method": record.get(
                            "identification_method", ""
                        ),
                        "life_stage": record.get("life_stage", ""),
                        "sex": record.get("sex", ""),
                        "reproduction": record.get("reproduction", ""),
                        "habitat": record.get("habitat", ""),
                        "sampling_protocol": record.get("sampling_protocol", ""),
                        "notes": record.get("notes", ""),
                    },
                )

                sequences.append(sequence_record)

            except Exception as e:
                logger.warning(
                    f"Error parsing BOLD record {record.get('processid', 'unknown')}: {e}"
                )
                continue

        return sequences

    def _build_title(self, record: Dict) -> str:
        """Build sequence title from BOLD record."""
        parts = []

        if record.get("species_name"):
            parts.append(record["species_name"])

        if record.get("marker_code"):
            parts.append(f"gene, {record['marker_code']}")

        if record.get("processid"):
            parts.append(f"specimen {record['processid']}")

        return " ".join(parts) if parts else "BOLD sequence"

    def _get_organism_name(self, record: Dict) -> Optional[str]:
        """Extract organism name from BOLD record."""
        return record.get("species_name") or record.get("genus_name")

    def _parse_length(self, length_str: str) -> int:
        """Parse sequence length from string."""
        try:
            return int(length_str) if length_str else 0
        except ValueError:
            return 0

    def _format_coordinates(self, record: Dict) -> Optional[str]:
        """Format latitude/longitude coordinates."""
        lat = record.get("lat")
        lon = record.get("lon")

        if lat and lon:
            try:
                return f"{float(lat):.6f},{float(lon):.6f}"
            except ValueError:
                pass

        return None

    def _add_sequence_data(self, sequences: List[SequenceRecord]):
        """Fetch actual sequence data for BOLD records."""
        # BOLD sequence data endpoint
        for sequence in sequences:
            if not sequence.id:
                continue

            try:
                # Rate limiting
                self._enforce_rate_limit()

                # Fetch sequence using BOLD's sequence API
                url = f"{self.base_url}sequence"
                params = {"ids": sequence.id, "format": "fasta"}

                response = self.session.get(url, params=params, timeout=30)
                response.raise_for_status()

                # Parse FASTA response
                if response.text.strip():
                    fasta_io = StringIO(response.text)
                    seq_records = list(SeqIO.parse(fasta_io, "fasta"))

                    if seq_records:
                        sequence.sequence = str(seq_records[0].seq)
                        # Update length with actual sequence length
                        sequence.length = len(sequence.sequence)

            except Exception as e:
                logger.warning(
                    f"Error fetching sequence data for BOLD ID {sequence.id}: {e}"
                )
                continue

    def _enforce_rate_limit(self):
        """Enforce rate limiting for API requests."""
        current_time = time.time()
        time_since_last = current_time - self.last_request_time
        min_interval = 1.0 / self.rate_limit

        if time_since_last < min_interval:
            sleep_time = min_interval - time_since_last
            time.sleep(sleep_time)

        self.last_request_time = time.time()

    async def fetch_sequences_async(
        self, taxa: List[str], max_per_taxon: int = 100
    ) -> Dict[str, SearchResult]:
        """
        Asynchronously fetch sequences for multiple taxa.

        Parameters
        ----------
        taxa : List[str]
            List of taxon names
        max_per_taxon : int
            Maximum sequences per taxon

        Returns
        -------
        Dict[str, SearchResult]
            Results keyed by taxon name
        """
        async with aiohttp.ClientSession() as session:
            tasks = []

            for taxon in taxa:
                task = self._fetch_taxon_async(session, taxon, max_per_taxon)
                tasks.append(task)

            results = await asyncio.gather(*tasks, return_exceptions=True)

            # Compile results
            taxon_results = {}
            for i, result in enumerate(results):
                if isinstance(result, Exception):
                    logger.error(f"Error fetching {taxa[i]}: {result}")
                    taxon_results[taxa[i]] = SearchResult(
                        query=SearchParameters(
                            search_term=taxa[i],
                            database=DatabaseType.BOLD,
                            max_results=max_per_taxon,
                        ),
                        total_found=0,
                        sequences=[],
                        search_time=0.0,
                        warnings=[str(result)],
                    )
                else:
                    taxon_results[taxa[i]] = result

            return taxon_results

    async def _fetch_taxon_async(
        self, session: aiohttp.ClientSession, taxon: str, max_results: int
    ) -> SearchResult:
        """Asynchronously fetch sequences for a single taxon."""
        start_time = time.time()

        try:
            # Build query
            params = {"taxon": taxon, "marker": "COI-5P", "format": "tsv"}

            url = f"{self.base_url}specimen"

            # Add rate limiting
            await asyncio.sleep(1.0 / self.rate_limit)

            async with session.get(url, params=params, timeout=30) as response:
                response.raise_for_status()
                text = await response.text()

                # Parse response
                records = self._parse_tsv_response(text)

                if len(records) > max_results:
                    records = records[:max_results]

                sequences = self._parse_bold_records(records)

                return SearchResult(
                    query=SearchParameters(
                        search_term=taxon,
                        database=DatabaseType.BOLD,
                        max_results=max_results,
                        taxon=taxon,
                        marker=MarkerType.COI,
                    ),
                    total_found=len(sequences),
                    sequences=sequences,
                    search_time=time.time() - start_time,
                    database_version="BOLD Public API",
                )

        except Exception as e:
            logger.error(f"Error in async BOLD fetch for {taxon}: {e}")
            raise
