"""
SILVA Database Client

This module provides integration with the SILVA rRNA database for 16S/18S/23S/28S 
sequences as specified in Phase 1.2 of the roadmap.
"""

import re
import time
from io import StringIO
from typing import Dict, List, Optional
from urllib.parse import urljoin

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


class SILVAClient:
    """Client for accessing SILVA rRNA database sequences."""

    def __init__(self, config: Optional[DatabaseConfig] = None):
        """
        Initialize SILVA client.

        Parameters
        ----------
        config : DatabaseConfig, optional
            Database configuration
        """
        self.config = config or DatabaseConfig()
        self.base_url = self.config.silva_base_url
        self.version = self.config.silva_version

        # Setup caching if enabled
        if self.config.cache_enabled:
            self.session = requests_cache.CachedSession(
                "silva_cache", expire_after=self.config.cache_expire_hours * 3600
            )
        else:
            self.session = requests.Session()

        # SILVA-specific URLs (these may need updating based on SILVA's current API)
        self.download_urls = {
            "SSU": f"{self.base_url}fileadmin/silva_databases/release_{self.version}/Exports/",
            "LSU": f"{self.base_url}fileadmin/silva_databases/release_{self.version}/Exports/",
        }

        # Marker to SILVA database mapping
        self.marker_mapping = {
            MarkerType.SSU_16S: "SSU",
            MarkerType.SSU_18S: "SSU",
            MarkerType.LSU_23S: "LSU",
            MarkerType.LSU_28S: "LSU",
        }

    def search_sequences(self, params: SearchParameters) -> SearchResult:
        """
        Search SILVA database for rRNA sequences.

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
            # Determine which SILVA dataset to use
            silva_type = self._get_silva_type(params.marker)

            if not silva_type:
                return SearchResult(
                    query=params,
                    total_found=0,
                    sequences=[],
                    search_time=time.time() - start_time,
                    warnings=["Marker not supported by SILVA database"],
                )

            # Download and parse SILVA data
            sequences = self._fetch_silva_sequences(silva_type, params)

            # Filter sequences based on parameters
            filtered_sequences = self._filter_sequences(sequences, params)

            search_time = time.time() - start_time

            return SearchResult(
                query=params,
                total_found=len(sequences),
                sequences=filtered_sequences[: params.max_results],
                search_time=search_time,
                database_version=f"SILVA {self.version}",
            )

        except Exception as e:
            logger.error(f"Error searching SILVA database: {e}", exc_info=True)
            return SearchResult(
                query=params,
                total_found=0,
                sequences=[],
                search_time=time.time() - start_time,
                warnings=[f"SILVA search failed: {str(e)}"],
            )

    def fetch_16s_sequences(
        self, domain: str = "Bacteria", max_results: int = 1000
    ) -> SearchResult:
        """
        Fetch 16S rRNA sequences for a specific domain.

        Parameters
        ----------
        domain : str
            Taxonomic domain (Bacteria, Archaea, Eukaryota)
        max_results : int
            Maximum number of sequences

        Returns
        -------
        SearchResult
            16S rRNA sequences
        """
        params = SearchParameters(
            search_term=domain,
            database=DatabaseType.SILVA,
            max_results=max_results,
            marker=MarkerType.SSU_16S,
            taxon=domain,
        )

        return self.search_sequences(params)

    def fetch_18s_sequences(
        self, taxon: str = "Eukaryota", max_results: int = 1000
    ) -> SearchResult:
        """
        Fetch 18S rRNA sequences for eukaryotes.

        Parameters
        ----------
        taxon : str
            Taxonomic group
        max_results : int
            Maximum number of sequences

        Returns
        -------
        SearchResult
            18S rRNA sequences
        """
        params = SearchParameters(
            search_term=taxon,
            database=DatabaseType.SILVA,
            max_results=max_results,
            marker=MarkerType.SSU_18S,
            taxon=taxon,
        )

        return self.search_sequences(params)

    def _get_silva_type(self, marker: Optional[MarkerType]) -> Optional[str]:
        """Determine SILVA database type from marker."""
        if not marker:
            return None
        return self.marker_mapping.get(marker)

    def _fetch_silva_sequences(
        self, silva_type: str, params: SearchParameters
    ) -> List[SequenceRecord]:
        """Fetch sequences from SILVA database files."""
        # This is a simplified implementation
        # In practice, you might need to download specific SILVA files
        # or use their web services if available

        sequences = []

        try:
            # Construct download URL for SILVA FASTA file
            # Note: Actual URLs may vary based on SILVA's current structure
            filename = self._get_silva_filename(silva_type, params)
            download_url = urljoin(self.download_urls[silva_type], filename)

            logger.info(f"Attempting to download SILVA data from: {download_url}")

            # Download SILVA data
            response = self.session.get(download_url, timeout=300)  # 5 minute timeout

            if response.status_code == 200:
                # Parse FASTA content
                sequences = self._parse_silva_fasta(response.text, params.database)
            else:
                logger.warning(
                    f"Could not download SILVA data: HTTP {response.status_code}"
                )

        except Exception as e:
            logger.error(f"Error fetching SILVA sequences: {e}")
            # Fall back to mock data or cached results

        return sequences

    def _get_silva_filename(self, silva_type: str, params: SearchParameters) -> str:
        """Generate SILVA filename based on parameters."""
        # This is a simplified mapping - actual filenames vary
        if silva_type == "SSU":
            if params.marker == MarkerType.SSU_16S:
                return f"SILVA_{self.version}_SSURef_Nr99_tax_silva.fasta.gz"
            else:  # 18S
                return f"SILVA_{self.version}_SSURef_Nr99_tax_silva.fasta.gz"
        else:  # LSU
            return f"SILVA_{self.version}_LSURef_tax_silva.fasta.gz"

    def _parse_silva_fasta(
        self, fasta_content: str, database: DatabaseType
    ) -> List[SequenceRecord]:
        """Parse SILVA FASTA content into SequenceRecord objects."""
        sequences = []

        try:
            fasta_io = StringIO(fasta_content)
            seq_records = SeqIO.parse(fasta_io, "fasta")

            for seq_record in seq_records:
                # Parse SILVA header format
                # Example: >AACY020068177.1.1441 Bacteria;Proteobacteria;Gammaproteobacteria;...
                header_parts = seq_record.description.split(" ", 1)
                accession = header_parts[0]
                taxonomy = header_parts[1] if len(header_parts) > 1 else ""

                # Extract organism name from taxonomy
                organism = self._extract_organism_from_taxonomy(taxonomy)

                # Determine marker type
                marker = self._infer_marker_from_length(len(seq_record.seq))

                sequence_record = SequenceRecord(
                    id=accession,
                    accession=accession,
                    title=seq_record.description,
                    organism=organism,
                    sequence=str(seq_record.seq),
                    length=len(seq_record.seq),
                    database=database,
                    marker=marker,
                    metadata={"taxonomy": taxonomy, "silva_version": self.version},
                )

                sequences.append(sequence_record)

        except Exception as e:
            logger.error(f"Error parsing SILVA FASTA: {e}")

        return sequences

    def _extract_organism_from_taxonomy(self, taxonomy: str) -> Optional[str]:
        """Extract organism name from SILVA taxonomy string."""
        if not taxonomy:
            return None

        # SILVA taxonomy format: Domain;Phylum;Class;Order;Family;Genus;Species
        parts = taxonomy.split(";")

        # Try to get genus and species
        if len(parts) >= 7:
            genus = parts[-2].strip()
            species = parts[-1].strip()

            if genus and species and not genus.startswith("uncultured"):
                return f"{genus} {species}"

        # Fall back to genus
        if len(parts) >= 6:
            genus = parts[-2].strip()
            if genus and not genus.startswith("uncultured"):
                return genus

        return None

    def _infer_marker_from_length(self, length: int) -> Optional[MarkerType]:
        """Infer marker type from sequence length."""
        if 1200 <= length <= 1800:
            return MarkerType.SSU_16S
        elif 1600 <= length <= 2000:
            return MarkerType.SSU_18S
        elif 2500 <= length <= 3500:
            return MarkerType.LSU_23S
        elif 3000 <= length <= 4500:
            return MarkerType.LSU_28S

        return None

    def _filter_sequences(
        self, sequences: List[SequenceRecord], params: SearchParameters
    ) -> List[SequenceRecord]:
        """Filter sequences based on search parameters."""
        filtered = []

        for seq in sequences:
            # Taxonomic filter
            if params.taxon:
                if not self._matches_taxon(seq, params.taxon):
                    continue

            # Length filters
            if params.min_length and seq.length < params.min_length:
                continue
            if params.max_length and seq.length > params.max_length:
                continue

            # Marker filter
            if params.marker and seq.marker != params.marker:
                continue

            filtered.append(seq)

        return filtered

    def _matches_taxon(self, sequence: SequenceRecord, taxon: str) -> bool:
        """Check if sequence matches taxonomic filter."""
        if not sequence.metadata.get("taxonomy"):
            return False

        taxonomy = sequence.metadata["taxonomy"].lower()
        taxon_lower = taxon.lower()

        return taxon_lower in taxonomy

    def get_available_taxa(self, marker: MarkerType) -> List[str]:
        """
        Get list of available taxa for a specific marker.

        Parameters
        ----------
        marker : MarkerType
            Target marker

        Returns
        -------
        List[str]
            Available taxonomic groups
        """
        # This would require parsing SILVA taxonomy files
        # For now, return common groups
        if marker in [MarkerType.SSU_16S, MarkerType.LSU_23S]:
            return [
                "Bacteria",
                "Archaea",
                "Proteobacteria",
                "Firmicutes",
                "Bacteroidetes",
                "Actinobacteria",
                "Cyanobacteria",
            ]
        elif marker in [MarkerType.SSU_18S, MarkerType.LSU_28S]:
            return [
                "Eukaryota",
                "Metazoa",
                "Fungi",
                "Plantae",
                "Protista",
                "Chordata",
                "Arthropoda",
            ]
        else:
            return []
