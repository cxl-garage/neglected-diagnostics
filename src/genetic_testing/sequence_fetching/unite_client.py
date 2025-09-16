"""
UNITE Database Client

This module provides integration with the UNITE database for fungal ITS sequences
as specified in Phase 1.2 of the roadmap.
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


class UNITEClient:
    """Client for accessing UNITE fungal ITS database sequences."""

    def __init__(self, config: Optional[DatabaseConfig] = None):
        """
        Initialize UNITE client.

        Parameters
        ----------
        config : DatabaseConfig, optional
            Database configuration
        """
        self.config = config or DatabaseConfig()
        self.base_url = self.config.unite_base_url
        self.version = self.config.unite_version

        # Setup caching if enabled
        if self.config.cache_enabled:
            self.session = requests_cache.CachedSession(
                "unite_cache", expire_after=self.config.cache_expire_hours * 3600
            )
        else:
            self.session = requests.Session()

        # UNITE download URLs (may need updating based on current UNITE structure)
        self.download_base = f"{self.base_url}repository.php"

        # Common UNITE datasets
        self.datasets = {
            "general": f"UNITE_general_release_{self.version}",
            "dynamic": f"UNITE_dynamic_{self.version}",
            "all_eukaryotes": f"UNITE_all_eukaryotes_{self.version}",
        }

    def search_sequences(self, params: SearchParameters) -> SearchResult:
        """
        Search UNITE database for ITS sequences.

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
            # Check if marker is supported
            if params.marker and params.marker not in [
                MarkerType.ITS,
                MarkerType.ITS1,
                MarkerType.ITS2,
            ]:
                return SearchResult(
                    query=params,
                    total_found=0,
                    sequences=[],
                    search_time=time.time() - start_time,
                    warnings=["Only ITS markers are supported by UNITE database"],
                )

            # Fetch UNITE sequences
            sequences = self._fetch_unite_sequences(params)

            # Filter sequences
            filtered_sequences = self._filter_sequences(sequences, params)

            search_time = time.time() - start_time

            return SearchResult(
                query=params,
                total_found=len(sequences),
                sequences=filtered_sequences[: params.max_results],
                search_time=search_time,
                database_version=f"UNITE {self.version}",
            )

        except Exception as e:
            logger.error(f"Error searching UNITE database: {e}", exc_info=True)
            return SearchResult(
                query=params,
                total_found=0,
                sequences=[],
                search_time=time.time() - start_time,
                warnings=[f"UNITE search failed: {str(e)}"],
            )

    def fetch_its_sequences(
        self, taxon: str = "Fungi", its_region: str = "ITS", max_results: int = 1000
    ) -> SearchResult:
        """
        Fetch ITS sequences for fungi.

        Parameters
        ----------
        taxon : str
            Target taxonomic group
        its_region : str
            ITS region (ITS, ITS1, ITS2)
        max_results : int
            Maximum number of sequences

        Returns
        -------
        SearchResult
            ITS sequences
        """
        marker = MarkerType.ITS
        if its_region.upper() == "ITS1":
            marker = MarkerType.ITS1
        elif its_region.upper() == "ITS2":
            marker = MarkerType.ITS2

        params = SearchParameters(
            search_term=taxon,
            database=DatabaseType.UNITE,
            max_results=max_results,
            marker=marker,
            taxon=taxon,
        )

        return self.search_sequences(params)

    def _fetch_unite_sequences(self, params: SearchParameters) -> List[SequenceRecord]:
        """Fetch sequences from UNITE database."""
        sequences = []

        try:
            # For this implementation, we'll simulate downloading from UNITE
            # In practice, you would download the actual UNITE FASTA files

            # Determine which dataset to use
            dataset = self._select_dataset(params)

            # Download UNITE data (this is a placeholder)
            fasta_content = self._download_unite_data(dataset)

            if fasta_content:
                sequences = self._parse_unite_fasta(fasta_content)

        except Exception as e:
            logger.error(f"Error fetching UNITE sequences: {e}")

        return sequences

    def _select_dataset(self, params: SearchParameters) -> str:
        """Select appropriate UNITE dataset based on parameters."""
        # Use general release by default
        return self.datasets["general"]

    def _download_unite_data(self, dataset: str) -> Optional[str]:
        """Download UNITE FASTA data."""
        try:
            # This is a placeholder - actual UNITE download would be more complex
            # You might need to construct specific URLs or use their API

            # For now, return None to indicate no data downloaded
            # In a real implementation, you would:
            # 1. Construct the proper download URL
            # 2. Handle authentication if required
            # 3. Download and potentially decompress the data

            logger.info(f"Would download UNITE dataset: {dataset}")
            return None

        except Exception as e:
            logger.error(f"Error downloading UNITE data: {e}")
            return None

    def _parse_unite_fasta(self, fasta_content: str) -> List[SequenceRecord]:
        """Parse UNITE FASTA content into SequenceRecord objects."""
        sequences = []

        try:
            fasta_io = StringIO(fasta_content)
            seq_records = SeqIO.parse(fasta_io, "fasta")

            for seq_record in seq_records:
                # Parse UNITE header format
                # Example: >UDB000001|SH000001.07FU|refs|k__Fungi;p__Ascomycota;...

                header_parts = seq_record.description.split("|")
                unite_id = header_parts[0] if header_parts else seq_record.id
                sh_id = header_parts[1] if len(header_parts) > 1 else ""

                # Extract taxonomy if present
                taxonomy = ""
                organism = None

                if len(header_parts) >= 4:
                    taxonomy = header_parts[3]
                    organism = self._extract_organism_from_taxonomy(taxonomy)

                # Determine ITS region from sequence length and content
                its_marker = self._determine_its_region(
                    seq_record.seq, seq_record.description
                )

                sequence_record = SequenceRecord(
                    id=unite_id,
                    accession=sh_id,
                    title=seq_record.description,
                    organism=organism,
                    sequence=str(seq_record.seq),
                    length=len(seq_record.seq),
                    database=DatabaseType.UNITE,
                    marker=its_marker,
                    metadata={
                        "unite_id": unite_id,
                        "sh_id": sh_id,
                        "taxonomy": taxonomy,
                        "unite_version": self.version,
                    },
                )

                sequences.append(sequence_record)

        except Exception as e:
            logger.error(f"Error parsing UNITE FASTA: {e}")

        return sequences

    def _extract_organism_from_taxonomy(self, taxonomy: str) -> Optional[str]:
        """Extract organism name from UNITE taxonomy string."""
        if not taxonomy:
            return None

        # UNITE taxonomy format: k__Kingdom;p__Phylum;...;g__Genus;s__Species
        parts = taxonomy.split(";")

        genus = None
        species = None

        for part in reversed(parts):
            if part.startswith("s__") and len(part) > 3:
                species = part[3:].strip().replace("_", " ")
            elif part.startswith("g__") and len(part) > 3:
                genus = part[3:].strip()

        if genus and species and not species.startswith("unidentified"):
            # Handle cases where species includes genus
            if " " in species:
                return species
            else:
                return f"{genus} {species}"
        elif genus:
            return genus

        return None

    def _determine_its_region(self, sequence, description: str) -> MarkerType:
        """Determine ITS region from sequence and description."""
        desc_lower = description.lower()

        if "its1" in desc_lower:
            return MarkerType.ITS1
        elif "its2" in desc_lower:
            return MarkerType.ITS2
        else:
            # Try to infer from sequence length
            seq_len = len(sequence)
            if seq_len < 300:
                # Shorter sequences are more likely ITS1 or ITS2
                return MarkerType.ITS1 if seq_len < 250 else MarkerType.ITS2
            else:
                # Longer sequences likely include full ITS region
                return MarkerType.ITS

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

            # Exclude taxa filter
            if params.exclude_taxa:
                if any(
                    self._matches_taxon(seq, exclude) for exclude in params.exclude_taxa
                ):
                    continue

            filtered.append(seq)

        return filtered

    def _matches_taxon(self, sequence: SequenceRecord, taxon: str) -> bool:
        """Check if sequence matches taxonomic filter."""
        taxon_lower = taxon.lower()

        # Check organism name
        if sequence.organism and taxon_lower in sequence.organism.lower():
            return True

        # Check taxonomy metadata
        if sequence.metadata.get("taxonomy"):
            taxonomy = sequence.metadata["taxonomy"].lower()
            return taxon_lower in taxonomy

        # Check title
        if taxon_lower in sequence.title.lower():
            return True

        return False

    def get_available_taxa(self) -> List[str]:
        """
        Get list of major fungal taxonomic groups available in UNITE.

        Returns
        -------
        List[str]
            Major fungal taxonomic groups
        """
        return [
            "Fungi",
            "Ascomycota",
            "Basidiomycota",
            "Chytridiomycota",
            "Glomeromycota",
            "Mucoromycota",
            "Zoopagomycota",
            "Saccharomycetes",
            "Eurotiomycetes",
            "Sordariomycetes",
            "Leotiomycetes",
            "Dothideomycetes",
            "Agaricomycetes",
            "Tremellomycetes",
            "Ustilaginomycetes",
            "Pucciniomycetes",
        ]

    def get_its_statistics(self) -> Dict[str, int]:
        """
        Get statistics about ITS sequences in UNITE.

        Returns
        -------
        Dict[str, int]
            Statistics about sequence counts
        """
        # This would require analyzing the actual UNITE database
        # For now, return placeholder statistics
        return {
            "total_sequences": 0,
            "its1_sequences": 0,
            "its2_sequences": 0,
            "full_its_sequences": 0,
            "species_represented": 0,
            "genera_represented": 0,
        }
