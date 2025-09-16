"""
Enhanced NCBI Client

This module provides an enhanced NCBI client with advanced search capabilities,
taxonomy integration, and marker-specific queries as specified in Phase 1.1 of the roadmap.
"""

import asyncio
import os
import time
from datetime import datetime
from io import StringIO
from typing import Dict, List, Optional, Set, Tuple, Union
from xml.etree import ElementTree as ET

import pandas as pd
import streamlit as st
from Bio import Entrez, SeqIO
from Bio.Entrez.Parser import DictionaryElement, ListElement, StringElement

from utils.log import _init_logger

from .database_models import (
    DatabaseConfig,
    DatabaseType,
    MarkerType,
    MoleculeType,
    SearchParameters,
    SearchResult,
    SequenceRecord,
)
from .sequence_filters import SequenceFilter

logger = _init_logger(__name__)


class NCBIClient:
    """Enhanced NCBI client with advanced search and filtering capabilities."""

    def __init__(self, config: Optional[DatabaseConfig] = None):
        """
        Initialize NCBI client with configuration.

        Parameters
        ----------
        config : DatabaseConfig, optional
            Database configuration. If not provided, uses defaults and environment variables.
        """
        self.config = config or DatabaseConfig()
        self._setup_entrez()

        # Marker-specific search terms
        self.marker_terms = {
            MarkerType.COI: [
                "COI",
                "COX1",
                "cytochrome c oxidase subunit 1",
                "cytochrome oxidase 1",
            ],
            MarkerType.COX1: ["COX1", "COI", "cytochrome c oxidase subunit 1"],
            MarkerType.RBCL: [
                "rbcL",
                "ribulose-1,5-bisphosphate carboxylase/oxygenase large subunit",
            ],
            MarkerType.MATK: ["matK", "maturase K"],
            MarkerType.ITS1: ["ITS1", "internal transcribed spacer 1"],
            MarkerType.ITS2: ["ITS2", "internal transcribed spacer 2"],
            MarkerType.ITS: ["ITS", "internal transcribed spacer"],
            MarkerType.SSU_16S: ["16S", "16S ribosomal RNA", "16S rRNA"],
            MarkerType.SSU_18S: ["18S", "18S ribosomal RNA", "18S rRNA"],
            MarkerType.LSU_23S: ["23S", "23S ribosomal RNA", "23S rRNA"],
            MarkerType.LSU_28S: ["28S", "28S ribosomal RNA", "28S rRNA"],
        }

        # Common sequence length ranges for quality filtering
        self.marker_length_ranges = {
            MarkerType.COI: (600, 700),
            MarkerType.COX1: (600, 700),
            MarkerType.RBCL: (1300, 1500),
            MarkerType.MATK: (800, 900),
            MarkerType.ITS1: (200, 400),
            MarkerType.ITS2: (200, 400),
            MarkerType.ITS: (400, 800),
            MarkerType.SSU_16S: (1400, 1600),
            MarkerType.SSU_18S: (1700, 1900),
            MarkerType.LSU_23S: (2800, 3200),
            MarkerType.LSU_28S: (3000, 4000),
        }

    def _setup_entrez(self):
        """Setup Entrez configuration."""
        # Try to get email and API key from Streamlit secrets first, then environment
        email = ""
        api_key = ""

        try:
            email = st.secrets.get("ENTREZ_EMAIL", "")
            api_key = st.secrets.get("NCBI_API_KEY", "")
            logger.debug("Successfully accessed Streamlit secrets")
        except (FileNotFoundError, AttributeError, KeyError) as e:
            logger.debug(
                f"Streamlit secrets not available ({e}), trying environment variables"
            )
            email = os.environ.get("ENTREZ_EMAIL", "")
            api_key = os.environ.get("NCBI_API_KEY", "")

        # Fallback to config object
        if not email and self.config.ncbi_email:
            email = self.config.ncbi_email
        if not api_key and self.config.ncbi_api_key:
            api_key = self.config.ncbi_api_key

        # Set email (required)
        if email:
            Entrez.email = email
            logger.info(f"NCBI Entrez email configured: {email}")
        else:
            logger.warning(
                "No email configured for NCBI Entrez. This may limit API access."
            )

        # Set API key (optional but recommended for better rate limits)
        if api_key:
            Entrez.api_key = api_key
            logger.info(
                "NCBI API key configured successfully - increased rate limits enabled"
            )
        else:
            logger.info(
                "No NCBI API key configured. Using standard rate limits (3 requests/second)"
            )

    def search_sequences(self, params: SearchParameters) -> SearchResult:
        """
        Search for sequences using enhanced parameters.

        Parameters
        ----------
        params : SearchParameters
            Search parameters including taxonomic and marker filters

        Returns
        -------
        SearchResult
            Search results with metadata
        """
        start_time = time.time()

        # Build search query
        query = self._build_search_query(params)
        logger.info(f"Built NCBI search query: {query}")

        try:
            # Perform search
            uids = self._search_database(
                params.database.value, query, params.max_results
            )
            logger.info(f"Found {len(uids)} sequences for query: {query}")

            if not uids:
                return SearchResult(
                    query=params,
                    total_found=0,
                    sequences=[],
                    search_time=time.time() - start_time,
                    warnings=["No sequences found for the given search criteria"],
                )

            # Get sequence summaries
            summaries = self._get_summaries(params.database.value, uids)

            # Convert to SequenceRecord objects
            sequences = self._parse_summaries_to_records(summaries, params.database)

            # Apply sequence filters if specified
            if any(
                [
                    params.min_length,
                    params.max_length,
                    params.exclude_partial,
                    params.exclude_predicted,
                ]
            ):
                seq_filter = SequenceFilter(
                    min_length=params.min_length,
                    max_length=params.max_length,
                    exclude_partial=params.exclude_partial,
                    exclude_predicted=params.exclude_predicted,
                )
                sequences = seq_filter.filter_sequences(sequences)

            search_time = time.time() - start_time

            return SearchResult(
                query=params,
                total_found=len(uids),
                sequences=sequences,
                search_time=search_time,
            )

        except Exception as e:
            logger.error(f"Error searching NCBI: {e}", exc_info=True)
            return SearchResult(
                query=params,
                total_found=0,
                sequences=[],
                search_time=time.time() - start_time,
                warnings=[f"Search failed: {str(e)}"],
            )

    def fetch_target_sequences(
        self, taxon: str, marker: MarkerType, max_results: int = 1000
    ) -> SearchResult:
        """
        Fetch target sequences for a specific taxon and marker.

        Parameters
        ----------
        taxon : str
            Target taxon (e.g., "Salmo salar" or "Salmonidae")
        marker : MarkerType
            Molecular marker to search for
        max_results : int
            Maximum number of results to return

        Returns
        -------
        SearchResult
            Target sequences
        """
        params = SearchParameters(
            search_term="",  # Will be built from taxon and marker
            database=DatabaseType.NCBI_NUCLEOTIDE,
            max_results=max_results,
            taxon=taxon,
            marker=marker,
            exclude_partial=True,
            exclude_predicted=True,
        )

        return self.search_sequences(params)

    def fetch_off_target_sequences(
        self,
        genus: str,
        exclude_species: List[str],
        marker: MarkerType,
        max_results: int = 1000,
    ) -> SearchResult:
        """
        Fetch off-target sequences from related taxa.

        Parameters
        ----------
        genus : str
            Genus to search within
        exclude_species : List[str]
            Species to exclude from results
        marker : MarkerType
            Molecular marker to search for
        max_results : int
            Maximum number of results to return

        Returns
        -------
        SearchResult
            Off-target sequences
        """
        params = SearchParameters(
            search_term="",  # Will be built from genus and marker
            database=DatabaseType.NCBI_NUCLEOTIDE,
            max_results=max_results,
            taxon=genus,
            marker=marker,
            exclude_taxa=exclude_species,
            exclude_partial=True,
            exclude_predicted=True,
        )

        return self.search_sequences(params)

    def download_mitogenomes(self, taxonomy_ids: List[int]) -> SearchResult:
        """
        Download complete mitochondrial genomes by taxonomy ID.

        Parameters
        ----------
        taxonomy_ids : List[int]
            NCBI taxonomy IDs

        Returns
        -------
        SearchResult
            Mitochondrial genome sequences
        """
        # Build query for mitochondrial genomes
        taxid_query = " OR ".join([f"txid{tid}[Organism]" for tid in taxonomy_ids])
        search_term = f"({taxid_query}) AND mitochondrion[filter] AND complete[title]"

        params = SearchParameters(
            search_term=search_term,
            database=DatabaseType.NCBI_NUCLEOTIDE,
            max_results=len(taxonomy_ids) * 10,  # Allow multiple genomes per species
            molecule_type=MoleculeType.GENOMIC,
            min_length=10000,  # Mitogenomes are typically 15-20kb
        )

        return self.search_sequences(params)

    def fetch_sequences_with_data(
        self, params: SearchParameters, include_sequences: bool = True
    ) -> SearchResult:
        """
        Fetch sequences with actual sequence data.

        Parameters
        ----------
        params : SearchParameters
            Search parameters
        include_sequences : bool
            Whether to download actual sequence data

        Returns
        -------
        SearchResult
            Sequences with sequence data included
        """
        # First get the basic search results
        result = self.search_sequences(params)

        if include_sequences and result.sequences:
            # Extract IDs for sequence download
            ids = [seq.id for seq in result.sequences]

            # Download sequences in batches
            batch_size = 100
            for i in range(0, len(ids), batch_size):
                batch_ids = ids[i : i + batch_size]
                sequences_data = self._fetch_sequence_data(
                    params.database.value, batch_ids
                )

                # Update sequence records with actual sequence data
                for seq_record in result.sequences[i : i + len(batch_ids)]:
                    if seq_record.id in sequences_data:
                        seq_record.sequence = sequences_data[seq_record.id]

        return result

    def _build_search_query(self, params: SearchParameters) -> str:
        """Build NCBI search query from parameters."""
        query_parts = []

        # Base search term
        if params.search_term:
            query_parts.append(params.search_term)

        # Taxonomic filters
        if params.taxon:
            query_parts.append(f"{params.taxon}[Organism]")
        elif params.taxonomy_id:
            query_parts.append(f"txid{params.taxonomy_id}[Organism]")

        # Exclude taxa
        for exclude_taxon in params.exclude_taxa:
            query_parts.append(f"NOT {exclude_taxon}[Organism]")

        # Marker-specific terms
        if params.marker and params.marker in self.marker_terms:
            marker_terms = self.marker_terms[params.marker]
            marker_query = " OR ".join([f"{term}[All Fields]" for term in marker_terms])
            query_parts.append(f"({marker_query})")

        # Gene name filter
        if params.gene_name:
            query_parts.append(f"{params.gene_name}[Gene Name]")

        # Molecule type filter
        if params.molecule_type:
            if params.molecule_type == MoleculeType.GENOMIC:
                query_parts.append("biomol_genomic[Properties]")
            elif params.molecule_type == MoleculeType.MRNA:
                query_parts.append("biomol_mrna[Properties]")

        # Length filters
        if params.min_length or params.max_length:
            min_len = params.min_length or 1
            max_len = params.max_length or 100000
            query_parts.append(f"{min_len}:{max_len}[Sequence Length]")

        # Quality filters
        if params.exclude_partial:
            query_parts.append("NOT partial[Title]")

        if params.exclude_predicted:
            query_parts.append("NOT predicted[Title]")
            query_parts.append("NOT hypothetical[Title]")

        # Geographic filter
        if params.country:
            query_parts.append(f"{params.country}[Country]")

        # Date filters
        if params.date_from or params.date_to:
            date_from = (
                params.date_from.strftime("%Y/%m/%d")
                if params.date_from
                else "1900/01/01"
            )
            date_to = (
                params.date_to.strftime("%Y/%m/%d") if params.date_to else "2030/12/31"
            )
            query_parts.append(f"{date_from}:{date_to}[Publication Date]")

        # Join all parts
        query = " AND ".join(query_parts)

        # If no specific query was built, use a default
        if not query:
            query = "all[filter]"

        return query

    def _search_database(
        self, database: str, query: str, max_results: int
    ) -> List[str]:
        """Search NCBI database and return UIDs."""
        try:
            handle = Entrez.esearch(
                db=database, term=query, retmax=max_results, usehistory="y"
            )

            if handle:
                record = Entrez.read(handle)
                uids = record.get("IdList", [])
                handle.close()
                return uids

        except Exception as e:
            logger.error(f"Error searching NCBI database {database}: {e}")
            raise

        return []

    def _get_summaries(self, database: str, uids: List[str]) -> List[DictionaryElement]:
        """Get document summaries for UIDs."""
        if not uids:
            return []

        try:
            # Process in batches to avoid timeout
            batch_size = 200
            all_summaries = []

            for i in range(0, len(uids), batch_size):
                batch_uids = uids[i : i + batch_size]

                handle = Entrez.esummary(db=database, id=batch_uids)
                if handle:
                    summaries = Entrez.read(handle)
                    if isinstance(summaries, list):
                        all_summaries.extend(summaries)
                    else:
                        all_summaries.append(summaries)
                    handle.close()

                # Add small delay between batches to be respectful
                if i + batch_size < len(uids):
                    time.sleep(0.1)

            return all_summaries

        except Exception as e:
            logger.error(f"Error getting summaries from NCBI: {e}")
            raise

    def _parse_summaries_to_records(
        self, summaries: List[DictionaryElement], database: DatabaseType
    ) -> List[SequenceRecord]:
        """Convert NCBI summaries to SequenceRecord objects."""
        records = []

        for summary in summaries:
            try:
                # Extract basic information
                record = SequenceRecord(
                    id=str(summary.get("Id", "")),
                    accession=summary.get("AccessionVersion", ""),
                    title=summary.get("Title", ""),
                    organism=self._extract_organism(summary.get("Title", "")),
                    taxonomy_id=summary.get("TaxId"),
                    length=int(summary.get("Length", 0)),
                    database=database,
                    create_date=self._parse_date(summary.get("CreateDate")),
                    update_date=self._parse_date(summary.get("UpdateDate")),
                )

                # Try to infer marker type from title
                record.marker = self._infer_marker_type(record.title)

                # Add to metadata
                record.metadata = {
                    "gi": summary.get("Gi", ""),
                    "status": summary.get("Status", ""),
                    "extra_info": summary.get("Extra", ""),
                }

                records.append(record)

            except Exception as e:
                logger.warning(
                    f"Error parsing summary {summary.get('Id', 'unknown')}: {e}"
                )
                continue

        return records

    def _extract_organism(self, title: str) -> Optional[str]:
        """Extract organism name from sequence title."""
        import re

        # Look for binomial nomenclature pattern
        pattern = r"\b([A-Z][a-z]+ [a-z]+)\b"
        match = re.search(pattern, title)

        if match:
            return match.group(1)

        return None

    def _infer_marker_type(self, title: str) -> Optional[MarkerType]:
        """Infer marker type from sequence title."""
        title_lower = title.lower()

        for marker, terms in self.marker_terms.items():
            if any(term.lower() in title_lower for term in terms):
                return marker

        return None

    def _parse_date(self, date_str: Optional[str]) -> Optional[datetime]:
        """Parse date string to datetime object."""
        if not date_str:
            return None

        try:
            # NCBI typically uses YYYY/MM/DD format
            return datetime.strptime(date_str, "%Y/%m/%d")
        except ValueError:
            try:
                # Try alternative format
                return datetime.strptime(date_str, "%Y-%m-%d")
            except ValueError:
                logger.warning(f"Could not parse date: {date_str}")
                return None

    def _fetch_sequence_data(self, database: str, ids: List[str]) -> Dict[str, str]:
        """Fetch actual sequence data for given IDs."""
        if not ids:
            return {}

        try:
            handle = Entrez.efetch(db=database, id=ids, rettype="fasta", retmode="text")

            if handle:
                # Parse FASTA sequences
                sequences = {}
                sequences_iterator = SeqIO.parse(handle, "fasta")

                for seq_record in sequences_iterator:
                    # Extract ID from description (NCBI format: gi|ID|...)
                    seq_id = (
                        seq_record.id.split("|")[-1]
                        if "|" in seq_record.id
                        else seq_record.id
                    )
                    sequences[seq_id] = str(seq_record.seq)

                handle.close()
                return sequences

        except Exception as e:
            logger.error(f"Error fetching sequence data: {e}")

        return {}

    def get_taxonomy_info(self, taxonomy_id: int) -> Optional[Dict]:
        """
        Get detailed taxonomy information for a taxonomy ID.

        Parameters
        ----------
        taxonomy_id : int
            NCBI taxonomy ID

        Returns
        -------
        Dict, optional
            Taxonomy information
        """
        try:
            handle = Entrez.efetch(db="taxonomy", id=str(taxonomy_id), retmode="xml")
            if handle:
                tree = ET.parse(handle)
                root = tree.getroot()
                handle.close()

                # Extract taxonomy information
                taxon = root.find(".//Taxon")
                if taxon is not None:
                    return {
                        "taxid": taxonomy_id,
                        "scientific_name": taxon.findtext("ScientificName"),
                        "common_name": taxon.findtext("CommonName"),
                        "rank": taxon.findtext("Rank"),
                        "lineage": self._extract_lineage(taxon),
                    }

        except Exception as e:
            logger.error(f"Error fetching taxonomy info for ID {taxonomy_id}: {e}")

        return None

    def _extract_lineage(self, taxon_element) -> List[Dict[str, str]]:
        """Extract taxonomic lineage from XML element."""
        lineage = []
        lineage_element = taxon_element.find("LineageEx")

        if lineage_element is not None:
            for taxon in lineage_element.findall("Taxon"):
                lineage.append(
                    {
                        "taxid": taxon.findtext("TaxId"),
                        "scientific_name": taxon.findtext("ScientificName"),
                        "rank": taxon.findtext("Rank"),
                    }
                )

        return lineage
