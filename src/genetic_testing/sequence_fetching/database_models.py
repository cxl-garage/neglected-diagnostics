"""
Database Models for Sequence Fetching

Pydantic models for representing sequences, search parameters, and database configurations.
"""

from datetime import datetime
from enum import Enum
from typing import Dict, List, Optional, Union

from pydantic import BaseModel, Field, validator


class DatabaseType(str, Enum):
    """Supported database types."""

    NCBI_NUCLEOTIDE = "nucleotide"
    NCBI_GENE = "gene"
    NCBI_GENOME = "genome"
    BOLD = "bold"
    SILVA = "silva"
    UNITE = "unite"


class MarkerType(str, Enum):
    """Common molecular markers for taxonomic identification."""

    COI = "COI"  # Cytochrome c oxidase subunit I
    COX1 = "COX1"  # Alternative name for COI
    RBCL = "rbcL"  # RuBisCO large subunit
    MATK = "matK"  # Maturase K
    ITS1 = "ITS1"  # Internal transcribed spacer 1
    ITS2 = "ITS2"  # Internal transcribed spacer 2
    ITS = "ITS"  # Full ITS region
    SSU_16S = "16S"  # 16S ribosomal RNA
    SSU_18S = "18S"  # 18S ribosomal RNA
    LSU_23S = "23S"  # 23S ribosomal RNA
    LSU_28S = "28S"  # 28S ribosomal RNA
    CUSTOM = "custom"


class MoleculeType(str, Enum):
    """Molecule types for sequence searches."""

    DNA = "DNA"
    RNA = "RNA"
    PROTEIN = "protein"
    GENOMIC = "genomic"
    MRNA = "mRNA"


class SequenceRecord(BaseModel):
    """Represents a sequence record from any database."""

    id: str = Field(..., description="Unique identifier for the sequence")
    accession: Optional[str] = Field(None, description="Accession number")
    title: str = Field(..., description="Sequence title/description")
    organism: Optional[str] = Field(None, description="Source organism")
    taxonomy_id: Optional[int] = Field(None, description="NCBI taxonomy ID")
    sequence: Optional[str] = Field(None, description="Nucleotide/protein sequence")
    length: int = Field(..., description="Sequence length")
    molecule_type: Optional[MoleculeType] = Field(None, description="Type of molecule")
    marker: Optional[MarkerType] = Field(None, description="Molecular marker type")
    create_date: Optional[datetime] = Field(None, description="Creation date")
    update_date: Optional[datetime] = Field(None, description="Last update date")
    database: DatabaseType = Field(..., description="Source database")
    quality_score: Optional[float] = Field(None, description="Quality assessment score")
    country: Optional[str] = Field(None, description="Country of origin")
    lat_lon: Optional[str] = Field(None, description="Latitude/longitude coordinates")
    specimen_voucher: Optional[str] = Field(
        None, description="Specimen voucher information"
    )

    # Additional metadata
    metadata: Dict[str, Union[str, int, float]] = Field(default_factory=dict)

    @validator("sequence")
    def validate_sequence(cls, v):
        """Validate sequence contains only valid nucleotide/amino acid characters."""
        if v is not None:
            # Remove whitespace and convert to uppercase
            v = "".join(v.split()).upper()
            # Basic validation - could be enhanced
            valid_chars = set("ACGTUNRYSWKMBDHV-")  # IUPAC nucleotide codes + gap
            if not all(c in valid_chars for c in v):
                # Check if it might be a protein sequence
                protein_chars = set("ACDEFGHIKLMNPQRSTVWY*-")
                if not all(c in protein_chars for c in v):
                    raise ValueError("Sequence contains invalid characters")
        return v


class SearchParameters(BaseModel):
    """Parameters for database searches."""

    # Basic search parameters
    search_term: str = Field(..., description="Search query term")
    database: DatabaseType = Field(..., description="Target database")
    max_results: int = Field(default=10000, description="Maximum number of results")

    # Taxonomic filters
    taxon: Optional[str] = Field(
        None, description="Taxonomic filter (e.g., 'Salmo salar')"
    )
    taxonomy_id: Optional[int] = Field(None, description="NCBI taxonomy ID")
    exclude_taxa: List[str] = Field(default_factory=list, description="Taxa to exclude")

    # Molecular marker filters
    marker: Optional[MarkerType] = Field(None, description="Target molecular marker")
    gene_name: Optional[str] = Field(None, description="Specific gene name")

    # Sequence property filters
    min_length: Optional[int] = Field(None, description="Minimum sequence length")
    max_length: Optional[int] = Field(None, description="Maximum sequence length")
    molecule_type: Optional[MoleculeType] = Field(
        None, description="Molecule type filter"
    )

    # Quality filters
    exclude_partial: bool = Field(
        default=False, description="Exclude partial sequences"
    )
    exclude_predicted: bool = Field(
        default=False, description="Exclude predicted sequences"
    )
    require_voucher: bool = Field(default=False, description="Require specimen voucher")

    # Geographic filters
    country: Optional[str] = Field(None, description="Country filter")

    # Date filters
    date_from: Optional[datetime] = Field(None, description="Earliest publication date")
    date_to: Optional[datetime] = Field(None, description="Latest publication date")


class DatabaseConfig(BaseModel):
    """Configuration for database connections."""

    # NCBI configuration
    ncbi_email: Optional[str] = Field(None, description="Email for NCBI Entrez")
    ncbi_api_key: Optional[str] = Field(
        None, description="NCBI API key for higher rate limits"
    )

    # BOLD configuration
    bold_base_url: str = Field(
        default="http://www.boldsystems.org/index.php/API_Public/",
        description="BOLD API base URL",
    )
    bold_rate_limit: int = Field(default=1, description="Requests per second for BOLD")

    # SILVA configuration
    silva_base_url: str = Field(
        default="https://www.arb-silva.de/", description="SILVA database URL"
    )
    silva_version: str = Field(default="138", description="SILVA database version")

    # UNITE configuration
    unite_base_url: str = Field(
        default="https://unite.ut.ee/", description="UNITE database URL"
    )
    unite_version: str = Field(default="8.3", description="UNITE database version")

    # Caching configuration
    cache_enabled: bool = Field(default=True, description="Enable request caching")
    cache_expire_hours: int = Field(
        default=24, description="Cache expiration time in hours"
    )

    # Rate limiting
    default_rate_limit: int = Field(
        default=3, description="Default requests per second"
    )
    max_concurrent_requests: int = Field(
        default=10, description="Max concurrent async requests"
    )


class SearchResult(BaseModel):
    """Results from a database search."""

    query: SearchParameters = Field(..., description="Original search parameters")
    total_found: int = Field(..., description="Total number of sequences found")
    sequences: List[SequenceRecord] = Field(
        ..., description="Retrieved sequence records"
    )
    search_time: float = Field(..., description="Search execution time in seconds")
    database_version: Optional[str] = Field(None, description="Database version used")
    warnings: List[str] = Field(default_factory=list, description="Search warnings")

    @property
    def retrieved_count(self) -> int:
        """Number of sequences actually retrieved."""
        return len(self.sequences)

    @property
    def truncated(self) -> bool:
        """Whether results were truncated due to max_results limit."""
        return self.total_found > self.retrieved_count


class SequenceAlignment(BaseModel):
    """Represents a multiple sequence alignment."""

    sequences: List[SequenceRecord] = Field(..., description="Aligned sequences")
    alignment_method: str = Field(..., description="Alignment algorithm used")
    alignment_score: Optional[float] = Field(
        None, description="Overall alignment quality score"
    )
    consensus_sequence: Optional[str] = Field(None, description="Consensus sequence")
    conserved_regions: List[Dict[str, Union[int, float]]] = Field(
        default_factory=list, description="Identified conserved regions"
    )
