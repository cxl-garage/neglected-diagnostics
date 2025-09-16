import asyncio
from datetime import datetime
from typing import Dict, List, Optional

import pandas as pd
import streamlit as st
from common.render_method import render_markdown

from utils.log import _init_logger

logger = _init_logger(__name__)

from app.common import data_processing, setup
from app.common.constants import (
    MAIN_PAGE_COLS_GAP,
    MAIN_PAGE_COLS_SIZE,
    MAX_SEQUENCE_LENGTH_DOWNLOAD,
    NAVIGATE_WARNING_MD,
    NCBI_DF,
    NCBI_DF_FILTER,
    NCBI_SUMMARY_FORM,
    TOP_N_ORGANISMS,
)
from app.frontend.aggrid_table import aggrid_table
from app.frontend.download_data import download_data
from app.frontend.filter_dataframe import filter_dataframe
from app.frontend.filters_info import show_filters_info
from genetic_testing.routers import ncbi

# Import new Phase 1 modules
from genetic_testing.sequence_fetching import (
    BatchSequenceFetcher,
    BOLDClient,
    DatabaseConfig,
    DatabaseType,
    MarkerType,
    MoleculeType,
    NCBIClient,
    SearchParameters,
    SILVAClient,
    UNITEClient,
)

# Initialize the Streamlit session state keys
setup.initialize()

# Initialize session state for new features
if "search_results" not in st.session_state:
    st.session_state.search_results = None
if "selected_databases" not in st.session_state:
    st.session_state.selected_databases = ["nucleotide"]
if "search_mode" not in st.session_state:
    st.session_state.search_mode = "basic"

# Define constants for download options
CSV_DOWNLOAD = "Download table as CSV"
FASTA_DOWNLOAD = "Download sequences as a fasta file"


# Helper functions
@st.cache_data(show_spinner=False)
def convert_search_result_to_df(search_results: Dict) -> pd.DataFrame:
    """Convert search results to pandas DataFrame for display."""
    all_sequences = []

    logger.info(f"Converting search results from {len(search_results)} databases")

    for db_name, result in search_results.items():
        logger.info(f"Processing database: {db_name}, result type: {type(result)}")

        if result and hasattr(result, "sequences") and result.sequences:
            logger.info(f"Found {len(result.sequences)} sequences in {db_name}")

            for seq in result.sequences:
                try:
                    seq_dict = {
                        "Id": str(seq.id) if seq.id else "",
                        "Accession": str(seq.accession) if seq.accession else "",
                        "Title": str(seq.title) if seq.title else "No title",
                        "Organism": str(seq.organism) if seq.organism else "",
                        "Length": int(seq.length) if seq.length else 0,
                        "Database": str(db_name),
                        "Marker": seq.marker.value if seq.marker else "",
                        "Quality Score": f"{seq.quality_score:.3f}"
                        if seq.quality_score
                        else "",
                        "Country": str(seq.country) if seq.country else "",
                        "Create Date": seq.create_date.strftime("%Y-%m-%d")
                        if seq.create_date
                        else "",
                    }
                    all_sequences.append(seq_dict)

                except Exception as e:
                    logger.warning(f"Error processing sequence {seq.id}: {e}")
                    continue
        else:
            if result:
                logger.warning(f"No sequences found in {db_name} result: {result}")
            else:
                logger.warning(f"Empty result for {db_name}")

    logger.info(f"Converted {len(all_sequences)} sequences total")

    if all_sequences:
        df = pd.DataFrame(all_sequences)
        logger.info(f"Created DataFrame with columns: {df.columns.tolist()}")
        return df
    else:
        # Return empty DataFrame with expected columns
        logger.warning("No sequences to convert, returning empty DataFrame")
        return pd.DataFrame(
            columns=[
                "Id",
                "Accession",
                "Title",
                "Organism",
                "Length",
                "Database",
                "Marker",
                "Quality Score",
                "Country",
                "Create Date",
            ]
        )


def perform_enhanced_search(params: SearchParameters) -> Dict:
    """Perform search using enhanced Phase 1 clients."""
    results = {}

    try:
        logger.info(f"Starting search for database: {params.database}")

        if params.database == DatabaseType.NCBI_NUCLEOTIDE:
            logger.info("Using enhanced NCBI client")
            client = NCBIClient()
            search_result = client.search_sequences(params)
            results["NCBI Nucleotide"] = search_result
            logger.info(
                f"NCBI search completed: {len(search_result.sequences) if search_result.sequences else 0} sequences"
            )

        elif params.database == DatabaseType.BOLD:
            logger.info("Using BOLD client")
            client = BOLDClient()
            search_result = client.search_sequences(params)
            results["BOLD"] = search_result
            logger.info(
                f"BOLD search completed: {len(search_result.sequences) if search_result.sequences else 0} sequences"
            )

        elif params.database == DatabaseType.SILVA:
            logger.info("Using SILVA client")
            client = SILVAClient()
            search_result = client.search_sequences(params)
            results["SILVA"] = search_result
            logger.info(
                f"SILVA search completed: {len(search_result.sequences) if search_result.sequences else 0} sequences"
            )

        elif params.database == DatabaseType.UNITE:
            logger.info("Using UNITE client")
            client = UNITEClient()
            search_result = client.search_sequences(params)
            results["UNITE"] = search_result
            logger.info(
                f"UNITE search completed: {len(search_result.sequences) if search_result.sequences else 0} sequences"
            )

        else:
            # Fall back to original NCBI client for compatibility
            logger.info(
                f"Using legacy NCBI client for database: {params.database.value}"
            )
            try:
                df = ncbi.get_data(params.database.value, params.search_term)
                # Convert legacy format to our enhanced format
                legacy_sequences = []
                for _, row in df.iterrows():
                    from genetic_testing.sequence_fetching.database_models import (
                        SequenceRecord,
                    )

                    seq_record = SequenceRecord(
                        id=str(row.get("Id", "")),
                        title=str(row.get("Title", "")),
                        length=int(row.get("Length", 0)),
                        database=params.database,
                    )
                    legacy_sequences.append(seq_record)

                # Create a mock SearchResult
                from genetic_testing.sequence_fetching.database_models import (
                    SearchResult,
                )

                legacy_result = SearchResult(
                    query=params,
                    total_found=len(legacy_sequences),
                    sequences=legacy_sequences,
                    search_time=0.0,
                )
                results["NCBI Legacy"] = legacy_result
                logger.info(
                    f"Legacy search completed: {len(legacy_sequences)} sequences"
                )

            except Exception as legacy_error:
                logger.error(f"Legacy search also failed: {legacy_error}")
                raise

    except Exception as e:
        logger.error(f"Search failed: {e}", exc_info=True)
        st.error(f"Search failed: {str(e)}")
        results = {}

    return results


# Main UI Layout
# Divides the page layout into two columns of relative width
query_col, summary_col = st.columns(MAIN_PAGE_COLS_SIZE, gap=MAIN_PAGE_COLS_GAP)

with query_col:
    st.header("Enhanced Sequence Search & Retrieval")
    st.markdown(
        "🧬 **Phase 1 Enhanced Features**: Multi-database support, advanced filtering, and marker-specific searches"
    )
    st.markdown(NAVIGATE_WARNING_MD)

    # Show NCBI configuration status
    try:
        email_configured = bool(st.secrets.get("ENTREZ_EMAIL", ""))
        api_key_configured = bool(st.secrets.get("NCBI_API_KEY", ""))

        if email_configured or api_key_configured:
            with st.expander("⚙️ NCBI Configuration Status", expanded=False):
                col1, col2 = st.columns(2)
                with col1:
                    if email_configured:
                        st.success(f"✅ Email: {st.secrets.get('ENTREZ_EMAIL')}")
                    else:
                        st.warning("⚠️ Email not configured")

                with col2:
                    if api_key_configured:
                        st.success("✅ API Key: Configured")
                        st.info("🚀 Enhanced rate limits enabled (10 req/sec)")
                    else:
                        st.info("ℹ️ API Key: Not configured (3 req/sec)")

        # Debug option
        if st.checkbox("🔧 Enable Debug Mode", help="Show detailed logging information"):
            st.info("Debug mode enabled - check logs for detailed information")

    except (FileNotFoundError, AttributeError, KeyError):
        pass  # Secrets not available

    # Search mode selection
    search_mode = st.radio(
        "Search Mode",
        ["Basic Search", "Advanced Search", "Marker-Specific Search"],
        key="search_mode_radio",
        help="Choose your search approach",
    )

    # Database selection
    st.subheader("Database Selection")

    # Add database guidance
    with st.expander("🗄️ Database Guide - Choose the Right Database", expanded=False):
        col1, col2 = st.columns(2)
        with col1:
            st.markdown("### **NCBI Databases**")
            st.markdown(
                """
            **🧬 NCBI Nucleotide**: 
            - Comprehensive DNA/RNA sequences
            - All organisms, all markers
            - Best for: General searches, any marker
            
            **🧪 NCBI Gene**: 
            - Gene-centric records
            - Curated gene information
            - Best for: Specific gene searches
            """
            )

        with col2:
            st.markdown("### **Specialized Databases**")
            st.markdown(
                """
            **🐟 BOLD (Barcode of Life)**:
            - COI barcode sequences
            - Animal identification focus
            - Best for: COI, animal barcoding
            
            **🦠 SILVA (rRNA)**:
            - Curated ribosomal RNA
            - High-quality 16S/18S sequences
            - Best for: 16S, 18S, microbiome
            
            **🍄 UNITE (Fungal ITS)**:
            - Fungal identification
            - Curated ITS sequences
            - Best for: ITS, fungal studies
            """
            )

        st.info(
            "💡 **Tip**: Select multiple databases for comprehensive results. Each database has unique strengths!"
        )

    database_options = {
        "NCBI Nucleotide": "nucleotide",
        "NCBI Gene": "gene",
        "BOLD (COI Barcodes)": "bold",
        "SILVA (rRNA)": "silva",
        "UNITE (Fungal ITS)": "unite",
    }

    selected_db_names = st.multiselect(
        "Select databases to search",
        options=list(database_options.keys()),
        default=["NCBI Nucleotide"],
        help="💡 Click 'Database Guide' above to learn about each database's strengths",
    )

    if not selected_db_names:
        st.warning("Please select at least one database")
        st.stop()

    # Convert to database types
    selected_databases = [
        DatabaseType(database_options[name]) for name in selected_db_names
    ]

    # Search form based on mode
    with st.form("enhanced_query"):
        if search_mode == "Basic Search":
            # Basic search interface (similar to original)
            st.subheader("Basic Search")

            # Enhanced guidance for search terms
            with st.expander("📖 Search Term Examples & Tips", expanded=False):
                st.markdown("### 🧬 **Species & Organism Searches:**")
                st.code(
                    """
# Single species
Salmo salar
Atlantic salmon

# Multiple species (OR)
"Salmo salar" OR "Oncorhynchus mykiss"

# Genus level
Salmo[organism]
                """
                )

                st.markdown("### 🧪 **Gene & Marker Searches:**")
                st.code(
                    """
# Specific genes
COI[gene]
cytochrome oxidase subunit 1

# Gene + organism
human[organism] AND COI[gene]
Salmo salar[organism] AND 16S[gene]

# Multiple markers
COI[gene] OR cytochrome[gene]
                """
                )

                st.markdown("### 🔬 **Advanced NCBI Syntax:**")
                st.code(
                    """
# Exclude terms
Salmo[organism] NOT farmed
COI[gene] NOT partial

# Date ranges
Salmo salar AND 2020:2024[pdat]

# Sequence length
Salmo salar AND 500:800[slen]

# Complete sequences only
Salmo salar AND complete[title]
                """
                )

                st.markdown("### 💡 **Pro Tips:**")
                st.info(
                    """
                • Use quotes for exact phrases: "Atlantic salmon"
                • Use [organism] to specify taxonomic searches
                • Use [gene] to search for specific genes
                • Combine with AND, OR, NOT for complex queries
                • Check different databases - BOLD is great for COI, SILVA for rRNA
                """
                )

            search_term = st.text_input(
                "Search Term",
                placeholder="Try: Salmo salar, COI[gene], or human[organism] AND 16S[gene]",
                help="💡 Click 'Search Term Examples & Tips' above for detailed guidance and examples",
            )

            # Basic filters
            col1, col2 = st.columns(2)
            with col1:
                min_length = st.number_input("Min Length", min_value=0, value=0)
            with col2:
                max_length = st.number_input("Max Length", min_value=0, value=10000)

        elif search_mode == "Advanced Search":
            # Advanced search interface
            st.subheader("Advanced Search Parameters")

            st.info(
                "💡 **Advanced Search**: Use specific fields below for precise filtering. Leave fields empty to ignore them."
            )

            # Taxonomic search
            col1, col2 = st.columns(2)
            with col1:
                search_term = st.text_input(
                    "General Search Term",
                    placeholder="Optional: COI, mitochondrial, ribosomal",
                    help="General keywords to search for across all fields",
                )
                taxon = st.text_input(
                    "Target Taxon",
                    placeholder="Salmo salar, Salmonidae, Chordata",
                    help="Species name, genus, family, or higher taxonomic group",
                )
            with col2:
                country = st.text_input(
                    "Country/Region",
                    placeholder="Canada, Norway, North Atlantic",
                    help="Geographic origin of specimens",
                )
                exclude_taxa_input = st.text_area(
                    "Exclude Taxa",
                    placeholder="Salmo trutta\nOncorhynchus\nfarmed",
                    help="Taxa to exclude (one per line)",
                )
                exclude_taxa = (
                    [t.strip() for t in exclude_taxa_input.split("\n") if t.strip()]
                    if exclude_taxa_input
                    else []
                )

            # Sequence properties
            st.subheader("Sequence Filters")
            col1, col2, col3 = st.columns(3)
            with col1:
                min_length = st.number_input("Min Length", min_value=0, value=0)
                max_length = st.number_input("Max Length", min_value=0, value=50000)
            with col2:
                molecule_type = st.selectbox(
                    "Molecule Type",
                    ["Any", "DNA", "RNA", "Genomic", "mRNA"],
                    help="Filter by molecule type",
                )
            with col3:
                exclude_partial = st.checkbox("Exclude Partial Sequences", value=True)
                exclude_predicted = st.checkbox(
                    "Exclude Predicted Sequences", value=True
                )

        else:  # Marker-Specific Search
            st.subheader("Marker-Specific Search")

            st.success(
                "🎯 **Optimized for Molecular Markers**: This mode provides targeted searches with appropriate databases and length filters for each marker type."
            )

            # Marker selection with detailed descriptions
            marker_options = {
                "COI (Cytochrome Oxidase I) - Animal Barcoding": MarkerType.COI,
                "16S rRNA - Bacterial/Archaeal Identification": MarkerType.SSU_16S,
                "18S rRNA - Eukaryotic Identification": MarkerType.SSU_18S,
                "ITS (Internal Transcribed Spacer) - Fungal ID": MarkerType.ITS,
                "ITS1 - Fungal Identification": MarkerType.ITS1,
                "ITS2 - Fungal/Plant Identification": MarkerType.ITS2,
                "rbcL - Plant Barcoding": MarkerType.RBCL,
                "matK - Plant Identification": MarkerType.MATK,
            }

            selected_marker_name = st.selectbox(
                "Select Molecular Marker",
                options=list(marker_options.keys()),
                help="Choose the molecular marker - each has optimized search parameters and database recommendations",
            )
            selected_marker = marker_options[selected_marker_name]

            # Show marker-specific information
            marker_info = {
                MarkerType.COI: {
                    "description": "🐟 **COI (Cytochrome Oxidase I)**: Standard animal DNA barcode. ~650bp region used for species identification.",
                    "databases": "NCBI Nucleotide + BOLD",
                    "typical_length": "600-700 bp",
                    "use_cases": "Animal species identification, biodiversity studies, food authentication",
                },
                MarkerType.SSU_16S: {
                    "description": "🦠 **16S rRNA**: Bacterial and archaeal identification marker. Highly conserved with variable regions.",
                    "databases": "NCBI Nucleotide + SILVA",
                    "typical_length": "1400-1600 bp",
                    "use_cases": "Bacterial identification, microbiome studies, phylogenetics",
                },
                MarkerType.SSU_18S: {
                    "description": "🧬 **18S rRNA**: Eukaryotic small subunit ribosomal RNA for broad taxonomic identification.",
                    "databases": "NCBI Nucleotide + SILVA",
                    "typical_length": "1700-1900 bp",
                    "use_cases": "Eukaryotic identification, environmental DNA, phylogenetics",
                },
                MarkerType.ITS: {
                    "description": "🍄 **ITS (Internal Transcribed Spacer)**: Primary fungal identification marker with high variability.",
                    "databases": "NCBI Nucleotide + UNITE",
                    "typical_length": "400-800 bp",
                    "use_cases": "Fungal identification, mycological studies, environmental surveys",
                },
                MarkerType.RBCL: {
                    "description": "🌱 **rbcL**: Plant DNA barcode. Large subunit of RuBisCO enzyme, highly conserved.",
                    "databases": "NCBI Nucleotide",
                    "typical_length": "1300-1500 bp",
                    "use_cases": "Plant identification, botanical studies, phylogenetics",
                },
                MarkerType.MATK: {
                    "description": "🌿 **matK**: Plant identification marker. Maturase K gene with good species resolution.",
                    "databases": "NCBI Nucleotide",
                    "typical_length": "800-900 bp",
                    "use_cases": "Plant species identification, complementary to rbcL",
                },
            }

            if selected_marker in marker_info:
                info = marker_info[selected_marker]
                with st.expander(f"ℹ️ About {selected_marker.value}", expanded=False):
                    st.markdown(info["description"])
                    col1, col2, col3 = st.columns(3)
                    with col1:
                        st.metric("Recommended DBs", info["databases"])
                    with col2:
                        st.metric("Typical Length", info["typical_length"])
                    with col3:
                        st.metric("Primary Use", "ID & Phylo")
                    st.markdown(f"**Common Applications**: {info['use_cases']}")

            # Taxonomic input
            col1, col2 = st.columns(2)
            with col1:
                taxon = st.text_input("Target Taxon", placeholder="e.g., Salmo salar")
                search_term = taxon  # Use taxon as search term for marker-specific
            with col2:
                country = st.text_input("Country Filter", placeholder="Optional")

            # Sequence quality filters
            col1, col2 = st.columns(2)
            with col1:
                min_length = st.number_input("Min Length", min_value=0, value=200)
                exclude_partial = st.checkbox("Exclude Partial Sequences", value=True)
            with col2:
                max_length = st.number_input("Max Length", min_value=0, value=2000)
                exclude_predicted = st.checkbox(
                    "Exclude Predicted Sequences", value=True
                )

        # Common parameters
        max_results = st.slider(
            "Maximum Results per Database", min_value=10, max_value=5000, value=1000
        )

        # Submit button
        submitted = st.form_submit_button("🔍 Search Sequences", type="primary")

        if submitted:
            if not search_term.strip():
                st.error("Please enter a search term")
            else:
                with st.spinner("Searching databases..."):
                    # Prepare search parameters
                    search_results = {}

                    for db_type in selected_databases:
                        # Build parameters based on search mode
                        params = SearchParameters(
                            search_term=search_term,
                            database=db_type,
                            max_results=max_results,
                        )

                        # Add advanced parameters if available
                        if search_mode != "Basic Search":
                            if "taxon" in locals() and taxon:
                                params.taxon = taxon
                            if "country" in locals() and country:
                                params.country = country
                            if "exclude_taxa" in locals():
                                params.exclude_taxa = [
                                    t.strip() for t in exclude_taxa if t.strip()
                                ]
                            if "exclude_partial" in locals():
                                params.exclude_partial = exclude_partial
                            if "exclude_predicted" in locals():
                                params.exclude_predicted = exclude_predicted

                        # Set length filters
                        if min_length > 0:
                            params.min_length = min_length
                        if (
                            max_length > 0 and max_length != 50000
                        ):  # Don't set if default
                            params.max_length = max_length

                        # Set molecule type
                        if search_mode == "Advanced Search" and molecule_type != "Any":
                            params.molecule_type = MoleculeType(molecule_type.upper())

                        # Set marker for marker-specific search
                        if search_mode == "Marker-Specific Search":
                            params.marker = selected_marker

                        # Perform search
                        try:
                            db_result = perform_enhanced_search(params)
                            search_results.update(db_result)
                        except Exception as e:
                            logger.error(f"Error searching {db_type}: {e}")
                            st.error(f"Error searching {db_type.value}: {str(e)}")

                    # Store results in session state
                    st.session_state.search_results = search_results
                    st.session_state[NCBI_SUMMARY_FORM] = bool(search_results)

                    # Convert to DataFrame for compatibility with existing code
                    if search_results:
                        try:
                            df_combined = convert_search_result_to_df(search_results)

                            if not df_combined.empty:
                                # Debug: Log DataFrame info
                                logger.info(
                                    f"DataFrame columns: {df_combined.columns.tolist()}"
                                )
                                logger.info(f"DataFrame shape: {df_combined.shape}")

                                st.session_state[NCBI_DF] = df_combined

                                # Only try to format if we have the expected columns
                                if "Title" in df_combined.columns:
                                    try:
                                        data_processing.format_ncbi_summary()
                                        logger.info(
                                            "Successfully formatted NCBI summary"
                                        )
                                    except Exception as format_error:
                                        logger.warning(
                                            f"Failed to format NCBI summary: {format_error}"
                                        )
                                        st.warning(
                                            "⚠️ Data formatting failed, but results are still available"
                                        )
                                else:
                                    logger.warning(
                                        f"Missing 'Title' column. Available columns: {df_combined.columns.tolist()}"
                                    )
                                    st.warning(
                                        "⚠️ Enhanced search format detected - some legacy features may not work"
                                    )

                                st.session_state[NCBI_SUMMARY_FORM] = True
                                st.success(
                                    f"Found {len(df_combined)} sequences across {len(search_results)} databases"
                                )
                            else:
                                st.warning(
                                    "Search completed but no sequences were retrieved"
                                )
                                st.session_state[NCBI_SUMMARY_FORM] = False

                        except Exception as df_error:
                            logger.error(
                                f"Error processing search results: {df_error}",
                                exc_info=True,
                            )
                            st.error(
                                f"Error processing search results: {str(df_error)}"
                            )
                            st.session_state[NCBI_SUMMARY_FORM] = False
                    else:
                        st.error("No sequences found matching your criteria")
                    st.session_state[NCBI_SUMMARY_FORM] = False
    df_aggrid = pd.DataFrame()
    if st.session_state[NCBI_SUMMARY_FORM]:
        aggrid_table = aggrid_table()  # Initialize and Render Aggrid Table
        df_aggrid = aggrid_table["data"]
    # Streamlit UI to filter the dataset
    df_filtered, applied_filters = filter_dataframe(df_aggrid)
    if st.session_state[NCBI_DF_FILTER]:
        st.dataframe(df_filtered)
        # Streamlit UI to show the applied filters
        show_filters_info(applied_filters)
    # User has made a search query
    if not df_aggrid.empty:
        df_download = df_filtered if applied_filters else df_aggrid
        total_length = df_download["Length"].sum()
        st.write(
            f"Size of final data: {len(df_download)} rows, {total_length:,} total base pairs"
        )
        if total_length > MAX_SEQUENCE_LENGTH_DOWNLOAD:
            st.warning(
                f"The total length of the sequences is {total_length:,}, which is more than the allowed max "
                f"({MAX_SEQUENCE_LENGTH_DOWNLOAD:,}). Please apply filters to reduce the size of the dataset."
            )
        else:
            # Dropdown menu for selecting download format
            download_format = st.selectbox(
                "Select an option:", [FASTA_DOWNLOAD, CSV_DOWNLOAD]
            )
            # Option 1: Download the correct DataFrame as CSV
            if download_format == CSV_DOWNLOAD:
                st.download_button(
                    label=CSV_DOWNLOAD,
                    key="download_csv",
                    data=df_download.to_csv(index=False).encode("utf-8"),
                    file_name="sequence_search_table.csv",
                    mime="text/csv",
                    help="Click button to Download this table as a CSV file",
                    args={"as_attachment": True},
                )
            # Option 2: Download the correct DataFrame as FASTA
            elif download_format == FASTA_DOWNLOAD:
                # Use nucleotide database as default for FASTA download
                # In the enhanced version, we'll need to determine the appropriate database
                default_database = "nucleotide"
                default_search_term = "sequence_search"

                # Try to get search term from session state if available
                if (
                    hasattr(st.session_state, "search_results")
                    and st.session_state.search_results
                ):
                    # Use the first database found in search results
                    db_names = list(st.session_state.search_results.keys())
                    if db_names:
                        if "NCBI Nucleotide" in db_names:
                            default_database = "nucleotide"
                        elif "NCBI Gene" in db_names:
                            default_database = "gene"
                        default_search_term = "enhanced_search"

                download_data(
                    default_database, df_download["Id"].unique(), default_search_term
                )
with summary_col:
    # Show the summary only if user has submitted the query form
    if st.session_state[NCBI_SUMMARY_FORM]:
        st.header("Search Results Summary")

        # Database-specific results
        if (
            hasattr(st.session_state, "search_results")
            and st.session_state.search_results
        ):
            st.subheader("Results by Database")

            for db_name, result in st.session_state.search_results.items():
                if hasattr(result, "sequences"):  # New format
                    seq_count = len(result.sequences)
                    search_time = (
                        f"{result.search_time:.2f}s"
                        if hasattr(result, "search_time")
                        else "N/A"
                    )

                    with st.expander(
                        f"📊 {db_name} ({seq_count} sequences)", expanded=True
                    ):
                        st.metric("Sequences Found", seq_count)
                        st.metric("Search Time", search_time)

                        if hasattr(result, "warnings") and result.warnings:
                            st.warning("⚠️ " + "; ".join(result.warnings))

                        # Show marker distribution if available
                        if result.sequences:
                            markers = [
                                seq.marker.value
                                for seq in result.sequences
                                if seq.marker
                            ]
                            if markers:
                                marker_counts = pd.Series(markers).value_counts()
                                st.write("**Markers Found:**")
                                for marker, count in marker_counts.items():
                                    st.write(f"• {marker}: {count}")

                            # Show organism distribution
                            organisms = [
                                seq.organism for seq in result.sequences if seq.organism
                            ]
                            if organisms:
                                org_counts = pd.Series(organisms).value_counts().head(5)
                                st.write("**Top Organisms:**")
                                for org, count in org_counts.items():
                                    st.write(f"• {org}: {count}")

        # Overall statistics
        if hasattr(st.session_state, NCBI_DF) and not st.session_state[NCBI_DF].empty:
            df = st.session_state[NCBI_DF]

            st.subheader("Overall Statistics")

            # Total sequences and databases
            col1, col2 = st.columns(2)
            with col1:
                st.metric("Total Sequences", len(df))
            with col2:
                if "Database" in df.columns:
                    st.metric("Databases Used", df["Database"].nunique())

            # Length statistics
            if "Length" in df.columns:
                st.subheader("Sequence Length Distribution")
                length_stats = df["Length"].describe()

                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Min Length", f"{int(length_stats['min']):,} bp")
                with col2:
                    st.metric("Mean Length", f"{int(length_stats['mean']):,} bp")
                with col3:
                    st.metric("Max Length", f"{int(length_stats['max']):,} bp")

                # Length histogram
                st.bar_chart(df["Length"].value_counts().sort_index())

            # Quality scores if available
            if "Quality Score" in df.columns:
                quality_scores = df["Quality Score"].replace("", None).dropna()
                if not quality_scores.empty:
                    quality_numeric = pd.to_numeric(
                        quality_scores, errors="coerce"
                    ).dropna()
                    if not quality_numeric.empty:
                        st.subheader("Quality Score Distribution")
                        st.metric("Average Quality", f"{quality_numeric.mean():.3f}")
                        st.bar_chart(quality_numeric.value_counts().sort_index())

            # Top organisms (legacy compatibility)
            st.subheader("Top Organisms")
            if "Organism" in df.columns:
                top_orgs = df["Organism"].value_counts().head(TOP_N_ORGANISMS)
                for org, count in top_orgs.items():
                    if org:  # Skip empty organism names
                        st.write(f"• **{org}**: {count} sequences")
            else:
                # Fall back to original method
                st.write(data_processing.get_top_organisms_counts(TOP_N_ORGANISMS))

st.sidebar.image("src/app/Conservation X Labs CXL logo.png", width="stretch")
