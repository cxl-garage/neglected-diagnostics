"""
Page Context Manager

This module manages context information for different pages to provide
tailored assistance based on the current page and user state.
"""

from datetime import datetime
from typing import Any, Dict, List, Optional

try:
    import streamlit as st
except ImportError:
    # Mock streamlit for non-web environments
    class MockSessionState:
        def __init__(self):
            self._state = {}

        def get(self, key, default=None):
            return self._state.get(key, default)

        def __contains__(self, key):
            return key in self._state

        def __getitem__(self, key):
            return self._state[key]

        def __setitem__(self, key, value):
            self._state[key] = value

    class MockStreamlit:
        session_state = MockSessionState()

    st = MockStreamlit()

from utils.log import _init_logger

logger = _init_logger(__name__)


class PageContextManager:
    """Manages page-specific context for the assistant."""

    def __init__(self):
        """Initialize the context manager."""
        self.page_configs = {
            "Sequence_Search": {
                "name": "Sequence Search",
                "description": "Search and retrieve sequences from multiple databases",
                "key_features": [
                    "Multi-database search (NCBI, BOLD, SILVA, UNITE)",
                    "Advanced filtering options",
                    "Marker-specific searches",
                    "Quality control filters",
                ],
                "common_tasks": [
                    "Find sequences for specific species",
                    "Search by molecular markers",
                    "Filter by sequence quality",
                    "Download sequences in FASTA format",
                ],
                "databases": ["NCBI Nucleotide", "NCBI Gene", "BOLD", "SILVA", "UNITE"],
                "markers": ["COI", "16S", "18S", "ITS", "rbcL", "matK"],
            },
            "Multisequences_Alignment": {
                "name": "Multiple Sequence Alignment",
                "description": "Align multiple sequences and analyze conservation",
                "key_features": [
                    "Multiple sequence alignment",
                    "Consensus sequence generation",
                    "Conservation analysis",
                    "Alignment visualization",
                ],
                "common_tasks": [
                    "Align sequences from search results",
                    "Generate consensus sequences",
                    "Identify conserved regions",
                    "Export alignment data",
                ],
            },
            "Sequence_Haplotypes": {
                "name": "Sequence Haplotypes",
                "description": "Analyze sequence haplotypes and variants",
                "key_features": [
                    "Haplotype identification",
                    "Variant analysis",
                    "Population genetics metrics",
                ],
                "common_tasks": [
                    "Identify unique haplotypes",
                    "Analyze sequence variants",
                    "Calculate diversity metrics",
                ],
            },
            "Sequence_Variability": {
                "name": "Sequence Variability",
                "description": "Analyze sequence variability and conservation",
                "key_features": [
                    "Variability analysis",
                    "Conservation scoring",
                    "Position-specific analysis",
                ],
                "common_tasks": [
                    "Calculate sequence variability",
                    "Identify conserved regions",
                    "Generate variability plots",
                ],
            },
            "Find_Assay_Design_Area": {
                "name": "Assay Design Area",
                "description": "Find optimal regions for primer/probe design",
                "key_features": [
                    "Target region identification",
                    "Primer design optimization",
                    "Specificity analysis",
                ],
                "common_tasks": [
                    "Identify target regions",
                    "Design primers and probes",
                    "Validate specificity",
                ],
            },
            "Galaxy_Workflow_Assistant": {
                "name": "Galaxy Workflow Assistant",
                "description": "Manage Galaxy workflows and analyses",
                "key_features": [
                    "Workflow management",
                    "Galaxy integration",
                    "Analysis automation",
                ],
                "common_tasks": ["Create workflows", "Run analyses", "Manage data"],
            },
        }

    def get_current_context(self) -> Dict[str, Any]:
        """Get current page context from Streamlit session state."""

        # Determine current page
        current_page = self._get_current_page()
        page_config = self.page_configs.get(current_page, {})

        # Build context dictionary
        context = {
            "page_name": page_config.get("name", current_page),
            "page_description": page_config.get("description", ""),
            "key_features": page_config.get("key_features", []),
            "common_tasks": page_config.get("common_tasks", []),
            "databases": page_config.get("databases", []),
            "markers": page_config.get("markers", []),
            "timestamp": datetime.now().isoformat(),
            "session_state": self._get_relevant_session_state(),
            "user_inputs": self._get_current_user_inputs(),
        }

        return context

    def _get_current_page(self) -> str:
        """Determine the current page from Streamlit context."""
        try:
            # Try to get page from Streamlit's page context
            if hasattr(st, "_get_page_name"):
                return st._get_page_name()

            # Try to get from session state
            if "current_page" in st.session_state:
                return st.session_state.current_page

            # Try to infer from URL or other context
            # This is a fallback - in practice, you might set this explicitly
            return "Sequence_Search"  # Default to main search page

        except Exception as e:
            logger.warning(f"Could not determine current page: {e}")
            return "Unknown"

    def _get_relevant_session_state(self) -> Dict[str, Any]:
        """Extract relevant session state information."""
        relevant_keys = [
            "search_results",
            "selected_databases",
            "search_mode",
            "NCBI_SUMMARY_FORM",
            "NCBI_DF",
        ]

        session_info = {}
        for key in relevant_keys:
            if key in st.session_state:
                value = st.session_state[key]
                # Serialize complex objects safely
                if isinstance(value, (str, int, float, bool, list)):
                    session_info[key] = value
                elif hasattr(value, "__len__"):
                    session_info[
                        key
                    ] = f"<{type(value).__name__} with {len(value)} items>"
                else:
                    session_info[key] = f"<{type(value).__name__}>"

        return session_info

    def _get_current_user_inputs(self) -> Dict[str, Any]:
        """Get current user inputs from the page (if available)."""
        # This would be populated by the UI components
        # For now, return empty dict as it requires integration with UI
        return {}

    def get_search_context(self) -> Dict[str, Any]:
        """Get specific context for sequence search page."""
        context = self.get_current_context()

        # Add search-specific information
        if "search_results" in st.session_state and st.session_state.search_results:
            results = st.session_state.search_results
            context["recent_search"] = {
                "databases_searched": list(results.keys()),
                "total_sequences": sum(
                    len(result.sequences)
                    if hasattr(result, "sequences") and result.sequences
                    else 0
                    for result in results.values()
                ),
                "search_successful": bool(results),
            }

        return context

    def get_alignment_context(self) -> Dict[str, Any]:
        """Get specific context for alignment page."""
        context = self.get_current_context()

        # Add alignment-specific information
        # This would be populated based on alignment page state
        context["alignment_status"] = {
            "has_sequences": "NCBI_DF" in st.session_state
            and not st.session_state.NCBI_DF.empty,
            "ready_for_alignment": False,  # Would be determined by actual state
        }

        return context

    def update_context(self, key: str, value: Any) -> None:
        """Update context information."""
        if "assistant_context" not in st.session_state:
            st.session_state.assistant_context = {}

        st.session_state.assistant_context[key] = value
        logger.debug(f"Updated assistant context: {key} = {value}")

    def get_help_topics_for_page(self, page_name: str) -> List[str]:
        """Get relevant help topics for a specific page."""
        page_config = self.page_configs.get(page_name, {})

        base_topics = [
            "Getting started",
            "Common workflows",
            "Troubleshooting",
            "Best practices",
        ]

        # Add page-specific topics
        if page_name == "Sequence_Search":
            base_topics.extend(
                [
                    "Search term syntax",
                    "Database selection",
                    "Filtering options",
                    "Molecular markers",
                    "Download formats",
                ]
            )
        elif page_name == "Multisequences_Alignment":
            base_topics.extend(
                [
                    "Alignment parameters",
                    "Quality control",
                    "Consensus sequences",
                    "Export options",
                ]
            )

        return base_topics

    def get_suggested_queries(self, page_name: str) -> List[str]:
        """Get suggested queries for the assistant based on current page."""
        suggestions = {
            "Sequence_Search": [
                "How do I search for COI sequences?",
                "What databases should I use for bacterial identification?",
                "Help me craft a search term for salmon sequences",
                "How do I filter out partial sequences?",
                "What's the difference between NCBI and BOLD databases?",
            ],
            "Multisequences_Alignment": [
                "How do I prepare sequences for alignment?",
                "What alignment parameters should I use?",
                "How do I interpret the alignment results?",
                "Can you help me generate a consensus sequence?",
            ],
            "default": [
                "How do I get started?",
                "What can this tool do?",
                "Help me with my research workflow",
                "What are the best practices?",
            ],
        }

        return suggestions.get(page_name, suggestions["default"])
