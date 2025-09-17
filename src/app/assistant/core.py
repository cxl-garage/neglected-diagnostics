"""
Core Agentic Assistant Implementation

This module contains the main assistant logic, OpenAI integration, and
context-aware response generation.
"""

import json
import os
import time
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

try:
    import streamlit as st
except ImportError:
    # Mock streamlit for non-web environments
    class MockStreamlit:
        class secrets:
            @staticmethod
            def get(key, default=""):
                return os.environ.get(key, default)

    st = MockStreamlit()

try:
    from openai import OpenAI
except ImportError:
    # Mock OpenAI for environments without the package
    class MockOpenAI:
        def __init__(self, api_key=None):
            pass

        class chat:
            class completions:
                @staticmethod
                def create(**kwargs):
                    class MockResponse:
                        class choices:
                            class message:
                                content = (
                                    "OpenAI not available - using fallback response"
                                )
                                function_call = None

                        choices = [choices()]

                    return MockResponse()

    OpenAI = MockOpenAI

from utils.log import _init_logger

from .context_manager import PageContextManager
from .prompts import PromptTemplates

logger = _init_logger(__name__)


class AgenticAssistant:
    """Main agentic assistant class with OpenAI integration."""

    def __init__(self):
        """Initialize the assistant with OpenAI client and context manager."""
        self.client = self._initialize_openai_client()
        self.context_manager = PageContextManager()
        self.prompts = PromptTemplates()
        self.conversation_history = []

        # Assistant configuration
        self.model = "gpt-4o-mini"  # Cost-effective model for assistance
        self.max_tokens = 1000
        self.temperature = 0.7

    def _initialize_openai_client(self) -> Optional[OpenAI]:
        """Initialize OpenAI client with API key."""
        try:
            # Try to get API key from Streamlit secrets first
            api_key = None
            try:
                api_key = st.secrets.get("OPENAI_API_KEY", "")
            except (FileNotFoundError, AttributeError):
                pass

            # Fall back to environment variable
            if not api_key:
                api_key = os.environ.get("OPENAI_API_KEY", "")

            if api_key:
                client = OpenAI(api_key=api_key)
                logger.info("OpenAI client initialized successfully")
                return client
            else:
                logger.warning(
                    "No OpenAI API key found. Assistant will use fallback responses."
                )
                return None

        except Exception as e:
            logger.error(f"Failed to initialize OpenAI client: {e}")
            return None

    def get_assistance(
        self,
        user_query: str,
        page_context: Dict[str, Any],
        conversation_history: Optional[List[Dict]] = None,
    ) -> Tuple[str, Dict[str, Any]]:
        """
        Get assistance response for user query with page context.

        Parameters
        ----------
        user_query : str
            User's question or request for help
        page_context : Dict[str, Any]
            Current page context and state
        conversation_history : List[Dict], optional
            Previous conversation messages

        Returns
        -------
        Tuple[str, Dict[str, Any]]
            Assistant response and any suggested actions
        """
        try:
            if self.client:
                return self._get_ai_response(
                    user_query, page_context, conversation_history
                )
            else:
                return self._get_fallback_response(user_query, page_context)

        except Exception as e:
            logger.error(f"Error getting assistance: {e}")
            return self._get_error_response(str(e))

    def _get_ai_response(
        self,
        user_query: str,
        page_context: Dict[str, Any],
        conversation_history: Optional[List[Dict]] = None,
    ) -> Tuple[str, Dict[str, Any]]:
        """Get response from OpenAI API."""

        # Build system prompt with context
        system_prompt = self.prompts.get_system_prompt(page_context)

        # Build conversation messages
        messages = [{"role": "system", "content": system_prompt}]

        # Add conversation history if provided
        if conversation_history:
            messages.extend(
                conversation_history[-6:]
            )  # Keep last 6 messages for context

        # Add current user query
        messages.append({"role": "user", "content": user_query})

        try:
            response = self.client.chat.completions.create(
                model=self.model,
                messages=messages,
                max_tokens=self.max_tokens,
                temperature=self.temperature,
                functions=[
                    {
                        "name": "suggest_search_parameters",
                        "description": "Suggest search parameters and database selections",
                        "parameters": {
                            "type": "object",
                            "properties": {
                                "search_term": {
                                    "type": "string",
                                    "description": "Suggested search term",
                                },
                                "databases": {
                                    "type": "array",
                                    "items": {"type": "string"},
                                    "description": "Recommended databases",
                                },
                                "filters": {
                                    "type": "object",
                                    "description": "Suggested filter parameters",
                                },
                                "explanation": {
                                    "type": "string",
                                    "description": "Explanation of the suggestions",
                                },
                            },
                        },
                    }
                ],
                function_call="auto",
            )

            message = response.choices[0].message

            # Check if function was called
            if message.function_call:
                function_response = json.loads(message.function_call.arguments)
                return self._format_function_response(
                    message.content or "", function_response
                )
            else:
                return message.content or "I'm here to help!", {}

        except Exception as e:
            logger.error(f"OpenAI API error: {e}")
            return self._get_fallback_response(user_query, page_context)

    def _get_fallback_response(
        self, user_query: str, page_context: Dict[str, Any]
    ) -> Tuple[str, Dict[str, Any]]:
        """Provide fallback response when OpenAI is not available."""

        query_lower = user_query.lower()
        page_name = page_context.get("page_name", "unknown")

        # Pattern matching for common queries
        if any(word in query_lower for word in ["search", "query", "term"]):
            return self._get_search_help(page_context)
        elif any(word in query_lower for word in ["database", "db", "ncbi", "bold"]):
            return self._get_database_help(page_context)
        elif any(word in query_lower for word in ["filter", "length", "quality"]):
            return self._get_filter_help(page_context)
        elif any(
            word in query_lower for word in ["marker", "gene", "coi", "16s", "its"]
        ):
            return self._get_marker_help(page_context)
        else:
            return self._get_general_help(page_context)

    def _get_search_help(
        self, page_context: Dict[str, Any]
    ) -> Tuple[str, Dict[str, Any]]:
        """Provide search term assistance."""
        suggestions = {
            "search_term": "Salmo salar",
            "databases": ["NCBI Nucleotide"],
            "explanation": "Try starting with a species name like 'Salmo salar' for Atlantic salmon",
        }

        response = """🔍 **Search Term Help**

For effective searches, try these approaches:

**Species Names**: Use scientific names like "Salmo salar" or "Homo sapiens"
**Gene Names**: Add [gene] tag like "COI[gene]" or "16S[gene]"
**Combined**: Mix species and genes: "Salmo salar AND COI[gene]"

**Examples**:
• `Salmo salar` - All sequences for Atlantic salmon
• `COI[gene]` - All COI gene sequences
• `human[organism] AND mitochondrial` - Human mitochondrial sequences

Would you like me to suggest specific search terms for your research?"""

        return response, suggestions

    def _get_database_help(
        self, page_context: Dict[str, Any]
    ) -> Tuple[str, Dict[str, Any]]:
        """Provide database selection assistance."""
        response = """🗄️ **Database Selection Guide**

Choose databases based on your research focus:

**NCBI Nucleotide**: 📊 General purpose, all organisms, comprehensive
**BOLD**: 🐟 Animal barcoding, COI sequences, specimen data
**SILVA**: 🦠 rRNA sequences, microbial identification
**UNITE**: 🍄 Fungal ITS sequences, mycology research

**Quick Recommendations**:
• Animal identification → NCBI + BOLD
• Microbial studies → NCBI + SILVA  
• Fungal research → NCBI + UNITE
• Plant studies → NCBI Nucleotide

What type of organism are you studying?"""

        return response, {}

    def _get_filter_help(
        self, page_context: Dict[str, Any]
    ) -> Tuple[str, Dict[str, Any]]:
        """Provide filtering assistance."""
        response = """⚙️ **Filtering Options**

Optimize your results with smart filtering:

**Length Filters**:
• COI: 600-700 bp (animal barcoding)
• 16S: 1400-1600 bp (bacterial ID)
• ITS: 400-800 bp (fungal ID)

**Quality Filters**:
• ✅ Exclude partial sequences (recommended)
• ✅ Exclude predicted sequences (recommended)
• Consider geographic filters for regional studies

**Pro Tips**:
• Start broad, then narrow down
• Use quality filters to improve results
• Length filters help focus on complete sequences

What type of sequences are you looking for?"""

        return response, {}

    def _get_marker_help(
        self, page_context: Dict[str, Any]
    ) -> Tuple[str, Dict[str, Any]]:
        """Provide molecular marker assistance."""
        response = """🧬 **Molecular Marker Guide**

Choose the right marker for your study:

**COI (Cytochrome Oxidase I)**: 🐟
• Animal species identification
• DNA barcoding standard
• ~650 bp, mitochondrial

**16S rRNA**: 🦠  
• Bacterial/archaeal identification
• Phylogenetic studies
• ~1500 bp, highly conserved

**ITS (Internal Transcribed Spacer)**: 🍄
• Fungal identification
• High species resolution
• Variable length, nuclear

**18S rRNA**: 🧬
• Eukaryotic identification
• Broad taxonomic coverage
• ~1800 bp

What organism group are you studying?"""

        return response, {}

    def _get_general_help(
        self, page_context: Dict[str, Any]
    ) -> Tuple[str, Dict[str, Any]]:
        """Provide general assistance."""
        page_name = page_context.get("page_name", "this page")

        response = f"""👋 **Hello! I'm your research assistant.**

I'm here to help you with {page_name}. I can assist with:

🔍 **Search Terms**: Crafting effective queries
🗄️ **Database Selection**: Choosing the right databases  
⚙️ **Filtering**: Optimizing your results
🧬 **Markers**: Understanding molecular markers
🛠️ **Troubleshooting**: Solving search issues

**Quick Start**:
• Ask me "How do I search for salmon sequences?"
• Say "Help me find bacterial 16S sequences"
• Request "Suggest databases for fungal research"

What would you like help with today?"""

        return response, {}

    def _get_error_response(self, error_msg: str) -> Tuple[str, Dict[str, Any]]:
        """Provide error response."""
        response = f"""⚠️ **Assistant Temporarily Unavailable**

I encountered an issue: {error_msg}

**In the meantime, here are some quick tips**:
• Use scientific names for species searches
• Try the "Search Term Examples & Tips" sections
• Check the database guide for recommendations
• Start with basic searches and add filters gradually

I'll be back online shortly to provide personalized assistance!"""

        return response, {}

    def _format_function_response(
        self, content: str, function_response: Dict[str, Any]
    ) -> Tuple[str, Dict[str, Any]]:
        """Format response when function calling is used."""

        # Extract suggestions
        search_term = function_response.get("search_term", "")
        databases = function_response.get("databases", [])
        filters = function_response.get("filters", {})
        explanation = function_response.get("explanation", "")

        # Build formatted response
        response_parts = []

        if content:
            response_parts.append(content)

        if explanation:
            response_parts.append(f"💡 **Recommendation**: {explanation}")

        if search_term:
            response_parts.append(f"🔍 **Suggested Search**: `{search_term}`")

        if databases:
            db_list = ", ".join(databases)
            response_parts.append(f"🗄️ **Recommended Databases**: {db_list}")

        if filters:
            response_parts.append("⚙️ **Suggested Filters**:")
            for key, value in filters.items():
                response_parts.append(f"• {key}: {value}")

        formatted_response = "\n\n".join(response_parts)

        return formatted_response, function_response

    def get_page_specific_help(self, page_name: str) -> str:
        """Get page-specific help content."""
        help_content = {
            "Sequence_Search": """
🧬 **Sequence Search Assistant**

I can help you with:
• Crafting effective search queries
• Selecting appropriate databases  
• Setting up filters for quality results
• Understanding molecular markers
• Troubleshooting search issues

Try asking: "Help me find COI sequences for salmon"
""",
            "Multisequences_Alignment": """
🧬 **Alignment Assistant** 

I can help you with:
• Preparing sequences for alignment
• Choosing alignment parameters
• Interpreting alignment results
• Quality control of alignments

Try asking: "How do I align my sequences?"
""",
            "default": """
👋 **Research Assistant**

I'm here to help with your bioinformatics research!
Ask me about search strategies, databases, or analysis methods.
""",
        }

        return help_content.get(page_name, help_content["default"])
