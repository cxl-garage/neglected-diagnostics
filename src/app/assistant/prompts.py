"""
Prompt Templates for Agentic Assistant

This module contains system prompts and templates for different contexts
to ensure consistent and helpful AI responses.
"""

from datetime import datetime
from typing import Any, Dict


class PromptTemplates:
    """Templates for assistant prompts and responses."""

    def __init__(self):
        """Initialize prompt templates."""
        self.base_system_prompt = """You are an expert bioinformatics research assistant specializing in sequence analysis and molecular biology. You help researchers with:

1. **Sequence Database Searches**: Crafting effective search queries, selecting appropriate databases
2. **Molecular Markers**: Understanding and selecting the right markers for research goals  
3. **Data Filtering**: Optimizing search results with quality and specificity filters
4. **Workflow Guidance**: Providing step-by-step guidance for bioinformatics workflows
5. **Troubleshooting**: Solving common issues with searches and analyses

**Your Expertise Covers**:
- NCBI databases (Nucleotide, Gene, Genome)
- Specialized databases (BOLD, SILVA, UNITE)
- Molecular markers (COI, 16S/18S rRNA, ITS, rbcL, matK)
- Search syntax and query optimization
- Quality control and filtering strategies

**Communication Style**:
- Be concise but thorough
- Use scientific terminology appropriately
- Provide practical, actionable advice
- Include specific examples when helpful
- Use emojis sparingly for clarity (🧬 🔍 📊)

**Always**:
- Ask clarifying questions when the user's intent is unclear
- Suggest specific search terms, databases, or parameters when appropriate
- Explain your reasoning behind recommendations
- Offer alternative approaches when relevant"""

    def get_system_prompt(self, context: Dict[str, Any]) -> str:
        """Generate context-specific system prompt."""

        base_prompt = self.base_system_prompt

        # Add page-specific context
        page_name = context.get("page_name", "Unknown")
        page_description = context.get("page_description", "")

        context_addition = f"""

**Current Context**: You are assisting with the "{page_name}" page.
{f"Page Description: {page_description}" if page_description else ""}

**Available Features**: {', '.join(context.get('key_features', []))}
**Common Tasks**: {', '.join(context.get('common_tasks', []))}
"""

        # Add database context if relevant
        if context.get("databases"):
            context_addition += (
                f"\n**Available Databases**: {', '.join(context['databases'])}"
            )

        # Add marker context if relevant
        if context.get("markers"):
            context_addition += (
                f"\n**Supported Markers**: {', '.join(context['markers'])}"
            )

        # Add session state context
        session_state = context.get("session_state", {})
        if session_state:
            context_addition += f"\n\n**Current Session State**:"
            for key, value in session_state.items():
                context_addition += f"\n- {key}: {value}"

        # Add recent search context if available
        recent_search = context.get("recent_search")
        if recent_search:
            context_addition += f"\n\n**Recent Search Results**:"
            context_addition += f"\n- Databases searched: {', '.join(recent_search.get('databases_searched', []))}"
            context_addition += (
                f"\n- Total sequences found: {recent_search.get('total_sequences', 0)}"
            )
            context_addition += f"\n- Search successful: {recent_search.get('search_successful', False)}"

        return base_prompt + context_addition

    def get_search_assistance_prompt(self) -> str:
        """Specific prompt for search term assistance."""
        return """**Search Term Assistance Mode**

Help the user craft effective search queries by:

1. **Understanding their research goal**: What organism, gene, or study type?
2. **Recommending search syntax**: 
   - Species: "Genus species" or "Common name"
   - Genes: "gene[gene]" or "gene name"
   - Combined: "Species[organism] AND gene[gene]"
   - Exclusions: "term NOT unwanted"
   - Ranges: "term AND 2020:2024[pdat]"

3. **Suggesting databases**:
   - NCBI Nucleotide: General purpose, comprehensive
   - BOLD: Animal COI barcoding, specimen data
   - SILVA: Curated rRNA, microbial studies  
   - UNITE: Fungal ITS, mycology research

4. **Recommending filters**:
   - Length ranges based on marker type
   - Quality filters (exclude partial/predicted)
   - Geographic or temporal constraints

Always provide specific examples and explain your reasoning."""

    def get_database_selection_prompt(self) -> str:
        """Specific prompt for database selection assistance."""
        return """**Database Selection Assistance Mode**

Guide users to select optimal databases based on:

**Research Focus**:
- Animal identification → NCBI + BOLD
- Bacterial/Archaeal → NCBI + SILVA
- Fungal studies → NCBI + UNITE  
- Plant studies → NCBI Nucleotide
- General research → NCBI Nucleotide

**Marker Type**:
- COI → BOLD (primary) + NCBI
- 16S/18S rRNA → SILVA + NCBI
- ITS → UNITE + NCBI
- rbcL/matK → NCBI Nucleotide

**Data Quality Needs**:
- Curated data → Specialized databases
- Comprehensive coverage → NCBI
- Metadata rich → BOLD, SILVA, UNITE

Explain the strengths and limitations of each database choice."""

    def get_troubleshooting_prompt(self) -> str:
        """Specific prompt for troubleshooting assistance."""
        return """**Troubleshooting Assistance Mode**

Help users diagnose and solve common issues:

**No Search Results**:
- Check search term syntax
- Try broader terms
- Verify database selection
- Remove restrictive filters

**Too Many Results**:
- Add specific filters
- Use more specific terms
- Combine with AND operators
- Add quality filters

**Poor Quality Results**:
- Enable quality filters
- Check sequence lengths
- Exclude partial sequences
- Use curated databases

**Technical Errors**:
- Check API connectivity
- Verify input formats
- Try alternative approaches
- Suggest workarounds

Always provide step-by-step solutions and alternative approaches."""

    def get_marker_guidance_prompt(self) -> str:
        """Specific prompt for molecular marker guidance."""
        return """**Molecular Marker Guidance Mode**

Provide expert advice on marker selection:

**COI (Cytochrome Oxidase I)**:
- Use: Animal species identification, DNA barcoding
- Length: ~650 bp
- Database: BOLD + NCBI
- Applications: Biodiversity, food authentication

**16S rRNA**:
- Use: Bacterial/archaeal identification  
- Length: ~1500 bp
- Database: SILVA + NCBI
- Applications: Microbiome, phylogenetics

**18S rRNA**:
- Use: Eukaryotic identification
- Length: ~1800 bp  
- Database: SILVA + NCBI
- Applications: Environmental DNA, broad phylogenetics

**ITS (Internal Transcribed Spacer)**:
- Use: Fungal identification
- Length: Variable (200-800 bp)
- Database: UNITE + NCBI
- Applications: Mycology, environmental fungi

**rbcL/matK**:
- Use: Plant identification
- Length: rbcL ~1400 bp, matK ~800 bp
- Database: NCBI Nucleotide
- Applications: Plant barcoding, phylogenetics

Match marker recommendations to user's research goals and organism groups."""

    def format_response_with_suggestions(
        self, response: str, suggestions: Dict[str, Any]
    ) -> str:
        """Format response with actionable suggestions."""

        formatted_response = response

        if suggestions:
            formatted_response += "\n\n**💡 My Suggestions:**\n"

            if "search_term" in suggestions:
                formatted_response += (
                    f"🔍 **Search Term**: `{suggestions['search_term']}`\n"
                )

            if "databases" in suggestions:
                db_list = ", ".join(suggestions["databases"])
                formatted_response += f"🗄️ **Databases**: {db_list}\n"

            if "filters" in suggestions:
                formatted_response += "⚙️ **Filters**:\n"
                for key, value in suggestions["filters"].items():
                    formatted_response += f"   • {key}: {value}\n"

            if "explanation" in suggestions:
                formatted_response += f"📝 **Why**: {suggestions['explanation']}\n"

        return formatted_response

    def get_welcome_message(self, page_name: str) -> str:
        """Get welcome message for specific page."""

        messages = {
            "Sequence_Search": """👋 **Welcome to Sequence Search!**

I'm your bioinformatics assistant. I can help you:
• Craft effective search queries
• Choose the right databases  
• Set up optimal filters
• Understand molecular markers

Try asking: *"Help me find COI sequences for salmon"*""",
            "Multisequences_Alignment": """👋 **Welcome to Sequence Alignment!**

I can assist you with:
• Preparing sequences for alignment
• Choosing alignment parameters
• Interpreting results
• Generating consensus sequences

Try asking: *"How do I align my sequences?"*""",
            "default": """👋 **Hello! I'm your research assistant.**

I'm here to help with your bioinformatics workflows.
Ask me about search strategies, databases, or analysis methods!""",
        }

        return messages.get(page_name, messages["default"])
