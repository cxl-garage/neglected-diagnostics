"""
UI Components for Agentic Assistant

This module contains Streamlit UI components for the assistant interface,
including chat interface, sidebar integration, and suggestion displays.
"""

import json
import time
from datetime import datetime
from typing import Any, Dict, List, Optional, Tuple

try:
    import streamlit as st
except ImportError:
    # Mock streamlit for non-web environments
    class MockStreamlit:
        class session_state:
            _state = {}

            def get(self, key, default=None):
                return self._state.get(key, default)

            def __contains__(self, key):
                return key in self._state

            def __getitem__(self, key):
                return self._state[key]

            def __setitem__(self, key, value):
                self._state[key] = value

        @staticmethod
        def subheader(text):
            print(f"SUBHEADER: {text}")

        @staticmethod
        def info(text):
            print(f"INFO: {text}")

        @staticmethod
        def chat_input(prompt):
            return None

        @staticmethod
        def chat_message(role):
            return MockContext()

        @staticmethod
        def write(text):
            print(f"WRITE: {text}")

        @staticmethod
        def button(text, **kwargs):
            return False

        @staticmethod
        def container():
            return MockContext()

        @staticmethod
        def columns(n):
            return [MockContext() for _ in range(n)]

        @staticmethod
        def expander(title, expanded=False):
            return MockContext()

        @staticmethod
        def text_input(label, **kwargs):
            return ""

        @staticmethod
        def success(text):
            print(f"SUCCESS: {text}")

        @staticmethod
        def rerun():
            pass

        session_state = session_state()

    class MockContext:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def write(self, text):
            print(f"CONTEXT: {text}")

        def button(self, text, **kwargs):
            return False

    st = MockStreamlit()

from utils.log import _init_logger

from .context_manager import PageContextManager
from .core import AgenticAssistant

logger = _init_logger(__name__)


class AssistantChat:
    """Main chat interface for the assistant."""

    def __init__(self):
        """Initialize the chat interface."""
        self.assistant = AgenticAssistant()
        self.context_manager = PageContextManager()

        # Initialize session state
        if "assistant_messages" not in st.session_state:
            st.session_state.assistant_messages = []
        if "assistant_suggestions" not in st.session_state:
            st.session_state.assistant_suggestions = {}

    def render_chat_interface(self, container_key: str = "main") -> None:
        """Render the main chat interface."""

        st.subheader("🤖 Research Assistant")

        # Show welcome message if no conversation yet
        if not st.session_state.assistant_messages:
            page_name = self._get_current_page_name()
            welcome_msg = self.assistant.prompts.get_welcome_message(page_name)
            st.info(welcome_msg)

        # Chat container
        chat_container = st.container()

        with chat_container:
            # Display conversation history
            for message in st.session_state.assistant_messages:
                with st.chat_message(message["role"]):
                    st.write(message["content"])

                    # Show suggestions if available
                    if message["role"] == "assistant" and "suggestions" in message:
                        self._render_suggestions(message["suggestions"])

        # Chat input
        if prompt := st.chat_input("Ask me anything about your research..."):
            self._handle_user_input(prompt)

    def render_sidebar_assistant(self) -> None:
        """Render compact assistant in sidebar."""

        with st.sidebar:
            st.markdown("---")
            st.subheader("🤖 Assistant")

            # Quick help button
            if st.button("💡 Quick Help", use_container_width=True):
                self._show_quick_help()

            # Suggested queries
            page_name = self._get_current_page_name()
            suggestions = self.context_manager.get_suggested_queries(page_name)

            if suggestions:
                st.markdown("**Quick Questions:**")
                for i, suggestion in enumerate(suggestions[:3]):  # Show top 3
                    if st.button(f"❓ {suggestion[:30]}...", key=f"suggestion_{i}"):
                        self._handle_user_input(suggestion)

            # Chat toggle
            if st.button("💬 Open Chat", use_container_width=True):
                st.session_state.show_assistant_chat = True
                st.rerun()

    def render_floating_assistant(self) -> None:
        """Render floating assistant widget."""

        # This would create a floating chat widget
        # For now, we'll use an expander as a simpler implementation
        with st.expander("🤖 Research Assistant", expanded=False):
            self.render_compact_chat()

    def render_compact_chat(self) -> None:
        """Render a compact version of the chat interface."""

        # Show last few messages
        recent_messages = (
            st.session_state.assistant_messages[-4:]
            if st.session_state.assistant_messages
            else []
        )

        for message in recent_messages:
            role_icon = "🧑" if message["role"] == "user" else "🤖"
            st.write(f"{role_icon} {message['content'][:100]}...")

        # Input field
        user_input = st.text_input("Ask assistant:", key="compact_chat_input")
        if st.button("Send", key="compact_send"):
            if user_input:
                self._handle_user_input(user_input)
                st.rerun()

    def _handle_user_input(self, user_input: str) -> None:
        """Process user input and get assistant response."""

        # Add user message to history
        st.session_state.assistant_messages.append(
            {
                "role": "user",
                "content": user_input,
                "timestamp": datetime.now().isoformat(),
            }
        )

        # Get page context
        context = self.context_manager.get_current_context()

        # Get assistant response
        try:
            response, suggestions = self.assistant.get_assistance(
                user_input,
                context,
                st.session_state.assistant_messages[
                    -10:
                ],  # Last 10 messages for context
            )

            # Add assistant response to history
            assistant_message = {
                "role": "assistant",
                "content": response,
                "timestamp": datetime.now().isoformat(),
            }

            if suggestions:
                assistant_message["suggestions"] = suggestions

            st.session_state.assistant_messages.append(assistant_message)

            # Apply suggestions if they include actionable items
            self._apply_suggestions(suggestions)

        except Exception as e:
            logger.error(f"Error getting assistant response: {e}")
            st.session_state.assistant_messages.append(
                {
                    "role": "assistant",
                    "content": "I apologize, but I encountered an error. Please try again.",
                    "timestamp": datetime.now().isoformat(),
                }
            )

    def _render_suggestions(self, suggestions: Dict[str, Any]) -> None:
        """Render actionable suggestions from assistant."""

        if not suggestions:
            return

        st.markdown("**💡 Suggestions:**")

        col1, col2 = st.columns(2)

        with col1:
            if "search_term" in suggestions:
                if st.button(
                    f"📝 Use: {suggestions['search_term']}",
                    key=f"use_search_{hash(suggestions['search_term'])}",
                ):
                    self._apply_search_term(suggestions["search_term"])

        with col2:
            if "databases" in suggestions:
                db_list = ", ".join(suggestions["databases"])
                if st.button(
                    f"🗄️ Select: {db_list[:20]}...",
                    key=f"use_db_{hash(str(suggestions['databases']))}",
                ):
                    self._apply_database_selection(suggestions["databases"])

        if "filters" in suggestions:
            st.json(suggestions["filters"])

    def _apply_suggestions(self, suggestions: Dict[str, Any]) -> None:
        """Apply suggestions to the current page state."""

        if not suggestions:
            return

        # Store suggestions for potential application
        st.session_state.assistant_suggestions = suggestions

        # Auto-apply certain suggestions if appropriate
        # This would depend on the specific page and user preferences
        logger.info(f"Stored suggestions: {suggestions}")

    def _apply_search_term(self, search_term: str) -> None:
        """Apply suggested search term to the search interface."""
        # This would interact with the search form
        # For now, we'll store it in session state
        st.session_state.suggested_search_term = search_term
        st.success(f"Applied search term: {search_term}")
        st.rerun()

    def _apply_database_selection(self, databases: List[str]) -> None:
        """Apply suggested database selection."""
        # This would interact with the database selection interface
        st.session_state.suggested_databases = databases
        st.success(f"Suggested databases: {', '.join(databases)}")
        st.rerun()

    def _show_quick_help(self) -> None:
        """Show quick help for the current page."""
        page_name = self._get_current_page_name()
        help_content = self.assistant.get_page_specific_help(page_name)
        st.info(help_content)

    def _get_current_page_name(self) -> str:
        """Get current page name for context."""
        try:
            # This would be implemented based on your page detection logic
            return "Sequence_Search"  # Default for now
        except:
            return "Unknown"


class AssistantSidebar:
    """Sidebar integration for the assistant."""

    def __init__(self):
        """Initialize sidebar assistant."""
        self.chat = AssistantChat()

    def render(self) -> None:
        """Render assistant sidebar components."""
        self.chat.render_sidebar_assistant()


class SuggestionApplicator:
    """Handles applying assistant suggestions to page interfaces."""

    def __init__(self):
        """Initialize the suggestion applicator."""
        pass

    def apply_to_search_form(self, suggestions: Dict[str, Any]) -> None:
        """Apply suggestions to sequence search form."""

        if not suggestions:
            return

        applied_any = False

        # Apply search term
        if "search_term" in suggestions:
            # This would set the search term in the form
            # Implementation depends on how forms are structured
            logger.info(f"Would apply search term: {suggestions['search_term']}")
            applied_any = True

        # Apply database selection
        if "databases" in suggestions:
            logger.info(f"Would apply databases: {suggestions['databases']}")
            applied_any = True

        # Apply filters
        if "filters" in suggestions:
            logger.info(f"Would apply filters: {suggestions['filters']}")
            applied_any = True

        if applied_any:
            st.success("✅ Applied assistant suggestions to your search!")

    def get_applicable_suggestions(self, page_name: str) -> Dict[str, Any]:
        """Get suggestions that can be applied to current page."""

        if "assistant_suggestions" not in st.session_state:
            return {}

        suggestions = st.session_state.assistant_suggestions

        # Filter suggestions based on page capabilities
        applicable = {}

        if page_name == "Sequence_Search":
            for key in ["search_term", "databases", "filters"]:
                if key in suggestions:
                    applicable[key] = suggestions[key]

        return applicable


class AssistantIntegration:
    """Main class for integrating assistant across the application."""

    def __init__(self):
        """Initialize assistant integration."""
        self.chat = AssistantChat()
        self.sidebar = AssistantSidebar()
        self.applicator = SuggestionApplicator()

    def render_page_assistant(self, page_name: str, location: str = "sidebar") -> None:
        """Render assistant for a specific page."""

        if location == "sidebar":
            self.sidebar.render()
        elif location == "main":
            self.chat.render_chat_interface()
        elif location == "floating":
            self.chat.render_floating_assistant()
        # bubble location removed - use popup_chat module instead

    def check_and_apply_suggestions(self, page_name: str) -> None:
        """Check for and apply any pending suggestions."""

        suggestions = self.applicator.get_applicable_suggestions(page_name)

        if suggestions:
            with st.expander("🤖 Assistant Suggestions Available", expanded=False):
                st.write("I have some suggestions for your search:")

                for key, value in suggestions.items():
                    st.write(f"**{key.title()}**: {value}")

                col1, col2 = st.columns(2)
                with col1:
                    if st.button("✅ Apply Suggestions"):
                        self.applicator.apply_to_search_form(suggestions)
                        # Clear applied suggestions
                        st.session_state.assistant_suggestions = {}
                        st.rerun()

                with col2:
                    if st.button("❌ Dismiss"):
                        st.session_state.assistant_suggestions = {}
                        st.rerun()


# Convenience function for easy integration
def render_assistant(page_name: str = "Unknown", location: str = "sidebar") -> None:
    """Convenience function to render assistant on any page."""

    integration = AssistantIntegration()
    integration.render_page_assistant(page_name, location)
    integration.check_and_apply_suggestions(page_name)


# Note: Bubble assistant functionality moved to popup_chat module
