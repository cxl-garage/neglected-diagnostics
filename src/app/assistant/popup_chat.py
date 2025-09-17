"""
Popup Chat Interface

This module provides a popup chatbot using native Streamlit components
for reliable chat functionality without HTML overlay issues.
"""

from datetime import datetime

try:
    import streamlit as st
except ImportError:
    # Mock for testing
    class MockStreamlit:
        @staticmethod
        def markdown(text, **kwargs):
            pass

        @staticmethod
        def button(text, **kwargs):
            return False

        @staticmethod
        def empty():
            return MockContext()

        @staticmethod
        def columns(cols):
            return [MockContext() for _ in range(cols)]

        @staticmethod
        def chat_message(role):
            return MockContext()

        @staticmethod
        def chat_input(prompt):
            return ""

        @staticmethod
        def container():
            return MockContext()

        @staticmethod
        def expander(title, **kwargs):
            return MockContext()

        @staticmethod
        def rerun():
            pass

        experimental_rerun = rerun

        class session_state:
            _state = {}

            @classmethod
            def get(cls, key, default=None):
                return cls._state.get(key, default)

            @classmethod
            def __contains__(cls, key):
                return key in cls._state

            @classmethod
            def __getitem__(cls, key):
                return cls._state[key]

            @classmethod
            def __setitem__(cls, key, value):
                cls._state[key] = value

        session_state = session_state()

    class MockContext:
        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def markdown(self, text, **kwargs):
            pass

        def button(self, text, **kwargs):
            return False

        def container(self):
            return self

    st = MockStreamlit()

from utils.log import _init_logger

from .context_manager import PageContextManager
from .conversation_manager import get_conversation_manager
from .core import AgenticAssistant

logger = _init_logger(__name__)


# Expander-based popup chat (working implementation)
class ExpanderPopupChat:
    """Popup chat using Streamlit expander (guaranteed to work)."""

    def __init__(self):
        self.assistant = AgenticAssistant()
        self.context_manager = PageContextManager()
        self.conversation_manager = get_conversation_manager()
        self._initialize_session_state()

    def _initialize_session_state(self):
        if "chat_messages" not in st.session_state:
            st.session_state.chat_messages = [
                {
                    "role": "assistant",
                    "content": "Hi! I'm your research assistant. I can help you with sequence searching, database selection, filtering criteria, and research workflows. Click any suggestion below or ask me anything!",
                }
            ]
        if "show_suggestions" not in st.session_state:
            st.session_state.show_suggestions = True
        if "suggestion_generation" not in st.session_state:
            st.session_state.suggestion_generation = 0  # Track suggestion refreshes
        if "show_chat" not in st.session_state:
            st.session_state.show_chat = False
        if "chat_conversation_id" not in st.session_state:
            st.session_state.chat_conversation_id = (
                f"chat_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
            )

    def render_expander_chat(self, page_name: str = "Unknown"):
        """Render chat using expander (always works)."""

        # Add floating bubble style
        st.markdown(
            """
        <style>
        .chat-expander {
            position: fixed;
            right: 24px;
            bottom: 24px;
            z-index: 1000;
            width: 400px;
            background: white;
            border-radius: 16px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.3);
        }
        </style>
        """,
            unsafe_allow_html=True,
        )

        # Floating button
        col1, col2, col3 = st.columns([8, 1, 1])
        with col3:
            if st.button("🤖", key="expander_bubble", help="Research Assistant"):
                st.session_state.show_chat = not st.session_state.get(
                    "show_chat", False
                )
                st.rerun()

        # Show chat in expander
        if st.session_state.get("show_chat", False):
            with st.expander("🤖 Research Assistant", expanded=True):
                # Close button
                col1, col2, col3 = st.columns([1, 6, 1])
                with col3:
                    if st.button("✕", key="close_chat", help="Close chat"):
                        st.session_state.show_chat = False
                        st.rerun()

                # Display messages
                for message in st.session_state.chat_messages:
                    with st.chat_message(message["role"]):
                        st.markdown(message["content"])

                # Show suggested prompts if it's the first interaction or user requested them
                if (
                    len(st.session_state.chat_messages) <= 1
                    and st.session_state.get("show_suggestions", True)
                ) or st.session_state.get("force_show_suggestions", False):
                    self._render_adaptive_suggested_prompts(page_name)
                    if st.session_state.get("force_show_suggestions", False):
                        st.session_state.force_show_suggestions = False

                # Show "Get suggestions" button for ongoing conversations
                elif len(st.session_state.chat_messages) > 1:
                    col1, col2, col3 = st.columns([1, 2, 1])
                    with col2:
                        suggestions_btn_key = f"show_suggestions_btn_{st.session_state.suggestion_generation}"
                        if st.button(
                            "💡 Get smart suggestions",
                            key=suggestions_btn_key,
                            help="Get contextual suggestions based on our conversation",
                        ):
                            try:
                                st.session_state.force_show_suggestions = True
                                st.session_state.suggestion_generation += (
                                    1  # Refresh suggestions
                                )
                                logger.info(
                                    f"Smart suggestions button clicked, generation: {st.session_state.suggestion_generation}"
                                )
                                st.rerun()
                            except Exception as e:
                                logger.error(f"Error showing suggestions: {e}")
                                st.error(f"Error loading suggestions: {e}")

                # Debug info (remove in production)
                if st.session_state.get("force_show_suggestions", False):
                    st.info(
                        f"Debug: force_show_suggestions = {st.session_state.force_show_suggestions}, generation = {st.session_state.suggestion_generation}"
                    )

                # Chat input
                if user_input := st.chat_input(
                    "Ask me about sequence analysis, databases, or research workflows...",
                    key="expander_input",
                ):
                    st.session_state.show_suggestions = (
                        False  # Hide suggestions after first message
                    )
                    self._handle_message(user_input, page_name)
                    st.rerun()

    def _render_adaptive_suggested_prompts(self, page_name: str):
        """Render adaptive suggested prompts based on conversation context."""

        # Get context-aware suggestions
        suggestions = self._get_adaptive_suggestions(page_name)

        # Debug: Log suggestions being rendered
        logger.info(
            f"Rendering {len(suggestions)} adaptive suggestions for page: {page_name}"
        )

        # Add custom CSS for suggestion buttons
        st.markdown(
            """
        <style>
        .suggestion-container {
            background: #f8f9fa;
            border-radius: 8px;
            padding: 12px;
            margin: 8px 0;
            border-left: 3px solid #667eea;
        }
        .adaptive-badge {
            background: #667eea;
            color: white;
            padding: 2px 8px;
            border-radius: 12px;
            font-size: 0.8em;
            margin-left: 8px;
        }
        </style>
        """,
            unsafe_allow_html=True,
        )

        st.markdown('<div class="suggestion-container">', unsafe_allow_html=True)

        # Show different headers based on conversation stage
        conversation_length = len(st.session_state.chat_messages)
        if conversation_length <= 1:
            st.markdown("**💡 Try asking:**")
        else:
            st.markdown(
                "**🧠 Smart suggestions based on our conversation:** <span class='adaptive-badge'>ADAPTIVE</span>",
                unsafe_allow_html=True,
            )

        # Create columns for better layout
        cols = st.columns(2)

        # Use generation number to make keys unique and allow multiple clicks
        generation = st.session_state.get("suggestion_generation", 0)

        for i, suggestion in enumerate(suggestions[:6]):  # Show up to 6 suggestions
            col_idx = i % 2
            with cols[col_idx]:
                button_key = f"adaptive_suggestion_{generation}_{i}"
                if st.button(
                    f"💬 {suggestion['text']}",
                    key=button_key,
                    help=suggestion.get("help", ""),
                    use_container_width=True,
                ):
                    try:
                        # Debug info (remove in production)
                        logger.info(
                            f"Adaptive suggestion clicked: {suggestion['text']}"
                        )

                        # Hide suggestions and process message
                        st.session_state.show_suggestions = False
                        self._handle_message(suggestion["prompt"], page_name)
                        st.rerun()
                    except Exception as e:
                        logger.error(f"Error processing suggestion click: {e}")
                        st.error(f"Error processing suggestion: {e}")

        st.markdown("</div>", unsafe_allow_html=True)

    def _get_page_suggestions(self, page_name: str):
        """Get page-specific suggestion prompts."""

        # General suggestions that work on any page
        general_suggestions = [
            {
                "text": "Help me search for sequences",
                "prompt": "I want to search for DNA sequences. Can you guide me through the process and explain what search terms work best?",
                "help": "Get guidance on sequence searching",
            },
            {
                "text": "Explain database differences",
                "prompt": "What are the differences between NCBI, BOLD, SILVA, and UNITE databases? Which should I use for my research?",
                "help": "Learn about available databases",
            },
            {
                "text": "Filter sequences by quality",
                "prompt": "How do I filter sequences by length, quality, and other criteria? What are good filtering parameters?",
                "help": "Learn about sequence filtering",
            },
            {
                "text": "Research workflow tips",
                "prompt": "Can you suggest a complete workflow for analyzing genetic sequences from search to analysis?",
                "help": "Get workflow recommendations",
            },
        ]

        # Page-specific suggestions
        page_specific = {
            "Sequence_Search": [
                {
                    "text": "Search term examples",
                    "prompt": "Can you give me examples of effective search terms for different types of genetic research (taxonomy, genes, proteins)?",
                    "help": "Get search term examples",
                },
                {
                    "text": "Batch download tips",
                    "prompt": "How do I efficiently download multiple sequences? What are the best practices for batch processing?",
                    "help": "Learn about batch operations",
                },
            ],
            "Home": [
                {
                    "text": "Getting started guide",
                    "prompt": "I'm new to sequence analysis. Can you give me a beginner's guide to using this platform?",
                    "help": "Get started with the platform",
                },
                {
                    "text": "Platform features",
                    "prompt": "What are the main features of this neglected diagnostics platform? How can it help my research?",
                    "help": "Learn about platform capabilities",
                },
            ],
        }

        # Combine general and page-specific suggestions
        suggestions = general_suggestions.copy()
        if page_name in page_specific:
            suggestions.extend(page_specific[page_name])

        return suggestions

    def _get_adaptive_suggestions(self, page_name: str):
        """Get adaptive suggestions based on conversation context and user interests."""

        conversation_length = len(st.session_state.chat_messages)

        # For first interaction, use initial suggestions
        if conversation_length <= 1:
            return self._get_initial_suggestions(page_name)

        # Analyze conversation context for adaptive suggestions
        return self._analyze_and_suggest(page_name)

    def _get_initial_suggestions(self, page_name: str):
        """Get initial suggestions for new conversations."""

        base_suggestions = [
            {
                "text": "Help me search sequences",
                "prompt": "I want to search for DNA sequences. Can you guide me through the process and explain what search terms work best?",
                "help": "Get guidance on sequence searching",
            },
            {
                "text": "Explain databases",
                "prompt": "What are the differences between NCBI, BOLD, SILVA, and UNITE databases? Which should I use for my research?",
                "help": "Learn about available databases",
            },
            {
                "text": "Filter by quality",
                "prompt": "How do I filter sequences by length, quality, and other criteria? What are good filtering parameters?",
                "help": "Learn about sequence filtering",
            },
            {
                "text": "Research workflow",
                "prompt": "Can you suggest a complete workflow for analyzing genetic sequences from search to analysis?",
                "help": "Get workflow recommendations",
            },
        ]

        # Add page-specific suggestions
        page_specific = {
            "Sequence_Search": [
                {
                    "text": "Search term examples",
                    "prompt": "Can you give me examples of effective search terms for different types of genetic research?",
                    "help": "Get search term examples",
                },
                {
                    "text": "Batch operations",
                    "prompt": "How do I efficiently download multiple sequences? What are the best practices?",
                    "help": "Learn about batch operations",
                },
            ],
            "Home": [
                {
                    "text": "Getting started",
                    "prompt": "I'm new to sequence analysis. Can you give me a beginner's guide?",
                    "help": "Get started with the platform",
                },
                {
                    "text": "Platform features",
                    "prompt": "What are the main features of this platform? How can it help my research?",
                    "help": "Learn about platform capabilities",
                },
            ],
        }

        if page_name in page_specific:
            base_suggestions.extend(page_specific[page_name])

        return base_suggestions

    def _analyze_and_suggest(self, page_name: str):
        """Analyze conversation and provide contextually relevant suggestions."""

        # Get recent messages for context analysis
        recent_messages = st.session_state.chat_messages[-6:]  # Last 6 messages
        user_messages = [msg for msg in recent_messages if msg["role"] == "user"]

        # Extract topics and interests from user messages
        conversation_text = " ".join([msg["content"].lower() for msg in user_messages])

        # Define topic-based adaptive suggestions
        adaptive_suggestions = []

        # Database-related follow-ups
        if any(
            db in conversation_text
            for db in ["ncbi", "bold", "silva", "unite", "database"]
        ):
            adaptive_suggestions.extend(
                [
                    {
                        "text": "Compare database results",
                        "prompt": "How do results differ between databases? Can you help me understand which database gives the most relevant results for my research?",
                        "help": "Compare database effectiveness",
                    },
                    {
                        "text": "Database-specific tips",
                        "prompt": "What are the best practices for searching in the specific database I mentioned? Any advanced search features?",
                        "help": "Get database-specific guidance",
                    },
                ]
            )

        # Search-related follow-ups
        if any(
            term in conversation_text for term in ["search", "query", "term", "find"]
        ):
            adaptive_suggestions.extend(
                [
                    {
                        "text": "Refine my search",
                        "prompt": "Based on what we discussed, can you help me refine my search strategy? What specific terms should I try?",
                        "help": "Improve search strategy",
                    },
                    {
                        "text": "Search troubleshooting",
                        "prompt": "I'm not getting good results with my searches. Can you help troubleshoot and suggest alternative approaches?",
                        "help": "Fix search issues",
                    },
                ]
            )

        # Filtering-related follow-ups
        if any(
            term in conversation_text
            for term in ["filter", "quality", "length", "criteria"]
        ):
            adaptive_suggestions.extend(
                [
                    {
                        "text": "Optimize my filters",
                        "prompt": "Based on our discussion, what specific filtering parameters would work best for my research goals?",
                        "help": "Customize filtering strategy",
                    },
                    {
                        "text": "Advanced filtering",
                        "prompt": "Are there advanced filtering techniques I should know about? How do I handle edge cases?",
                        "help": "Learn advanced filtering",
                    },
                ]
            )

        # Analysis-related follow-ups
        if any(
            term in conversation_text
            for term in ["analyze", "analysis", "workflow", "process"]
        ):
            adaptive_suggestions.extend(
                [
                    {
                        "text": "Next analysis steps",
                        "prompt": "Now that I have my sequences, what are the recommended next steps for analysis? What tools should I use?",
                        "help": "Continue analysis workflow",
                    },
                    {
                        "text": "Analysis best practices",
                        "prompt": "What are common pitfalls in sequence analysis? How can I ensure my analysis is robust?",
                        "help": "Avoid analysis mistakes",
                    },
                ]
            )

        # Technical issues follow-ups
        if any(
            term in conversation_text
            for term in ["error", "problem", "issue", "not working", "help"]
        ):
            adaptive_suggestions.extend(
                [
                    {
                        "text": "Troubleshoot together",
                        "prompt": "Let's work through this issue step by step. Can you walk me through exactly what's happening?",
                        "help": "Get detailed troubleshooting",
                    },
                    {
                        "text": "Alternative approaches",
                        "prompt": "Since we're having issues, what are some alternative approaches or workarounds I could try?",
                        "help": "Find alternative solutions",
                    },
                ]
            )

        # Research context follow-ups
        if any(
            term in conversation_text
            for term in ["research", "project", "study", "paper"]
        ):
            adaptive_suggestions.extend(
                [
                    {
                        "text": "Research design help",
                        "prompt": "Can you help me think through my research design? Am I approaching this in the most effective way?",
                        "help": "Optimize research approach",
                    },
                    {
                        "text": "Publication tips",
                        "prompt": "What should I document about my sequence analysis methods for publication? Any reporting standards?",
                        "help": "Prepare for publication",
                    },
                ]
            )

        # Beginner-friendly follow-ups
        if any(
            term in conversation_text
            for term in ["new", "beginner", "start", "learn", "don't know"]
        ):
            adaptive_suggestions.extend(
                [
                    {
                        "text": "Learn more basics",
                        "prompt": "I'd like to understand more fundamentals. Can you explain the key concepts I should know?",
                        "help": "Build foundational knowledge",
                    },
                    {
                        "text": "Practice exercises",
                        "prompt": "Do you have any practice exercises or examples I can work through to build my skills?",
                        "help": "Get hands-on practice",
                    },
                ]
            )

        # If no specific context, provide general follow-ups
        if not adaptive_suggestions:
            adaptive_suggestions = [
                {
                    "text": "Dive deeper",
                    "prompt": "Can you elaborate more on what we just discussed? I'd like to understand this topic better.",
                    "help": "Get more detailed information",
                },
                {
                    "text": "Related topics",
                    "prompt": "What related topics should I know about? What else is important in this area?",
                    "help": "Explore related concepts",
                },
                {
                    "text": "Practical examples",
                    "prompt": "Can you give me some concrete examples of how to apply what we've discussed?",
                    "help": "See practical applications",
                },
                {
                    "text": "Common questions",
                    "prompt": "What are other common questions people have about this topic? What else should I be thinking about?",
                    "help": "Discover related questions",
                },
            ]

        # Always include some general helpful options
        adaptive_suggestions.extend(
            [
                {
                    "text": "Summarize our chat",
                    "prompt": "Can you summarize the key points from our conversation and give me actionable next steps?",
                    "help": "Get conversation summary",
                },
                {
                    "text": "Start something new",
                    "prompt": "I'd like to explore a different topic. What other areas can you help me with?",
                    "help": "Explore new topics",
                },
            ]
        )

        # Return up to 6 most relevant suggestions
        return adaptive_suggestions[:6]

    def _handle_message(self, user_input: str, page_name: str):
        """Handle user message."""
        try:
            # Add user message
            st.session_state.chat_messages.append(
                {
                    "role": "user",
                    "content": user_input,
                    "timestamp": datetime.now().isoformat(),
                }
            )

            # Get assistant response
            context = self.context_manager.get_current_context()
            context["page_name"] = page_name

            # Get recent conversation history for context
            recent_messages = st.session_state.chat_messages[-6:]

            # Try to get assistant response
            try:
                response, suggestions = self.assistant.get_assistance(
                    user_input, context, recent_messages
                )
            except Exception as assistant_error:
                logger.warning(f"Assistant error, using fallback: {assistant_error}")
                # Fallback response if assistant fails
                response = f"I understand you're asking about: '{user_input[:100]}...'\n\nI'm having trouble accessing my full capabilities right now, but I can still help! Could you try rephrasing your question or ask me something specific about sequence analysis, databases, or research workflows?"
                suggestions = None

            # Add assistant response
            assistant_message = {
                "role": "assistant",
                "content": response,
                "timestamp": datetime.now().isoformat(),
            }
            if suggestions:
                assistant_message["suggestions"] = suggestions

            st.session_state.chat_messages.append(assistant_message)

            # Update conversation metadata
            self.conversation_manager.update_session_metadata(
                page_name, user_input[:50]
            )

            # Auto-save periodically
            if len(st.session_state.chat_messages) % 10 == 0:
                self.conversation_manager.auto_save_current_conversation()

        except Exception as e:
            logger.error(f"Error handling message: {e}")
            st.session_state.chat_messages.append(
                {
                    "role": "assistant",
                    "content": f"I apologize, but I encountered an error processing your request. Please try again or rephrase your question.",
                    "timestamp": datetime.now().isoformat(),
                }
            )


# Convenience function


def render_expander_popup_chat(page_name: str = "Unknown") -> None:
    """Render expander-based popup chat (guaranteed to work)."""
    try:
        popup = ExpanderPopupChat()
        popup.render_expander_chat(page_name)
    except Exception as e:
        logger.error(f"Expander popup failed: {e}")
        # Show a simple fallback chat
        st.error(f"Chat assistant error: {e}")

        # Simple fallback interface
        if st.button("🤖 Simple Chat (Fallback)", key="fallback_chat"):
            if "fallback_messages" not in st.session_state:
                st.session_state.fallback_messages = []

            with st.expander("Simple Chat", expanded=True):
                for msg in st.session_state.fallback_messages:
                    st.write(f"{msg['role']}: {msg['content']}")

                if user_input := st.text_input("Ask something:", key="fallback_input"):
                    st.session_state.fallback_messages.append(
                        {"role": "user", "content": user_input}
                    )
                    st.session_state.fallback_messages.append(
                        {"role": "assistant", "content": f"You asked: {user_input}"}
                    )
                    st.rerun()
