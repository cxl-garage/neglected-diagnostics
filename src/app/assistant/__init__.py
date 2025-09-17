"""
Agentic Assistant Module

This module provides a site-wide intelligent assistant that helps users with:
- Crafting search terms and queries
- Selecting appropriate databases and filters
- Understanding molecular markers and their applications
- Troubleshooting search issues
- Providing context-specific guidance

The assistant is context-aware and provides tailored help based on the current page
and user's specific needs.
"""

from .context_manager import PageContextManager
from .conversation_manager import ConversationManager, get_conversation_manager
from .core import AgenticAssistant
from .popup_chat import ExpanderPopupChat, render_expander_popup_chat
from .prompts import PromptTemplates
from .ui_components import AssistantChat, AssistantSidebar, render_assistant

__all__ = [
    "AgenticAssistant",
    "AssistantChat",
    "AssistantSidebar",
    "PageContextManager",
    "PromptTemplates",
    "ExpanderPopupChat",
    "ConversationManager",
    "render_assistant",
    "render_expander_popup_chat",
    "get_conversation_manager",
]
