"""
Conversation Manager for Persistent Chat History

This module manages conversation persistence, history tracking, and 
conversation analytics for the agentic assistant.
"""

import json
import os
import pickle
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    import streamlit as st
except ImportError:
    # Mock for testing
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


class ConversationManager:
    """Manages conversation persistence and history tracking."""

    def __init__(self, storage_dir: Optional[str] = None):
        """Initialize conversation manager."""
        self.storage_dir = storage_dir or self._get_default_storage_dir()
        self._ensure_storage_dir()

        # Initialize session-based tracking
        self._initialize_session_tracking()

    def _get_default_storage_dir(self) -> str:
        """Get default storage directory for conversation history."""
        # Use a hidden directory in user's home or temp directory
        home_dir = Path.home()
        storage_path = home_dir / ".neglected_diagnostics" / "conversations"
        return str(storage_path)

    def _ensure_storage_dir(self):
        """Ensure storage directory exists."""
        try:
            Path(self.storage_dir).mkdir(parents=True, exist_ok=True)
            logger.info(f"Conversation storage directory: {self.storage_dir}")
        except Exception as e:
            logger.warning(f"Could not create storage directory: {e}")
            # Fall back to session-only storage
            self.storage_dir = None

    def _initialize_session_tracking(self):
        """Initialize session-based conversation tracking."""
        if "conversation_sessions" not in st.session_state:
            st.session_state.conversation_sessions = {}

        if "current_session_id" not in st.session_state:
            st.session_state.current_session_id = self._generate_session_id()

        if "conversation_metadata" not in st.session_state:
            st.session_state.conversation_metadata = {
                "session_start": datetime.now().isoformat(),
                "total_messages": 0,
                "pages_visited": set(),
                "topics_discussed": set(),
            }

    def _generate_session_id(self) -> str:
        """Generate unique session ID."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        return f"session_{timestamp}"

    def save_conversation(
        self, conversation_id: str, messages: List[Dict], metadata: Dict = None
    ):
        """Save conversation to persistent storage."""
        if not self.storage_dir:
            logger.warning("No storage directory available, conversation not saved")
            return False

        try:
            conversation_data = {
                "id": conversation_id,
                "messages": messages,
                "metadata": metadata or {},
                "created_at": datetime.now().isoformat(),
                "message_count": len(messages),
                "last_updated": datetime.now().isoformat(),
            }

            file_path = Path(self.storage_dir) / f"{conversation_id}.json"
            with open(file_path, "w", encoding="utf-8") as f:
                json.dump(conversation_data, f, indent=2, ensure_ascii=False)

            logger.info(
                f"Conversation saved: {conversation_id} ({len(messages)} messages)"
            )
            return True

        except Exception as e:
            logger.error(f"Failed to save conversation {conversation_id}: {e}")
            return False

    def load_conversation(self, conversation_id: str) -> Optional[Dict]:
        """Load conversation from persistent storage."""
        if not self.storage_dir:
            return None

        try:
            file_path = Path(self.storage_dir) / f"{conversation_id}.json"
            if not file_path.exists():
                return None

            with open(file_path, "r", encoding="utf-8") as f:
                conversation_data = json.load(f)

            logger.info(f"Conversation loaded: {conversation_id}")
            return conversation_data

        except Exception as e:
            logger.error(f"Failed to load conversation {conversation_id}: {e}")
            return None

    def list_conversations(self, limit: int = 50) -> List[Dict]:
        """List recent conversations with metadata."""
        if not self.storage_dir:
            return []

        try:
            conversations = []
            storage_path = Path(self.storage_dir)

            for file_path in storage_path.glob("*.json"):
                try:
                    with open(file_path, "r", encoding="utf-8") as f:
                        data = json.load(f)

                    # Extract summary info
                    summary = {
                        "id": data["id"],
                        "created_at": data["created_at"],
                        "last_updated": data.get("last_updated", data["created_at"]),
                        "message_count": data.get(
                            "message_count", len(data.get("messages", []))
                        ),
                        "preview": self._generate_conversation_preview(
                            data.get("messages", [])
                        ),
                        "topics": data.get("metadata", {}).get("topics_discussed", []),
                    }
                    conversations.append(summary)

                except Exception as e:
                    logger.warning(f"Could not read conversation file {file_path}: {e}")
                    continue

            # Sort by last updated, most recent first
            conversations.sort(key=lambda x: x["last_updated"], reverse=True)
            return conversations[:limit]

        except Exception as e:
            logger.error(f"Failed to list conversations: {e}")
            return []

    def _generate_conversation_preview(self, messages: List[Dict]) -> str:
        """Generate a preview of the conversation."""
        if not messages:
            return "Empty conversation"

        # Find first user message
        for msg in messages:
            if msg.get("role") == "user":
                content = msg.get("content", "")
                # Truncate and clean
                preview = content.replace("\n", " ").strip()
                return preview[:60] + "..." if len(preview) > 60 else preview

        return "Assistant conversation"

    def delete_conversation(self, conversation_id: str) -> bool:
        """Delete a conversation from storage."""
        if not self.storage_dir:
            return False

        try:
            file_path = Path(self.storage_dir) / f"{conversation_id}.json"
            if file_path.exists():
                file_path.unlink()
                logger.info(f"Conversation deleted: {conversation_id}")
                return True
            return False

        except Exception as e:
            logger.error(f"Failed to delete conversation {conversation_id}: {e}")
            return False

    def cleanup_old_conversations(self, days_old: int = 30):
        """Clean up conversations older than specified days."""
        if not self.storage_dir:
            return

        try:
            cutoff_date = datetime.now() - timedelta(days=days_old)
            deleted_count = 0

            storage_path = Path(self.storage_dir)
            for file_path in storage_path.glob("*.json"):
                try:
                    with open(file_path, "r", encoding="utf-8") as f:
                        data = json.load(f)

                    created_at = datetime.fromisoformat(data["created_at"])
                    if created_at < cutoff_date:
                        file_path.unlink()
                        deleted_count += 1

                except Exception as e:
                    logger.warning(
                        f"Could not process file {file_path} for cleanup: {e}"
                    )
                    continue

            if deleted_count > 0:
                logger.info(f"Cleaned up {deleted_count} old conversations")

        except Exception as e:
            logger.error(f"Failed to cleanup old conversations: {e}")

    def get_conversation_analytics(self) -> Dict[str, Any]:
        """Get analytics about conversation history."""
        analytics = {
            "total_conversations": 0,
            "total_messages": 0,
            "average_messages_per_conversation": 0,
            "most_active_day": None,
            "common_topics": [],
            "conversation_lengths": [],
            "recent_activity": [],
        }

        if not self.storage_dir:
            return analytics

        try:
            conversations = []
            daily_counts = {}
            all_topics = []

            storage_path = Path(self.storage_dir)
            for file_path in storage_path.glob("*.json"):
                try:
                    with open(file_path, "r", encoding="utf-8") as f:
                        data = json.load(f)

                    conversations.append(data)

                    # Count by day
                    created_date = datetime.fromisoformat(data["created_at"]).date()
                    daily_counts[created_date] = daily_counts.get(created_date, 0) + 1

                    # Collect topics
                    topics = data.get("metadata", {}).get("topics_discussed", [])
                    all_topics.extend(topics)

                except Exception as e:
                    continue

            # Calculate analytics
            analytics["total_conversations"] = len(conversations)
            analytics["total_messages"] = sum(
                len(c.get("messages", [])) for c in conversations
            )

            if conversations:
                analytics["average_messages_per_conversation"] = analytics[
                    "total_messages"
                ] / len(conversations)

            # Most active day
            if daily_counts:
                most_active_date = max(daily_counts, key=daily_counts.get)
                analytics["most_active_day"] = {
                    "date": most_active_date.isoformat(),
                    "conversations": daily_counts[most_active_date],
                }

            # Common topics (simplified - would need NLP for better topic extraction)
            from collections import Counter

            topic_counts = Counter(all_topics)
            analytics["common_topics"] = topic_counts.most_common(5)

            # Conversation lengths
            lengths = [len(c.get("messages", [])) for c in conversations]
            analytics["conversation_lengths"] = {
                "min": min(lengths) if lengths else 0,
                "max": max(lengths) if lengths else 0,
                "average": sum(lengths) / len(lengths) if lengths else 0,
            }

            # Recent activity (last 7 days)
            recent_date = datetime.now().date() - timedelta(days=7)
            recent_conversations = [
                c
                for c in conversations
                if datetime.fromisoformat(c["created_at"]).date() >= recent_date
            ]
            analytics["recent_activity"] = len(recent_conversations)

        except Exception as e:
            logger.error(f"Failed to generate conversation analytics: {e}")

        return analytics

    def update_session_metadata(self, page_name: str, topic: str = None):
        """Update current session metadata."""
        if "conversation_metadata" in st.session_state:
            metadata = st.session_state.conversation_metadata
            metadata["pages_visited"].add(page_name)

            if topic:
                metadata["topics_discussed"].add(topic)

            st.session_state.conversation_metadata = metadata

    def get_current_session_summary(self) -> Dict[str, Any]:
        """Get summary of current session."""
        if "conversation_metadata" not in st.session_state:
            return {}

        metadata = st.session_state.conversation_metadata
        session_duration = datetime.now() - datetime.fromisoformat(
            metadata["session_start"]
        )

        return {
            "session_id": st.session_state.get("current_session_id", "unknown"),
            "duration_minutes": int(session_duration.total_seconds() / 60),
            "total_messages": len(st.session_state.get("bubble_chat_history", [])),
            "pages_visited": list(metadata.get("pages_visited", [])),
            "topics_discussed": list(metadata.get("topics_discussed", [])),
            "session_start": metadata["session_start"],
        }

    def export_conversation_data(self, format: str = "json") -> str:
        """Export conversation data in specified format."""
        try:
            conversations = self.list_conversations(limit=1000)  # Export all

            if format.lower() == "json":
                return json.dumps(conversations, indent=2, ensure_ascii=False)
            elif format.lower() == "csv":
                # Simple CSV export
                import csv
                import io

                output = io.StringIO()
                writer = csv.writer(output)

                # Header
                writer.writerow(["ID", "Created", "Messages", "Preview", "Topics"])

                # Data
                for conv in conversations:
                    writer.writerow(
                        [
                            conv["id"],
                            conv["created_at"],
                            conv["message_count"],
                            conv["preview"],
                            ", ".join(conv.get("topics", [])),
                        ]
                    )

                return output.getvalue()
            else:
                return json.dumps(conversations, indent=2)

        except Exception as e:
            logger.error(f"Failed to export conversation data: {e}")
            return ""

    def auto_save_current_conversation(self):
        """Auto-save current bubble conversation."""
        if (
            "bubble_chat_history" in st.session_state
            and st.session_state.bubble_chat_history
        ):
            session_id = st.session_state.get("current_session_id", "unknown")
            messages = st.session_state.bubble_chat_history
            metadata = st.session_state.get("conversation_metadata", {})

            # Only save if there are meaningful messages
            if len(messages) >= 2:  # At least one exchange
                self.save_conversation(session_id, messages, metadata)


# Global conversation manager instance
_conversation_manager = None


def get_conversation_manager() -> ConversationManager:
    """Get global conversation manager instance."""
    global _conversation_manager
    if _conversation_manager is None:
        _conversation_manager = ConversationManager()
    return _conversation_manager
