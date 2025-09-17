# Agentic Research Assistant

## Overview

The Agentic Research Assistant is a site-wide AI-powered assistant that helps users with bioinformatics workflows, search term crafting, database selection, and analysis guidance. It provides context-aware assistance tailored to each page and user's specific needs.

## Features

### 🤖 **Core Capabilities**
- **Search Term Assistance**: Helps craft effective database queries
- **Database Selection**: Recommends appropriate databases based on research goals
- **Filter Optimization**: Suggests optimal filtering parameters
- **Molecular Marker Guidance**: Provides expert advice on marker selection
- **Troubleshooting**: Helps solve common search and analysis issues

### 🎯 **Context Awareness**
- **Page-Specific Help**: Tailored assistance based on current page
- **Session State Integration**: Understands current user state and preferences
- **Search History**: Learns from previous searches to provide better suggestions
- **Workflow Continuity**: Maintains context across different analysis steps

### 💬 **Multiple Interface Options**
- **Bubble Assistant**: 🎈 Floating bubble at bottom-right corner with expandable chat window
- **Sidebar Assistant**: Quick help and suggestions in sidebar
- **Main Chat Interface**: Full conversational interface
- **Floating Widget**: Compact assistant overlay
- **Suggestion Applicator**: Direct application of AI suggestions to forms

### 🎈 **NEW: Bubble Assistant Interface**
- **Floating Bubble**: Always-present bubble at bottom-right corner
- **Expandable Window**: 350x500px conversation window with smooth animations
- **Persistent History**: Chat history saved across sessions with timestamps
- **Unread Notifications**: Badge showing unread message count
- **Quick Suggestions**: Smart suggestion chips for new conversations
- **Auto-Save**: Conversations automatically saved every 10 messages
- **Responsive Design**: Hover effects and smooth transitions

## Architecture

```
src/app/assistant/
├── __init__.py              # Module exports
├── core.py                  # Main assistant logic & OpenAI integration
├── context_manager.py       # Page context and state management
├── prompts.py              # AI prompt templates and formatting
├── ui_components.py        # Traditional Streamlit UI components
├── popup_chat.py           # Popup chat interface (working implementation)
├── conversation_manager.py # Persistent conversation tracking
└── README.md              # This documentation
```

### Core Components

#### **AgenticAssistant (core.py)**
- Main assistant class with OpenAI integration
- Handles conversation management and response generation
- Provides fallback responses when AI is unavailable
- Supports function calling for structured suggestions

#### **PageContextManager (context_manager.py)**
- Manages page-specific context and state
- Extracts relevant session state information
- Provides page-specific help topics and suggestions
- Maintains workflow continuity

#### **PromptTemplates (prompts.py)**
- System prompts for different assistance modes
- Context-specific prompt generation
- Response formatting and suggestion structuring
- Welcome messages and help content

#### **UI Components (ui_components.py)**
- `AssistantChat`: Main chat interface
- `AssistantSidebar`: Sidebar integration
- `SuggestionApplicator`: Applies AI suggestions to forms
- `AssistantIntegration`: Main integration class

## Configuration

### OpenAI API Setup

1. **Streamlit Secrets** (Recommended):
```toml
# ~/.streamlit/secrets.toml
OPENAI_API_KEY = "your_openai_api_key_here"
```

2. **Environment Variable**:
```bash
export OPENAI_API_KEY="your_openai_api_key_here"
```

### Assistant Configuration

The assistant uses GPT-4o-mini by default for cost-effectiveness. Configuration can be adjusted in `core.py`:

```python
self.model = "gpt-4o-mini"  # Model selection
self.max_tokens = 1000      # Response length
self.temperature = 0.7      # Response creativity
```

## Integration

### Basic Integration

Add to any Streamlit page:

```python
from app.assistant.ui_components import render_assistant
from app.assistant import render_assistant_bubble

# Traditional interfaces
render_assistant("Page_Name", "sidebar")  # Sidebar
render_assistant("Page_Name", "main")     # Main chat
render_assistant("Page_Name", "floating") # Floating widget

# NEW: Bubble assistant (recommended)
render_assistant_bubble("Page_Name")
```

### Easy Popup Chat Integration

```python
# Simple integration
from app.assistant.popup_chat import render_expander_popup_chat
render_expander_popup_chat("My Page Name")

# Or use the module import
from app.assistant import render_expander_popup_chat
render_expander_popup_chat("My Page Name")
```

### Advanced Integration

For pages with forms and interactive elements:

```python
from app.assistant.ui_components import AssistantIntegration

# Initialize
assistant = AssistantIntegration()

# Check for and apply suggestions
assistant.check_and_apply_suggestions("Page_Name")

# Render assistant
assistant.render_page_assistant("Page_Name", "sidebar")
```

## Usage Examples

### Sequence Search Assistance

**User**: "Help me find COI sequences for salmon"

**Assistant**: 
```
🔍 I can help you craft an effective search for salmon COI sequences!

💡 My Suggestions:
🔍 Search Term: `Salmo salar AND COI[gene]`
🗄️ Databases: NCBI Nucleotide, BOLD
⚙️ Filters:
   • Min Length: 600
   • Max Length: 700
   • Exclude Partial: true

📝 Why: COI is the standard animal barcode, ~650bp. BOLD specializes in COI sequences with rich metadata, while NCBI provides comprehensive coverage.
```

### Database Selection Guidance

**User**: "What databases should I use for bacterial identification?"

**Assistant**:
```
🗄️ For bacterial identification, I recommend:

**Primary**: NCBI Nucleotide + SILVA
- SILVA provides curated, high-quality 16S rRNA sequences
- NCBI offers comprehensive coverage and recent submissions

**Marker**: 16S rRNA (bacterial universal marker)
**Length**: ~1500 bp for full-length sequences

Would you like help crafting a specific search term?
```

## Fallback Mode

When OpenAI API is unavailable, the assistant provides intelligent fallback responses based on pattern matching:

- Search term queries → Search syntax help
- Database questions → Database selection guide  
- Filter questions → Filtering best practices
- Marker questions → Molecular marker guide

## Page-Specific Features

### Sequence Search Page
- Search term crafting assistance
- Database selection recommendations
- Filter optimization suggestions
- Real-time search result analysis

### Alignment Page
- Sequence preparation guidance
- Parameter optimization
- Result interpretation help
- Consensus sequence advice

### Analysis Pages
- Workflow guidance
- Parameter recommendations
- Result interpretation
- Best practices

## Best Practices

### For Users
1. **Be Specific**: "Help me find 16S sequences for gut bacteria" vs "help with search"
2. **Provide Context**: Mention your research goals and organism types
3. **Use Suggestions**: Apply AI recommendations when appropriate
4. **Ask Follow-ups**: Refine suggestions based on your specific needs

### For Developers
1. **Context Integration**: Always provide rich page context
2. **Graceful Degradation**: Ensure fallback responses work well
3. **Session State**: Update context as user interacts with forms
4. **Error Handling**: Handle API failures gracefully

## Troubleshooting

### Common Issues

**Assistant not responding**:
- Check OpenAI API key configuration
- Verify internet connectivity
- Check for API quota limits

**Suggestions not applying**:
- Ensure page integration is properly configured
- Check session state management
- Verify form field mapping

**Context not working**:
- Confirm page name is correctly set
- Check session state keys
- Verify context manager integration

### Debug Mode

Enable debug mode in the UI to see:
- Current page context
- Session state information
- API request/response details
- Error messages and warnings

## Future Enhancements

- **Multi-language Support**: Support for different languages
- **Custom Models**: Integration with domain-specific models
- **Advanced Function Calling**: More sophisticated form interactions
- **Learning Capabilities**: Personalized assistance based on user history
- **Voice Interface**: Speech-to-text and text-to-speech capabilities
- **Workflow Automation**: Automated execution of common workflows

## Contributing

When extending the assistant:

1. **Add New Contexts**: Update `PageContextManager` for new pages
2. **Enhance Prompts**: Add specialized prompts in `PromptTemplates`
3. **UI Components**: Create new UI components for specific use cases
4. **Function Calling**: Add new function definitions for structured responses
5. **Testing**: Test both AI and fallback modes thoroughly

## License

This assistant module is part of the neglected-diagnostics project and follows the same license terms.
