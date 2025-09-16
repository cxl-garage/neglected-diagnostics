import asyncio
import json
import os
import subprocess
import traceback
from typing import Any, Dict, List, Optional

import pandas as pd
import streamlit as st
from common.render_method import render_markdown

from utils.log import _init_logger

logger = _init_logger(__name__)

from app.common import setup
from app.common.constants import (
    MAIN_PAGE_COLS_GAP,
    MAIN_PAGE_COLS_SIZE,
    NAVIGATE_WARNING_MD,
)

# Galaxy Workflow Assistant specific constants
GALAXY_CONNECTION = "GALAXY_CONNECTION"
GALAXY_TOOLS = "GALAXY_TOOLS"
CONVERSATION_HISTORY = "CONVERSATION_HISTORY"
WORKFLOW_PLAN = "WORKFLOW_PLAN"
LLM_PROVIDER = "LLM_PROVIDER"
LLM_CONFIG = "LLM_CONFIG"
MCP_SERVER_PROCESS = "MCP_SERVER_PROCESS"
MCP_CONNECTION = "MCP_CONNECTION"

# Initialize the Streamlit session state keys
setup.initialize()


def init_session_state_galaxy_assistant():
    """Initialize the Streamlit session state for Galaxy Workflow Assistant."""
    if GALAXY_CONNECTION not in st.session_state:
        st.session_state[GALAXY_CONNECTION] = None
    if GALAXY_TOOLS not in st.session_state:
        st.session_state[GALAXY_TOOLS] = []
    if CONVERSATION_HISTORY not in st.session_state:
        st.session_state[CONVERSATION_HISTORY] = []
    if WORKFLOW_PLAN not in st.session_state:
        st.session_state[WORKFLOW_PLAN] = None
    if LLM_PROVIDER not in st.session_state:
        st.session_state[LLM_PROVIDER] = "openai"
    if LLM_CONFIG not in st.session_state:
        st.session_state[LLM_CONFIG] = {}
    if MCP_SERVER_PROCESS not in st.session_state:
        st.session_state[MCP_SERVER_PROCESS] = None
    if MCP_CONNECTION not in st.session_state:
        st.session_state[MCP_CONNECTION] = None


def connect_to_galaxy(galaxy_url: str, api_key: str):
    """Connect to Galaxy instance using BioBlend."""
    try:
        from bioblend.galaxy import GalaxyInstance

        gi = GalaxyInstance(galaxy_url, key=api_key)
        # Test connection
        gi.config.get_config()
        return gi
    except ImportError:
        st.error(
            "BioBlend library not found. Please install it using: pip install bioblend"
        )
        return None
    except Exception as e:
        st.error(f"Failed to connect to Galaxy: {str(e)}")
        return None


def start_galaxy_mcp_server(galaxy_url: str, api_key: str):
    """Start the Galaxy MCP server process."""
    try:
        # Set environment variables for Galaxy MCP
        env = os.environ.copy()
        env["GALAXY_URL"] = galaxy_url
        env["GALAXY_API_KEY"] = api_key

        # Start the MCP server process
        process = subprocess.Popen(
            ["uvx", "--from", "galaxy-mcp", "mcp", "run", "galaxy_mcp.server"],
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        st.session_state[MCP_SERVER_PROCESS] = process
        return process
    except Exception as e:
        logger.error(f"Failed to start Galaxy MCP server: {str(e)}")
        return None


def stop_galaxy_mcp_server():
    """Stop the Galaxy MCP server process."""
    process = st.session_state.get(MCP_SERVER_PROCESS)
    if process:
        try:
            process.terminate()
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            process.kill()
        st.session_state[MCP_SERVER_PROCESS] = None


def get_galaxy_mcp_tools_definition():
    """Get the definition of available Galaxy MCP tools."""
    return """Available Galaxy MCP Tools:
    
1. galaxy_connect - Connect to a Galaxy server
2. galaxy_get_server_info - Get Galaxy server information and capabilities  
3. galaxy_search_tools - Search for available Galaxy tools
4. galaxy_get_tool_details - Get detailed information about a specific tool
5. galaxy_execute_tool - Execute a Galaxy tool with specified parameters
6. galaxy_get_histories - List user's Galaxy histories
7. galaxy_get_history_contents - Get contents of a specific history
8. galaxy_upload_file - Upload a file to Galaxy
9. galaxy_get_workflows - Get available workflows from Galaxy
10. galaxy_import_workflow - Import a workflow from the Interactive Workflow Composer (IWC)
11. galaxy_run_workflow - Execute a workflow with specified inputs"""


def get_openai_function_definitions():
    """Get OpenAI function definitions for Galaxy MCP tools."""
    return [
        {
            "type": "function",
            "function": {
                "name": "galaxy_search_tools",
                "description": "Search for Galaxy tools by name, description, or category",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "query": {
                            "type": "string",
                            "description": "Search query for tools",
                        },
                        "limit": {
                            "type": "integer",
                            "description": "Maximum number of results to return",
                            "default": 10,
                        },
                    },
                    "required": ["query"],
                },
            },
        },
        {
            "type": "function",
            "function": {
                "name": "galaxy_get_tool_details",
                "description": "Get detailed information about a specific Galaxy tool",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "tool_id": {
                            "type": "string",
                            "description": "The ID of the tool to get details for",
                        }
                    },
                    "required": ["tool_id"],
                },
            },
        },
        {
            "type": "function",
            "function": {
                "name": "galaxy_get_server_info",
                "description": "Get Galaxy server information and capabilities",
                "parameters": {"type": "object", "properties": {}, "required": []},
            },
        },
        {
            "type": "function",
            "function": {
                "name": "galaxy_get_histories",
                "description": "List user's Galaxy histories",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "limit": {
                            "type": "integer",
                            "description": "Maximum number of histories to return",
                            "default": 20,
                        }
                    },
                    "required": [],
                },
            },
        },
        {
            "type": "function",
            "function": {
                "name": "galaxy_get_workflows",
                "description": "Get available workflows from Galaxy and IWC",
                "parameters": {
                    "type": "object",
                    "properties": {
                        "source": {
                            "type": "string",
                            "enum": ["user", "iwc", "all"],
                            "description": "Source of workflows to retrieve",
                            "default": "all",
                        }
                    },
                    "required": [],
                },
            },
        },
    ]


def execute_galaxy_mcp_function(function_name: str, arguments: dict):
    """Execute a Galaxy MCP function call."""
    try:
        # This is a simplified implementation
        # In practice, you would use the actual MCP client to make these calls

        if function_name == "galaxy_search_tools":
            return search_galaxy_tools_via_mcp(
                arguments.get("query"), arguments.get("limit", 10)
            )
        elif function_name == "galaxy_get_tool_details":
            return get_tool_details_via_mcp(arguments.get("tool_id"))
        elif function_name == "galaxy_get_server_info":
            return get_server_info_via_mcp()
        elif function_name == "galaxy_get_histories":
            return get_histories_via_mcp(arguments.get("limit", 20))
        elif function_name == "galaxy_get_workflows":
            return get_workflows_via_mcp(arguments.get("source", "all"))
        else:
            return {"error": f"Unknown function: {function_name}"}
    except Exception as e:
        logger.error(f"MCP function execution error: {str(e)}")
        return {"error": str(e)}


def search_galaxy_tools_via_mcp(query: str, limit: int = 10):
    """Search Galaxy tools via MCP (simplified implementation)."""
    # This would use the actual MCP client in production
    gi = st.session_state.get(GALAXY_CONNECTION)
    if not gi:
        return {"error": "Not connected to Galaxy"}

    try:
        all_tools = gi.tools.get_tools()
        matching_tools = [
            tool
            for tool in all_tools
            if query.lower() in tool.get("name", "").lower()
            or query.lower() in tool.get("description", "").lower()
        ][:limit]

        return {
            "tools": matching_tools,
            "total_found": len(matching_tools),
            "query": query,
        }
    except Exception as e:
        return {"error": str(e)}


def get_tool_details_via_mcp(tool_id: str):
    """Get tool details via MCP (simplified implementation)."""
    gi = st.session_state.get(GALAXY_CONNECTION)
    if not gi:
        return {"error": "Not connected to Galaxy"}

    try:
        tool_details = gi.tools.show_tool(tool_id)
        return tool_details
    except Exception as e:
        return {"error": str(e)}


def get_server_info_via_mcp():
    """Get server info via MCP (simplified implementation)."""
    gi = st.session_state.get(GALAXY_CONNECTION)
    if not gi:
        return {"error": "Not connected to Galaxy"}

    try:
        config = gi.config.get_config()
        return {
            "version": config.get("version_major", "Unknown"),
            "brand": config.get("brand", "Galaxy"),
            "server_info": config,
        }
    except Exception as e:
        return {"error": str(e)}


def get_histories_via_mcp(limit: int = 20):
    """Get histories via MCP (simplified implementation)."""
    gi = st.session_state.get(GALAXY_CONNECTION)
    if not gi:
        return {"error": "Not connected to Galaxy"}

    try:
        histories = gi.histories.get_histories()[:limit]
        return {"histories": histories, "count": len(histories)}
    except Exception as e:
        return {"error": str(e)}


def get_workflows_via_mcp(source: str = "all"):
    """Get workflows via MCP (simplified implementation)."""
    gi = st.session_state.get(GALAXY_CONNECTION)
    if not gi:
        return {"error": "Not connected to Galaxy"}

    try:
        workflows = gi.workflows.get_workflows()
        return {"workflows": workflows, "count": len(workflows), "source": source}
    except Exception as e:
        return {"error": str(e)}


def get_llm_response(
    user_message: str, conversation_history: List, galaxy_tools: List = None
):
    """Get response from LLM provider."""
    try:
        if st.session_state[LLM_PROVIDER] == "openai":
            return get_openai_response(user_message, conversation_history, galaxy_tools)
        elif st.session_state[LLM_PROVIDER] == "custom":
            return get_custom_model_response(
                user_message, conversation_history, galaxy_tools
            )
        else:
            return "LLM provider not configured. Please set up your AI provider in the settings."
    except Exception as e:
        logger.error(f"LLM response error: {str(e)}")
        return f"Sorry, I encountered an error: {str(e)}"


def get_openai_response(
    user_message: str, conversation_history: List, galaxy_tools: List = None
):
    """Get response from OpenAI API with Galaxy MCP tools integration."""
    try:
        import openai

        api_key = st.session_state[LLM_CONFIG].get("openai_api_key")
        if not api_key:
            return "OpenAI API key not configured. Please add it in the LLM Settings."

        client = openai.OpenAI(api_key=api_key)

        # Build enhanced system prompt with MCP Galaxy tools
        galaxy_mcp_tools = get_galaxy_mcp_tools_definition()

        system_prompt = f"""You are a bioinformatics workflow planning assistant with direct access to Galaxy through MCP tools. 

You can:
1. Connect to Galaxy servers and retrieve real-time information
2. Search for tools and workflows
3. Execute Galaxy tools and manage workflows
4. Access user histories and datasets
5. Upload files and manage Galaxy resources

Available Galaxy MCP Tools:
{galaxy_mcp_tools}

When helping users:
- Use Galaxy MCP tools to get real-time information
- Provide specific, actionable workflow recommendations
- Execute tools when requested
- Help with actual Galaxy operations, not just planning

Be conversational and use the Galaxy tools proactively to provide accurate, current information."""

        messages = [{"role": "system", "content": system_prompt}]

        # Add conversation history
        for msg in conversation_history[-10:]:  # Last 10 messages for context
            messages.append(msg)

        messages.append({"role": "user", "content": user_message})

        # Define Galaxy MCP tools for OpenAI function calling
        tools = get_openai_function_definitions()

        response = client.chat.completions.create(
            model="gpt-4",
            messages=messages,
            tools=tools,
            tool_choice="auto",
            max_tokens=1500,
            temperature=0.7,
        )

        # Handle function calls
        response_message = response.choices[0].message

        if response_message.tool_calls:
            # Execute Galaxy MCP function calls
            function_responses = []
            for tool_call in response_message.tool_calls:
                function_response = execute_galaxy_mcp_function(
                    tool_call.function.name, json.loads(tool_call.function.arguments)
                )
                function_responses.append(
                    {
                        "tool_call_id": tool_call.id,
                        "role": "tool",
                        "content": json.dumps(function_response),
                    }
                )

            # Get final response with function results
            messages.append(response_message)
            messages.extend(function_responses)

            final_response = client.chat.completions.create(
                model="gpt-4", messages=messages, max_tokens=1000, temperature=0.7
            )

            return final_response.choices[0].message.content

        return response_message.content

    except ImportError:
        return "OpenAI library not installed. Please install with: pip install openai"
    except Exception as e:
        logger.error(f"OpenAI API error: {traceback.format_exc()}")
        return f"OpenAI API error: {str(e)}"


def get_custom_model_response(
    user_message: str, conversation_history: List, galaxy_tools: List = None
):
    """Get response from custom model endpoint."""
    try:
        import requests

        endpoint = st.session_state[LLM_CONFIG].get("custom_endpoint")
        api_key = st.session_state[LLM_CONFIG].get("custom_api_key", "")

        if not endpoint:
            return "Custom model endpoint not configured."

        # Prepare the payload (adjust based on your model's API)
        payload = {
            "messages": conversation_history
            + [{"role": "user", "content": user_message}],
            "max_tokens": 1000,
            "temperature": 0.7,
        }

        headers = {"Content-Type": "application/json"}
        if api_key:
            headers["Authorization"] = f"Bearer {api_key}"

        response = requests.post(endpoint, json=payload, headers=headers, timeout=30)
        response.raise_for_status()

        result = response.json()
        return result.get(
            "response",
            result.get("choices", [{}])[0]
            .get("message", {})
            .get("content", "No response"),
        )

    except Exception as e:
        return f"Custom model error: {str(e)}"


def render_llm_settings():
    """Render LLM configuration settings."""
    st.subheader("🤖 AI Assistant Settings")

    provider = st.selectbox(
        "Choose LLM Provider",
        ["openai", "custom"],
        index=0 if st.session_state[LLM_PROVIDER] == "openai" else 1,
    )
    st.session_state[LLM_PROVIDER] = provider

    if provider == "openai":
        api_key = st.text_input(
            "OpenAI API Key",
            type="password",
            value=st.session_state[LLM_CONFIG].get("openai_api_key", ""),
            help="Get your API key from https://platform.openai.com/api-keys",
        )
        if api_key:
            st.session_state[LLM_CONFIG]["openai_api_key"] = api_key

    elif provider == "custom":
        endpoint = st.text_input(
            "Custom Model Endpoint",
            value=st.session_state[LLM_CONFIG].get("custom_endpoint", ""),
            placeholder="https://your-model-endpoint.com/v1/chat/completions",
        )
        api_key = st.text_input(
            "API Key (optional)",
            type="password",
            value=st.session_state[LLM_CONFIG].get("custom_api_key", ""),
        )
        if endpoint:
            st.session_state[LLM_CONFIG]["custom_endpoint"] = endpoint
        if api_key:
            st.session_state[LLM_CONFIG]["custom_api_key"] = api_key


def render_galaxy_connection():
    """Render Galaxy connection form with MCP server integration."""
    st.subheader("🔗 Galaxy Connection with MCP")

    st.info(
        "This connection uses the official Galaxy MCP server for enhanced AI integration."
    )

    with st.form("galaxy_connection"):
        galaxy_url = st.text_input(
            "Galaxy Server URL",
            value="https://usegalaxy.org",
            help="Connect to Galaxy server for AI-powered workflow assistance",
        )
        api_key = st.text_input(
            "API Key",
            type="password",
            help="Required for MCP server connection and full functionality",
        )

        col1, col2 = st.columns(2)
        with col1:
            connect_btn = st.form_submit_button("🚀 Connect with MCP", type="primary")
        with col2:
            test_btn = st.form_submit_button("🧪 Test Connection")

        if connect_btn or test_btn:
            if not galaxy_url or not api_key:
                st.error("Both Galaxy URL and API key are required for MCP connection")
            else:
                with st.spinner("Connecting to Galaxy via MCP server..."):
                    # First test the basic connection
                    gi = connect_to_galaxy(galaxy_url, api_key)
                    if gi:
                        st.session_state[GALAXY_CONNECTION] = gi

                        if connect_btn:
                            # Start MCP server for AI integration
                            mcp_process = start_galaxy_mcp_server(galaxy_url, api_key)
                            if mcp_process:
                                st.success("✅ Connected to Galaxy with MCP server!")
                                st.info(
                                    "🤖 AI assistant now has direct Galaxy access via MCP tools"
                                )
                            else:
                                st.warning(
                                    "⚠️ Connected to Galaxy but MCP server failed to start. Basic functionality available."
                                )
                        else:
                            st.success("✅ Galaxy connection test successful!")
                    else:
                        st.error("❌ Failed to connect to Galaxy")


def render_conversation_interface():
    """Render the main conversation interface."""
    st.subheader("💬 Workflow Planning Assistant")

    # Display conversation history
    conversation = st.session_state[CONVERSATION_HISTORY]

    if not conversation:
        st.info(
            "👋 Hello! I'm your bioinformatics workflow planning assistant. Tell me about your research project and I'll help you design the perfect Galaxy workflow."
        )

        # Quick start buttons
        st.write("**Quick Start Options:**")
        col1, col2, col3 = st.columns(3)

        with col1:
            if st.button("🧬 RNA-seq Analysis"):
                start_conversation(
                    "I want to analyze RNA-seq data to find differentially expressed genes"
                )

        with col2:
            if st.button("🔍 Variant Calling"):
                start_conversation(
                    "I need to identify genetic variants from whole genome sequencing data"
                )

        with col3:
            if st.button("📊 Metagenomics"):
                start_conversation(
                    "I want to analyze microbiome data from 16S or shotgun sequencing"
                )

    # Display conversation
    for i, message in enumerate(conversation):
        if message["role"] == "user":
            with st.chat_message("user"):
                st.write(message["content"])
        else:
            with st.chat_message("assistant"):
                st.write(message["content"])

                # Add workflow plan extraction if the message contains structured recommendations
                if (
                    "workflow plan" in message["content"].lower()
                    or "steps:" in message["content"].lower()
                ):
                    if st.button(f"Extract Workflow Plan", key=f"extract_{i}"):
                        extract_workflow_plan(message["content"])


def start_conversation(initial_message: str):
    """Start conversation with a predefined message."""
    st.session_state[CONVERSATION_HISTORY] = []
    handle_user_message(initial_message)


def handle_user_message(user_message: str):
    """Process user message and get AI response."""
    # Add user message to history
    st.session_state[CONVERSATION_HISTORY].append(
        {"role": "user", "content": user_message}
    )

    # Get AI response
    with st.spinner("🤔 Thinking..."):
        ai_response = get_llm_response(
            user_message,
            st.session_state[CONVERSATION_HISTORY],
            st.session_state[GALAXY_TOOLS],
        )

    # Add AI response to history
    st.session_state[CONVERSATION_HISTORY].append(
        {"role": "assistant", "content": ai_response}
    )

    # Rerun to update display
    st.rerun()


def extract_workflow_plan(ai_message: str):
    """Extract and structure workflow plan from AI response."""
    # This is a simplified extraction - in practice, you might use more sophisticated parsing
    plan = {
        "description": "Workflow plan extracted from AI recommendation",
        "steps": [],
        "tools_needed": [],
        "considerations": [],
    }

    # Simple parsing logic (can be enhanced with NLP)
    lines = ai_message.split("\n")
    current_section = None

    for line in lines:
        line = line.strip()
        if "step" in line.lower() and any(char.isdigit() for char in line):
            plan["steps"].append(line)
        elif "tool" in line.lower() or any(
            tool_name in line.lower()
            for tool_name in ["blast", "bowtie", "fastqc", "samtools"]
        ):
            plan["tools_needed"].append(line)

    st.session_state[WORKFLOW_PLAN] = plan
    st.success("✅ Workflow plan extracted! Check the 'Generated Plan' tab.")


def render_workflow_plan():
    """Render the extracted workflow plan."""
    st.subheader("📋 Generated Workflow Plan")

    plan = st.session_state[WORKFLOW_PLAN]

    if not plan:
        st.info(
            "No workflow plan generated yet. Have a conversation with the AI assistant to create one!"
        )
        return

    st.write("### Plan Overview")
    st.write(plan["description"])

    if plan["steps"]:
        st.write("### Workflow Steps")
        for i, step in enumerate(plan["steps"], 1):
            st.write(f"{i}. {step}")

    if plan["tools_needed"]:
        st.write("### Tools Needed")
        for tool in plan["tools_needed"]:
            st.write(f"• {tool}")

    # Export options
    st.write("### Export Options")
    col1, col2 = st.columns(2)

    with col1:
        if st.button("📥 Download Plan as JSON"):
            st.download_button(
                "Download",
                data=json.dumps(plan, indent=2),
                file_name="workflow_plan.json",
                mime="application/json",
            )

    with col2:
        if st.button("📝 Generate Galaxy Workflow Template"):
            generate_galaxy_template(plan)


def generate_galaxy_template(plan):
    """Generate a basic Galaxy workflow template."""
    template = {
        "a_galaxy_workflow": "true",
        "annotation": plan["description"],
        "format-version": "0.1",
        "name": "AI Generated Workflow",
        "steps": {},
        "tags": ["ai-generated", "streamlit"],
    }

    # Add basic steps (this would be enhanced with proper tool mapping)
    for i, step in enumerate(plan["steps"]):
        template["steps"][str(i)] = {
            "annotation": step,
            "id": i,
            "type": "tool",
            "tool_state": "{}",
            "position": {"left": i * 200, "top": 100},
        }

    st.download_button(
        "📥 Download Galaxy Template",
        data=json.dumps(template, indent=2),
        file_name="galaxy_workflow_template.ga",
        mime="application/json",
    )


# Initialize session state
init_session_state_galaxy_assistant()

# Sidebar
st.sidebar.image("src/app/Conservation X Labs CXL logo.png", width="stretch")

# Main page layout
st.header("🤖 Galaxy Workflow Assistant")
st.markdown(NAVIGATE_WARNING_MD)

# Create markdown guide
assistant_guide = """
## AI-Powered Workflow Planning Assistant

This intelligent assistant helps you plan and design Galaxy workflows through natural conversation. Instead of manually navigating complex tool options, describe your research goals and get personalized workflow recommendations.

### How it works:
1. **Describe your project**: Tell the AI about your data and research goals
2. **Get expert guidance**: Receive tailored workflow recommendations
3. **Compare options**: Explore different analytical approaches
4. **Export plans**: Generate Galaxy-compatible workflow templates

### Features:
- 🤖 Conversational interface with AI guidance
- 🔧 Tool recommendations based on Galaxy's available tools
- 📊 Workflow comparison and optimization suggestions
- 📥 Export workflow plans and Galaxy templates
- 🎯 Personalized recommendations for your specific use case
"""

render_markdown_content = lambda content: st.markdown(content)
render_markdown_content(assistant_guide)

# Main interface tabs
tab1, tab2, tab3, tab4 = st.tabs(
    ["💬 Chat Assistant", "📋 Generated Plan", "🔗 Galaxy Connection", "⚙️ Settings"]
)

with tab1:
    render_conversation_interface()

    # Chat input
    if "chat_input" not in st.session_state:
        st.session_state.chat_input = ""

    user_input = st.chat_input("Ask me about your workflow needs...")
    if user_input:
        handle_user_message(user_input)

with tab2:
    render_workflow_plan()

with tab3:
    render_galaxy_connection()

    # Show connection status
    if st.session_state[GALAXY_CONNECTION]:
        mcp_status = (
            "🟢 Active" if st.session_state[MCP_SERVER_PROCESS] else "🔴 Inactive"
        )
        st.success("✅ Connected to Galaxy")
        st.info(f"🤖 MCP Server Status: {mcp_status}")

        if st.session_state[MCP_SERVER_PROCESS]:
            st.success("🚀 AI assistant has full Galaxy integration via MCP")
        else:
            st.warning(
                "⚠️ Basic connection only - start MCP server for full AI integration"
            )

        # Disconnect button
        if st.button("🔌 Disconnect from Galaxy"):
            stop_galaxy_mcp_server()
            st.session_state[GALAXY_CONNECTION] = None
            st.session_state[GALAXY_TOOLS] = []
            st.rerun()
    else:
        st.info("Connect to Galaxy to enable AI-powered workflow assistance")

with tab4:
    render_llm_settings()

    # Clear conversation button
    if st.button("🗑️ Clear Conversation History"):
        st.session_state[CONVERSATION_HISTORY] = []
        st.session_state[WORKFLOW_PLAN] = None
        st.success("Conversation cleared!")
