import streamlit as st

from app.common import setup

setup.initialize()
st.sidebar.image("src/app/Conservation X Labs CXL logo.png", width="stretch")

# Add expander popup chat assistant
try:
    from app.assistant.popup_chat import render_expander_popup_chat

    render_expander_popup_chat("Home")
except ImportError:
    pass  # Assistant not available

from common.render_method import render_markdown

render_markdown("src/app/home_markdown.md")
