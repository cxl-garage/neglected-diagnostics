import streamlit as st

from common import setup

setup.initialize()
st.sidebar.image("Conservation X Labs CXL logo.png", use_column_width=True)

from common.render_method import render_markdown

render_markdown("home_markdown.md")
