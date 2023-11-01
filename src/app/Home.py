import streamlit as st

st.sidebar.image("Conservation X Labs CXL logo.png", use_column_width=True)

from common.render_method import render_markdown

# st.sidebar.image('Conservation X Labs CXL logo.png', use_column_width=True)
render_markdown("home_markdown.md")
