import streamlit as st

st.sidebar.image("src/app/Conservation X Labs CXL logo.png", use_column_width=True)

from common.render_method import render_markdown

render_markdown("src/app/home_markdown.md")
