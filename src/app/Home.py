import streamlit as st

from app.common import setup

setup.initialize()
st.sidebar.image("src/app/Conservation X Labs CXL logo.png", width="stretch")

from common.render_method import render_markdown

render_markdown("src/app/home_markdown.md")
