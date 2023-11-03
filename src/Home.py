import streamlit as st

st.sidebar.image("assets//Conservation X Labs CXL logo.png", use_column_width=True)

st.markdown("""
    <style>
        .reportview-container {
            margin-top: -2em;
        }
        #MainMenu {visibility: hidden;}
        .stDeployButton {display:none;}
        footer {visibility: hidden;}
        #stDecoration {display:none;}
    </style>
""", unsafe_allow_html=True)

from app.common.render_method import render_markdown

render_markdown("assets/home_markdown.md")
