import streamlit as st


def render_markdown(file_name):
    # Read and display the content from the Markdown file
    with open(file_name, "r") as f:
        content = f.read()
    st.markdown(content)
