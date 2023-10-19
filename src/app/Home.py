import streamlit as st

# Read and display the content from the Markdown file
with open("src/app/home_markdown.md", "r") as f:
    content = f.read()
st.markdown(content)
