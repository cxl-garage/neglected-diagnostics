import streamlit as st
import streamlit.components.v1 as components
from common.render_method import render_markdown

from app.common import setup
from app.common.constants import NAVIGATE_WARNING_MD, SEQVAR_DF
from app.common.setup import init_session_state_seq_var
from genetic_testing.sequence_tool.msaViewer import start_sequence_viewer

setup.initialize()
st.sidebar.image("src/app/Conservation X Labs CXL logo.png", use_column_width=True)

# Initialize the Streamlit session state for this page
init_session_state_seq_var()

# Streamlit app header
st.header("Sequence data cleaning and MSA")

st.markdown(NAVIGATE_WARNING_MD)
render_markdown("src/app/multisequence_alignment.md")


# File uploader
uploaded_file = st.file_uploader("Upload a FASTA file", type=["fasta"])

# Display the uploaded file's content
if uploaded_file is not None:
    with uploaded_file as file:
        sequences = file.read()

with st.form("seq_viewer"):
    # Button to trigger viewer
    submitted = st.form_submit_button("Show sequence viewer")
    if submitted:
        # Make sure kill other viewer thr

        # Launch seq viewer
        st.session_state[SEQVAR_DF] = start_sequence_viewer(sequences)
        # Display the sequence viewer at full width
        components.iframe(
            "http://127.0.0.1:8050/", width=1200, height=1800, scrolling=True
        )
