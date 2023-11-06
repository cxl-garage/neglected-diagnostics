from io import BytesIO

import streamlit as st
from common.render_method import render_markdown

from app.common import setup
from app.common.constants import NAVIGATE_WARNING_MD, SEQVAR_DF
from app.common.data_processing import read_fasta
from app.common.setup import init_session_state_seq_var
from app.frontend.download_data import download_seq_var_data
from genetic_testing.sequence_analysis.sequence_variability import (
    calculate_sequence_variability,
)

setup.initialize()
st.sidebar.image("src/app/Conservation X Labs CXL logo.png", use_column_width=True)

# Initialize the Streamlit session state for this page
init_session_state_seq_var()

# Streamlit app header
st.header("Determine Sequence Haplotypes")

st.markdown(NAVIGATE_WARNING_MD)
render_markdown("src/app/sequence_haplotypes_guide.md")


# File uploader
uploaded_file = st.file_uploader("Upload a FASTA file", type=["fasta"])

# Display the uploaded file's content
if uploaded_file is not None:
    with uploaded_file as file:
        sequences = read_fasta(BytesIO(file.read()))
        st.write("Uploaded sequences:")
        st.write(sequences)

with st.form("seq_haplotypes"):
    # Display a form with an input box
    base_sequence = st.text_input("Enter base sequence")

    # Button to trigger sequence haplotypes calculation
    submitted = st.form_submit_button("Find Sequence Haplotypes")
    if submitted:
        if not base_sequence:
            st.warning("Please enter a base sequence before calculating.")
        else:
            # Calculate sequence haplotypes
            st.session_state[SEQVAR_DF] = calculate_sequence_variability(
                sequences, base_sequence
            )


# Display the sequence variability data
if st.session_state[SEQVAR_DF] is not None:
    st.write("Sequence Variability Data:")
    st.dataframe(st.session_state[SEQVAR_DF])

    # Enter the species name for downloading the data
    species = st.text_input("Enter the species name:")

    if st.button("Prepare for download"):
        download_seq_var_data(st.session_state[SEQVAR_DF], species)
else:
    st.write("No sequence variability data to download")
