from io import BytesIO

import streamlit as st

from app.common.constants import SEQVAR_CALC_BTN, SEQVAR_CALCULATED, SEQVAR_DF
from app.common.data_processing import read_fasta
from app.common.setup import init_session_state_seq_var
from app.frontend.download_data import download_seq_var_data
from genetic_testing.sequence_analysis.sequence_variability import (
    calculate_sequence_variability,
)

# Initialize the Streamlit session state for this page
init_session_state_seq_var()

# Streamlit app header
st.header("Calculate Sequence Variability")

# File uploader
uploaded_file = st.file_uploader("Upload a FASTA file", type=["fasta"])

# Display the uploaded file's content
if uploaded_file is not None:
    with uploaded_file as file:
        sequences = read_fasta(BytesIO(file.read()))
        st.write("Uploaded sequences:")
        st.write(sequences)

    # Display a form with an input box
    base_sequence = st.text_input("Enter base sequence")

    # Display a button to calculate sequence variability
    if st.button("Calculate sequence variability"):
        st.session_state[SEQVAR_CALC_BTN] = True

# Button to calculate sequence variability is clicked
if st.session_state[SEQVAR_CALC_BTN]:
    if not base_sequence:
        st.warning("Please enter a base sequence before calculating.")
    else:
        # Calculate sequence variability
        if not st.session_state[SEQVAR_CALCULATED]:
            st.session_state[SEQVAR_DF] = calculate_sequence_variability(
                sequences, base_sequence
            )
            # Set the flag for calculation of sequence variability to True
            st.session_state[SEQVAR_CALCULATED] = True

        # Display the sequence variability data
        if not st.session_state[SEQVAR_DF].empty:
            st.write("Sequence Variability Data:")
            st.dataframe(st.session_state[SEQVAR_DF])

            # Enter the species name for downloading the data
            species = st.text_input("Enter the species name:")

            if st.button("Prepare for download"):
                download_seq_var_data(st.session_state[SEQVAR_DF], species)
        else:
            st.write("No sequence variability data to download")
