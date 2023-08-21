from io import BytesIO

import streamlit as st

from app.common.data_processing import read_fasta
from genetic_testing.sequence_variability import calculate_sequence_variability

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
    base_sequence = st.text_input("Enter Base Sequence")

    # Display a button to calculate sequence variability
    if st.button("Calculate Sequence Variability"):
        if not base_sequence:
            st.warning("Please enter a base sequence before calculating.")
        else:
            df = calculate_sequence_variability(sequences, base_sequence)
            st.write("Sequence Variability Calculation Results:")
            st.dataframe(df)
