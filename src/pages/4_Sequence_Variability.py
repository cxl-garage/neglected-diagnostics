from io import BytesIO

import streamlit as st
from app.common.render_method import render_markdown

from app.common.constants import SEQVAR_TABLE, SEQVAR_TABLE_BTN, SEQVAR_TABLE_CALC
from app.common.data_processing import read_fasta
from app.common.setup import init_session_state_seq_var_table
from genetic_testing.sequence_analysis.sequence_variability import (
    calculate_seq_variability_table,
)

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

st.sidebar.image("assets//Conservation X Labs CXL logo.png", use_column_width=True)


# Initialize the Streamlit session state for this page
init_session_state_seq_var_table()

# Streamlit app header
st.header("Calculate Sequence Variability")

render_markdown("assets//sequence_variability_guide.md")

# File uploader
uploaded_file = st.file_uploader("Upload a FASTA file", type=["fasta"])

# Display the uploaded file's content
if uploaded_file is not None:
    with uploaded_file as file:
        sequences = read_fasta(BytesIO(file.read()))
        st.write("Uploaded sequences:")
        st.write(sequences)

# Display a form with an input box for filter threshold
seqvar_threshold = st.text_input("Enter the maximum sequence variability threshold")

# Display a button to display the sequence variability table
if st.button("Display sequence variability table"):
    st.session_state[SEQVAR_TABLE_BTN] = True


# Button to display the sequence variability table is clicked
if st.session_state[SEQVAR_TABLE_BTN]:
    # Calculate sequence variability table
    if not st.session_state[SEQVAR_TABLE_CALC]:
        st.session_state[SEQVAR_TABLE] = calculate_seq_variability_table(
            list(sequences.values())
        )
        # Set the flag for calculation of sequence variability table to True
        st.session_state[SEQVAR_TABLE_CALC] = True

# Filter the sequence variability table based on the threshold
if st.session_state[SEQVAR_TABLE_CALC]:
    # Check if the threshold is entered
    if not seqvar_threshold:
        seqvar_threshold = 100.00

    # Filter all the rows to be less than threshold
    df = st.session_state[SEQVAR_TABLE].copy()
    df = df[df <= float(seqvar_threshold)].dropna()
    st.dataframe(df)

    # Download the transpose of sequence variability table as a CSV file
    st.download_button(
        label="Download sequence variability table as CSV",
        data=df.transpose().to_csv().encode("utf-8"),
        file_name="sequence_variability_table.csv",
        mime="text/csv",
        help="Download the transpose of sequence variability table as a CSV file",
    )
