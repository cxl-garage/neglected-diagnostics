import streamlit as st
import streamlit.components.v1 as components
from common.render_method import render_markdown

from app.common import setup
from app.common.constants import MSA_DF, MSA_VIEWER, NAVIGATE_WARNING_MD
from app.common.setup import init_session_state_msa
from app.frontend.download_data import download_msa_data
from genetic_testing.sequence_tool.consensusSeq import ConsensusSeq, alignMs
from genetic_testing.sequence_tool.msaViewer import MsaViewer

setup.initialize()
st.sidebar.image("src/app/Conservation X Labs CXL logo.png", use_column_width=True)

# Initialize the Streamlit session state for this page
init_session_state_msa()

# Streamlit app header
st.header("Sequence data cleaning and MSA")

st.markdown(NAVIGATE_WARNING_MD)
render_markdown("src/app/multisequence_alignment.md")

# File uploader
multifastaFile = st.file_uploader("Upload a FASTA file", type=["fasta"])

with st.form("seq_viewer"):
    # Button to trigger viewer
    submitted = st.form_submit_button("Show sequence viewer")
    if submitted:
        # Make sure kill other viewer threads
        if st.session_state[MSA_VIEWER] is None:
            st.session_state[MSA_VIEWER] = MsaViewer(multifastaFile)
            st.session_state[MSA_VIEWER].start()
        else:
            st.session_state[MSA_VIEWER].updateData(multifastaFile)

        components.iframe(
            "http://127.0.0.1:8050/", width=1200, height=1200, scrolling=True
        )
if multifastaFile is not None:
    consensus = ConsensusSeq(multifastaFile)

with st.form("multi_seq_alignment"):
    # Button to trigger alignment
    submitted = st.form_submit_button("Align sequences")

    if submitted:
        if not multifastaFile:
            st.warning("Please upload a FASTA file before calculating.")
        # Execute the alignment
        st.session_state[MSA_DF] = alignMs(multifastaFile)

if st.session_state[MSA_DF] is not None:
    st.write("Multisequence alignment Data:")
    mda_df = st.session_state[MSA_DF]
    st.dataframe(mda_df)

    # Enter the species name for downloading the data
    species = st.text_input("Enter the species name:")
    if st.button("Prepare for download"):
        download_msa_data(mda_df, species)
else:
    st.warning("No multisequence alignment data to download.")
