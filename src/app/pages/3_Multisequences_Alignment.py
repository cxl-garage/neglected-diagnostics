import configparser
import os

import streamlit as st
import streamlit.components.v1 as components
from common.render_method import render_markdown

from app.common import setup
from app.common.constants import (
    MSA_CLEAN_ZIP,
    MSA_DF,
    MSA_VIEWER,
    MSA_VIEWER_SERVICE_URL,
    NAVIGATE_WARNING_MD,
)
from app.common.setup import init_session_state_msa
from app.frontend.download_data import download_msa_data
from genetic_testing.sequence_tool.consensusSeq import ConsensusSeq, alignMs
from genetic_testing.sequence_tool.msaViewer import MsaViewer, msa_cleaner

setup.initialize()
st.sidebar.image("src/app/Conservation X Labs CXL logo.png", use_column_width=True)

# Initialize the Streamlit session state for this page
init_session_state_msa()

# Streamlit app header
st.header("Sequence data alignment and cleanup")

st.markdown(NAVIGATE_WARNING_MD)
render_markdown("src/app/multisequence_alignment.md")

st.markdown("### Align multisequences")

# File uploader
multifastaFile = st.file_uploader("Upload a FASTA file", type=["fasta"])

with st.form("seq_viewer"):
    # Button to trigger viewer
    submitted = st.form_submit_button("Show sequence viewer")
    if submitted:
        if not multifastaFile:
            st.warning("Please upload a FASTA file for sequence viewer")
        else:
            # Make sure kill other viewer threads
            if st.session_state[MSA_VIEWER] is None:
                st.session_state[MSA_VIEWER] = MsaViewer(multifastaFile)
                st.session_state[MSA_VIEWER].run()
            else:
                st.session_state[MSA_VIEWER].updateData(multifastaFile)

            components.iframe(
                MSA_VIEWER_SERVICE_URL, width=1200, height=1200, scrolling=True
            )
if multifastaFile is not None:
    consensus = ConsensusSeq(multifastaFile)

with st.form("multi_seq_alignment"):
    # Button to trigger alignment
    submitted = st.form_submit_button("Align sequences")

    if submitted:
        if not multifastaFile:
            st.warning("Please upload a FASTA file before calculating.")
        else:
            # Execute the alignment
            st.session_state[MSA_DF] = alignMs(multifastaFile)

if not st.session_state[MSA_DF].empty:
    st.markdown("### Alignned multisequences")
    mda_df = st.session_state[MSA_DF]
    st.dataframe(mda_df)

    # Enter the species name for downloading the data
    species = st.text_input("Enter the species name:")
    if st.button("Prepare for download"):
        download_msa_data(mda_df, species)
else:
    st.write("No multisequence alignment data to download.")


st.markdown("### Multisequence alignment clean up")
with st.form("msa_cleanup"):
    # File uploader
    msafastaFile = st.file_uploader("Upload MSA FASTA for cleanup", type=["fasta"])

    # CIAlign user options
    st.markdown("#### MSA cleanup options")

    col1, col2, col3, col4 = st.columns(4)
    config = configparser.ConfigParser()
    config["remove_divergent"] = {}
    config["remove_insertions"] = {}
    config["remove_short"] = {}
    config["crop_ends"] = {}
    config["crop_divergent"] = {}
    config["unalign"] = {}
    config["replace"] = {}
    config["section"] = {}

    with col1:
        config["remove_divergent"]["remove_divergent"] = str(
            st.toggle("Enable Divergent cleaning")
        )
        config["remove_divergent"]["remove_divergent_minperc"] = st.text_input(
            "Remove Divergent - min_perc (0~1)",
            value="0.65",
            placeholder="Default: 0.65",
        )
        config["remove_short"]["remove_short"] = str(st.toggle("Enable Short cleaning"))
        config["remove_short"]["remove_min_length"] = st.text_input(
            "Remove Short - min_length",
            value="50",
            placeholder="Default: 50",
        )

    with col2:
        config["remove_insertions"]["remove_insertions"] = str(
            st.toggle("Enable Insertion cleaning")
        )
        config["remove_insertions"]["insertion_min_perc"] = st.text_input(
            "Remove Insertion - min_perc (0~1)",
            value="0.5",
            placeholder="Default: 0.5",
        )
        config["remove_insertions"]["insertion_min_size"] = st.text_input(
            "Remove Insertion - min_size",
            value="3",
            placeholder="Default: 3",
        )
        config["remove_insertions"]["insertion_max_size"] = st.text_input(
            "Remove Insertion - max_size",
            value="200",
            placeholder="Default: 200",
        )
        config["remove_insertions"]["insertion_min_flank"] = st.text_input(
            "Remove Insertion - min_flank",
            value="5",
            placeholder="Default: 5",
        )

    with col3:
        config["crop_ends"]["crop_ends"] = str(st.toggle("Enable Crop Ends cleaning"))
        config["crop_ends"]["crop_ends_mingap_perc"] = st.text_input(
            "Crop Ends - min_gap_perc (0~0.06)",
            value="0.05",
            placeholder="Default: 0.05",
        )
        config["crop_ends"]["crop_ends_redefine_perc"] = st.text_input(
            "Crop Ends - redefine_perc (0~0.5)",
            value="0.1",
            placeholder="Default: 0.1",
        )

    with col4:
        config["crop_divergent"]["crop_divergent"] = str(
            st.toggle("Enable Crop Divergent cleaning")
        )
        config["crop_divergent"]["crop_divergent_min_prop_ident"] = st.text_input(
            "Crop Divergent - min_prop_ident",
            value="0.5",
            placeholder="Default: 0.5",
        )
        config["crop_divergent"]["crop_divergent_min_prop_nongap"] = st.text_input(
            "Crop Divergent - min_prop_nongap",
            value="0.5",
            placeholder="Default: 0.5",
        )
        config["crop_divergent"]["crop_divergent_buffer_size"] = st.text_input(
            "Crop Divergent - buffer_size",
            value="5",
            placeholder="Default: 5",
        )

    colA, colB = st.columns(2)

    with colA:
        config["unalign"]["unalign_input"] = str(
            st.toggle("Generate a copy of the input alignment with no gaps")
        )
        config["unalign"]["unalign_output"] = str(
            st.toggle("Generate a copy of the cleaned alignment with no gaps")
        )
        config["replace"]["replace_input_ut"] = str(
            st.toggle("Replaces all Us by Ts in input alignment")
        )
        config["replace"]["replace_output_ut"] = str(
            st.toggle("Replaces all Us by Ts in output alignment")
        )

    with colB:
        config["replace"]["replace_input_tu"] = str(
            st.toggle("Replaces all Ts by Us in input alignment")
        )
        config["replace"]["replace_output_tu"] = str(
            st.toggle("Replaces all Ts by Us in output alignment")
        )
        config["section"]["get_section"] = str(
            st.toggle("Retrieve a section of the alignment")
        )

    # Button to trigger alignment
    submitted = st.form_submit_button("MSA cleanup")
    if submitted:
        if not msafastaFile:
            st.warning("Please upload a MSA FASTA file before cleanup.")
        else:
            # Execute the MSA cleanup
            zipFilePath = msa_cleaner(msafastaFile, config)
            if st.session_state[MSA_CLEAN_ZIP] is not None:
                preZipFile = st.session_state[MSA_CLEAN_ZIP]
                os.remove(preZipFile)
            st.session_state[MSA_CLEAN_ZIP] = zipFilePath

if st.session_state[MSA_CLEAN_ZIP] is not None:
    zipFilePath = st.session_state[MSA_CLEAN_ZIP]
    zipFileName = os.path.basename(zipFilePath)
    # Display a download button
    with open(zipFilePath, "rb") as fp:
        btn = st.download_button(
            label="Download ZIP", data=fp, file_name=zipFileName, mime="application/zip"
        )
        if btn:
            os.remove(zipFilePath)
            st.session_state[MSA_CLEAN_ZIP] = None
