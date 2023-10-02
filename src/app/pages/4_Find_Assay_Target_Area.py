from io import StringIO

import streamlit as st
from Bio import SeqIO

from app.common.constants import TGT_AREA_BTN
from app.common.setup import init_session_state_tgt_area
from genetic_testing.assay_target.primer_design import find_target_area

# Initialize the Streamlit session state for this page
init_session_state_tgt_area()

# Streamlit app header
st.header("Find Assay Target Area")

target_container = st.container()

with target_container:
    st.subheader("Upload Target Sequences")
    # Create a place to upload FASTA files
    target_files = st.file_uploader(
        "Upload FASTA files",
        key="target_files_uploader",
        type=["fasta", "fas"],
        accept_multiple_files=True,
        help="Only single FASTA files are allowed",
    )

    # # Create a checkbox to enable multiselect
    # tgt_split_flag = st.checkbox(
    #     "Split Target Files",
    #     help="Select this option to split a target file into single FASTA files",
    # )

    # if tgt_split_flag:
    #     # Create checkboxes to allow user to select files for preprocessing
    #     tgt_split_files = st.multiselect(
    #         "Select Target Files for Splitting", [file.name for file in target_files]
    #     )

    #     st.write(tgt_split_files)

st.divider()

off_target_container = st.container()
with off_target_container:
    st.subheader("Off-Target Sequences")

    # Create a place to upload FASTA files
    off_target_files = st.file_uploader(
        "Upload FASTA files",
        key="off_target_files_uploader",
        type=["fasta", "fas"],
        accept_multiple_files=True,
        help="Only single FASTA files are allowed",
    )

    # # Create a checkbox to enable multiselect
    # off_tgt_split_flag = st.checkbox(
    #     "Split Off-Target Files",
    #     help="Select this option to split a off-target file into single FASTA files",
    # )

    # if off_tgt_split_flag:
    #     # Create checkboxes to allow user to select files for preprocessing
    #     off_tgt_split_files = st.multiselect(
    #         "Select Off-Target Files for Splitting",
    #         [file.name for file in off_target_files],
    #     )

    #     st.write(off_tgt_split_files)

st.divider()

target_area_container = st.container()
with target_area_container:
    if st.button("Find Target Area", help="Click this button to find the target area"):
        st.session_state[TGT_AREA_BTN] = True

if st.session_state[TGT_AREA_BTN]:
    st.write("Target Area Button clicked")
    df = find_target_area(target_files, off_target_files)
    st.write(df)

    # Download the transpose of sequence variability table as a CSV file
    st.download_button(
        label="Download",
        data=df.to_csv().encode("utf-8"),
        file_name="sequence_variability.csv",
        mime="text/csv",
        help="Download the transpose of sequence variability table as a CSV file",
    )
