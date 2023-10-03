from io import StringIO

import streamlit as st
from Bio import SeqIO

from app.common.constants import TGT_AREA_BTN
from app.common.data_processing import parse_fasta_files
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

    target_sequences = parse_fasta_files(target_files)

    # Sort the files by name
    # target_files = sorted(target_files, key=lambda x: x.name)

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

    off_target_sequences = parse_fasta_files(off_target_files)

    # Sort the files by name
    # off_target_files = sorted(off_target_files, key=lambda x: x.name)

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
    st.subheader("Target Area Parameters")

    # Reference Sequence
    reference_sequence = st.selectbox(
        "Select Reference Sequence",
        [file.name for file in target_files],
        help="If no selection is made, the first file will be used as the reference sequence",
    )

    # Target region size
    tgt_region_size = st.number_input(
        label="Enter the target region size (# of base pairs)",
        min_value=0,
        value=300,
        step=1,
        help="If no value is entered, default value of 300 will be used as the target region size",
    )

    # Target region slide size
    slide_size = st.number_input(
        label="Enter the target region slide size (# of base pairs)",
        min_value=0,
        value=20,
        step=1,
        help="If no value is entered, default value of 20 will be used as the target region slide size",
    )

    # Maximum allowed differences between primers and targets
    max_difference = st.number_input(
        label="Enter the maximum differences allowed between target area and target sequences",
        min_value=0,
        value=5,
        step=1,
        help="If no value is entered, default value of 5 will be used as the maximum differences allowed between target area and target sequences",
    )

    # Create a button to find the target area
    if st.button("Find Target Area", help="Click this button to find the target area"):
        st.session_state[TGT_AREA_BTN] = True

if st.session_state[TGT_AREA_BTN]:
    st.write(reference_sequence)
    st.write(tgt_region_size)
    st.write(slide_size)
    st.write(max_difference)

    # Find the target area
    df = find_target_area(
        target_files=target_files,
        off_target_files=off_target_files,
        reference_sequence=reference_sequence,
        window=tgt_region_size,
        slide=slide_size,
        max_dif=max_difference,
    )
    st.write(df)

    # Download the transpose of sequence variability table as a CSV file
    st.download_button(
        label="Download",
        data=df.to_csv().encode("utf-8"),
        file_name="sequence_variability.csv",
        mime="text/csv",
        help="Download the transpose of sequence variability table as a CSV file",
    )
