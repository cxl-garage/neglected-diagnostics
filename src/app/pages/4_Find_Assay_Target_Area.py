from io import StringIO

import streamlit as st
from Bio import SeqIO

from app.common.constants import TGT_AREA_DF, TGT_AREA_FORM

# from app.common.data_processing import parse_fasta_files
from app.common.setup import init_session_state_tgt_area
from genetic_testing.assay_target.datatypes import AssayTargetColumns
from genetic_testing.assay_target.primer_design import find_target_area

# Initialize the Streamlit session state for this page
init_session_state_tgt_area()

CSV_DOWNLOAD = "Download as CSV"
FASTA_DOWNLOAD = "Download sequences as a fasta file"

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

st.divider()

target_area_container = st.container()
with target_area_container:
    st.subheader("Target Area Parameters")
    with st.form("find_target_area"):
        # Reference Sequence
        reference_sequence = st.text_input(
            "Select Reference Sequence",
            # [file.name for file in target_files],
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

        submitted = st.form_submit_button("Find Target Area")
        if submitted:
            # Find the target area
            try:
                st.session_state[TGT_AREA_DF] = find_target_area(
                    target_files=target_files,
                    off_target_files=off_target_files,
                    reference_sequence=reference_sequence,
                    window=tgt_region_size,
                    slide=slide_size,
                    max_dif=max_difference,
                )
                st.session_state[TGT_AREA_FORM] = True
            except ValueError as e:
                st.error(e)
                st.session_state[TGT_AREA_FORM] = False

if st.session_state[TGT_AREA_FORM]:
    st.write(st.session_state[TGT_AREA_DF])

    # Create a dropdown menu
    download_format = st.selectbox("Select an option:", [CSV_DOWNLOAD, FASTA_DOWNLOAD])

    if download_format == CSV_DOWNLOAD:
        # Option 1: Download the entire DataFrame as CSV
        st.download_button(
            label="Download",
            data=st.session_state[TGT_AREA_DF].to_csv(index=False).encode("utf-8"),
            file_name="potential_target_area.csv",
            mime="text/csv",
            key="csv_button",
            help="Download this table as a CSV file",
            args={"as_attachment": True},
        )

    elif download_format == FASTA_DOWNLOAD:
        cols = AssayTargetColumns()
        fasta_data = "\n".join(st.session_state[TGT_AREA_DF][cols.assay_design_area])
        st.download_button(
            label="Download Fasta",
            data=fasta_data,
            file_name="potential_target_area.fasta",
            key="fasta_button",
            help="Download sequences as a FASTA file",
            args={"as_attachment": True},
        )
