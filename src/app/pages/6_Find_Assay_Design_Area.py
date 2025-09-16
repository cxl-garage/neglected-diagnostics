import time

import streamlit as st
from common.render_method import render_markdown

from app.common import setup
from app.common.constants import NAVIGATE_WARNING_MD, TGT_AREA_DF, TGT_AREA_FORM
from app.common.data_processing import get_headers
from app.common.setup import init_session_state_tgt_area
from genetic_testing.assay_target.datatypes import AssayTargetColumns
from genetic_testing.assay_target.primer_design import find_target_area

setup.initialize()
st.sidebar.image("src/app/Conservation X Labs CXL logo.png", width="stretch")


# Initialize the Streamlit session state for this page
init_session_state_tgt_area()

# Define constants for download options
CSV_DOWNLOAD = "Download as CSV"
FASTA_DOWNLOAD = "Download sequences as a fasta file"
TARGET_AREA_PREFIX = "Assay_Design_Area_"

# Streamlit app header
st.header("Find Assay Design Area")

st.markdown(NAVIGATE_WARNING_MD)
render_markdown("src/app/find_assay_design_guide.md")

# Create a container for uploading target sequences
target_container = st.container()
with target_container:
    st.subheader("Upload Target Sequences")
    # Create a place to upload FASTA files
    target_files = st.file_uploader(
        "Upload FASTA files",
        key="target_files_uploader",
        type=["fasta", "fas"],
        accept_multiple_files=True,
        help="Only single FASTA or non-aligned multi-FASTA files are allowed",
    )

st.divider()

# Create a container for uploading off-target sequences
off_target_container = st.container()
with off_target_container:
    st.subheader("Upload Off-Target Sequences")

    # Create a place to upload FASTA files
    off_target_files = st.file_uploader(
        "Upload FASTA files",
        key="off_target_files_uploader",
        type=["fasta", "fas"],
        accept_multiple_files=True,
        help="Only single FASTA or non-aligned multi-FASTA files are allowed",
    )

st.divider()
cols = AssayTargetColumns()
# Create a container for target area parameters
target_area_container = st.container()
with target_area_container:
    st.subheader("Design Area Parameters")
    with st.form("find_target_area"):
        # Input field for selecting Reference sequence, Target region Size, Target region Slide size, and Maximum
        # allowed differences between primers and targets
        reference_sequence = st.selectbox(
            "Select Reference Sequence",
            get_headers(target_files),
            help="Select from the dropdown or type in the sequence header for the desired reference sequence",
            key="sequence",
        )

        tgt_region_size = st.number_input(
            label="Enter the assay design region size (# of base pairs)",
            min_value=60,
            value=300,
            step=1,
            help="If no value is entered, default value of 300 will be used as the target region size",
        )

        slide_size = st.number_input(
            label="Enter the assay design region slide size (# of base pairs)",
            min_value=0,
            value=20,
            step=1,
            help="If no value is entered, default value of 20 will be used as the target region slide size",
        )

        max_differenc_o = st.number_input(
            label="Enter the maximum differences allowed between assay area and target sequences",
            min_value=0,
            value=5,
            step=1,
            help="If no value is entered, default value of 5 will be used as the maximum differences allowed between "
            "primer and target sequences",
        )

        max_difference_ot = st.number_input(
            label="Enter the maximum differences allowed between assay area and off-target sequences",
            min_value=10,
            value=30,
            step=1,
            help="If no value is entered, default value of 15 will be used as the maximum differences allowed between "
            "primer and off-target sequences",
        )

        # Button to trigger target area calculation
        submitted = st.form_submit_button("Find Assay Design Area")
        if submitted:
            if len(target_files) < 2:
                st.error("A minimum of 2 Target Sequence files are required.")
                st.session_state[TGT_AREA_FORM] = False
            else:
                progress_bar = st.progress(0, text="Beginning process...")
                # Find the target area
                try:
                    df = find_target_area(
                        target_files=target_files,
                        off_target_files=off_target_files,
                        reference_sequence=reference_sequence,
                        progress_bar=progress_bar,
                        window=tgt_region_size,
                        slide=slide_size,
                        maxDif_t=max_differenc_o,
                        maxDif_ot=max_difference_ot,
                    )
                    target_id_col = df.index.map(lambda i: f"{TARGET_AREA_PREFIX}{i+1}")
                    ref_col = df.index.map(lambda i: reference_sequence)
                    df.insert(
                        loc=0, column=cols.assay_design_area_reference, value=ref_col
                    )
                    df.insert(loc=0, column=cols.target_id, value=target_id_col)
                    st.session_state[TGT_AREA_DF] = df
                    st.session_state[TGT_AREA_FORM] = True
                except ValueError as e:
                    # Handle the case where an error occurs during target area calculation
                    st.error(e)
                    st.session_state[TGT_AREA_FORM] = False

# If target area has been successfully calculated, display it and provide download options
if st.session_state[TGT_AREA_FORM]:
    st.write(st.session_state[TGT_AREA_DF])

    # Dropdown menu for selecting download format
    download_format = st.selectbox("Select an option:", [CSV_DOWNLOAD, FASTA_DOWNLOAD])

    df = st.session_state[TGT_AREA_DF]
    if download_format == CSV_DOWNLOAD:
        # Option 1: Download the entire DataFrame as CSV
        st.download_button(
            label="Download",
            data=df.to_csv(index=False).encode("utf-8"),
            file_name="potential_target_area.csv",
            mime="text/csv",
            key="csv_button",
            help="Download this table as a CSV file",
            args={"as_attachment": True},
        )

    elif download_format == FASTA_DOWNLOAD:
        # Option 2: Download sequences as a FASTA file
        seqs = df[cols.assay_design_area].tolist()
        target_ids = df[cols.target_id].tolist()
        fasta_data = "\n".join([f">{id}\n{seq}" for id, seq in zip(target_ids, seqs)])
        st.download_button(
            label="Download Fasta",
            data=fasta_data,
            file_name="potential_target_area.fasta",
            key="fasta_button",
            help="Download sequences as a FASTA file",
            args={"as_attachment": True},
        )
