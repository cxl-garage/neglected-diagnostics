import streamlit as st
from common.render_method import render_markdown

from app.common import setup
from app.common.constants import NAVIGATE_WARNING_MD, TGT_AREA_DF, TGT_AREA_FORM
from app.common.data_processing import get_first_header
from app.common.setup import init_session_state_tgt_area
from genetic_testing.assay_target.datatypes import AssayTargetColumns
from genetic_testing.assay_target.primer_design import find_target_area

setup.initialize()
st.sidebar.image("src/app/Conservation X Labs CXL logo.png", use_column_width=True)


# Initialize the Streamlit session state for this page
init_session_state_tgt_area()

# Define constants for download options
CSV_DOWNLOAD = "Download as CSV"
FASTA_DOWNLOAD = "Download sequences as a fasta file"
TARGET_AREA_PREFIX = "Targetarea"

# Streamlit app header
st.header("Find Assay Target Area")

st.markdown(NAVIGATE_WARNING_MD)
render_markdown("src/app/find_assay_target_guide.md")

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

# Create a container for target area parameters
target_area_container = st.container()
with target_area_container:
    st.subheader("Target Area Parameters")
    with st.form("find_target_area"):
        # Input field for selecting Reference sequence, Target region Size, Target region Slide size, and Maximum
        # allowed differences between primers and targets
        reference_sequence = st.text_input(
            "Select Reference Sequence",
            get_first_header(target_files),
            help="If no selection is made, the first sequence header in the first file will be used as the reference "
            "sequence",

        )

        tgt_region_size = st.number_input(
            label="Enter the target region size (# of base pairs)",
            min_value=60,
            value=300,
            step=1,
            help="If no value is entered, default value of 300 will be used as the target region size",
        )

        slide_size = st.number_input(
            label="Enter the target region slide size (# of base pairs)",
            min_value=0,
            value=20,
            step=1,
            help="If no value is entered, default value of 20 will be used as the target region slide size",
        )

        max_differenc_o = st.number_input(
            label="Enter the maximum differences allowed between primer and target sequences",
            min_value=0,
            value=5,
            step=1,
            help="If no value is entered, default value of 5 will be used as the maximum differences allowed between "
            "primer and target sequences",
        )

        max_difference_ot = st.number_input(
            label="Enter the maximum differences allowed between primer and off-target sequences",
            min_value=10,
            value=30,
            step=1,
            help="If no value is entered, default value of 15 will be used as the maximum differences allowed between "
            "primer and off-target sequences",
        )

        # Button to trigger target area calculation
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
                    maxDif_t=max_differenc_o,
                    maxDif_ot=max_difference_ot,
                )
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
        # Option 2: Download sequences as a FASTA file
        cols = AssayTargetColumns()
        df = st.session_state[TGT_AREA_DF]
        seqs = df[cols.assay_design_area].tolist()
        fasta_data = "\n".join(
            [f"{TARGET_AREA_PREFIX}{i+1}\n{seq}" for i, seq in enumerate(seqs)]
        )
        st.download_button(
            label="Download Fasta",
            data=fasta_data,
            file_name="potential_target_area.fasta",
            key="fasta_button",
            help="Download sequences as a FASTA file",
            args={"as_attachment": True},
        )
