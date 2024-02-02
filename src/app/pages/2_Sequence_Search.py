import pandas as pd
import streamlit as st
from common.render_method import render_markdown

from app.common import data_processing, setup
from app.common.constants import (
    MAIN_PAGE_COLS_GAP,
    MAIN_PAGE_COLS_SIZE,
    MAX_SEQUENCE_LENGTH_DOWNLOAD,
    NAVIGATE_WARNING_MD,
    NCBI_DF,
    NCBI_DF_FILTER,
    NCBI_SUMMARY_FORM,
    TOP_N_ORGANISMS,
)
from app.frontend.aggrid_table import aggrid_table
from app.frontend.download_data import download_data
from app.frontend.filter_dataframe import filter_dataframe
from app.frontend.filters_info import show_filters_info
from genetic_testing.routers import ncbi

# Initialize the Streamlit session state keys
setup.initialize()
# Define constants for download options
CSV_DOWNLOAD = "Download table as CSV"
FASTA_DOWNLOAD = "Download sequences as a fasta file"
# Divides the page layout into two columns of relative width
query_col, summary_col = st.columns(MAIN_PAGE_COLS_SIZE, gap=MAIN_PAGE_COLS_GAP)
with query_col:
    st.header("Search and Retrieve Sequences")
    st.markdown(NAVIGATE_WARNING_MD)
    render_markdown("src/app/sequence_search_quick_guide.md")
    # Streamlit form to capture search conditions
    with st.form("query"):
        database = st.selectbox("Select the database to search", ("nucleotide", "gene"))
        search_term = st.text_input(
            label="Enter the search term",
            placeholder="Example: human[organism] AND topoisomerase[protein name]",
        )
        submitted = st.form_submit_button("Submit")
        if submitted:
            # search term is empty
            if search_term == "":
                st.error("Please input a search term before you submit")
            else:
                st.session_state[NCBI_SUMMARY_FORM] = True
                st.session_state[NCBI_DF] = ncbi.get_data(database, search_term)
                data_processing.format_ncbi_summary()
    df_aggrid = pd.DataFrame()
    if st.session_state[NCBI_SUMMARY_FORM]:
        aggrid_table = aggrid_table()  # Initialize and Render Aggrid Table
        df_aggrid = aggrid_table["data"]
    # Streamlit UI to filter the dataset
    df_filtered, applied_filters = filter_dataframe(df_aggrid)
    if st.session_state[NCBI_DF_FILTER]:
        st.dataframe(df_filtered)
        # Streamlit UI to show the applied filters
        show_filters_info(applied_filters)
    # User has made a search query
    if not df_aggrid.empty:
        df_download = df_filtered if applied_filters else df_aggrid
        total_length = df_download["Length"].sum()
        st.write(
            f"Size of final data: {len(df_download)} rows, {total_length:,} total base pairs"
        )
        if total_length > MAX_SEQUENCE_LENGTH_DOWNLOAD:
            st.warning(
                f"The total length of the sequences is {total_length:,}, which is more than the allowed max "
                f"({MAX_SEQUENCE_LENGTH_DOWNLOAD:,}). Please apply filters to reduce the size of the dataset."
            )
        else:
            # Dropdown menu for selecting download format
            download_format = st.selectbox(
                "Select an option:", [FASTA_DOWNLOAD, CSV_DOWNLOAD]
            )
            # Option 1: Download the correct DataFrame as CSV
            if download_format == CSV_DOWNLOAD:
                st.download_button(
                    label=CSV_DOWNLOAD,
                    key="download_csv",
                    data=df_download.to_csv(index=False).encode("utf-8"),
                    file_name="sequence_search_table.csv",
                    mime="text/csv",
                    help="Click button to Download this table as a CSV file",
                    args={"as_attachment": True},
                )
            # Option 2: Download the correct DataFrame as FASTA
            elif download_format == FASTA_DOWNLOAD:
                download_data(database, df_download["Id"].unique(), search_term)
with summary_col:
    # Show the summary only if user has submitted the query form
    if st.session_state[NCBI_SUMMARY_FORM]:
        st.header("Top Organisms")
        st.write(data_processing.get_top_organisms_counts(TOP_N_ORGANISMS))

st.sidebar.image("src/app/Conservation X Labs CXL logo.png", use_column_width=True)
