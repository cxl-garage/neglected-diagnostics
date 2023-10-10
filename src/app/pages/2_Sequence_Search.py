"""Demo App for Neglected Diagnostics.
This is a demo app for neglected diagnostics built using Streamlit and Biopython. The PIs want me to build a quick prototype. So, I am first focusing on building a quick prototype that implements things end to end and then work towards improving individual modules.
"""


import pandas as pd
import streamlit as st

from app.common import data_processing, setup
from app.common.constants import (
    MAIN_PAGE_COLS_GAP,
    MAIN_PAGE_COLS_SIZE,
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
    st.header("Neglected Diagnostics: Democratizing Genetic Testing!")
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
        st.write(f"Size of final data: {len(df_filtered)}")
        # Streamlit UI to show the applied filters
        show_filters_info(applied_filters)
    # User has made a search query
    if not df_aggrid.empty:
        # Dropdown menu for selecting download format
        download_format = st.selectbox(
            "Select an option:", [CSV_DOWNLOAD, FASTA_DOWNLOAD]
        )
        # Option 1: Download the correct DataFrame as CSV
        if download_format == CSV_DOWNLOAD:
            # User has applied filters using `Add filters` checkbox
            if applied_filters:
                st.download_button(
                    label="Download table as CSV",
                    data=df_filtered.to_csv(index=False).encode("utf-8"),
                    file_name="sequence_search_table.csv",
                    mime="text/csv",
                    key="df_filtered_csv_button",
                    help="Click button to Download this table as a CSV file",
                    args={"as_attachment": True},
                )
            # User has not applied filters using `Add filters` checkbox
            else:
                st.download_button(
                    label="Download table as CSV",
                    data=df_aggrid.to_csv(index=False).encode("utf-8"),
                    file_name="sequence_search_table.csv",
                    mime="text/csv",
                    key="df_aggrid_csv_button",
                    help="Click button to Download this table as a CSV file",
                    args={"as_attachment": True},
                )
        # Option 2: Download the correct DataFrame as FASTA
        elif download_format == FASTA_DOWNLOAD:
            # User has applied filters using `Add filters` checkbox
            if applied_filters:
                download_data(database, df_filtered["Id"].unique(), search_term)
            # User has not applied filters using `Add filters` checkbox
            else:
                download_data(database, df_aggrid["Id"].unique(), search_term)
with summary_col:
    # Show the summary only if user has submitted the query form
    if st.session_state[NCBI_SUMMARY_FORM]:
        st.header("Top Organisms")
        st.write(data_processing.get_top_organisms_counts(TOP_N_ORGANISMS))
