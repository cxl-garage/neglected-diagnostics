from typing import List

import pandas as pd
import streamlit as st

from app.common.data_processing import msa_df_to_fasta, seqvar_df_to_fasta
from genetic_testing.routers import ncbi


def download_data(database: str, ids: List[int], search_term: str) -> None:
    """
    Download sequence data from NCBI's database using the given UIDs.
    Parameters
    ----------
    database : str
        The name of the NCBI database from which data will be fetched.
    ids : List[int]
        A list of UIDs representing the entries to be fetched from the database.
    search_term : str
        The search term used to filter the data to be fetched.
    Returns
    -------
    None
        This function does not return any value. The fetched data is directly written to the Streamlit app.
    Notes
    -----
    This function uses the `fetch_data` function from the `genetic_testing.routers.ncbi` module to retrieve data
    from the specified NCBI database based on the provided list of IDs. If the `ids` list is empty, the function writes "Filtered Data is empty" to the Streamlit app. Otherwise, it fetches the data, converts it to a string, and writes "Fetched Data Successfully!" to the Streamlit app, along with a download button for the data. The downloaded file will be named "data.fasta" by default.
    """
    if len(ids) == 0:
        st.write("Filtered Data is empty")
    else:
        if st.button("Prepare Download"):
            file_buffer = ncbi.fetch_data(database, ids)
            file_str = file_buffer.getvalue()
            if len(file_str) == 0:
                st.write(
                    f"NCBI's '{database}' database has no data for '{search_term}' search term"
                )
            else:
                st.write("Fetched Data Successfully!")
                st.download_button(
                    label="Download",
                    data=file_str,
                    file_name=f"{database}_{search_term}.fasta",
                )


def download_seq_var_data(df: pd.DataFrame, species: str) -> None:
    """Downloads sequence variability data as a FASTA file.
    This function converts a DataFrame containing sequence variability data to
    FASTA format and provides the option to
    download the resulting FASTA file.
    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing sequence variability data.
    species : str
        The name of the species for which the data is being processed.
    """
    # Convert the DataFrame to FASTA format string
    fasta_str = seqvar_df_to_fasta(df, species)
    st.write(
        "FASTA file is ready for download. Please click the button below to proceed."
    )
    # Display a download button
    st.download_button(
        label="Download FASTA file",
        data=fasta_str,
        file_name=f"sequence_haplotypes_{species}.fasta",
    )


def download_msa_data(df: pd.DataFrame, species: str) -> None:
    """Downloads multisequence alignment data as a FASTA file.
    This function converts a DataFrame containing multisequence alignment data to
    FASTA format and provides the option to
    download the resulting FASTA file.
    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing sequence variability data.
    species : str
        The name of the species for which the data is being processed.
    """

    # Convert the DataFrame to FASTA format string
    fasta_str = msa_df_to_fasta(df, species)
    st.write(
        "FASTA file is ready for download. Please click the button below to proceed."
    )
    # Display a download button
    st.download_button(
        label="Download FASTA file",
        data=fasta_str,
        file_name=f"{species}_aligned.fasta",
    )
