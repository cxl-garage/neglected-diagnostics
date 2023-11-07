from io import BytesIO
from typing import Dict, List, Optional

import pandas as pd
import streamlit as st
from streamlit.runtime.uploaded_file_manager import UploadedFile

from app.common.constants import NCBI_DF
from genetic_testing.sequence_analysis.datatypes import GroupSequenceColumns


def format_ncbi_summary() -> None:
    """Format and restructure the NCBI summary DataFrame.

    This function performs formatting and restructuring tasks on the NCBI summary DataFrame, which is assumed to be
    stored in the Streamlit session state under the key specified by the constant `NCBI_DF`. The function performs
    the following actions:

    1. Extract the "Species" information from the "Title" column and store it in a new "Species" column.
    2. Convert the "Species" column to a categorical data type.
    3. Reorder the columns of the DataFrame to a predefined order.
    """
    # Extract the species
    st.session_state[NCBI_DF]["Species"] = st.session_state[NCBI_DF]["Title"].apply(
        lambda x: " ".join(x.split()[:2])
    )
    # Convert the species to categorical
    st.session_state[NCBI_DF]["Species"] = st.session_state[NCBI_DF]["Species"].astype(
        "category"
    )

    # Reorder the columns
    column_order = [
        "Species",
        "Title",
        "TaxId",
        "Id",
        "Length",
        "Gi",
        "CreateDate",
        "UpdateDate",
        "Status",
    ]
    st.session_state[NCBI_DF] = st.session_state[NCBI_DF].reindex(columns=column_order)


def get_top_organisms_counts(top_n=15) -> str:
    """Get the top N organisms and their counts in the session's NCBI DataFrame.

    Parameters
    ----------
    top_n : int, optional
        Number of top organisms to show, by default 15

    Returns
    -------
    str
        A formatted string containing the top N organisms and their counts, separated by two newlines ("\n\n").
    """
    # Calculate counts of species and sort by counts in descending order
    species_counts = (
        st.session_state[NCBI_DF]["Species"]
        .value_counts()
        .sort_values(ascending=False)
        .head(top_n)
    )

    # Create the desired output format
    species_counts_list = [
        f"{value} ({count})"
        for value, count in zip(species_counts.index, species_counts)
    ]

    return "\n\n".join(species_counts_list)


def read_fasta(file: BytesIO) -> Dict[str, str]:
    """
    Read a FASTA file and return a dictionary of sequences and headers.

    Parameters
    ----------
    file : BytesIO
        A BytesIO object containing the FASTA file data.

    Returns
    -------
    Dict[str, str]
        A dictionary containing the sequences and headers from the FASTA file.
    """
    sequences = {}  # Dictionary to store sequences and headers
    current_header = None
    current_sequence = []

    with file as f:
        for line in f:
            line = line.strip()  # Remove leading/trailing whitespace

            if line.startswith(b">"):  # Header line (as bytes)
                if current_header is not None:  # Store previous sequence
                    sequences[current_header] = "".join(current_sequence)
                current_header = line[1:].decode("utf-8")  # Decode bytes to string
                current_sequence = []  # Reset sequence list
            else:
                current_sequence.append(line.decode("utf-8"))  # Decode bytes to str

        if current_header is not None:  # Store the last sequence
            sequences[current_header] = "".join(current_sequence)

    return sequences


def get_headers(files: List[UploadedFile]) -> List[str]:
    """
    Read the FASTA files in the given list and return a list of headers.

    Parameters
    ----------
    files : List[str]
        A list of file paths to FASTA files.

    Returns
    -------
    List[str]
        A list of headers from the given FASTA files.
    """
    headers = []
    for file in files:
        sequences = read_fasta(BytesIO(file.read()))
        headers.extend(sequences.keys())
    return sorted(headers)


def seqvar_df_to_fasta(df: pd.DataFrame, species: str = "species") -> str:
    """Convert a sequence variability DataFrame to FASTA format string.

    Parameters
    ----------
    df : pd.DataFrame
        A DataFrame containing sequence variability data. It is expected to have `Sequence`, `Group_Number`,
        and `Count` columns.

    species : str, optional
        The name of the species, by default `species`.

    Returns
    -------
    str
        A string containing the sequence variability data in FASTA format.
    """
    cols = GroupSequenceColumns()
    fasta_str = []
    for _, row in df.iterrows():
        group_number = int(row[cols.group_number])
        count = int(row[cols.count])
        sequence = row[cols.seq]

        fasta_str.append(f">Group{group_number}_{species}_{count}\n{sequence}")

    return "\n".join(fasta_str)
    return "\n".join(fasta_str)
