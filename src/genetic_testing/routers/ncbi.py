"""Module to search and download sequence data from the NCBI datasource. 

This module is still work in progress. Currently, I am working on implementing the search 
functionality
"""

import os
from typing import List

import streamlit as st
from Bio import Entrez

ENTREZ_EMAIL = "ENTREZ_EMAIL"

if ENTREZ_EMAIL in st.secrets:
    Entrez.email = st.secrets[ENTREZ_EMAIL]
else:
    Entrez.email = os.environ.get(ENTREZ_EMAIL, "")


def search(database: str, term: str) -> List[int]:
    """Get the List of UIDs matching the Entrez query.

    This function searches the provided NCBI database and returns the list of UIDs matching the
    Entrez query.

    Args:
        database (str): Name of the Entrez database.
        term (str): Entrez text query.

    Returns:
        List[int]: A list of UIDs.
    """
    handle = Entrez.esearch(db=database, term=term)
    uids = []
    if handle:
        record = Entrez.read(handle)
        uids = record["IdList"]

    return uids
