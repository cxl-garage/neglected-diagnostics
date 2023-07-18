import streamlit as st
from app.utils.constants import NCBI_DF


def format_ncbi_summary():
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
