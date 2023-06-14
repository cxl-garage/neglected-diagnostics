"""Demo App for Neglected Diagnostics.

This is a demo app for neglected diagnostics built using Streamlit and Biopython. The PIs want me 
to build a quick prototype. So, I am first focusing on building a quick prototype that implements 
things end to end and then work towards improving individual modules.
"""

import streamlit as st

from genetic_testing.routers import ncbi

st.title("Neglected Diagnostics: Perform Genetic Testing At Scale!")

database = st.selectbox("Select the database to search", ("gene", "nucleotide"))
term = st.text_input(
    label="Enter the search term",
    placeholder="Example: human[organism] AND topoisomerase[protein name]")

if st.button('Perform Operation'):
    data = ncbi.search(database, term)
    st.write(data)
