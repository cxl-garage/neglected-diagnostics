import streamlit as st
from routers import ncbi

st.title("Neglected Diagnostics: Perform Genetic Testing At Scale!")

database = st.selectbox("Select the database to search", ("gene", "nucleotide"))
operation = st.selectbox("Select the operation to perform", ("Search", "Fetch"))
term = st.text_input(
    label="Enter the search term",
    placeholder="Example: human[organism] AND topoisomerase[protein name]",
)

if st.button('Perform Operation'):
    data = ncbi.get(operation, database, term)
    st.write(data)

