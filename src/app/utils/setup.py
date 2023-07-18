import pandas as pd
import streamlit as st
from app.utils.constants import NCBI_SUMMARY_FORM, NCBI_DF, NCBI_DF_FILTER 

def initialize():
    if NCBI_SUMMARY_FORM not in st.session_state:
        st.session_state[NCBI_SUMMARY_FORM] = False
        st.session_state[NCBI_DF_FILTER] = False

    if NCBI_DF not in st.session_state:
        st.session_state[NCBI_DF] = pd.DataFrame()
        
