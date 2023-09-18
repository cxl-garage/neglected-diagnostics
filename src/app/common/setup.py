import pandas as pd
import streamlit as st

from app.common.constants import (
    NCBI_DF,
    NCBI_DF_FILTER,
    NCBI_SUMMARY_FORM,
    PAGE_LAYOUT,
    SEQVAR_CALC_BTN,
    SEQVAR_CALCULATED,
    SEQVAR_DF,
)


def _initialize_UI() -> None:
    """Initialize the user interface (UI) settings for the Streamlit app."""
    st.set_page_config(layout=PAGE_LAYOUT)  # "wide" layout uses the entire screen


def _initialize_session_state() -> None:
    """Initialize the Streamlit session state.

    This function initializes the Streamlit session state by checking the existence of specific keys for storing data and state variables. If any of these keys are not present in the session state, the function creates them and sets their initial values.
    """
    if NCBI_SUMMARY_FORM not in st.session_state:
        st.session_state[NCBI_SUMMARY_FORM] = False
        st.session_state[NCBI_DF_FILTER] = False

    if NCBI_DF not in st.session_state:
        st.session_state[NCBI_DF] = pd.DataFrame()


def initialize_session_state_seq_var() -> None:
    """Initialize the Streamlit session state for the Sequence Variability page."""
    if SEQVAR_CALC_BTN not in st.session_state:
        st.session_state[SEQVAR_CALC_BTN] = False
        st.session_state[SEQVAR_CALCULATED] = False
        st.session_state[SEQVAR_DF] = pd.DataFrame()


def initialize() -> None:
    """Initialize the Streamlit application

    This function sets up the necessary UI configurations and session states for the application.
    """
    _initialize_UI()
    _initialize_session_state()
