import pandas as pd
import streamlit as st

from app.common.constants import (
    CONVERSATION_HISTORY,
    GALAXY_CONNECTION,
    GALAXY_TOOLS,
    LLM_CONFIG,
    LLM_PROVIDER,
    MCP_CONNECTION,
    MCP_SERVER_PROCESS,
    MSA_CLEAN_ZIP,
    MSA_DF,
    MSA_FORM,
    MSA_VIEWER,
    NCBI_DF,
    NCBI_DF_FILTER,
    NCBI_SUMMARY_FORM,
    PAGE_LAYOUT,
    SEQVAR_DF,
    SEQVAR_TABLE,
    SEQVAR_TABLE_BTN,
    SEQVAR_TABLE_CALC,
    TGT_AREA_DF,
    TGT_AREA_FORM,
    WORKFLOW_PLAN,
)


def _initialize_UI() -> None:
    """Initialize the user interface (UI) settings for the Streamlit app."""
    st.set_page_config(layout=PAGE_LAYOUT)  # "wide" layout uses the entire screen


def _initialize_session_state() -> None:
    """Initialize the Streamlit session state.

    This function initializes the Streamlit session state by checking the existence of specific keys for storing data
    and state variables. If any of these keys are not present in the session state, the function creates them and sets
    their initial values.
    """
    if NCBI_SUMMARY_FORM not in st.session_state:
        st.session_state[NCBI_SUMMARY_FORM] = False
        st.session_state[NCBI_DF_FILTER] = False

    if NCBI_DF not in st.session_state:
        st.session_state[NCBI_DF] = pd.DataFrame()


def init_session_state_seq_var() -> None:
    """Initialize the Streamlit session state for the Sequence Variability page."""
    if SEQVAR_DF not in st.session_state:
        st.session_state[SEQVAR_DF] = pd.DataFrame()


def init_session_state_seq_var_table() -> None:
    """Initialize the Streamlit session state for the Sequence Variability Table."""
    if SEQVAR_TABLE_BTN not in st.session_state:
        st.session_state[SEQVAR_TABLE_BTN] = False

    if SEQVAR_TABLE_CALC not in st.session_state:
        st.session_state[SEQVAR_TABLE_CALC] = False

    if SEQVAR_TABLE not in st.session_state:
        st.session_state[SEQVAR_TABLE] = pd.DataFrame()


def init_session_state_msa() -> None:
    """Initialize the Streamlit session state for the Multisequence Alignment page."""
    if MSA_FORM not in st.session_state:
        st.session_state[MSA_FORM] = False

    if MSA_DF not in st.session_state:
        st.session_state[MSA_DF] = pd.DataFrame()

    if MSA_VIEWER not in st.session_state:
        st.session_state[MSA_VIEWER] = None

    if MSA_CLEAN_ZIP not in st.session_state:
        st.session_state[MSA_CLEAN_ZIP] = None


def init_session_state_tgt_area() -> None:
    """Initialize the Streamlit session state for the Target Area page."""
    if TGT_AREA_FORM not in st.session_state:
        st.session_state[TGT_AREA_FORM] = False

    if TGT_AREA_DF not in st.session_state:
        st.session_state[TGT_AREA_DF] = pd.DataFrame()


def init_session_state_galaxy_assistant() -> None:
    """Initialize the Streamlit session state for Galaxy Workflow Assistant."""
    if GALAXY_CONNECTION not in st.session_state:
        st.session_state[GALAXY_CONNECTION] = None
    if GALAXY_TOOLS not in st.session_state:
        st.session_state[GALAXY_TOOLS] = []
    if CONVERSATION_HISTORY not in st.session_state:
        st.session_state[CONVERSATION_HISTORY] = []
    if WORKFLOW_PLAN not in st.session_state:
        st.session_state[WORKFLOW_PLAN] = None
    if LLM_PROVIDER not in st.session_state:
        st.session_state[LLM_PROVIDER] = "openai"
    if LLM_CONFIG not in st.session_state:
        st.session_state[LLM_CONFIG] = {}
    if MCP_SERVER_PROCESS not in st.session_state:
        st.session_state[MCP_SERVER_PROCESS] = None
    if MCP_CONNECTION not in st.session_state:
        st.session_state[MCP_CONNECTION] = None


def initialize() -> None:
    """Initialize the Streamlit app by setting up the UI and session state."""
    _initialize_UI()
    _initialize_session_state()
