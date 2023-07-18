import streamlit as st
from typing import Dict, DefaultDict
from st_aggrid import (
    AgGrid,
    GridOptionsBuilder,
    GridUpdateMode,
    DataReturnMode,
    ColumnsAutoSizeMode,
)
from app.utils.constants import NCBI_DF

DEFAULT_COL_DEF = {
    "filter": True,
    "resizable": True,
    "sortable": True,
    "editable": True,
}


def _initialize_grid_options() -> DefaultDict:
    """_summary_

    Returns
    -------
    DefaultDict
        _description_
    """
    options_builder = GridOptionsBuilder.from_dataframe(st.session_state[NCBI_DF])
    # options_builder.configure_default_column(
    #     resizable=True, filterable=True, sortable=True, editable=True
    # )
    options_builder.configure_pagination(
        paginationAutoPageSize=False, paginationPageSize=10
    )
    options_builder.configure_default_column(**DEFAULT_COL_DEF)
    grid_options = options_builder.build()

    return grid_options


def aggrid_table() -> Dict:
    """_summary_

    Returns
    -------
    Dict
        _description_
    """
    grid_options = _initialize_grid_options()
    grid_table = AgGrid(
        st.session_state[NCBI_DF],
        gridOptions=grid_options,
        columns_auto_size_mode=ColumnsAutoSizeMode.FIT_CONTENTS,
        update_mode=GridUpdateMode.MODEL_CHANGED,
        data_return_mode=DataReturnMode.FILTERED,
        enable_enterprise_modules=False,
    )

    return grid_table
