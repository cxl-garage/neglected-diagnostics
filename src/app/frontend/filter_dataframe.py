from pandas.api.types import (
    is_categorical_dtype,
    is_datetime64_any_dtype,
    is_numeric_dtype,
)
import pandas as pd
import streamlit as st

from app.utils.constants import NCBI_DF_FILTER

def _preprocessing(df: pd.DataFrame):
    """
    Preprocesses the NCBI dataframe before filtering

    Parameters
    ----------
    df : pd.DataFrame
        Original dataframe

    Returns
    -------
    pd.DataFrame
        Processed dataframe
    """    
    # Convert columns to categorical
    df["Species"] = df["Species"].astype("category")
    df["TaxId"] = df["TaxId"].astype("category")
    df["Id"] = df["Id"].astype("category")
    df["Gi"] = df["Gi"].astype("category")
    df["Status"] = df["Status"].astype("category")
    
    # Convert Date objects to Datetime format
    date_strings = [str(date_element) for date_element in df["CreateDate"]]
    df["CreateDate"] = pd.to_datetime(date_strings, format="%Y/%m/%d")

    date_strings = [str(date_element) for date_element in df["UpdateDate"]]
    df["UpdateDate"] = pd.to_datetime(date_strings, format="%Y/%m/%d")

    # Try to convert datetimes into a standard format (datetime, no timezone)
    for col in df.columns:
        if is_datetime64_any_dtype(df[col]):
            df[col] = df[col].dt.tz_localize(None)

def filter_dataframe(df: pd.DataFrame) -> pd.DataFrame:
    """
    Adds a UI on top of a dataframe to let viewers filter columns

    Parameters
    ----------
    df : pd.DataFrame
        Original dataframe

    Returns
    -------
    pd.DataFrame
        Filtered dataframe
    """    
    st.session_state[NCBI_DF_FILTER] = st.checkbox("Add filters")

    if st.session_state[NCBI_DF_FILTER] is False:
        return df

    df = df.copy()
    _preprocessing(df)

    filter_container = st.container()
    with filter_container:
        to_filter_columns = st.multiselect("Filter dataframe on", df.columns)
        for column in to_filter_columns:
            left, right = st.columns((1, 20))
            left.write("↳")
            
            if is_categorical_dtype(df[column]):
                user_cat_input = right.multiselect(
                    f"Values for {column}",
                    df[column].unique(),
                )
                df = df[df[column].isin(user_cat_input)]
            elif is_numeric_dtype(df[column]):
                with right:
                    # Create a sub column layout
                    sub_columns = st.columns(2)
                    # Add input boxes for minimum and maximum values in the row
                    with sub_columns[0]:
                        min_length = st.number_input("Enter Minimum Length")

                    with sub_columns[1]:
                        max_length = st.number_input("Enter Maximum Length")
                df = df[df[column].between(min_length, max_length)]
            elif is_datetime64_any_dtype(df[column]):
                user_date_input = right.date_input(
                    f"Values for {column}",
                    value=(
                        df[column].min(),
                        df[column].max(),
                    ),
                )
                if len(user_date_input) == 2:
                    user_date_input = tuple(map(pd.to_datetime, user_date_input))
                    start_date, end_date = user_date_input
                    df = df.loc[df[column].between(start_date, end_date)]
            else:
                # TODO: Add support for Contains, Not Contains
                user_text_input = right.text_input(
                    f"Substring or regex in {column}",
                )
                if user_text_input:
                    df = df[df[column].astype(str).str.contains(user_text_input)]

    return df
