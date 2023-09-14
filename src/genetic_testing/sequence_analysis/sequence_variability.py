from collections import defaultdict
from typing import DefaultDict, Dict, List, Optional

import pandas as pd
from Levenshtein import distance as levenshtein_distance

from .datatypes import GroupSequenceColumns


def _preprocess(sequences: Dict[str, str], base_sequence: str) -> Dict[str, str]:
    """
    Preprocesses the given sequences by comparing them with the base sequence

    Replace the characters of the input sequences that match with characters of the base sequence with '-' character. If the characters of the input sequences do not match with the characters of the base sequence, then the characters of the input sequences are kept as it is.

    Parameters
    ----------
    sequences : Dict[str, str]
        A dictionary containing sequence IDs as keys and their corresponding sequences as values.
    base_sequence : str
        The base sequence to compare the given sequences with.

    Returns
    -------
    Dict[str, str]
        A dictionary containing sequence IDs as keys and their corresponding processed sequences as values.
    """
    processed_sequences = {}  # Dictionary to store the processed sequences
    for seq_id, seq in sequences.items():
        new_seq = "".join(
            "-" if seq_char != "-" and seq_char == base_char else seq_char
            for seq_char, base_char in zip(seq, base_sequence)
        )
        processed_sequences[seq_id] = new_seq

    return processed_sequences


def _group_sequences(
    sequences: Dict[str, str], base_sequence: str
) -> Optional[pd.DataFrame]:
    """
    Groups sequences by their similarity to a base sequence and returns a pandas DataFrame with the count of each group.

    Parameters
    ----------
    sequences : Dict[str, str]
        A dictionary containing the sequences to group, with the keys being the sequence IDs and the values being the sequences.
    base_sequence : str
        The base sequence to compare the other sequences to.

    Returns
    -------
    Optional[pd.DataFrame]
        If there are sequences to group, returns a pandas DataFrame with the count of each group, sorted by count in descending order.
        The DataFrame has four columns:
        1. SequenceID: the ID of the sequence.
        2. Sequence: the sequence that belongs to the group.
        3. Count: the number of sequences that belong to the group.
        4. Group_Number: the number of the group, starting from 1 for the most common group.
        If there are no sequences to group, returns None.
    """
    try:
        processed_sequences = _preprocess(sequences, base_sequence)

        cols = GroupSequenceColumns()
        df_sequence = pd.DataFrame(
            processed_sequences.items(), columns=[cols.id, cols.seq]
        )

        # Group by Sequence, get the count for each sequence, and sort by counts descending order
        df_sequence_counts = (
            df_sequence.groupby(cols.seq)
            .size()
            .reset_index(name=cols.count)
            .sort_values(by=cols.count, ascending=False)
        )

        # Add Group_Number column
        df_sequence_counts[cols.group_number] = df_sequence_counts.index + 1

        return df_sequence_counts

    except ValueError:
        return None


def _calculate_distances(
    curr_pattern: str, patterns_dict: Dict[str, int]
) -> DefaultDict[int, List[int]]:
    """
    Calculate the Levenshtein distance between the current sequence pattern and all other patterns in the more frequent groups of sequence data.

    Parameters
    ----------
    curr_pattern : str
        The current pattern to compare with other patterns in the input sequence data.
    patterns_dict : Dict[str, int]
        A dictionary containing all the patterns belonging to more frequent groups and their corresponding group numbers.

    Returns
    -------
    DefaultDict[int, List[int]]
        A defaultdict containing the Levenshtein distances as keys and the corresponding group numbers of patterns with that distance as values.
    """
    # Initialize a defaultdict to store distances and their corresponding group numbers
    dist = defaultdict(list)
    for key, value in patterns_dict.items():
        curr_distance = levenshtein_distance(curr_pattern, key)
        dist[curr_distance].append(value)
    return dist


def _recreate_sequences(
    df_sequence_patterns: pd.DataFrame, base_sequence: str
) -> pd.DataFrame:
    """Recreate the sequences by replacing gaps with characters from the base sequence.

    This function takes a DataFrame containing sequence patterns and replaces gaps ('-') in each pattern with characters from the corresponding positions of the base sequence.

    Parameters
    ----------
    df_sequence_patterns : pd.DataFrame
        DataFrame containing sequence patterns with columns "Sequence", "Count", and "Group Number".
    base_sequence : str
        The base sequence used to fill gaps in the patterns.

    Returns
    -------
    pd.DataFrame
        DataFrame with recreated sequences.
    """
    cols = GroupSequenceColumns()
    # Create a copy of the input DataFrame to avoid modifying the original data
    df_final_sequences = df_sequence_patterns.copy()

    # Iterate through each row in the DataFrame
    for index, row in df_final_sequences.iterrows():
        new_pattern = ""  # Initialize an empty string to store the new sequence
        curr_pattern = row[cols.seq]  # Current sequence pattern from the DataFrame

        # Iterate through each character in the current pattern
        for i in range(len(curr_pattern)):
            character = curr_pattern[i]
            if character != "-":
                # If character is not a gap, append it in uppercase to the new pattern
                new_pattern += character.upper()
            else:
                # If character is a gap, replace it with the corresponding character from the base sequence
                new_pattern += base_sequence[i]

        # Update the "Sequence" column in the DataFrame with the new pattern
        df_final_sequences.at[index, cols.seq] = new_pattern

    return df_final_sequences


def _merge_sequence_data(
    df: pd.DataFrame, base_sequence: str, min_count: int
) -> pd.DataFrame:
    """
    Merge sequences having less than `min_count` with similar sequences having a count of atleast `min_count`.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing sequence data with columns "SequenceID", "Sequence", "Count", and "Group_Number".
    base_sequence : str
        The base sequence used to recreate the final sequences in the DataFrame.
    min_count : int
        The threshold count of a sequence to be considered for grouping.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing the merged sequences with updated counts.
    """
    cols = GroupSequenceColumns()

    # Filter sequences with counts less than minimum count
    df_min_count = df[df[cols.count] < min_count]

    # Filter sequences with counts greater than or equal to minimum count
    df_max_count = df[df[cols.count] >= min_count]

    if df_max_count.empty:
        # If there are no sequences with count greater than or equal to minimum count, return empty DataFrame
        return df_max_count

    # Create a dictionary of patterns and their corresponding group numbers
    patterns_dict = {
        row[cols.seq]: row[cols.group_number] for _, row in df_max_count.iterrows()
    }

    # Iterate over sequences with counts less than minimum count
    for _, row in df_min_count.iterrows():
        # Calculate the distances between the current sequence and all other sequences in the more frequent groups
        distances_groups = _calculate_distances(row[cols.seq], patterns_dict)
        sorted_distances_groups = sorted(distances_groups.items())

        # Get the group number of the closest sequence
        _, group_numbers = sorted_distances_groups[0]
        closest_group_number = group_numbers[0]

        # Add the count of the current sequence to the group with the closest sequence
        df_max_count.loc[
            df_max_count[cols.group_number] == closest_group_number, cols.count
        ] += row[cols.count]

    # Recreate the sequences
    df_final_sequences = _recreate_sequences(df_max_count, base_sequence)
    df_final_sequences.reset_index(drop=True, inplace=True)

    return df_final_sequences


def calculate_sequence_variability(
    sequences: Dict[str, str], base_sequence: str, min_count: int = 2
) -> Optional[pd.DataFrame]:
    """Calculate sequence variability based on similarity to a base sequence and sequence counts.

    Parameters
    ----------
    sequences : Dict[str, str]
        A dictionary containing sequence IDs as keys and their corresponding sequences as values.
    base_sequence : str
        The base sequence to compare the given sequences with.
    min_count : int, optional
        The minimum count of sequences in a group for it to be not considered for further merging with more frequent groups, by default 2

    Returns
    -------
    Optional[pd.DataFrame]
        If `sequences` and `base_sequence` are provided, returns a pandas DataFrame with grouped sequences and their counts.
        If either sequences or base_sequence is empty, returns None.
    """
    # Check if sequences or base_sequence is empty
    if not sequences or not base_sequence:
        return None

    # Group sequences by similarity and count
    df_sequence_counts = _group_sequences(sequences, base_sequence)

    # Merge sequences based on count and similarity
    df_final_sequences = _merge_sequence_data(
        df_sequence_counts, base_sequence, min_count
    )

    # Return None if there are no sequences to group
    if df_final_sequences.empty:
        return None

    return df_final_sequences
