import pandas as pd
from Levenshtein import distance as levenshtein_distance


from typing import Dict, Union
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
) -> Union[pd.DataFrame, None]:
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
    Union[pd.DataFrame, None]
        If there are sequences to group, returns a pandas DataFrame with the count of each group, sorted by count in descending order.
        The DataFrame has three columns:
        - Sequence: the sequence that belongs to the group.
        - Count: the number of sequences that belong to the group.
        - Group_Number: the number of the group, starting from 1 for the most common group.
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


def _calculate_distances(curr_pattern, patterns_dict):
    dist = {}
    for key, value in patterns_dict.items():
        curr_distance = levenshtein_distance(curr_pattern, key)
        if curr_distance not in dist:
            dist[curr_distance] = [value]
    return dist


def _recreate_sequences(df_sequence_patterns, base_sequence):
    df_final_sequences = df_sequence_patterns.copy()
    for index, row in df_final_sequences.iterrows():
        new_pattern = ""
        curr_pattern = row["Sequence"]
        for i in range(len(curr_pattern)):
            character = curr_pattern[i]
            if character != "-":
                new_pattern += character.upper()
            else:
                new_pattern += base_sequence[i]

        print(new_pattern)
        df_final_sequences.at[index, "Sequence"] = new_pattern

    return df_final_sequences


def _merge_sequence_data(df, base_sequence, min_count):
    # Get all the sequences less than minimum count
    df_min_count = df[df["Count"] < min_count]

    # Get all the sequences greater than equal to minimum count
    df_max_count = df[df["Count"] >= min_count]

    # Get all max count pattern
    patterns_dict = {
        row["Sequence"]: row["Group_Number"] for index, row in df_max_count.iterrows()
    }

    # Find the closest group for each of the sequences in the minimum count dataframe
    for _, row in df_min_count.iterrows():
        print("Current Sequence: ", row["Sequence"])
        distances_groups = _calculate_distances(row["Sequence"], patterns_dict)
        sorted_distances_groups = sorted(distances_groups.items())
        group_number = sorted_distances_groups[0][1][0]  # TODO: Improve this line
        print("Current Group Number: ", group_number)

        df_max_count.loc[df_max_count["Group_Number"] == group_number, "Count"] += row[
            "Count"
        ]

        print(df_max_count[df_max_count["Group_Number"] == group_number])
        print("\n\n")

    # Recreate the sequences
    df_final_sequences = _recreate_sequences(df_max_count, base_sequence)

    return df_final_sequences


def calculate_sequence_variability(sequences, base_sequence, min_count=2):
    if not sequences or not base_sequence:
        return None

    df_sequence_counts = _group_sequences(sequences, base_sequence)
    df_final_sequences = _merge_sequence_data(
        df_sequence_counts, base_sequence, min_count
    )

    return df_final_sequences
