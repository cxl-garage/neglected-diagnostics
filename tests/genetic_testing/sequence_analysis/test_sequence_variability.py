import pandas as pd
import pytest

from genetic_testing.sequence_analysis.datatypes import GroupSequenceColumns
from genetic_testing.sequence_analysis.sequence_variability import (
    calculate_sequence_variability,
)

cols = GroupSequenceColumns()

SEQUENCES = "sequences"
BASE_SEQUENCE = "base_sequence"
MIN_COUNT = "min_count"
EXPECTED_RESULT = "expected_result"

# Define expected results
EXPECTED_RESULTS = {
    "empty_sequences": ({}, "ACGT", 2, None),
    "empty_base_sequence": ({"seq1": "ACGT"}, "", 2, None),
    "single_sequence": (
        {"seq1": "ACGT"},
        "ACGT",
        1,
        pd.DataFrame(
            {"Sequence": ["ACGT"], "Count": [1], "Group_Number": [1]}
        ),  # Test case for single sequence
    ),
    "multiple_sequence_same_group": (
        {"seq1": "-CGT", "seq2": "ACGT"},
        "ACGT",
        1,
        pd.DataFrame(
            {"Sequence": ["ACGT"], "Count": [2], "Group_Number": [1]}
        ),  # Test case for multiple sequences belonging to same group
    ),
    "multiple_sequence_different_group": (
        {
            "seq1": "-CGT",
            "seq2": "ACGT",
            "seq3": "TGCA",
            "seq4": "TGCA",
            "seq5": "TGC-",
        },
        "ACGT",
        2,
        pd.DataFrame(
            {"Sequence": ["ACGT", "TGCA"], "Count": [2, 3], "Group_Number": [1, 2]},
        ),  # Test case for multiple sequences belonging to different groups
    ),
}


@pytest.mark.parametrize(
    [SEQUENCES, BASE_SEQUENCE, MIN_COUNT, EXPECTED_RESULT],
    list(EXPECTED_RESULTS.values()),
    ids=list(EXPECTED_RESULTS.keys()),
)
def test_calculate_sequence_variability(
    sequences, base_sequence, min_count, expected_result
):
    result = calculate_sequence_variability(sequences, base_sequence, min_count)
    if expected_result is None:
        assert result is None
    else:
        pd.testing.assert_frame_equal(
            result[["Sequence", "Count"]], expected_result[["Sequence", "Count"]]
        )
