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
    "empty_sequences": {
        SEQUENCES: {},
        BASE_SEQUENCE: "ACGT",
        MIN_COUNT: 2,
        EXPECTED_RESULT: None,
    },
    "empty_base_sequence": {
        SEQUENCES: {"seq1": "ACGT"},
        BASE_SEQUENCE: "",
        MIN_COUNT: 2,
        EXPECTED_RESULT: None,
    },
    "single_sequence": {
        SEQUENCES: {"seq1": "ACGT"},
        BASE_SEQUENCE: "ACGT",
        MIN_COUNT: 1,
        EXPECTED_RESULT: pd.DataFrame(
            {"Sequence": ["ACGT"], "Count": [1], "Group_Number": [1]}
        ),  # Test case for single sequence
    },
    "multiple_sequence_same_group": {
        SEQUENCES: {"seq1": "-CGT", "seq2": "ACGT"},
        BASE_SEQUENCE: "ACGT",
        MIN_COUNT: 1,
        EXPECTED_RESULT: pd.DataFrame(
            {"Sequence": ["ACGT"], "Count": [2], "Group_Number": [1]}
        ),  # Test case for multiple sequences belonging to same group
    },
    "multiple_sequence_different_group": {
        SEQUENCES: {
            "seq1": "-CGT",
            "seq2": "ACGT",
            "seq3": "TGCA",
            "seq4": "TGCA",
            "seq5": "TGC-",
        },
        BASE_SEQUENCE: "ACGT",
        MIN_COUNT: 2,
        EXPECTED_RESULT: pd.DataFrame(
            {"Sequence": ["ACGT", "TGCA"], "Count": [2, 3], "Group_Number": [1, 2]},
        ),  # Test case for multiple sequences belonging to different groups
    },
}


@pytest.mark.parametrize(
    "sequences, base_sequence, min_count, expected_result",
    [
        (
            data["sequences"],
            data["base_sequence"],
            data["min_count"],
            data["expected_result"],
        )
        for data in EXPECTED_RESULTS.values()
    ],
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
