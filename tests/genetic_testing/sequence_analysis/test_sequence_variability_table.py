import pandas as pd
import pytest

from genetic_testing.sequence_analysis.sequence_variability import (
    calculate_seq_variability_table,
)

SEQUENCES = "sequences"
EXPECTED_RESULT = "expected_result"

# Define expected results
EXPECTED_RESULTS = {
    "single_sequence": (
        ["AGCT"],
        pd.DataFrame(
            {
                "A": [100.0, 0.0, 0.0, 0.0],
                "G": [0.0, 100.0, 0.0, 0.0],
                "C": [0.0, 0.0, 100.0, 0.0],
                "T": [0.0, 0.0, 0.0, 100.0],
            },
            index=[1, 2, 3, 4],
        ).rename_axis("Position"),
    ),
    "single_sequence_with_gap": (
        ["AGC-"],
        pd.DataFrame(
            {
                "A": [100.0, 0.0, 0.0, 0.0],
                "G": [0.0, 100.0, 0.0, 0.0],
                "C": [0.0, 0.0, 100.0, 0.0],
            },
            index=[1, 2, 3, 4],
        ).rename_axis("Position"),
    ),
    "multiple_sequences": (
        ["AGCT", "ACGT"],
        pd.DataFrame(
            {
                "A": [100.0, 0.0, 0.0, 0.0],
                "G": [0.0, 50.0, 50.0, 0.0],
                "C": [0.0, 50.0, 50.0, 0.0],
                "T": [0.0, 0.0, 0.0, 100.0],
            },
            index=[1, 2, 3, 4],
        ).rename_axis("Position"),
    ),
    "multiple_sequence_with_gaps": (
        ["AGCT", "ACG-"],
        pd.DataFrame(
            {
                "A": [100.0, 0.0, 0.0, 0.0],
                "G": [0.0, 50.0, 50.0, 0.0],
                "C": [0.0, 50.0, 50.0, 0.0],
                "T": [0.0, 0.0, 0.0, 100.0],
            },
            index=[1, 2, 3, 4],
        ).rename_axis("Position"),
    ),
}


@pytest.mark.parametrize(
    [SEQUENCES, EXPECTED_RESULT],
    list(EXPECTED_RESULTS.values()),
    ids=list(EXPECTED_RESULTS.keys()),
)
def test_calculate_sequence_variability_table(sequences, expected_result):
    result = calculate_seq_variability_table(sequences)
    pd.testing.assert_frame_equal(result, expected_result)
