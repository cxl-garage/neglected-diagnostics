from io import BytesIO

import pandas as pd
import pytest

# test script to test the functionality of consensusSeq.py
from genetic_testing.sequence_tool.consensusSeq import ConsensusSeq, alignMs

MULTIFASTA = "multifasta"
EXPECTED_RESULT = "expected_result"


# define test data
inFile = "Sei whale_Dloop_408"

# test data for reading fasta file
consensusSeq_test_data = f"./tests/genetic_testing/sequence_tool/{inFile}.fasta"
file = open(consensusSeq_test_data, "rb")
sequences = BytesIO(file.read())
result = "TATATTGTACAATAACTACAAGGCCACAGTGTCATGTCCGTATCGAAAAATAACTTGTCCTATCACATATTATTATGTGATTTGTACATATATACACCTCCCACAACTTAATTAATAGTCTTCTCCTGTAGGTATGTATATAAATACATGCTATGTATAACTGTGCATTCAATTATTTTCACTACGAGCAGTTGAAGCTCGTATTAAATTTCATTAATTTTACATATTACATAATATTTATTAATAGTACATTAGCGCATGTTATTATGCATCCCCTGGTAAATTCTATTCAAATGATTCCTATGGCCGCTCCATTAGATCACGAGCTTAACCACCATGCCGCGTGAAACCAGCAACCCGCTTGGCAGGGATCCCTCTTCTCGCACCGGGCCCATCAATCGTGGGGGT"

# Define expected results
EXPECTED_RESULTS = {
    "test_multifasta": (sequences, result),
}


@pytest.mark.parametrize(
    [MULTIFASTA, EXPECTED_RESULT],
    list(EXPECTED_RESULTS.values()),
    ids=list(EXPECTED_RESULTS.keys()),
)
def test_consensusSeq(multifasta, expected_result):
    result = ConsensusSeq(multifasta).consensus_seq_gapless
    assert result == expected_result
