from io import BytesIO

import pandas as pd
import pytest

# test script to test the functionality of consensusSeq.py
from genetic_testing.sequence_tool.consensusSeq import ConsensusSeq, alignMs

MULTIFASTA = "multifasta"
EXPECTED_RESULT = "expected_result"

# Define expected results
EXPECTED_RESULTS = {}

# define test data
# inFile = "Sei whale_CytB_unaligned"
inFile = "Sei whale_CytB_unaligned"

# test data for reading fasta file
consensusSeq_test_data = f"./tests/genetic_testing/sequence_tool/{inFile}.fasta"
with open(consensusSeq_test_data, "rb") as file:
    sequences = BytesIO(file.read())

# testConsensusSeq = ConsensusSeq(sequences)

# test data for alignment calculation
print(alignMs(sequences))

# # test data for consensus score calculation
# testConsensusSeq.compute_consensus_sequence()
# print(testConsensusSeq.consensus_seq)

# # test data for gapless consensus sequence calculation
# testConsensusSeq.compute_consensus_seq_gapless()
# print(testConsensusSeq.consensus_seq_gapless)

# # test data for ideal consensus region calculation
# testConsensusSeq.compute_consensus_regions()
# print(testConsensusSeq.consensus_regions)
