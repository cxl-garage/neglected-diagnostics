import pandas as pd
import pytest

# test script to test the functionality of consensusSeq.py
from genetic_testing.sequence_tool.consensusSeq import ConsensusSeq

MULTIFASTA = "multifasta"
EXPECTED_RESULT = "expected_result"

# Define expected results
EXPECTED_RESULTS = {}

# define test data
# test data for reading fasta file
consensusSeq_test_data = "Sei whale_Dloop_408.fasta"
testConsensusSeq = ConsensusSeq(consensusSeq_test_data)

# test data for consensus score calculation
testConsensusSeq.compute_consensus_sequence()
print(testConsensusSeq.consensus_seq)

# test data for gapless consensus sequence calculation
testConsensusSeq.compute_consensus_seq_gapless()
print(testConsensusSeq.consensus_seq_gapless)

# test data for ideal consensus region calculation
testConsensusSeq.compute_consensus_regions()
print(testConsensusSeq.consensus_regions)
