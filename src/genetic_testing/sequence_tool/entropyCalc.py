# import jpype
# import jpype.imports
# from jpype.types import *

# from jpype import JImplements, JOverride, JImplementationFor


# # Launch the JVM
# jpype.startJVM(classpath=['./src/genetic_testing/sequence_tool/jars/*'])
import uuid
from functools import reduce

import snoop

from genetic_testing.sequence_tool.entropyCalcApi import *
from genetic_testing.sequence_tool.entropyCalcDataClass import *

uUid = str(uuid.uuid4())
multifasta = "./tests/genetic_testing/sequence_tool/COI_Atlantic Salmon_trimmed sequences_reduced.fasta"


# File input, compute entropy and group sequences
# @snoop
def seqInput(multifasta, isUseAll, mergedIntervals):
    try:
        config = Config()
        config.inFile = multifasta
        config.useDouble = True
        config.intervals = mergedIntervals
        config.wholeSequence = isUseAll
        config.threshold = 100

        apiResult = Api.group(config)
        apiResult.userSequence = config.useSequences

    except Exception as e:
        print(e)
        raise e

    seqResult = SeqBigResultModel()
    resultStore = ResultStore()

    seqResult.time = apiResult.time
    seqResult.seqLength = apiResult.positions
    seqResult.seqTotal = apiResult.sequencesTotal
    seqResult.groups = apiResult.groups
    seqResult.intervals = apiResult.intervals

    resultStore.add_result(uUid, seqResult)

    return resultStore


# Informatives, Compress and Stratify file output


def calcInformativeSeqs(result) -> list:
    groups = result.groups
    # Per result.groups, check group type equals NORMAL, if so, add 1 to groups_length
    groups_length = reduce(
        lambda acc, g: acc + 1 if g.type == GroupedSequenceType.NORMAL else acc,
        groups,
        0,
    )

    # Sort groups by priority, then by length of sequences
    result.groups.sort(
        key=lambda g: g.get_priority() * 1e10
        - g.get_priority() * 1e10
        + len(g.sequences)
        - len(g.sequences)
    )

    informative_seqs = []

    for i, g in enumerate(result.groups):
        g.number = g.name = i + 1
        g.freq = len(g.sequences)
        g.probability = g.freq / result.seqTotal * 100

        if g.get_priority() != 0:
            g.name = g.type.name.lower()
            g.color = "red"

        for s in g.sequences:
            s.group = g

        if g.type == GroupedSequenceType.NORMAL:
            informative_seqs += g.sequences

    return informative_seqs


def exportCompressSeqs(result):
    resultText = ""
    for g in result.groups:
        resultText += f">Group-{g.name}|frequency-{round(g.probability, 3)}%\r\n{g.sequences[0].data}\r\n"

    # export to 'compressed.txt'
    with open("compressed.txt", "w") as f:
        f.write(resultText)


def export(seqs, outFile="export.txt"):
    resultText = ""
    for s in seqs:
        resultText += f"{s.meta}{s.data}\n"

    # export to outFile
    with open(outFile, "w") as f:
        f.write(resultText)


if __name__ == "__main__":
    resultStore = seqInput(multifasta, True, None)

    infoSeqs = calcInformativeSeqs(resultStore.get_result(uUid))
    export(infoSeqs, "informatives.txt")
    exportCompressSeqs(resultStore.get_result(uUid))
