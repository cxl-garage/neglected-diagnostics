import datetime
import io
import logging
import math
import threading
import time
from abc import ABC, abstractmethod
from collections import defaultdict
from enum import Enum

import Apfloat

logger = logging.getLogger(__name__)


# Essentials data structures
class Config:
    def __init__(self):
        self.inFile = None
        self.outFile = None
        self.precision = None
        self.useDouble = None
        self.execServcie = None
        self.useSequences = None
        self.intervals = None
        self.wholeSequence = None


class State:
    def __init__(self, config, result):
        self.config = config
        self.result = result
        self.factory = Factory(self)
        self.start = None


class Factory:
    def __init__(self, state):
        self.state = state

    def getCalculations(self, pos, a, g, c, u, gap):
        if self.state.config.useDouble:
            return DoublePositionCalculation(pos, a, g, c, u, gap)
        return ApfloatPositionCalculation(
            pos, a, g, c, u, gap, self.state.config.precision
        )

    def getParser(self):
        if self.state.config.useSequences:
            return SequenceParser(self.state)
        return PositionFileParser(self.state)


class Interval:
    def __init__(self, start, end):
        self.start = start
        self.end = end

    def getEnd(self):
        return self.end

    def setEnd(self, end):
        self.end = end

    def getStart(self):
        return self.start

    def setStart(self, start):
        self.start = start


class Result:
    def __init__(self):
        self.time = None
        self.positions = None
        self.useSequence = None
        self.sequences = None
        self.sequencesTotal = None
        self.sequencesParcent = None
        self.intervals = None
        self.groups = None
        self.pnns = None

    def getTime(self):
        return self.time

    def getPossitions(self):
        return self.positions

    def isUseSequence(self):
        return self.useSequence

    def getSequences(self):
        return self.sequences

    def getSequencesTotal(self):
        return self.sequencesTotal

    def getSequencesParcent(self):
        return self.sequencesParcent

    def getIntervals(self):
        return self.intervals

    def getGroups(self):
        return self.groups


class Sequence:
    def __init__(self, number, meta, data):
        self.weight = 1
        self.number = number
        self.meta = meta
        self.data = data
        self.group = None

    def get_number(self):
        return self.number

    def set_number(self, number):
        self.number = number

    def get_weight(self):
        return self.weight

    def set_weight(self, weight):
        self.weight = weight

    def get_meta(self):
        return self.meta

    def set_meta(self, meta):
        self.meta = meta

    def get_data(self):
        return self.data

    def set_data(self, data):
        self.data = data

    def add_weight(self):
        self.weight += 1


class GroupedSequenceType(Enum):
    NORMAL = 0
    OUTGROUP1 = 1
    OUTGROUP2 = 2
    EXCLUDED = 3


class GroupedSequence:
    def __init__(self, type=GroupedSequenceType.NORMAL):
        self.sequences = []
        self.type = type
        self.number = 0
        self.name = None
        self.probability = 0.0
        self.fred = 0.0
        self.color = "red"

    def get_sequences(self):
        return self.sequences

    def get_type(self):
        return self.type

    def get_priority(self):
        return self.type.value

    def is_empty(self):
        return len(self.sequences) == 0

    def get_one_sample(self):
        if self.is_empty():
            return None
        return self.sequences[0].get_data()

    def compare_to(self, other):
        if self.get_priority() != other.get_priority():
            return self.get_priority() - other.get_priority()
        return len(other.get_sequences()) - len(self.get_sequences())


class DataParser(ABC):
    def __init__(self, state):
        self.state = state

    @abstractmethod
    def parse(self) -> list:
        pass


class SeqBigResultModel:
    def __init__(self):
        self.time = None
        self.seqLength = None
        self.seqTotal = None
        self.groups = None
        self.intervals = None

    def getTime(self):
        return self.time

    def getSeqLength(self):
        return self.seqLength

    def getSeqTotal(self):
        return self.seqTotal

    def getGroups(self):
        return self.groups

    def getIntervals(self):
        return self.intervals


MAX_SIZE = 100
CLEANUP_CHECK_PERIOD_MS = 20 * 1000
PERSISTENCE = datetime.timedelta(minutes=5)
PERSISTENCE_DELETED = datetime.timedelta(minutes=1)


class ResultStore(threading.Thread):
    class ValueItem:
        def __init__(self, uuid, value, persistence):
            self.uuid = uuid
            self.value = value
            self.delete_time = datetime.datetime.now() + persistence

        def touch(self):
            self.delete_time = datetime.datetime.now() + PERSISTENCE

        def mark_for_delete(self):
            self.delete_time = datetime.datetime.now() + PERSISTENCE_DELETED

    def __init__(self):
        super().__init__()
        self.results = {}
        self.start()

    def add_result(self, uuid, result):
        logger.debug(f"added: {uuid}")

        val = self.ValueItem(uuid, result, PERSISTENCE)

        if len(self.results) >= MAX_SIZE:
            oldest = min(self.results.values(), key=lambda v: v.delete_time)
            del self.results[oldest.uuid]
            logger.debug("too many results, removed: %s", oldest.uuid)

        self.results[uuid] = val

    def get_result(self, uuid):
        print(f"getting: {uuid}")

        v = self.results.get(uuid)

        if v is not None:
            v.touch()
            return v.value

        return None

    def mark_for_delete(self, uuid):
        logger.debug("marking for delete: %s", uuid)

        v = self.results.get(uuid)

        if v is not None:
            v.mark_for_delete()
            return True

        return False

    def run(self):
        # while True:
        # time.sleep(CLEANUP_CHECK_PERIOD_MS / 1000)

        # logger.trace("checking results, len: %s", len(self.results))

        to_delete = [
            v for v in self.results.values() if datetime.datetime.now() > v.delete_time
        ]

        for v in to_delete:
            logger.debug("deleting result: %s", v.uuid)
            del self.results[v.uuid]

        to_delete.clear()


# Calculations classes
class PositionCalculation(ABC):
    def __init__(self, pos, a, g, c, u, gap):
        self.ok = False
        self.position = pos
        self.a = a
        self.g = g
        self.c = c
        self.u = u
        self.gap = gap
        self.sum = 0
        self.total = 0
        self.coverage = 0.0

    def call(self):
        self.sum = self.a + self.g + self.c + self.u
        self.total = self.gap + self.sum
        self.coverage = (self.sum / self.total) * 100

        return True

    def is_ok(self):
        return self.ok

    def get_position(self):
        return self.position

    def get_a(self):
        return self.a

    def get_g(self):
        return self.g

    def get_c(self):
        return self.c

    def get_u(self):
        return self.u

    def get_gap(self):
        return self.gap

    def get_sum(self):
        return self.sum

    def get_total(self):
        return self.total

    def get_coverage(self):
        return self.coverage

    @abstractmethod
    def get_afreq(self):
        pass

    @abstractmethod
    def get_gfreq(self):
        pass

    @abstractmethod
    def get_cfreq(self):
        pass

    @abstractmethod
    def get_ufreq(self):
        pass

    @abstractmethod
    def get_lnafreq(self):
        pass

    @abstractmethod
    def get_lngfreq(self):
        pass

    @abstractmethod
    def get_lncfreq(self):
        pass

    @abstractmethod
    def get_lnufreq(self):
        pass

    @abstractmethod
    def get_aentropy(self):
        pass

    @abstractmethod
    def get_gentropy(self):
        pass

    @abstractmethod
    def get_centropy(self):
        pass

    @abstractmethod
    def get_uentropy(self):
        pass

    @abstractmethod
    def get_entropy(self):
        pass


class ApfloatPositionCalculation(threading.Thread):
    def __init__(self, pos, a, g, c, u, gap, precision):
        super().__init__()
        self.position = pos
        self.a = a
        self.g = g
        self.c = c
        self.u = u
        self.gap = gap
        self.precision = precision
        self.ok = False

        self.Afreq = None
        self.Gfreq = None
        self.Cfreq = None
        self.Ufreq = None
        self.lnAfreq = None
        self.lnGfreq = None
        self.lnCfreq = None
        self.lnUfreq = None
        self.AEntropy = None
        self.GEntropy = None
        self.CEntropy = None
        self.UEntropy = None
        self.entropy = None

    def run(self):
        try:
            sumApf = Apfloat(self.a + self.g + self.c + self.u, self.precision)

            self.Afreq = Apfloat(self.a, self.precision) / sumApf
            self.Gfreq = Apfloat(self.g, self.precision) / sumApf
            self.Cfreq = Apfloat(self.c, self.precision) / sumApf
            self.Ufreq = Apfloat(self.u, self.precision) / sumApf
            self.lnAfreq = self.log(self.Afreq)
            self.lnGfreq = self.log(self.Gfreq)
            self.lnCfreq = self.log(self.Cfreq)
            self.lnUfreq = self.log(self.Ufreq)
            self.AEntropy = self.Afreq * self.lnAfreq * -1
            self.GEntropy = self.Gfreq * self.lnGfreq * -1
            self.CEntropy = self.Cfreq * self.lnCfreq * -1
            self.UEntropy = self.Ufreq * self.lnUfreq * -1
            self.entropy = self.AEntropy + self.GEntropy + self.CEntropy + self.UEntropy
        except Exception as e:
            print(f"Position {self.position} calculation error, reason: {str(e)}")
            return False
        self.ok = True
        return True

    def log(self, x):
        if x == Apfloat(0, self.precision):
            return Apfloat(0, self.precision)

        return Apfloat.ln(x)

    def get_Afreq(self):
        return float(self.Afreq)

    def get_Gfreq(self):
        return float(self.Gfreq)

    def get_Cfreq(self):
        return float(self.Cfreq)

    def get_Ufreq(self):
        return float(self.Ufreq)

    def get_LnAfreq(self):
        return float(self.lnAfreq)

    def get_LnGfreq(self):
        return float(self.lnGfreq)

    def get_LnCfreq(self):
        return float(self.lnCfreq)

    def get_LnUfreq(self):
        return float(self.lnUfreq)

    def get_AEntropy(self):
        return float(self.AEntropy)

    def get_GEntropy(self):
        return float(self.GEntropy)

    def get_CEntropy(self):
        return float(self.CEntropy)

    def get_UEntropy(self):
        return float(self.UEntropy)

    def get_entropy(self):
        return float(self.entropy)


class DoublePositionCalculation(PositionCalculation):
    def __init__(self, pos, a, g, c, u, gap):
        super().__init__(pos, a, g, c, u, gap)
        self.aFreq = 0.0
        self.gFreq = 0.0
        self.cFreq = 0.0
        self.uFreq = 0.0
        self.lnAfreq = 0.0
        self.lnGfreq = 0.0
        self.lnCfreq = 0.0
        self.lnUfreq = 0.0
        self.aEntropy = 0.0
        self.gEntropy = 0.0
        self.cEntropy = 0.0
        self.uEntropy = 0.0
        self.entropy = 0.0
        self.ok = False

    def call(self):
        super().call()

        self.aFreq = self.a / self.sum
        self.gFreq = self.g / self.sum
        self.cFreq = self.c / self.sum
        self.uFreq = self.u / self.sum

        self.lnAfreq = self.log(self.aFreq)
        self.lnGfreq = self.log(self.gFreq)
        self.lnCfreq = self.log(self.cFreq)
        self.lnUfreq = self.log(self.uFreq)

        self.aEntropy = -self.aFreq * self.lnAfreq
        self.gEntropy = -self.gFreq * self.lnGfreq
        self.cEntropy = -self.cFreq * self.lnCfreq
        self.uEntropy = -self.uFreq * self.lnUfreq

        self.entropy = self.aEntropy + self.gEntropy + self.cEntropy + self.uEntropy

        self.ok = True
        return True

    def log(self, x):
        return 0 if x == 0 else math.log(x)

    def getAfreq(self):
        return self.aFreq

    def getGfreq(self):
        return self.gFreq

    def getCfreq(self):
        return self.cFreq

    def getUfreq(self):
        return self.uFreq

    def getLnAfreq(self):
        return self.lnAfreq

    def getLnGfreq(self):
        return self.lnGfreq

    def getLnCfreq(self):
        return self.lnCfreq

    def getLnUfreq(self):
        return self.lnUfreq

    def getAEntropy(self):
        return self.aEntropy

    def getGEntropy(self):
        return self.gEntropy

    def getCEntropy(self):
        return self.cEntropy

    def getUEntropy(self):
        return self.uEntropy

    def getEntropy(self):
        return self.entropy


class SequenceFileParser(DataParser):
    def __init__(self, state):
        super().__init__(state)
        self.sequences = None
        self.filtered_seqs = []
        self.grouped_sequences = []
        self.BASES1 = "ryswkmbdhvn"
        self.BASES2 = ".-_~?"

    def parse(self):
        self.sequences = self.parse_data()
        self.filter_data()
        return self.build_positions()

    def group(self):
        self.sequences = self.parse_data()
        self.filter_data()
        return self.grouped_sequences

    def pnns(self):
        seqs = self.parse_data()
        return self.build_pnns(seqs)

    def build_pnns(self, seqs):
        result = []

        for i in range(len(list(seqs.values())[0].get_data())):
            a = g = c = u = gap = misc = 0

            for entry in seqs.items():
                sq = entry[1]
                base = sq.get_data()[i]

                if base == "a":
                    a += 1
                elif base == "g":
                    g += 1
                elif base == "c":
                    c += 1
                elif base == "u" or base == "t":
                    u += 1
                elif self.BASES1.find(base) >= 0:
                    misc += 1
                elif self.BASES2.find(base) >= 0:
                    gap += 1

            result.append(PnnsCalculation(i + 1, a, g, c, u, gap, misc))

        return result

    def parse_data(self):
        self.state.start = time.time()
        seqs = {}

        with io.open(self.state.config.inFile, "r", encoding="utf-8") as reader:
            print("reading input file ...\n")

            i = 0
            while True:
                line = reader.readline()
                if not line:
                    break

                next_line = reader.readline()
                if not next_line:
                    break

                seqs[i] = Sequence(i + 1, line, next_line.lower())
                i += 1

            print(f"found {len(seqs)} sequences.\n")
            self.state.result.sequencesTotal = len(seqs)

            if seqs:
                self.state.result.positions = len(list(seqs.values())[0].get_data())

        return seqs

    def filter_data(self):
        if not self.state.config.wholeSequence and (
            not self.state.config.intervals or not self.state.config.intervals
        ):
            raise ValueError("Illegal argument")

        print(
            f"filtering sequences with {self.state.config.threshold}% threshold ...\n"
        )

        if self.state.config.wholeSequence:
            self.state.result.intervals = [Interval(1, self.state.result.positions)]
        else:
            self.make_substrings()

        self.filtered_seqs = []
        self.grouped_sequences = []

        outgroup1 = GroupedSequence(GroupedSequenceType.OUTGROUP1)
        outgroup2 = GroupedSequence(GroupedSequenceType.OUTGROUP2)
        excluded = GroupedSequence(GroupedSequenceType.EXCLUDED)

        to_delete = []

        for key, val in self.sequences.items():
            sq = val
            data = sq.get_data().lower()

            seq_type = self.tell_seq_type(data)

            if seq_type == GroupedSequenceType.EXCLUDED:
                excluded.sequences.append(sq)
                to_delete.append(key)
            elif seq_type == GroupedSequenceType.OUTGROUP2:
                outgroup2.sequences.append(sq)
                to_delete.append(key)
            elif seq_type == GroupedSequenceType.OUTGROUP1:
                outgroup1.sequences.append(sq)
                to_delete.append(key)

        for i in to_delete:
            self.sequences.pop(i)
        to_delete.clear()

        while self.sequences:
            key = next(iter(self.sequences))
            sq = self.sequences.pop(key)
            self.filtered_seqs.append(sq)

            gsq = GroupedSequence()
            gsq.sequences.append(sq)
            self.grouped_sequences.append(gsq)

            s1 = sq.get_data()

            for key, val in self.sequences.items():
                s2 = val.get_data()

                if s1 == s2:
                    sq.add_weight()
                    gsq.sequences.append(val)
                    to_delete.append(key)

            for i in to_delete:
                del self.sequences[i]
            to_delete.clear()

        if not outgroup1.is_empty():
            self.grouped_sequences.append(outgroup1)
        if not outgroup2.is_empty():
            self.grouped_sequences.append(outgroup2)
        if not excluded.is_empty():
            self.grouped_sequences.append(excluded)

        count = len(self.filtered_seqs)
        percent = 100 * count / self.state.result.sequencesTotal

        print(
            "done. "
            + str(count)
            + " "
            + "("
            + str(percent)
            + "%) unique sequences remaining."
        )
        self.state.result.sequences = count
        self.state.result.sequences_percent = percent
        self.state.result.time = int(round(time.time() * 1000)) - self.state.start

    def tell_seq_type(self, seq: str) -> GroupedSequenceType:
        inoutgr1 = False
        inoutgr2 = False
        inexcluded = True

        for i in range(len(seq)):
            base = seq[i]

            if not inoutgr2:
                if base in self.BASES2:
                    inoutgr2 = True

            if not inoutgr1 and not inoutgr2:
                if base in self.BASES1:
                    inoutgr1 = True

            if inexcluded:
                if base not in self.BASES2:
                    inexcluded = False

            if not inexcluded and inoutgr2:
                break

        if inexcluded:
            return GroupedSequenceType.EXCLUDED
        elif inoutgr2:
            return GroupedSequenceType.OUTGROUP2
        elif inoutgr1:
            return GroupedSequenceType.OUTGROUP1

        return GroupedSequenceType.NORMAL

    def make_substrings(self):
        len_ = self.state.result.positions
        intervals = [
            Interval(i.start, min(i.end, len_))
            for i in self.state.config.intervals
            if i.start < len_
        ]
        self.state.result.intervals = intervals

        for entry in self.sequences.items():
            seq_data = "".join([entry[1].data[i.start - 1 : i.end] for i in intervals])
            entry[1].data = seq_data

    def build_positions(self) -> list[PositionCalculation]:
        result = []

        for i in range(len(self.filtered_seqs[0].data)):
            a = g = c = u = 0

            for sq in self.filtered_seqs:
                if sq.data[i] == "a":
                    a += sq.weight
                elif sq.data[i] == "g":
                    g += sq.weight
                elif sq.data[i] == "c":
                    c += sq.weight
                elif sq.data[i] in ("u", "t"):
                    u += sq.weight

            result.append(self.state.factory.get_calculation(i + 1, a, g, c, u, 0))

        return result


class PnnsCalculation:
    def __init__(self, position, a, g, c, u, gap, misc):
        self.position = position
        self.a = a
        self.g = g
        self.c = c
        self.u = u
        self.gap = gap
        self.misc = misc

        self.sum = a + g + c + u
        self.total = gap + misc + self.sum

    def getPossition(self):
        return self.position

    def getA(self):
        return self.a

    def getG(self):
        return self.g

    def getC(self):
        return self.c

    def getU(self):
        return self.u

    def getGap(self):
        return self.gap

    def getSum(self):
        return self.sum

    def getTotal(self):
        return self.total
