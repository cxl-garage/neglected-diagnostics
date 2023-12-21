import os
import subprocess

import pandas as pd

outdir = os.getcwd()


# wrap consensus sequence to a class, with methods to identify consensus sequence, calculate consensus score and to get ideal consensus regions
class ConsensusSeq:
    def __init__(self, multifasta):
        self.infile = multifasta
        self.consensus_dict = {}
        self.consensus_seq = ""
        self.consensus_score = []
        self.consensus_seq_gapless = ""
        self.consensus_regions = []
        self.consensus_flag = []
        self.excluded_regions = []
        self.consensus_flag_inverted = []
        self.fasta = {}
        self.fastaheader = {}
        self.fasta, self.fastaheader = self.read_fasta(multifasta)
        self.number_of_sequences = len(self.fasta)
        self.fasta_pd = None
        self.fasta_to_panda()
        # self.compute_consensus_sequence()
        # self.compute_consensus_score()
        # self.compute_consensus_seq_gapless()
        # self.compute_consensus_regions()

    def read_fasta(self, multifasta, delim="\t", idpos=0):
        """reading input fasta file
        returns fasta dictionary with key=accessionNumber and value=Sequence
        returns fastaheader dictionary with key=accesstionNumber and value=originalFastaHeader
        """
        fasta = {}
        fastaheader = {}
        with open(multifasta, "r") as infile:
            acNumber = ""
            for line in infile:
                if line.startswith(">"):
                    if delim:
                        acNumber = line.split(delim)[idpos].strip().strip(">")
                        fastaheader[acNumber] = line.strip()
                    else:
                        acNumber = line.split()[idpos].strip().strip(">")
                        fastaheader[acNumber] = line.strip()
                else:
                    if acNumber in fasta:
                        fasta[acNumber] += line.strip().upper()
                    else:
                        fasta[acNumber] = line.strip().upper()
        return fasta, fastaheader

    def fasta_to_panda(self):
        """convert fasta to panda dataframe to calculate consensus scores"""
        fasta_for_panda = {}
        for gene in self.fasta:
            fasta_for_panda[gene] = list(self.fasta[gene])
        self.fasta_pd = pd.DataFrame.from_dict(fasta_for_panda, orient="index")

    def align(self, outfile):
        with open(outfile, "w") as outfile:
            subprocess.call(
                ["mafft", "--auto", "--adjustdirection", "--thread -1", self.infile],
                stdout=outfile,
                stderr=subprocess.DEVNULL,
            )

    def compute_consensus_sequence(self):
        self.compute_consensus_score()
        """get consensus sequence
        The score is the proportion of the most common amino acid in a column of the alignment."""
        for col in self.fasta_pd:
            consensus_nt = ""
            max_consensus_score = 0
            for nt in self.consensus_dict[col]:
                score = self.consensus_dict[col][nt]
                if score > max_consensus_score:
                    consensus_nt = nt
                    max_consensus_score = score
            self.consensus_seq += consensus_nt
            self.consensus_score.append(max_consensus_score)
            max_consensus_score = 0

    def compute_consensus_score(self):
        """calculate nucleotide ratio for every position in the alignment
        NOTE: non-standard nucleotides are counted as '-'"""
        for col in self.fasta_pd:
            nucl_score = {"A": 0, "T": 0, "C": 0, "G": 0, "-": 0}
            for nt in self.fasta_pd[col]:
                if nt in nucl_score:
                    nucl_score[nt] += 1
                else:
                    nucl_score["-"] += 1
            for key in nucl_score:
                nucl_score[key] = nucl_score[key] / self.number_of_sequences
            self.consensus_dict[col] = nucl_score

    def compute_consensus_seq_gapless(self):
        """remove gaps from the consensus sequence as primers can not be designed with gaps"""
        gap_positions = []
        for i in range(len(self.consensus_seq) - 1, -1, -1):
            if "-" == self.consensus_seq[i]:
                gap_positions.append(i)
        self.consensus_seq_gapless = list(self.consensus_seq)
        for pos in gap_positions:
            del self.consensus_seq_gapless[pos]
            del self.consensus_score[pos]
        self.consensus_seq_gapless = "".join(self.consensus_seq_gapless)

    def compute_consensus_regions(self, threshold=0.95):
        """get regions above the consensus threshold (ratio of most common amino acid per posisiton)
        a consensus region must be at least of the length 20 to be considered for primer design
        """

        print(f"Consensus threshold set to: {threshold}\n")
        minlength = 20
        start = 0
        stop = 0
        new = True
        old_start = -1
        for pos in range(len(self.consensus_score)):
            if self.consensus_score[pos] >= threshold:
                if new:
                    start = pos
                    stop = pos
                    new = False
                else:
                    stop = pos
            elif stop - start >= minlength and start != old_start:
                self.consensus_regions.append([start, stop])
                old_start = start
                new = True
            else:
                new = True
        if stop - start >= minlength and start != old_start:
            self.consensus_regions.append([start, stop])
        self.compute_consensus_flag()
        self.compute_excluded_regions()

    def compute_consensus_flag(self):
        """get regions that have to be excluded from the primer design
        Necessary because primer3 can only automatically search for primers over multiple excluded areas in multiple regions.
        """
        self.consensus_flag = [0] * len(self.consensus_score)
        for region in self.consensus_regions:
            for i in range(region[0], region[1] + 1):
                self.consensus_flag[i] = 1
        self.consensus_flag_inverted = [1 if x == 0 else 0 for x in self.consensus_flag]

    def compute_excluded_regions(self):
        """get regions that have to be excluded from the primer design
        Necessary because primer3 can only automatically search for primers over multiple excluded areas in multiple regions.
        """
        start = 0
        length = 0
        new = True
        for pos in range(len(self.consensus_flag)):
            if self.consensus_flag[pos] == 0:
                if new:
                    start = pos
                length += 1
                new = False
            elif not new:
                self.excluded_regions.append([start, length])
                new = True
                length = 0
        if self.excluded_regions:
            if not self.excluded_regions[-1] == [start, length] and length >= 1:
                self.excluded_regions.append([start, length])

    def export_consensus_regions(self, outfile):
        """write identified consensus regions to a csv table"""
        with open(outfile, "w") as outfile:
            outfile.write("Start\tStop\tLength\tConsensus_score\tSequence\n")
            outfile.write(
                f'1\t{len(self.consensus_seq_gapless)}\t{len(self.consensus_seq_gapless)}\t{"{:.3f}".format(avg(self.consensus_score))}\t{self.consensus_seq_gapless}\n'
            )
            if not self.consensus_regions:
                print(
                    "#####################################################################\n### No Consensus regions detected for current parameter settings. ###\n#####################################################################"
                )
                print("\nDone.")
                exit()
            for region in self.consensus_regions:
                scores = self.consensus_score[region[0] : region[1] + 1]
                score_avg = avg(scores)
                scores = "\t".join([str("{:.2f}".format(x)) for x in scores])
                print(
                    f'Start: {region[0]+1} Stop: {region[1]+1} Length: {region[1]+1 - region[0]} Consensus_score: {"{:.3f}".format(score_avg)} Sequence: {self.consensus_seq_gapless[region[0]:region[1]+1]}'
                )
                outfile.write(
                    f'{region[0]+1}\t{region[1]+1}\t{region[1]+1 - region[0]}\t{"{:.3f}".format(score_avg)}\t{self.consensus_seq_gapless[region[0]:region[1]+1]}\n'
                )
        #    print(scores)
