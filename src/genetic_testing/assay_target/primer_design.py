""" 
This script is taken from the Thylacine_Design repository (https://github.com/kmceres/Thylacine_Design/tree/main). It is used to to detect target sequences among whole genome sequences, or gene sequences and off-target genomes/genes.

Link to this script: https://github.com/kmceres/Thylacine_Design/blob/main/general_primer_design/scripts/Primer_opt.py
"""
import re
from copy import deepcopy
from io import StringIO

import numpy as np
import pandas as pd
from Bio import SeqIO
from fuzzysearch import find_near_matches

from genetic_testing.assay_target.datatypes import AssayTargetColumns


class primer_list:
    def __init__(self, primer_seq, mean_score, bpwise_score):
        """

        :rtype: object
        """
        self.primer_seq = primer_seq
        self.mean_score = mean_score
        self.bpwise_score = bpwise_score


class exact_match_primer:
    def __init__(self, primer_seq, window_l, n_match):
        self.primer_seq = primer_seq
        self.window_l = window_l
        self.n_match = n_match


def import_files(target_files, off_target_files, reference_sequence):
    target_files_sorted = sorted(target_files, key=lambda x: x.name)
    off_target_files_sorted = sorted(off_target_files, key=lambda x: x.name)
    # targ_fasta_paths = glob.glob(os.path.join(target_path, "*.fasta"))
    # off_fasta_paths = glob.glob(os.path.join(off_path, "*.fasta"))

    # create lists of sequences (target and off target)
    # and labels (0 = target, !0 = off target)
    aln = []
    labs = []
    index = -1  # default value for the index of the reference sequence
    for fasta_file in target_files_sorted:
        # To convert to a string based IO:
        stringio = StringIO(fasta_file.getvalue().decode("utf-8"))
        for i, seq_record in enumerate(SeqIO.parse(stringio, "fasta")):
            if seq_record.id == reference_sequence:
                index = i

            labs.append(0)
            aln.append(str(seq_record.seq))

    for fasta_file in off_target_files_sorted:
        # To convert to a string based IO:
        stringio = StringIO(fasta_file.getvalue().decode("utf-8"))
        for i, seq_record in enumerate(SeqIO.parse(stringio, "fasta")):
            labs.append(1)
            aln.append(str(seq_record.seq))

    target = [i for i, val in enumerate(labs) if val == 0]  # target seq has label 0
    off = [i for i, val in enumerate(labs) if val != 0]
    return (aln, labs, target, off, index)


# select reference genome
# remove reference from target list
def select_ref(aln, target, index=0):
    ref_ind = target[index]
    ref = aln[ref_ind]
    aln_l = deepcopy(target)
    # Remove reference from target list by index
    del aln_l[ref_ind]

    return (aln_l, ref)


def iter_window(ref, window, slide_ind, aln, aln_l):
    j = 0
    seq_l = []
    while j < len(ref) - window:
        primer_seq = ref[j : j + window]
        new_primer = detect_exact_seqs(primer_seq, aln_l, aln)
        seq_l.append(new_primer)
        j = j + slide_ind
    return seq_l


def detect_fuzzy_matches(seq_l, aln_l, off, aln, maxDif_t, maxDif_ot):
    # for low matches and no matches, get hamming dist
    primers_l = []
    matches = []  # number of exact matches
    target_snps = []  # avg hamming distance primer-target
    target_snps_bpwise = []  # bp wise error percentage
    off_snps = []  # avg hamming distance primer-off
    if len(seq_l) > 0:
        for i in range(len(seq_l)):
            primer_seq = seq_l[i].primer_seq
            new_target = calc_nuc_dif(primer_seq, aln_l, aln, maxDif_t)
            new_off_target = calc_nuc_dif(primer_seq, off, aln, maxDif_ot)
            primers_l.append(primer_seq)
            matches.append(seq_l[i].n_match / len(aln_l) * 100)
            target_snps.append(new_target.mean_score)
            target_snps_bpwise.append(new_target.bpwise_score)
            off_snps.append(new_off_target.mean_score)
    return (primers_l, target_snps, target_snps_bpwise, matches, off_snps)


def write_out(primers_l, target_snps, target_snps_bpwise, matches, off_snps):
    cols = AssayTargetColumns()
    df = pd.DataFrame(
        {
            cols.assay_design_area: primers_l,
            cols.perc_tgt_match: matches,
            cols.ratio_tgt_mismatch: target_snps,
            cols.bpwise_error_percentage_tgt_mismatch: target_snps_bpwise,
            cols.num_off_tgt_mismatch: off_snps,
        }
    )
    df.sort_values(
        by=[
            "% target match",
            "Average of target mismatches",
            "Bpwise error percentage of target mismatches",
            "Minimum # off-target mismatches",
        ],
        ascending=[False, True, True, False],
        inplace=True,
    )

    return df
    # df.to_csv(filename, index=False)


# finds exact matches of window sized sequence fragment created from reference genome
# compared first to all targets. If an exact match is found in > threshold fraction of
# targets, and the fragment is not found in any of the off target sequences,
# if not, None is returned


def detect_exact_seqs(primer_seq, aln_l, aln):
    window_l = []
    for i in range(len(aln_l)):
        seq2n = aln_l[i]
        s2 = aln[seq2n]
        if re.search(primer_seq, s2):
            window_l.append(seq2n)
    le = len(window_l)
    new_match = exact_match_primer(primer_seq, window_l, le)
    return new_match


# calculates pairwise nuc differences (score) for a given primer
# and all other target or off target sequences
def calc_nuc_dif(primer_seq, aln_l, aln, max_dif):
    mismatch_l = []
    for i in range(len(aln_l)):
        seq2n = aln_l[i]
        s2 = aln[seq2n]
        matches = find_near_matches(primer_seq, s2, max_l_dist=max_dif)
        dist = []
        if matches != []:
            for match in matches:
                dist.append(match.dist)
            mismatch_l.append(min(dist))
        else:
            mismatch_l.append(max_dif)
    mean_score = np.mean(mismatch_l)
    mean_seq_len = sum(map(len, aln)) / len(aln)
    bpwise_score = mean_score / mean_seq_len * 100
    primer_scores = primer_list(primer_seq, mean_score, bpwise_score)
    return primer_scores


def find_target_area(
    target_files,
    off_target_files,
    reference_sequence,
    window=300,
    slide=20,
    maxDif_t=5,
    maxDif_ot=15,
):
    aln, labs, target, off, index = import_files(
        target_files, off_target_files, reference_sequence
    )

    if index == -1:
        raise ValueError("Reference sequence not found in the target files")

    aln_l, ref = select_ref(aln, target, index=index)
    seq_l = iter_window(ref, window=window, slide_ind=slide, aln=aln, aln_l=aln_l)
    (
        primers_l,
        target_snps,
        target_snps_bpwise,
        matches,
        off_snps,
    ) = detect_fuzzy_matches(
        seq_l, aln_l, off, aln, maxDif_t=maxDif_t, maxDif_ot=maxDif_ot
    )
    df = write_out(primers_l, target_snps, target_snps_bpwise, matches, off_snps)

    return df
