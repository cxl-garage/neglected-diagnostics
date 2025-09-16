"""
Sequence Filtering and Quality Control

This module provides functions for filtering and validating sequences based on
various quality criteria.
"""

import re
from typing import List, Optional, Set, Tuple

try:
    from Bio.SeqUtils import GC
except ImportError:
    # Fallback for older BioPython versions or if GC is not available
    def GC(sequence):
        """Calculate GC content as percentage."""
        if not sequence:
            return 0.0
        sequence = sequence.upper()
        gc_count = sequence.count("G") + sequence.count("C")
        return (gc_count / len(sequence)) * 100.0


from .database_models import MoleculeType, SequenceRecord


class SequenceFilter:
    """Filters sequences based on quality and content criteria."""

    def __init__(
        self,
        min_length: Optional[int] = None,
        max_length: Optional[int] = None,
        max_n_content: float = 0.05,  # Max 5% N's
        min_gc_content: Optional[float] = None,
        max_gc_content: Optional[float] = None,
        exclude_partial: bool = True,
        exclude_predicted: bool = True,
        max_homopolymer_length: int = 10,
        require_start_codon: bool = False,
        require_stop_codon: bool = False,
    ):
        """
        Initialize sequence filter with quality criteria.

        Parameters
        ----------
        min_length : int, optional
            Minimum sequence length
        max_length : int, optional
            Maximum sequence length
        max_n_content : float
            Maximum proportion of N characters allowed (default 5%)
        min_gc_content : float, optional
            Minimum GC content proportion
        max_gc_content : float, optional
            Maximum GC content proportion
        exclude_partial : bool
            Exclude sequences marked as partial
        exclude_predicted : bool
            Exclude predicted/hypothetical sequences
        max_homopolymer_length : int
            Maximum allowed homopolymer run length
        require_start_codon : bool
            Require valid start codon (for coding sequences)
        require_stop_codon : bool
            Require valid stop codon (for coding sequences)
        """
        self.min_length = min_length
        self.max_length = max_length
        self.max_n_content = max_n_content
        self.min_gc_content = min_gc_content
        self.max_gc_content = max_gc_content
        self.exclude_partial = exclude_partial
        self.exclude_predicted = exclude_predicted
        self.max_homopolymer_length = max_homopolymer_length
        self.require_start_codon = require_start_codon
        self.require_stop_codon = require_stop_codon

        # Common partial/predicted sequence indicators
        self.partial_indicators = {
            "partial",
            "fragment",
            "incomplete",
            "truncated",
            "partial cds",
            "partial sequence",
        }
        self.predicted_indicators = {
            "predicted",
            "hypothetical",
            "putative",
            "probable",
            "uncharacterized",
            "similar to",
        }

    def filter_sequences(self, sequences: List[SequenceRecord]) -> List[SequenceRecord]:
        """
        Filter a list of sequences based on quality criteria.

        Parameters
        ----------
        sequences : List[SequenceRecord]
            Input sequences to filter

        Returns
        -------
        List[SequenceRecord]
            Filtered sequences that pass all criteria
        """
        filtered = []

        for seq in sequences:
            if self.passes_filters(seq):
                # Calculate quality score
                seq.quality_score = self.calculate_quality_score(seq)
                filtered.append(seq)

        return filtered

    def passes_filters(self, sequence: SequenceRecord) -> bool:
        """
        Check if a sequence passes all filter criteria.

        Parameters
        ----------
        sequence : SequenceRecord
            Sequence to evaluate

        Returns
        -------
        bool
            True if sequence passes all filters
        """
        # Length filters
        if self.min_length and sequence.length < self.min_length:
            return False
        if self.max_length and sequence.length > self.max_length:
            return False

        # Skip sequences without actual sequence data for content-based filters
        if not sequence.sequence:
            return True

        # N content filter
        if not self._check_n_content(sequence.sequence):
            return False

        # GC content filter
        if not self._check_gc_content(sequence.sequence):
            return False

        # Homopolymer filter
        if not self._check_homopolymers(sequence.sequence):
            return False

        # Partial sequence filter
        if self.exclude_partial and self._is_partial(sequence):
            return False

        # Predicted sequence filter
        if self.exclude_predicted and self._is_predicted(sequence):
            return False

        # Codon filters (for coding sequences)
        if sequence.molecule_type == MoleculeType.DNA:
            if self.require_start_codon and not self._has_start_codon(
                sequence.sequence
            ):
                return False
            if self.require_stop_codon and not self._has_stop_codon(sequence.sequence):
                return False

        return True

    def _check_n_content(self, sequence: str) -> bool:
        """Check if N content is within acceptable limits."""
        if not sequence:
            return True
        n_count = sequence.upper().count("N")
        n_proportion = n_count / len(sequence)
        return n_proportion <= self.max_n_content

    def _check_gc_content(self, sequence: str) -> bool:
        """Check if GC content is within acceptable limits."""
        if not sequence or len(sequence) < 10:  # Skip very short sequences
            return True

        try:
            gc_content = GC(sequence) / 100.0  # BioPython returns percentage

            if self.min_gc_content and gc_content < self.min_gc_content:
                return False
            if self.max_gc_content and gc_content > self.max_gc_content:
                return False

            return True
        except Exception:
            # If GC calculation fails, don't exclude the sequence
            return True

    def _check_homopolymers(self, sequence: str) -> bool:
        """Check for excessive homopolymer runs."""
        if not sequence:
            return True

        # Find longest homopolymer run
        max_run = 0
        current_run = 1
        prev_char = sequence[0] if sequence else ""

        for char in sequence[1:]:
            if char == prev_char:
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 1
            prev_char = char

        return max_run <= self.max_homopolymer_length

    def _is_partial(self, sequence: SequenceRecord) -> bool:
        """Check if sequence is marked as partial."""
        title_lower = sequence.title.lower()
        return any(indicator in title_lower for indicator in self.partial_indicators)

    def _is_predicted(self, sequence: SequenceRecord) -> bool:
        """Check if sequence is predicted/hypothetical."""
        title_lower = sequence.title.lower()
        return any(indicator in title_lower for indicator in self.predicted_indicators)

    def _has_start_codon(self, sequence: str) -> bool:
        """Check if sequence starts with a valid start codon."""
        if len(sequence) < 3:
            return False
        start_codon = sequence[:3].upper()
        valid_start_codons = {"ATG", "GTG", "TTG", "CTG"}
        return start_codon in valid_start_codons

    def _has_stop_codon(self, sequence: str) -> bool:
        """Check if sequence ends with a valid stop codon."""
        if len(sequence) < 3:
            return False
        stop_codon = sequence[-3:].upper()
        valid_stop_codons = {"TAA", "TAG", "TGA"}
        return stop_codon in valid_stop_codons

    def calculate_quality_score(self, sequence: SequenceRecord) -> float:
        """
        Calculate a quality score for the sequence (0-1 scale).

        Parameters
        ----------
        sequence : SequenceRecord
            Sequence to score

        Returns
        -------
        float
            Quality score between 0 and 1
        """
        if not sequence.sequence:
            return 0.5  # Neutral score for sequences without sequence data

        score = 1.0
        seq = sequence.sequence.upper()

        # Length score (optimal range varies by marker)
        length_score = self._calculate_length_score(len(seq))
        score *= length_score

        # N content penalty
        n_content = seq.count("N") / len(seq) if seq else 0
        n_score = max(0, 1 - (n_content / self.max_n_content))
        score *= n_score

        # GC content score (penalize extreme values)
        try:
            gc_content = GC(seq) / 100.0
            gc_score = self._calculate_gc_score(gc_content)
            score *= gc_score
        except Exception:
            pass  # Don't penalize if GC calculation fails

        # Homopolymer penalty
        max_homopolymer = self._get_max_homopolymer_length(seq)
        homopolymer_score = max(
            0, 1 - (max_homopolymer / (2 * self.max_homopolymer_length))
        )
        score *= homopolymer_score

        # Title quality indicators
        title_score = self._calculate_title_score(sequence.title)
        score *= title_score

        return max(0.0, min(1.0, score))

    def _calculate_length_score(self, length: int) -> float:
        """Calculate score based on sequence length."""
        if self.min_length and length < self.min_length:
            return 0.0
        if self.max_length and length > self.max_length:
            return 0.0

        # Optimal length ranges for common markers
        optimal_ranges = {
            "COI": (600, 700),
            "16S": (1400, 1600),
            "18S": (1700, 1900),
            "ITS": (400, 800),
            "rbcL": (1300, 1500),
        }

        # Use a general optimal range if no specific marker
        optimal_min, optimal_max = (500, 1500)

        if optimal_min <= length <= optimal_max:
            return 1.0
        elif length < optimal_min:
            return length / optimal_min
        else:
            return max(0.5, optimal_max / length)

    def _calculate_gc_score(self, gc_content: float) -> float:
        """Calculate score based on GC content."""
        # Optimal GC content is typically 40-60%
        if 0.4 <= gc_content <= 0.6:
            return 1.0
        elif gc_content < 0.4:
            return max(0.5, gc_content / 0.4)
        else:
            return max(0.5, (1.0 - gc_content) / 0.4)

    def _get_max_homopolymer_length(self, sequence: str) -> int:
        """Get the length of the longest homopolymer run."""
        if not sequence:
            return 0

        max_run = 1
        current_run = 1
        prev_char = sequence[0]

        for char in sequence[1:]:
            if char == prev_char:
                current_run += 1
                max_run = max(max_run, current_run)
            else:
                current_run = 1
            prev_char = char

        return max_run

    def _calculate_title_score(self, title: str) -> float:
        """Calculate score based on title quality indicators."""
        title_lower = title.lower()

        # Penalize partial/predicted sequences
        if any(indicator in title_lower for indicator in self.partial_indicators):
            return 0.7
        if any(indicator in title_lower for indicator in self.predicted_indicators):
            return 0.8

        # Reward sequences with species names
        if re.search(r"\b[A-Z][a-z]+ [a-z]+\b", title):  # Binomial nomenclature
            return 1.0

        return 0.9


def deduplicate_sequences(
    sequences: List[SequenceRecord],
    similarity_threshold: float = 0.95,
    keep_longest: bool = True,
) -> List[SequenceRecord]:
    """
    Remove highly similar sequences to reduce redundancy.

    Parameters
    ----------
    sequences : List[SequenceRecord]
        Input sequences
    similarity_threshold : float
        Similarity threshold for considering sequences duplicates
    keep_longest : bool
        If True, keep the longest sequence among duplicates

    Returns
    -------
    List[SequenceRecord]
        Deduplicated sequence list
    """
    if not sequences:
        return []

    # Simple deduplication based on exact matches for now
    # Could be enhanced with sequence alignment for similarity detection
    seen_sequences = set()
    deduplicated = []

    # Sort by length if keeping longest
    if keep_longest:
        sequences = sorted(sequences, key=lambda x: x.length, reverse=True)

    for seq in sequences:
        if seq.sequence and seq.sequence not in seen_sequences:
            seen_sequences.add(seq.sequence)
            deduplicated.append(seq)
        elif not seq.sequence:
            # Keep sequences without sequence data
            deduplicated.append(seq)

    return deduplicated


def mask_low_complexity_regions(sequence: str, window_size: int = 20) -> str:
    """
    Mask low-complexity regions in a sequence.

    Parameters
    ----------
    sequence : str
        Input DNA sequence
    window_size : int
        Window size for complexity analysis

    Returns
    -------
    str
        Sequence with low-complexity regions masked with 'N'
    """
    if not sequence or len(sequence) < window_size:
        return sequence

    masked_seq = list(sequence.upper())

    for i in range(len(sequence) - window_size + 1):
        window = sequence[i : i + window_size]
        complexity = calculate_sequence_complexity(window)

        # Mask if complexity is too low
        if complexity < 0.3:  # Threshold for low complexity
            for j in range(i, i + window_size):
                masked_seq[j] = "N"

    return "".join(masked_seq)


def calculate_sequence_complexity(sequence: str) -> float:
    """
    Calculate sequence complexity using Shannon entropy.

    Parameters
    ----------
    sequence : str
        Input sequence

    Returns
    -------
    float
        Complexity score (0-1, higher is more complex)
    """
    if not sequence:
        return 0.0

    # Count character frequencies
    char_counts = {}
    for char in sequence.upper():
        char_counts[char] = char_counts.get(char, 0) + 1

    # Calculate Shannon entropy
    length = len(sequence)
    entropy = 0.0

    for count in char_counts.values():
        if count > 0:
            probability = count / length
            entropy -= probability * (probability.bit_length() - 1)

    # Normalize by maximum possible entropy (log2 of alphabet size)
    max_entropy = (len(char_counts)).bit_length() - 1 if len(char_counts) > 1 else 1

    return entropy / max_entropy if max_entropy > 0 else 0.0
