"""Tests for the discovery module."""

import pytest
from motifqu.discovery import (
    kmer_to_index,
    index_to_kmer,
    enumerate_kmers,
    count_kmer_occurrences,
    build_significance_oracle,
)
from motifqu.bio_patterns import (
    expand_iupac,
    matches_consensus,
    get_reverse_complement,
)


class TestKmerEncoding:
    """Test k-mer encoding and decoding."""

    def test_kmer_to_index_aa(self):
        assert kmer_to_index("AA") == 0

    def test_kmer_to_index_ac(self):
        assert kmer_to_index("AC") == 1

    def test_kmer_to_index_tt(self):
        assert kmer_to_index("TT") == 15

    def test_index_to_kmer_roundtrip(self):
        for k in [3, 4, 5]:
            for idx in range(4**k):
                kmer = index_to_kmer(idx, k)
                assert kmer_to_index(kmer) == idx

    def test_enumerate_kmers_4(self):
        kmers = enumerate_kmers(2)
        assert len(kmers) == 16
        assert "AA" in kmers
        assert "TT" in kmers
        assert "GC" in kmers


class TestKmerCounting:
    """Test k-mer occurrence counting."""

    def test_count_simple(self):
        genome = "AAAAACGT"
        counts = count_kmer_occurrences(genome, 2, include_revcomp=False)
        assert "AA" in counts
        assert len(counts["AA"]) == 4  # positions 0,1,2,3

    def test_count_with_revcomp(self):
        genome = "ACGTACGT"
        counts = count_kmer_occurrences(genome, 4, include_revcomp=True)
        # ACGT is its own reverse complement
        assert "ACGT" in counts

    def test_build_significance_oracle(self):
        genome = "ATATAT"  # AT appears 3 times, TA appears 2 times
        marked, significant = build_significance_oracle(
            genome, k=2, min_count=2, include_revcomp=False
        )
        assert len(marked) > 0
        assert "AT" in significant or "TA" in significant


class TestBioPatterns:
    """Test biological pattern utilities."""

    def test_expand_iupac_simple(self):
        assert expand_iupac("A") == ["A"]
        assert expand_iupac("N") == ["A", "C", "G", "T"]

    def test_expand_iupac_pattern(self):
        result = expand_iupac("AN")
        assert len(result) == 4
        assert "AA" in result
        assert "AT" in result

    def test_matches_consensus_exact(self):
        assert matches_consensus("TATAAA", "TATAAA", 0) is True
        assert matches_consensus("TATAAG", "TATAAA", 0) is False

    def test_matches_consensus_with_iupac(self):
        # N matches any base
        assert matches_consensus("CACATG", "CANNTG", 0) is True
        assert matches_consensus("CAGGTG", "CANNTG", 0) is True

    def test_matches_consensus_with_mismatch(self):
        assert matches_consensus("TATAAG", "TATAAA", 1) is True
        assert matches_consensus("TATACG", "TATAAA", 1) is False

    def test_reverse_complement(self):
        assert get_reverse_complement("ACGT") == "ACGT"
        assert get_reverse_complement("AAAA") == "TTTT"
        assert get_reverse_complement("GCTA") == "TAGC"
