"""Biological motif patterns and IUPAC code handling.

This module provides:
- Known biological motif patterns (TFBS, regulatory elements)
- IUPAC ambiguity code expansion
- Consensus pattern matching
"""

from typing import Dict, List, Set

# IUPAC nucleotide ambiguity codes
IUPAC_CODES: Dict[str, Set[str]] = {
    "A": {"A"},
    "C": {"C"},
    "G": {"G"},
    "T": {"T"},
    "R": {"A", "G"},      # puRine
    "Y": {"C", "T"},      # pYrimidine
    "S": {"G", "C"},      # Strong
    "W": {"A", "T"},      # Weak
    "K": {"G", "T"},      # Keto
    "M": {"A", "C"},      # aMino
    "B": {"C", "G", "T"}, # not A
    "D": {"A", "G", "T"}, # not C
    "H": {"A", "C", "T"}, # not G
    "V": {"A", "C", "G"}, # not T
    "N": {"A", "C", "G", "T"},  # aNy
}

# Known biological motifs (consensus sequences)
KNOWN_MOTIFS: Dict[str, Dict[str, str]] = {
    # Core promoter elements
    "TATA_BOX": {
        "consensus": "TATAAA",
        "description": "TATA box, core promoter element (~25bp upstream of TSS)",
    },
    "GC_BOX": {
        "consensus": "GGGCGG",
        "description": "GC box, Sp1 transcription factor binding site",
    },
    "CAAT_BOX": {
        "consensus": "CCAAT",
        "description": "CAAT box, enhancer element (~80bp upstream of TSS)",
    },
    # E-box family
    "E_BOX": {
        "consensus": "CANNTG",
        "description": "E-box, helix-loop-helix TF binding site",
    },
    "E_BOX_CANONICAL": {
        "consensus": "CACGTG",
        "description": "Canonical E-box sequence",
    },
    # Restriction sites (common)
    "ECORI": {
        "consensus": "GAATTC",
        "description": "EcoRI restriction site",
    },
    "BAMHI": {
        "consensus": "GGATCC",
        "description": "BamHI restriction site",
    },
    "HINDIII": {
        "consensus": "AAGCTT",
        "description": "HindIII restriction site",
    },
    # Kozak sequence
    "KOZAK": {
        "consensus": "GCCGCCATGG",
        "description": "Kozak consensus sequence for translation initiation",
    },
    # Poly-A signal
    "POLYA_SIGNAL": {
        "consensus": "AATAAA",
        "description": "Polyadenylation signal",
    },
    # Splice sites
    "SPLICE_DONOR": {
        "consensus": "GTAAGT",
        "description": "Splice donor consensus (5' splice site)",
    },
    "SPLICE_ACCEPTOR": {
        "consensus": "NYAG",
        "description": "Splice acceptor consensus (3' splice site)",
    },
}


def expand_iupac(pattern: str) -> List[str]:
    """Expand an IUPAC pattern to all matching DNA sequences.

    Args:
        pattern: DNA pattern with IUPAC ambiguity codes (e.g., "CANNTG")

    Returns:
        List of all possible DNA sequences matching the pattern

    Example:
        >>> expand_iupac("CAT")
        ['CAT']
        >>> sorted(expand_iupac("CAN"))
        ['CAA', 'CAC', 'CAG', 'CAT']
    """
    pattern = pattern.upper()
    if not pattern:
        return [""]

    first = pattern[0]
    rest = pattern[1:]

    if first not in IUPAC_CODES:
        raise ValueError(f"Unknown IUPAC code: {first}")

    rest_expanded = expand_iupac(rest) if rest else [""]
    result = []
    for base in sorted(IUPAC_CODES[first]):
        for suffix in rest_expanded:
            result.append(base + suffix)

    return result


def hamming_distance(s1: str, s2: str) -> int:
    """Calculate Hamming distance between two equal-length strings."""
    if len(s1) != len(s2):
        raise ValueError("Strings must have equal length")
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def matches_consensus(kmer: str, consensus: str, max_mismatch: int = 0) -> bool:
    """Check if a k-mer matches a consensus pattern.

    The consensus can contain IUPAC ambiguity codes.

    Args:
        kmer: DNA sequence to check
        consensus: Consensus pattern (may contain IUPAC codes)
        max_mismatch: Maximum allowed mismatches

    Returns:
        True if kmer matches consensus within allowed mismatches
    """
    kmer = kmer.upper()
    consensus = consensus.upper()

    if len(kmer) != len(consensus):
        return False

    mismatches = 0
    for k_base, c_code in zip(kmer, consensus):
        if c_code not in IUPAC_CODES:
            raise ValueError(f"Unknown IUPAC code: {c_code}")
        if k_base not in IUPAC_CODES[c_code]:
            mismatches += 1
            if mismatches > max_mismatch:
                return False

    return True


def filter_by_consensus(
    kmers: List[str], consensus: str, max_mismatch: int = 1
) -> List[str]:
    """Filter k-mers that match a consensus pattern.

    Args:
        kmers: List of k-mers to filter
        consensus: IUPAC consensus pattern
        max_mismatch: Maximum allowed mismatches

    Returns:
        List of k-mers matching the consensus
    """
    return [k for k in kmers if matches_consensus(k, consensus, max_mismatch)]


def get_reverse_complement(seq: str) -> str:
    """Get the reverse complement of a DNA sequence.

    Args:
        seq: DNA sequence (ACGT)

    Returns:
        Reverse complement sequence
    """
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(complement.get(base, "N") for base in reversed(seq.upper()))


def list_known_motifs() -> List[str]:
    """List all known motif names."""
    return sorted(KNOWN_MOTIFS.keys())


def get_motif_info(name: str) -> Dict[str, str]:
    """Get information about a known motif.

    Args:
        name: Motif name (case-insensitive)

    Returns:
        Dict with 'consensus' and 'description'

    Raises:
        KeyError if motif not found
    """
    return KNOWN_MOTIFS[name.upper()]
