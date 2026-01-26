"""Quantum motif discovery using Grover's algorithm.

This module provides functions to discover significant motifs (k-mers)
in a genome sequence using quantum amplitude amplification.
"""

import math
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np
from qiskit import QuantumCircuit, transpile

from motifqu.bio_patterns import get_reverse_complement
from motifqu.grover import apply_diffuser, apply_mark_indices_phase_oracle
from motifqu.util import log

try:
    from qiskit_aer import AerSimulator
except ImportError as e:
    raise ImportError(
        "AerSimulator not available. Install qiskit-aer (pip install qiskit-aer)."
    ) from e

# DNA alphabet
DNA_BASES = "ACGT"
BASE_TO_INT = {"A": 0, "C": 1, "G": 2, "T": 3}
INT_TO_BASE = {0: "A", 1: "C", 2: "G", 3: "T"}


def kmer_to_index(kmer: str) -> int:
    """Encode a k-mer as an integer index.

    Uses base-4 encoding: A=0, C=1, G=2, T=3
    The first base is the most significant digit.

    Args:
        kmer: DNA k-mer (ACGT only)

    Returns:
        Integer index in range [0, 4^k)

    Example:
        >>> kmer_to_index("AA")
        0
        >>> kmer_to_index("AC")
        1
        >>> kmer_to_index("TT")
        15
    """
    kmer = kmer.upper()
    index = 0
    for base in kmer:
        if base not in BASE_TO_INT:
            raise ValueError(f"Invalid base in k-mer: {base}")
        index = index * 4 + BASE_TO_INT[base]
    return index


def index_to_kmer(idx: int, k: int) -> str:
    """Decode an integer index to a k-mer.

    Args:
        idx: Integer index
        k: Length of k-mer

    Returns:
        DNA k-mer string

    Example:
        >>> index_to_kmer(0, 2)
        'AA'
        >>> index_to_kmer(15, 2)
        'TT'
    """
    if idx < 0 or idx >= 4**k:
        raise ValueError(f"Index {idx} out of range for k={k}")

    bases = []
    for _ in range(k):
        bases.append(INT_TO_BASE[idx % 4])
        idx //= 4
    return "".join(reversed(bases))


def enumerate_kmers(k: int) -> List[str]:
    """Generate all possible k-mers of length k.

    Args:
        k: k-mer length

    Returns:
        List of all 4^k possible k-mers
    """
    return [index_to_kmer(i, k) for i in range(4**k)]


def count_kmer_occurrences(
    genome: str, k: int, include_revcomp: bool = True
) -> Dict[str, List[int]]:
    """Count occurrences of each k-mer in a genome.

    Args:
        genome: DNA sequence
        k: k-mer length
        include_revcomp: If True, count reverse complement as same k-mer

    Returns:
        Dict mapping k-mer -> list of 0-based positions
    """
    genome = genome.upper()
    occurrences: Dict[str, List[int]] = defaultdict(list)

    for i in range(len(genome) - k + 1):
        kmer = genome[i : i + k]
        # Skip k-mers with non-ACGT bases
        if all(b in DNA_BASES for b in kmer):
            occurrences[kmer].append(i)
            if include_revcomp:
                revcomp = get_reverse_complement(kmer)
                if revcomp != kmer:  # Avoid double-counting palindromes
                    occurrences[revcomp].append(i)

    return dict(occurrences)


def build_significance_oracle(
    genome: str, k: int, min_count: int, include_revcomp: bool = True
) -> Tuple[List[int], Dict[str, List[int]]]:
    """Build oracle marking k-mers that appear >= min_count times.

    Args:
        genome: DNA sequence
        k: k-mer length
        min_count: Minimum occurrence threshold
        include_revcomp: Count reverse complement as same k-mer

    Returns:
        Tuple of:
        - List of marked k-mer indices (for Grover oracle)
        - Dict of k-mer -> positions (for downstream analysis)
    """
    occurrences = count_kmer_occurrences(genome, k, include_revcomp)

    marked_indices = []
    significant_kmers = {}

    for kmer, positions in occurrences.items():
        if len(positions) >= min_count:
            idx = kmer_to_index(kmer)
            marked_indices.append(idx)
            significant_kmers[kmer] = positions

    return sorted(set(marked_indices)), significant_kmers


def grover_discover_motifs(
    contig: str,
    genome: str,
    k: int,
    min_count: int = 2,
    topk: int = 10,
    include_revcomp: bool = True,
    progress_every: int = 5,
    force_iters: Optional[int] = None,
    optimization_level: int = 1,
) -> List[Tuple[str, int, List[int]]]:
    """Discover significant motifs using Grover's algorithm.

    This function:
    1. Classically identifies k-mers with count >= min_count
    2. Builds a quantum oracle marking those k-mers
    3. Runs Grover iterations to amplify significant k-mers
    4. Returns top-k motifs with their counts and positions

    Args:
        contig: Contig/chromosome name
        genome: DNA sequence
        k: k-mer length (4-12 recommended)
        min_count: Minimum occurrences for significance
        topk: Number of top motifs to return
        include_revcomp: Count reverse complement as same k-mer
        progress_every: Print progress every N iterations
        force_iters: Override automatic iteration count
        optimization_level: Qiskit transpile optimization level

    Returns:
        List of (kmer, count, positions) tuples, sorted by probability
    """
    genome = genome.upper()

    log(f"Discovering {k}-mers with count >= {min_count}")
    log(f"Genome length: {len(genome)}, search space: 4^{k} = {4**k} k-mers")

    # Build oracle
    marked_indices, significant_kmers = build_significance_oracle(
        genome, k, min_count, include_revcomp
    )

    if not marked_indices:
        log("No k-mers found meeting significance threshold")
        return []

    M = len(marked_indices)
    log(f"Found {M} significant k-mers (meeting threshold)")

    # Calculate qubits needed: 2 bits per base
    n = 2 * k  # 4^k states = 2^(2k)
    N = 4**k

    if n > 25:
        log(f"WARNING: {n} qubits required, simulation may be slow")

    # Optimal Grover iterations
    if force_iters is None:
        iters = max(1, int(round((math.pi / 4.0) * math.sqrt(N / M))))
    else:
        iters = max(1, int(force_iters))

    log(f"Running Grover: qubits={n}, marked={M}, iterations={iters}")

    # Build circuit
    qc = QuantumCircuit(n)
    qc.h(range(n))

    for iteration in range(1, iters + 1):
        apply_mark_indices_phase_oracle(qc, n, marked_indices)
        apply_diffuser(qc, n)
        if progress_every and (iteration % progress_every == 0 or iteration == iters):
            log(f"  Grover iteration {iteration}/{iters} completed")

    qc.save_statevector()

    # Simulate
    sim = AerSimulator(method="statevector")
    tqc = transpile(qc, sim, optimization_level=optimization_level)
    result = sim.run(tqc).result()
    sv = result.get_statevector(tqc)

    amps = np.asarray(sv, dtype=complex)
    probs = (np.abs(amps) ** 2).real

    # Get top-k results
    top_indices = probs.argsort()[-topk:][::-1]

    results = []
    log(f"\nTop-{topk} discovered motifs:")
    log("-" * 60)

    for rank, idx in enumerate(top_indices, 1):
        idx = int(idx)
        prob = float(probs[idx])

        if idx >= N:
            continue  # Skip padding states

        kmer = index_to_kmer(idx, k)
        positions = significant_kmers.get(kmer, [])
        count = len(positions)

        # Check if this is a marked (significant) k-mer
        is_significant = idx in marked_indices
        tag = "SIG" if is_significant else "---"

        results.append((kmer, count, positions))

        # Format position preview
        if positions:
            pos_preview = ", ".join(str(p) for p in positions[:3])
            if len(positions) > 3:
                pos_preview += f", ... ({len(positions)} total)"
        else:
            pos_preview = "none"

        log(
            f"  #{rank:2d}: {kmer} | prob={prob:.6f} | count={count:4d} | {tag} | positions: {pos_preview}"
        )

    log("-" * 60)

    return results


def discover_and_report(
    contig: str,
    genome: str,
    k: int,
    min_count: int,
    topk: int = 10,
    include_revcomp: bool = True,
) -> None:
    """Run discovery and print a detailed report.

    Wrapper around grover_discover_motifs with formatted output.
    """
    log(f"=== Quantum Motif Discovery ===")
    log(f"Contig: {contig}")
    log(f"K-mer length: {k}")
    log(f"Minimum count threshold: {min_count}")
    log(f"Include reverse complement: {include_revcomp}")
    log("")

    results = grover_discover_motifs(
        contig=contig,
        genome=genome,
        k=k,
        min_count=min_count,
        topk=topk,
        include_revcomp=include_revcomp,
    )

    if not results:
        log("No significant motifs discovered.")
        return

    log("\n=== Discovery Summary ===")
    log(f"Total significant motifs found: {len(results)}")

    # Group by count
    by_count: Dict[int, List[str]] = defaultdict(list)
    for kmer, count, _ in results:
        by_count[count].append(kmer)

    log("\nMotifs by occurrence count:")
    for count in sorted(by_count.keys(), reverse=True):
        kmers = by_count[count]
        log(f"  Count {count}: {', '.join(kmers[:5])}" + (" ..." if len(kmers) > 5 else ""))
