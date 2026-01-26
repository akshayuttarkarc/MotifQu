import math
from typing import List, Optional, Tuple

import numpy as np
from qiskit import QuantumCircuit, transpile

from motifqu.util import log

try:
    from qiskit_aer import AerSimulator
except Exception as e:
    raise ImportError(
        "AerSimulator not available. Install qiskit-aer (pip install qiskit-aer)."
    ) from e


def hamming(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def find_positions(genome: str, motif: str, max_mismatches: int = 0) -> List[int]:
    L = len(motif)
    if L == 0 or L > len(genome):
        return []
    hits: List[int] = []
    for i in range(len(genome) - L + 1):
        if hamming(genome[i : i + L], motif) <= max_mismatches:
            hits.append(i)
    return hits


def idx_to_coords(contig: str, idx0: int, motif_len: int) -> str:
    """Return both 1-based inclusive and 0-based half-open coordinates."""
    start0 = idx0
    end0 = idx0 + motif_len
    start1 = start0 + 1
    end1 = end0
    return f"{contig}:{start1}-{end1} (1-based) | [{start0},{end0}) (0-based)"


def apply_mark_indices_phase_oracle(qc: QuantumCircuit, n: int, marked: List[int]) -> None:
    """Hard-coded phase oracle that flips phase of |idx> for idx in marked."""
    if not marked:
        raise ValueError("No marked states (M=0).")

    controls = list(range(n - 1))
    target = n - 1

    for idx in marked:
        # Map |idx> -> |11..1>
        for q in range(n):
            if ((idx >> q) & 1) == 0:
                qc.x(q)

        # Multi-controlled Z via H - MCX - H
        qc.h(target)
        qc.mcx(controls, target)
        qc.h(target)

        # Uncompute map
        for q in range(n):
            if ((idx >> q) & 1) == 0:
                qc.x(q)


def apply_diffuser(qc: QuantumCircuit, n: int) -> None:
    """Standard Grover diffuser over n qubits."""
    qc.h(range(n))
    qc.x(range(n))

    controls = list(range(n - 1))
    target = n - 1
    qc.h(target)
    qc.mcx(controls, target)
    qc.h(target)

    qc.x(range(n))
    qc.h(range(n))


def grover_run_aer_statevector(
    contig: str,
    genome: str,
    motif: str,
    mismatches: int = 0,
    topk: int = 5,
    progress_every: int = 5,
    force_iters: Optional[int] = None,
    optimization_level: int = 1,
) -> Tuple[int, float, List[Tuple[int, float]]]:
    """Run the demo Grover search using Aer statevector simulation.

    Returns:
      best_idx (0-based window start), best_prob, topk list of (idx, prob)

    Notes:
      - The oracle is defined by a classical scan of the genome.
      - Coordinates are relative to the FASTA sequence provided.
    """
    genome = genome.upper()
    motif = motif.upper().strip()
    L = len(motif)

    hits = find_positions(genome, motif, max_mismatches=mismatches)
    if not hits:
        raise RuntimeError("No hits found; cannot run Grover with M=0.")

    num_positions = len(genome) - L + 1
    n = math.ceil(math.log2(num_positions))
    N = 2**n
    M = len(hits)

    if force_iters is None:
        iters = max(1, int(round((math.pi / 4.0) * math.sqrt(N / M))))
    else:
        iters = max(1, int(force_iters))

    log(f"hits_count={M}; qubits={n}; padded_N={N}; real_positions={num_positions}; iters={iters}")

    qc = QuantumCircuit(n)
    qc.h(range(n))
    for k in range(1, iters + 1):
        apply_mark_indices_phase_oracle(qc, n, hits)
        apply_diffuser(qc, n)
        if progress_every and (k % progress_every == 0 or k == iters):
            log(f"  Grover iteration {k}/{iters} completed")

    qc.save_statevector()

    sim = AerSimulator(method="statevector")
    tqc = transpile(qc, sim, optimization_level=optimization_level)
    result = sim.run(tqc).result()
    sv = result.get_statevector(tqc)

    amps = np.asarray(sv, dtype=complex)
    probs = (np.abs(amps) ** 2).real

    top = probs.argsort()[-topk:][::-1]
    ranked = [(int(i), float(probs[int(i)])) for i in top]

    best_idx, best_prob = ranked[0]

    # Print best + top-k with coordinates when meaningful
    if best_idx < num_positions:
        log(f"Top outcome: idx={best_idx} prob≈{best_prob:.6f} coord={idx_to_coords(contig, best_idx, L)}")
        log(f"Recovered window: {genome[best_idx:best_idx+L]}")
        log(f"Hamming distance: {hamming(genome[best_idx:best_idx+L], motif)}")
    else:
        log(f"Top outcome: idx={best_idx} prob≈{best_prob:.6f} (PAD region; no coordinate)")

    log(f"Top-{topk} outcomes:")
    for r, (idx, p) in enumerate(ranked, start=1):
        tags = []
        if idx in hits:
            tags.append("HIT")
        if idx >= num_positions:
            tags.append("PAD")
        tag_str = ",".join(tags) if tags else "-"
        if idx < num_positions:
            log(f"  #{r}: idx={idx:6d} prob≈{p:.6f} {tag_str} | {idx_to_coords(contig, idx, L)}")
        else:
            log(f"  #{r}: idx={idx:6d} prob≈{p:.6f} {tag_str}")

    return best_idx, best_prob, ranked
