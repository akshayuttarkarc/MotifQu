"""Visualization and output utilities for MotifQu.

This module provides functions to:
- Draw and save quantum circuits
- Visualize motif positions on a genome
- Export results to CSV/JSON
"""

import csv
import json
import os
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from qiskit import QuantumCircuit

from motifqu.util import log


def ensure_output_dir(output_dir: str) -> Path:
    """Create output directory if it doesn't exist."""
    path = Path(output_dir)
    path.mkdir(parents=True, exist_ok=True)
    return path


def save_circuit_diagram(
    qc: QuantumCircuit,
    output_dir: str,
    filename: str = "quantum_circuit.png",
    style: str = "iqp",
) -> str:
    """Save quantum circuit diagram to file.

    Args:
        qc: Qiskit QuantumCircuit
        output_dir: Directory to save to
        filename: Output filename
        style: Circuit drawing style (iqp, bw, etc.)

    Returns:
        Path to saved file
    """
    outdir = ensure_output_dir(output_dir)
    filepath = outdir / filename

    fig = qc.draw(output="mpl", style=style, fold=-1)
    fig.savefig(filepath, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    log(f"Saved circuit diagram: {filepath}")
    return str(filepath)


def save_probability_histogram(
    probs: np.ndarray,
    top_indices: List[int],
    labels: List[str],
    output_dir: str,
    filename: str = "probability_distribution.png",
) -> str:
    """Save probability distribution histogram.

    Args:
        probs: Probability array from simulation
        top_indices: Indices of top results
        labels: Labels for top results (e.g., k-mer sequences)
        output_dir: Directory to save to
        filename: Output filename

    Returns:
        Path to saved file
    """
    outdir = ensure_output_dir(output_dir)
    filepath = outdir / filename

    fig, ax = plt.subplots(figsize=(12, 6))

    # Extract probabilities for top indices
    top_probs = [probs[i] for i in top_indices]

    # Create bar chart
    x = np.arange(len(labels))
    bars = ax.bar(x, top_probs, color="steelblue", edgecolor="navy", alpha=0.8)

    # Highlight significant hits
    for i, bar in enumerate(bars):
        if top_probs[i] > np.mean(top_probs):
            bar.set_color("darkgreen")
            bar.set_alpha(0.9)

    ax.set_xlabel("K-mer", fontsize=12)
    ax.set_ylabel("Probability", fontsize=12)
    ax.set_title("Grover Search - Top K-mer Probabilities", fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()
    fig.savefig(filepath, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    log(f"Saved probability histogram: {filepath}")
    return str(filepath)


def save_genome_visualization(
    genome: str,
    motif_positions: Dict[str, List[int]],
    output_dir: str,
    filename: str = "genome_motifs.png",
    max_positions: int = 100,
) -> str:
    """Visualize motif positions along the genome.

    Args:
        genome: Genome sequence
        motif_positions: Dict of motif -> list of positions
        output_dir: Directory to save to
        filename: Output filename
        max_positions: Maximum positions to show per motif

    Returns:
        Path to saved file
    """
    outdir = ensure_output_dir(output_dir)
    filepath = outdir / filename

    # Prepare data
    motifs = list(motif_positions.keys())[:20]  # Max 20 motifs
    genome_len = len(genome)

    fig, ax = plt.subplots(figsize=(14, max(4, len(motifs) * 0.5)))

    # Color palette
    colors = plt.cm.tab10(np.linspace(0, 1, len(motifs)))

    for i, motif in enumerate(motifs):
        positions = motif_positions[motif][:max_positions]
        y = [i] * len(positions)
        ax.scatter(positions, y, c=[colors[i]], s=30, alpha=0.7, label=motif)

    ax.set_xlim(0, genome_len)
    ax.set_ylim(-0.5, len(motifs) - 0.5)
    ax.set_yticks(range(len(motifs)))
    ax.set_yticklabels(motifs, fontsize=9)
    ax.set_xlabel("Genome Position (bp)", fontsize=12)
    ax.set_title("Motif Positions Across Genome", fontsize=14)
    ax.grid(axis="x", alpha=0.3)

    # Add genome length annotation
    ax.axvline(genome_len, color="red", linestyle="--", alpha=0.5)
    ax.text(genome_len * 0.98, len(motifs) - 0.3, f"{genome_len} bp",
            ha="right", fontsize=10, color="red")

    plt.tight_layout()
    fig.savefig(filepath, dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    log(f"Saved genome visualization: {filepath}")
    return str(filepath)


def save_results_csv(
    results: List[Tuple[str, int, List[int]]],
    output_dir: str,
    filename: str = "discovered_motifs.csv",
) -> str:
    """Save discovery results to CSV.

    Args:
        results: List of (kmer, count, positions) tuples
        output_dir: Directory to save to
        filename: Output filename

    Returns:
        Path to saved file
    """
    outdir = ensure_output_dir(output_dir)
    filepath = outdir / filename

    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["rank", "kmer", "count", "positions"])
        for rank, (kmer, count, positions) in enumerate(results, 1):
            pos_str = ";".join(str(p) for p in positions[:50])  # Limit positions
            if len(positions) > 50:
                pos_str += f";...({len(positions)} total)"
            writer.writerow([rank, kmer, count, pos_str])

    log(f"Saved results CSV: {filepath}")
    return str(filepath)


def save_results_json(
    results: List[Tuple[str, int, List[int]]],
    metadata: Dict,
    output_dir: str,
    filename: str = "results.json",
) -> str:
    """Save discovery results to JSON with metadata.

    Args:
        results: List of (kmer, count, positions) tuples
        metadata: Additional metadata (genome length, k, etc.)
        output_dir: Directory to save to
        filename: Output filename

    Returns:
        Path to saved file
    """
    outdir = ensure_output_dir(output_dir)
    filepath = outdir / filename

    output = {
        "metadata": metadata,
        "motifs": [
            {
                "rank": i + 1,
                "kmer": kmer,
                "count": count,
                "positions": positions,
            }
            for i, (kmer, count, positions) in enumerate(results)
        ],
    }

    with open(filepath, "w", encoding="utf-8") as f:
        json.dump(output, f, indent=2)

    log(f"Saved results JSON: {filepath}")
    return str(filepath)


def save_summary_report(
    contig: str,
    genome_len: int,
    k: int,
    min_count: int,
    results: List[Tuple[str, int, List[int]]],
    output_dir: str,
    filename: str = "summary_report.txt",
) -> str:
    """Save a human-readable summary report.

    Args:
        contig: Contig name
        genome_len: Genome length
        k: K-mer length
        min_count: Minimum count threshold
        results: Discovery results
        output_dir: Directory to save to
        filename: Output filename

    Returns:
        Path to saved file
    """
    outdir = ensure_output_dir(output_dir)
    filepath = outdir / filename

    with open(filepath, "w", encoding="utf-8") as f:
        f.write("=" * 60 + "\n")
        f.write("MotifQu - Quantum Motif Discovery Report\n")
        f.write("=" * 60 + "\n\n")

        f.write("PARAMETERS\n")
        f.write("-" * 40 + "\n")
        f.write(f"Contig:           {contig}\n")
        f.write(f"Genome length:    {genome_len:,} bp\n")
        f.write(f"K-mer length:     {k}\n")
        f.write(f"Min count:        {min_count}\n")
        f.write(f"Search space:     4^{k} = {4**k:,} k-mers\n")
        f.write(f"Qubits required:  {2*k}\n\n")

        f.write("DISCOVERED MOTIFS\n")
        f.write("-" * 40 + "\n")
        if results:
            f.write(f"{'Rank':<6} {'K-mer':<15} {'Count':<8} {'Positions'}\n")
            f.write("-" * 60 + "\n")
            for rank, (kmer, count, positions) in enumerate(results[:20], 1):
                pos_preview = ", ".join(str(p) for p in positions[:5])
                if len(positions) > 5:
                    pos_preview += f" ... ({len(positions)} total)"
                f.write(f"{rank:<6} {kmer:<15} {count:<8} {pos_preview}\n")
        else:
            f.write("No significant motifs discovered.\n")

        f.write("\n" + "=" * 60 + "\n")

    log(f"Saved summary report: {filepath}")
    return str(filepath)
