import argparse
import sys
import time

from motifqu.fasta import read_fasta
from motifqu.grover import grover_run_aer_statevector
from motifqu.discovery import grover_discover_motifs, index_to_kmer
from motifqu.bio_patterns import list_known_motifs, get_motif_info, expand_iupac
from motifqu.util import log


def cmd_search(args: argparse.Namespace) -> None:
    """Run motif search (original functionality)."""
    t0 = time.time()
    log(f"Loading FASTA: {args.fasta}")
    contig, genome = read_fasta(args.fasta)
    log(f"Loaded contig={contig} length={len(genome)}")

    motif = args.motif
    if motif is None:
        motif = input("Enter motif sequence (ACGTN): ").strip()

    if not motif:
        print("ERROR: motif is empty.", file=sys.stderr)
        raise SystemExit(1)

    if args.mismatches < 0:
        print("ERROR: mismatches must be >= 0.", file=sys.stderr)
        raise SystemExit(1)

    log("Running Grover (Aer statevector)...")
    grover_run_aer_statevector(
        contig=contig,
        genome=genome,
        motif=motif,
        mismatches=args.mismatches,
        topk=args.topk,
        progress_every=args.progress_every,
        force_iters=args.iters,
        optimization_level=args.opt_level,
    )

    log(f"TOTAL runtime: {time.time() - t0:.3f}s")


def cmd_discover(args: argparse.Namespace) -> None:
    """Run quantum motif discovery."""
    t0 = time.time()
    log(f"Loading FASTA: {args.fasta}")
    contig, genome = read_fasta(args.fasta)
    log(f"Loaded contig={contig} length={len(genome)}")

    k = args.kmer_length
    if k < 3 or k > 12:
        print(f"WARNING: k={k} may be suboptimal. Recommended range: 4-10bp", file=sys.stderr)

    log(f"=== Quantum Motif Discovery ===")
    log(f"K-mer length: {k}")
    log(f"Minimum count: {args.min_count}")
    log(f"Include reverse complement: {args.revcomp}")

    results, qc, probs, significant_kmers = grover_discover_motifs(
        contig=contig,
        genome=genome,
        k=k,
        min_count=args.min_count,
        topk=args.topk,
        include_revcomp=args.revcomp,
        progress_every=args.progress_every,
        force_iters=args.iters,
        optimization_level=args.opt_level,
        output_dir=args.output,
        backend=args.backend,
        ibm_token=args.ibm_token,
        ibm_backend=args.ibm_backend,
        shots=args.shots,
    )

    if not results:
        log("No significant motifs discovered.")
    else:
        log(f"\nDiscovered {len(results)} significant motifs")

    # Save outputs if output directory specified
    if args.output and results:
        from motifqu.visualization import (
            save_circuit_diagram,
            save_probability_histogram,
            save_genome_visualization,
            save_results_csv,
            save_results_json,
            save_summary_report,
        )

        log(f"\nSaving outputs to: {args.output}")

        # Save circuit diagram (simplified version for large circuits)
        if qc is not None:
            try:
                # For small circuits, save the actual circuit
                if qc.num_qubits <= 8:
                    save_circuit_diagram(qc, args.output, "grover_circuit.png")
                else:
                    # For larger circuits, create a simplified schematic
                    from qiskit import QuantumCircuit as QC
                    n = 2 * k
                    viz_qc = QC(n, name="Grover Discovery")
                    viz_qc.h(range(n))
                    viz_qc.barrier(label="Oracle")
                    viz_qc.barrier(label="Diffuser")
                    viz_qc.barrier(label=f"x{10} iters")
                    save_circuit_diagram(viz_qc, args.output, "grover_circuit.png")
            except Exception as e:
                log(f"Could not save circuit diagram: {e}")

        # Save probability histogram
        if probs is not None:
            top_indices = probs.argsort()[-args.topk:][::-1]
            labels = [index_to_kmer(int(i), k) for i in top_indices if i < 4**k]
            save_probability_histogram(
                probs, [int(i) for i in top_indices if i < 4**k],
                labels, args.output, "probability_distribution.png"
            )

        # Save genome visualization
        if significant_kmers:
            save_genome_visualization(
                genome, significant_kmers, args.output, "genome_motifs.png"
            )

        # Save results
        save_results_csv(results, args.output, "discovered_motifs.csv")
        save_results_json(
            results,
            {
                "contig": contig,
                "genome_length": len(genome),
                "k": k,
                "min_count": args.min_count,
                "include_revcomp": args.revcomp,
                "total_significant": len(significant_kmers),
            },
            args.output,
            "results.json",
        )
        save_summary_report(
            contig, len(genome), k, args.min_count, results, args.output
        )

        log(f"All outputs saved to: {args.output}")

    log(f"TOTAL runtime: {time.time() - t0:.3f}s")


def cmd_list_motifs(args: argparse.Namespace) -> None:
    """List known biological motifs."""
    motifs = list_known_motifs()
    print(f"\nKnown biological motifs ({len(motifs)} total):\n")
    print(f"{'Name':<20} {'Consensus':<15} {'Description'}")
    print("-" * 70)
    for name in motifs:
        info = get_motif_info(name)
        print(f"{name:<20} {info['consensus']:<15} {info['description']}")
    print()


def cmd_expand(args: argparse.Namespace) -> None:
    """Expand IUPAC pattern to all matching sequences."""
    pattern = args.pattern.upper()
    sequences = expand_iupac(pattern)
    print(f"\nPattern: {pattern}")
    print(f"Expansions ({len(sequences)} sequences):\n")
    for seq in sequences:
        print(f"  {seq}")
    print()


def main() -> None:
    ap = argparse.ArgumentParser(
        prog="motifqu",
        description="MotifQu: Quantum motif search and discovery using Grover's algorithm.",
    )
    subparsers = ap.add_subparsers(dest="command", help="Available commands")

    # === SEARCH subcommand ===
    search_parser = subparsers.add_parser(
        "search",
        help="Search for a specific motif in a genome (original functionality)",
    )
    search_parser.add_argument("--fasta", required=True, help="Path to FASTA file.")
    search_parser.add_argument("--motif", default=None, help="Motif sequence (ACGTN). If omitted, prompt.")
    search_parser.add_argument("--mismatches", type=int, default=0, help="Allowed Hamming mismatches (default 0).")
    search_parser.add_argument("--topk", type=int, default=5, help="Show top-k outcomes (default 5).")
    search_parser.add_argument("--progress-every", type=int, default=5, help="Progress print every N iterations.")
    search_parser.add_argument("--iters", type=int, default=None, help="Force Grover iteration count.")
    search_parser.add_argument("--opt-level", type=int, default=1, choices=[0, 1, 2, 3], help="Qiskit transpile optimization level.")
    search_parser.set_defaults(func=cmd_search)

    # === DISCOVER subcommand ===
    discover_parser = subparsers.add_parser(
        "discover",
        help="Discover all significant motifs in a genome using Grover's algorithm",
    )
    discover_parser.add_argument("--fasta", required=True, help="Path to FASTA file.")
    discover_parser.add_argument("-k", "--kmer-length", type=int, default=6, help="K-mer length (default 6, range 4-12).")
    discover_parser.add_argument("--min-count", type=int, default=2, help="Minimum occurrences for significance (default 2).")
    discover_parser.add_argument("--topk", type=int, default=10, help="Show top-k discovered motifs (default 10).")
    discover_parser.add_argument("--no-revcomp", dest="revcomp", action="store_false", help="Don't count reverse complement as same k-mer.")
    discover_parser.add_argument("--progress-every", type=int, default=5, help="Progress print every N iterations.")
    discover_parser.add_argument("--iters", type=int, default=None, help="Force Grover iteration count.")
    discover_parser.add_argument("--opt-level", type=int, default=1, choices=[0, 1, 2, 3], help="Qiskit transpile optimization level.")
    discover_parser.add_argument("-o", "--output", default=None, help="Output directory for results, plots, and circuit diagrams.")
    discover_parser.add_argument("--backend", default="aer", choices=["aer", "ibm"], help="Backend: 'aer' (simulator) or 'ibm' (real hardware).")
    discover_parser.add_argument("--ibm-token", default=None, help="IBM Quantum API token (or set IBMQ_TOKEN env var).")
    discover_parser.add_argument("--ibm-backend", default=None, help="Specific IBM backend name (e.g., 'ibm_brisbane').")
    discover_parser.add_argument("--shots", type=int, default=4096, help="Number of shots for hardware execution (default 4096).")
    discover_parser.set_defaults(func=cmd_discover, revcomp=True)

    # === LIST-MOTIFS subcommand ===
    list_parser = subparsers.add_parser(
        "list-motifs",
        help="List known biological motif patterns",
    )
    list_parser.set_defaults(func=cmd_list_motifs)

    # === EXPAND subcommand ===
    expand_parser = subparsers.add_parser(
        "expand",
        help="Expand IUPAC pattern to all matching DNA sequences",
    )
    expand_parser.add_argument("pattern", help="IUPAC pattern (e.g., CANNTG)")
    expand_parser.set_defaults(func=cmd_expand)

    args = ap.parse_args()

    # Handle no subcommand (backward compatibility)
    if args.command is None:
        # Check for legacy usage
        if len(sys.argv) > 1 and sys.argv[1].startswith("--"):
            print("Note: Direct usage is deprecated. Use 'motifqu search' or 'motifqu discover'.", file=sys.stderr)
            print("Run 'motifqu --help' for available commands.\n", file=sys.stderr)
        ap.print_help()
        raise SystemExit(0)

    args.func(args)


if __name__ == "__main__":
    main()
