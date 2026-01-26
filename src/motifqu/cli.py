import argparse
import sys
import time

from motifqu.fasta import read_fasta
from motifqu.grover import grover_run_aer_statevector
from motifqu.discovery import grover_discover_motifs
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

    results = grover_discover_motifs(
        contig=contig,
        genome=genome,
        k=k,
        min_count=args.min_count,
        topk=args.topk,
        include_revcomp=args.revcomp,
        progress_every=args.progress_every,
        force_iters=args.iters,
        optimization_level=args.opt_level,
    )

    if not results:
        log("No significant motifs discovered.")
    else:
        log(f"\nDiscovered {len(results)} significant motifs")

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
