# MotifQu

MotifQu is a Quantum motif search using Grover's algorithm. it reads a FASTA sequence, finds motif-hit indices classically, then constructs a Grover oracle that marks those indices and runs Grover iterations using Aer statevector simulation. This demonstrates amplitude amplification as a usage for Quantum motif search

## Install

```bash
pip install MotifQu
```

## Usage

Exact match:

```bash
motifqu --fasta genome.fa --motif GTTGTTGGAGAAG --mismatches 0
```

Interactive motif entry:

```bash
motifqu --fasta genome.fa --mismatches 1
```

## Coordinate output

MotifQu prints both:

- 1-based inclusive coordinates: contig:start-end
- 0-based half-open interval: [start,end)

These coordinates are relative to the FASTA sequence provided. If your FASTA is a slice of a reference genome, you must add the appropriate offset yourself.

## Notes

- If --mismatches > 0, multiple hits are common (M > 1). Grover iteration count depends on M.
- The oracle is hard-coded from classical hits. That is deliberate: building a coherent data-access oracle for a large genome is the real bottleneck.
