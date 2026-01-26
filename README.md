# MotifQu

MotifQu is a quantum motif search and discovery tool using Grover's algorithm. It provides two main functions:

1. **Motif Search**: Find occurrences of a known motif pattern in a genome
2. **Motif Discovery**: Discover all significant motifs (k-mers) in a genome using quantum amplitude amplification

## Install

```bash
pip install MotifQu
```

## Usage

### Discover Motifs (New!)

Find all significant k-mers in a genome:

```bash
# Discover 6-mers appearing at least 3 times
motifqu discover --fasta genome.fa -k 6 --min-count 3

# Discover 8-mers, show top 20 results
motifqu discover --fasta genome.fa -k 8 --min-count 2 --topk 20

# Ignore reverse complement
motifqu discover --fasta genome.fa -k 6 --min-count 3 --no-revcomp
```

### Search for Specific Motif

```bash
# Exact match
motifqu search --fasta genome.fa --motif GTTGTTGGAGAAG --mismatches 0

# Allow 1 mismatch
motifqu search --fasta genome.fa --motif TATAAA --mismatches 1
```

### List Known Biological Motifs

```bash
motifqu list-motifs
```

### Expand IUPAC Pattern

```bash
# Expand E-box pattern (CANNTG)
motifqu expand CANNTG
```

## Coordinate Output

MotifQu prints both:

- 1-based inclusive coordinates: contig:start-end
- 0-based half-open interval: [start,end)

These coordinates are relative to the FASTA sequence provided.

## Biological Context

The quantum motif discovery tool is designed for:

- **Transcription Factor Binding Sites (TFBS)** - identifying regulatory sequences
- **Repeat elements** - finding tandem repeats and microsatellites
- **Conserved sequences** - detecting evolutionarily preserved patterns

The algorithm uses Grover's search to amplify the probability of significant k-mers (those appearing >= threshold times), providing a quadratic speedup over classical enumeration for the 4^k k-mer search space.

## IBM Quantum Hardware (Optional)

Run on real quantum hardware with IBM Quantum access:

```bash
# Install IBM runtime support
pip install MotifQu[ibm]

# Run on IBM Quantum (requires API token)
export IBMQ_TOKEN="your-ibm-quantum-token"
motifqu discover --fasta genome.fa -k 8 --min-count 10 --backend ibm

# Or specify token directly
motifqu discover --fasta genome.fa -k 8 --min-count 10 --backend ibm --ibm-token YOUR_TOKEN

# Use specific backend
motifqu discover --fasta genome.fa -k 6 --backend ibm --ibm-backend ibm_brisbane --shots 8192
```

Get your IBM Quantum API token at: https://quantum.ibm.com/

## Output Files

Use `--output` to save results:

```bash
motifqu discover --fasta genome.fa -k 6 --min-count 5 --output ./results
```

Generated files:
- `grover_circuit.png` - Quantum circuit diagram
- `probability_distribution.png` - K-mer probabilities
- `genome_motifs.png` - Motif positions visualization
- `discovered_motifs.csv` - Results table
- `results.json` - Full results with metadata
- `summary_report.txt` - Human-readable report

## Notes

- For discovery, k-mer lengths 6-10bp are recommended (4^k states require 2k qubits)
- The oracle is built from classical pre-computation of k-mer counts
- Reverse complement is counted as the same motif by default
- Quantum advantage is greatest when few k-mers are significant (M << N)
