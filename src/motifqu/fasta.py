from typing import List, Tuple


def read_fasta(path: str) -> Tuple[str, str]:
    """Read a FASTA file and return (contig_name, sequence).

    contig_name is the first header token after '>' (or 'FASTA' if absent).
    The sequence is the concatenation of non-header, non-empty lines.
    """
    name = "FASTA"
    seq_chunks: List[str] = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name == "FASTA":
                    hdr = line[1:].strip()
                    name = hdr.split()[0] if hdr else "FASTA"
                continue
            seq_chunks.append(line)

    seq = "".join(seq_chunks).upper()
    if not seq:
        raise ValueError(f"No sequence found in FASTA: {path}")

    return name, seq
