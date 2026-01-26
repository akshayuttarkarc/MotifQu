from motifqu.fasta import read_fasta


def test_read_fasta(tmp_path):
    p = tmp_path / "t.fa"
    p.write_text(">chrTest\nACGT\nAC\n", encoding="utf-8")
    contig, seq = read_fasta(str(p))
    assert contig == "chrTest"
    assert seq == "ACGTAC"
