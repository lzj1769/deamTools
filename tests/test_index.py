"""Tests for deamtools.align.index (pure-Python FASTA conversion).

These do not require ``bwa`` or ``samtools``; they exercise the doubly-converted
reference writer directly, and the path logic of ``run_index`` with the external
tools monkeypatched out.
"""

import os

import deamtools.align.index as index_mod
from deamtools.align.index import LINE_WIDTH, _convert_fasta, run_index


def _parse_fasta(path: str) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:]
                seqs[name] = []
            else:
                assert name is not None
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def _write_fasta(path: str, records: dict[str, str]) -> None:
    with open(path, "w") as f:
        for name, seq in records.items():
            f.write(f">{name}\n{seq}\n")


class TestConvertFasta:
    def test_emits_f_and_r_entries_per_chrom(self, tmp_path):
        ref = str(tmp_path / "ref.fa")
        out = str(tmp_path / "ref.c2t")
        _write_fasta(ref, {"chr1": "ACGTCN"})
        _convert_fasta(ref, out)

        seqs = _parse_fasta(out)
        assert set(seqs) == {"fchr1", "rchr1"}
        # f = C->T conversion of the forward sequence
        assert seqs["fchr1"] == "ATGTTN"
        # r = G->A conversion of the forward sequence
        assert seqs["rchr1"] == "ACATCN"

    def test_multiple_chromosomes(self, tmp_path):
        ref = str(tmp_path / "ref.fa")
        out = str(tmp_path / "ref.c2t")
        _write_fasta(ref, {"chr1": "CCCC", "chr2": "GGGG"})
        _convert_fasta(ref, out)

        seqs = _parse_fasta(out)
        assert set(seqs) == {"fchr1", "rchr1", "fchr2", "rchr2"}
        assert seqs["fchr1"] == "TTTT"  # C->T
        assert seqs["rchr1"] == "CCCC"  # no G to convert
        assert seqs["fchr2"] == "GGGG"  # no C to convert
        assert seqs["rchr2"] == "AAAA"  # G->A

    def test_header_uses_first_token_only(self, tmp_path):
        ref = str(tmp_path / "ref.fa")
        out = str(tmp_path / "ref.c2t")
        # FASTA header with a trailing description.
        with open(ref, "w") as f:
            f.write(">chr1 some description here\nACGT\n")
        _convert_fasta(ref, out)
        seqs = _parse_fasta(out)
        assert set(seqs) == {"fchr1", "rchr1"}

    def test_lowercase_is_converted(self, tmp_path):
        ref = str(tmp_path / "ref.fa")
        out = str(tmp_path / "ref.c2t")
        _write_fasta(ref, {"chr1": "acgtACGT"})
        _convert_fasta(ref, out)
        seqs = _parse_fasta(out)
        assert seqs["fchr1"] == "atgtATGT"  # c->t and C->T
        assert seqs["rchr1"] == "acatACAT"  # g->a and G->A

    def test_run_index_default_writes_next_to_fasta(self, tmp_path, monkeypatch):
        # Monkeypatch the external tools so run_index works without bwa/samtools.
        monkeypatch.setattr(index_mod.shutil, "which", lambda name: "/usr/bin/" + name)

        def fake_run(cmd, check=False):
            if cmd[:2] == ["samtools", "faidx"]:
                open(cmd[2] + ".fai", "w").close()
            elif cmd[:2] == ["bwa", "index"]:
                open(cmd[2] + ".bwt", "w").close()

        monkeypatch.setattr(index_mod.subprocess, "run", fake_run)

        ref = str(tmp_path / "ref.fa")
        with open(ref, "w") as f:
            f.write(">chr1\nACGTACGT\n")
        run_index(ref)

        # Default: converted index sits right next to the FASTA.
        assert os.path.exists(ref + ".fai")
        assert os.path.exists(ref + ".deamtools.c2t")
        assert os.path.exists(ref + ".deamtools.c2t.bwt")

    def test_run_index_custom_out_dir_and_name(self, tmp_path, monkeypatch):
        monkeypatch.setattr(index_mod.shutil, "which", lambda name: "/usr/bin/" + name)
        bwa_paths = []

        def fake_run(cmd, check=False):
            if cmd[:2] == ["samtools", "faidx"]:
                open(cmd[2] + ".fai", "w").close()
            elif cmd[:2] == ["bwa", "index"]:
                bwa_paths.append(cmd[2])
                open(cmd[2] + ".bwt", "w").close()

        monkeypatch.setattr(index_mod.subprocess, "run", fake_run)

        ref = str(tmp_path / "ref.fa")
        with open(ref, "w") as f:
            f.write(">chr1\nACGTACGT\n")
        out_dir = str(tmp_path / "idx")
        run_index(ref, out_dir=out_dir, out_name="myref")

        # Converted reference + BWA index go to <out_dir>/<out_name>.*
        converted = os.path.join(out_dir, "myref.deamtools.c2t")
        assert os.path.exists(converted)
        assert os.path.exists(converted + ".bwt")
        assert bwa_paths == [converted]
        # The .fai still lives next to the original FASTA.
        assert os.path.exists(ref + ".fai")
        assert not os.path.exists(os.path.join(out_dir, "myref.fai"))

    def test_sequence_lines_are_wrapped(self, tmp_path):
        ref = str(tmp_path / "ref.fa")
        out = str(tmp_path / "ref.c2t")
        long_seq = "ACGT" * 50  # 200 bp, longer than LINE_WIDTH (80)
        _write_fasta(ref, {"chr1": long_seq})
        _convert_fasta(ref, out)

        with open(out) as f:
            seq_lines = [ln.rstrip("\n") for ln in f if not ln.startswith(">")]
        assert all(len(line) <= LINE_WIDTH for line in seq_lines)
        assert any(len(line) == LINE_WIDTH for line in seq_lines)
        # Round-trip: reconstructed length is preserved.
        seqs = _parse_fasta(out)
        assert len(seqs["fchr1"]) == len(long_seq)
