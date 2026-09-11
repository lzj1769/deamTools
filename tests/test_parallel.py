"""Worker processes: run_jobs, and each command's output under --threads > 1.

Every command used to fan out over a thread pool, which the GIL serialised.
They now use worker processes. The risks a process pool adds are pickling --
jobs and results cross a process boundary -- and ordering, since results come
back as they finish. So each command is run with one worker (in-process) and
with several (a real pool), over a reference with enough contigs that the pool
actually has more than one job, and the outputs must be identical.
"""

import json
import os
import random
from functools import partial

import numpy as np
import pyBigWig
import pysam
import pytest

from deamtools.footprint import run_footprint
from deamtools.preprocessing.bam2bw import run_bam2bw
from deamtools.preprocessing.bam2fragment import run_bam2fragment
from deamtools.qc import run_qc
from deamtools.utils import run_jobs

CHROMS = ["chr1", "chr2", "chr3"]
LENGTH = 400


class TestRunJobs:
    def test_serial_path_returns_every_result(self):
        jobs = [partial(pow, 2, i) for i in range(5)]
        assert sorted(run_jobs(jobs, 1)) == [1, 2, 4, 8, 16]

    def test_pool_path_returns_results_in_submission_order(self):
        # Builtins pickle cleanly, so this exercises the pool itself. Order is
        # part of the contract: callers fold floats, and addition order shows.
        jobs = [partial(pow, 2, i) for i in range(8)]
        assert list(run_jobs(jobs, 4)) == [2**i for i in range(8)]

    def test_single_job_runs_in_process(self):
        assert list(run_jobs([partial(pow, 3, 2)], 8)) == [9]

    def test_empty(self):
        assert list(run_jobs([], 4)) == []

    def test_worker_exception_propagates(self):
        with pytest.raises(ValueError):
            list(run_jobs([partial(int, "1"), partial(int, "not a number")], 2))


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    """A three-contig reference and a BAM of edited pairs and singles on each."""
    d = tmp_path_factory.mktemp("parallel")
    rng = random.Random(0)
    seqs = {c: "".join(rng.choice("ACGT") for _ in range(LENGTH)) for c in CHROMS}

    fasta = str(d / "ref.fa")
    with open(fasta, "w") as f:
        for c, s in seqs.items():
            f.write(f">{c}\n{s}\n")
    pysam.faidx(fasta)

    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": c, "LN": LENGTH} for c in CHROMS]}

    def edited(s):
        return "".join(
            ("T" if b == "C" else "A") if b in "CG" and rng.random() < 0.3 else b
            for b in s
        )

    def record(name, ref_id, pos, seq, flag, mate_pos=None, tlen=0):
        a = pysam.AlignedSegment()
        a.query_name, a.query_sequence, a.flag = name, seq, flag
        a.reference_id, a.reference_start, a.mapping_quality = ref_id, pos, 60
        a.cigar = [(0, len(seq))]
        a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
        if mate_pos is not None:
            a.next_reference_id, a.next_reference_start = ref_id, mate_pos
            a.template_length = tlen
        return a

    tmp = str(d / "u.bam")
    with pysam.AlignmentFile(tmp, "wb", header=header) as bam:
        for ref_id, c in enumerate(CHROMS):
            s = seqs[c]
            for i in range(40):
                p1 = rng.randint(0, LENGTH - 150)
                p2 = p1 + rng.randint(0, 60)  # mates often overlap
                frag = edited(s[p1 : p2 + 80])
                r1, r2 = frag[:80], frag[p2 - p1 : p2 - p1 + 80]
                tlen = p2 + 80 - p1
                bam.write(
                    record(
                        f"{c}p{i}", ref_id, p1, r1, 0x1 | 0x2 | 0x40 | 0x20, p2, tlen
                    )
                )
                bam.write(
                    record(
                        f"{c}p{i}", ref_id, p2, r2, 0x1 | 0x2 | 0x80 | 0x10, p1, -tlen
                    )
                )
            for i in range(20):
                p = rng.randint(0, LENGTH - 80)
                bam.write(record(f"{c}s{i}", ref_id, p, edited(s[p : p + 80]), 0))
    path = str(d / "x.bam")
    pysam.sort("-o", path, tmp)
    pysam.index(path)
    return {"dir": d, "bam": path, "fasta": fasta}


def test_qc_is_identical_with_worker_processes(dataset, tmp_path):
    outs = []
    for workers in (1, 3):
        m = run_qc(
            dataset["bam"],
            dataset["fasta"],
            str(tmp_path / f"w{workers}"),
            "q",
            min_mapq=0,
            min_baseq=0,
            threads=workers,
            plot=False,
        )
        outs.append(json.dumps(m, sort_keys=True))
    assert outs[0] == outs[1]
    # And the pool really had work on every contig.
    assert json.loads(outs[0])["fragments"]["total"] == 3 * (40 + 20)


def test_bam2bw_is_identical_with_worker_processes(dataset, tmp_path):
    tracks = []
    for workers in (1, 3):
        out = str(tmp_path / f"w{workers}")
        run_bam2bw(
            dataset["bam"],
            dataset["fasta"],
            out,
            "t",
            min_mapq=0,
            min_baseq=0,
            threads=workers,
        )
        with pyBigWig.open(os.path.join(out, "t.bw")) as bw:
            # Bases without an entry come back as NaN; NaN != NaN, so compare
            # with equal_nan rather than list equality.
            tracks.append({c: np.array(bw.values(c, 0, LENGTH)) for c in CHROMS})
    for c in CHROMS:
        assert np.array_equal(tracks[0][c], tracks[1][c], equal_nan=True)
    assert all(np.nansum(tracks[0][c]) > 0 for c in CHROMS)  # signal on every contig


def test_bam2fragment_is_identical_with_worker_processes(dataset, tmp_path):
    tables = []
    for workers in (1, 3):
        out = str(tmp_path / f"w{workers}")
        run_bam2fragment(
            dataset["bam"],
            dataset["fasta"],
            out,
            "f",
            min_mapq=0,
            min_baseq=0,
            threads=workers,
        )
        tables.append(open(os.path.join(out, "f.tsv")).read())
    assert tables[0] == tables[1]
    assert {line.split("\t")[0] for line in tables[0].splitlines()} == set(CHROMS)


def test_footprint_is_identical_with_worker_processes(dataset, tmp_path):
    bw_path = str(tmp_path / "sig.bw")
    rng = random.Random(1)
    with pyBigWig.open(bw_path, "w") as bw:
        bw.addHeader([(c, LENGTH) for c in CHROMS])
        for c in CHROMS:
            bw.addEntries(
                c,
                list(range(LENGTH)),
                values=[float(rng.randint(0, 5)) for _ in range(LENGTH)],
                span=1,
            )
    bed = str(tmp_path / "sites.bed")
    with open(bed, "w") as f:
        for c in CHROMS:
            for i in range(5):
                s = 50 + i * 60
                f.write(f"{c}\t{s}\t{s + 10}\tsite{i}\n")
    results = []
    for workers in (1, 3):
        out = str(tmp_path / f"w{workers}")
        run_footprint(bw_path, bed, out, "fp", n_shuffles=50, threads=workers, seed=7)
        results.append(open(os.path.join(out, "fp.bed")).read())
    assert results[0] == results[1]
    assert len(results[0].splitlines()) == 15


class TestUnimportableMain:
    """A program read from stdin cannot be re-imported by spawned workers."""

    def test_falls_back_to_in_process(self, monkeypatch, caplog):
        import multiprocessing
        import sys
        import types

        fake_main = types.ModuleType("__main__")
        fake_main.__file__ = "<stdin>"  # what Python records for `python -`
        fake_main.__spec__ = None
        monkeypatch.setitem(sys.modules, "__main__", fake_main)
        monkeypatch.setattr(
            multiprocessing, "get_start_method", lambda allow_none=False: "spawn"
        )
        jobs = [partial(pow, 2, i) for i in range(4)]
        with caplog.at_level("WARNING"):
            assert list(run_jobs(jobs, 4)) == [1, 2, 4, 8]
        assert "read from stdin" in caplog.text

    def test_real_main_still_uses_the_pool(self, monkeypatch):
        from deamtools.utils import parallel

        # Under pytest __main__ is a real file (or a -m module), so no fallback.
        assert parallel._workers_can_start()
