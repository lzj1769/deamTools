"""Fan independent jobs out to worker processes.

Every command that takes ``--threads`` splits its work by chromosome or region
and runs a per-read (or per-site) loop that is pure Python. Under a thread pool
the GIL serialises those loops, so ``--threads`` bought nothing: ``qc`` with
``--threads 12`` on a 3.2 M-record BAM ran at 105% CPU. Worker *processes*
each have their own interpreter, so the loops genuinely run in parallel.

Workers already open their own pysam/pyBigWig handles (they were never
thread-safe either), so nothing is shared between them; the only requirement a
process pool adds is that jobs and their results be picklable, which is why
jobs are passed as :func:`functools.partial` objects over module-level
functions.
"""

from __future__ import annotations

import logging
import multiprocessing
import os
import sys
from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import ProcessPoolExecutor

logger = logging.getLogger(__name__)


def _workers_can_start() -> bool:
    """Can a spawned worker re-import this program's ``__main__``?

    Under the ``spawn`` and ``forkserver`` start methods (``spawn`` is the
    default on macOS) each worker starts a fresh interpreter and re-imports the
    main module by file path. When the program was fed in on stdin -- ``python
    -`` or a heredoc -- that path is the literal ``<stdin>``, the import fails
    and every worker dies, which surfaces as a ``BrokenProcessPool``. Detect
    that up front instead. ``python -m``, console scripts, ordinary script
    files and notebooks (whose ``__main__`` has no file) all pass.
    """
    if multiprocessing.get_start_method(allow_none=False) == "fork":
        return True  # forked workers inherit __main__; nothing is re-imported
    main = sys.modules.get("__main__")
    if main is None or getattr(main, "__spec__", None) is not None:
        return True  # run as `python -m ...`: re-imported by name
    path = getattr(main, "__file__", None)
    return path is None or os.path.exists(path)


def run_jobs[T](jobs: Sequence[Callable[[], T]], workers: int) -> Iterator[T]:
    """Yield the result of each job, in the order the jobs were given.

    ``jobs`` are zero-argument callables -- typically
    ``functools.partial(module_level_function, ...)`` -- so that they pickle
    across the process boundary.

    Results come back in **submission** order, not completion order, even
    though the jobs run concurrently. That is what makes the output
    deterministic: callers fold results together (``qc`` sums per-chromosome
    floats), and floating-point addition is not associative, so merging in
    completion order made the last digit of ``qc``'s mean edit rate change from
    run to run (by 1.4e-17 on a real BAM). Holding a finished result until the
    ones before it are done costs nothing here, since every caller accumulates
    everything before writing.

    With one worker, or one job, everything runs in this process: there is
    nothing to parallelise, and skipping the pool avoids both the start-up cost
    of spawning an interpreter (which re-imports numpy and pysam) and the need
    to pickle anything. That is also the path the test suite exercises by
    default. The same in-process path is taken, with a warning, when worker
    processes could not start at all (see :func:`_workers_can_start`).

    Library callers on macOS (and anywhere ``spawn`` is the start method) must
    keep the usual ``if __name__ == "__main__":`` guard around the call in
    their script, or each worker re-runs the script on import; Python reports
    that with a clear ``RuntimeError`` about the bootstrapping phase.
    """
    if workers > 1 and len(jobs) > 1 and not _workers_can_start():
        logger.warning(
            "Running %d job(s) in this process instead of %d workers: the "
            "program was read from stdin, so worker processes cannot re-import "
            "it. Save it to a file to run in parallel.",
            len(jobs),
            workers,
        )
        workers = 1
    if workers <= 1 or len(jobs) <= 1:
        for job in jobs:
            yield job()
        return
    with ProcessPoolExecutor(max_workers=min(workers, len(jobs))) as pool:
        futures = [pool.submit(job) for job in jobs]
        for future in futures:
            yield future.result()
