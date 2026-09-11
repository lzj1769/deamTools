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

from collections.abc import Callable, Iterator, Sequence
from concurrent.futures import ProcessPoolExecutor


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
    default.
    """
    if workers <= 1 or len(jobs) <= 1:
        for job in jobs:
            yield job()
        return
    with ProcessPoolExecutor(max_workers=min(workers, len(jobs))) as pool:
        futures = [pool.submit(job) for job in jobs]
        for future in futures:
            yield future.result()
