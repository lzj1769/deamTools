"""BED region parsing utilities.

This module centralises the loader used by :mod:`deamtools.preprocessing.bam2bw`,
:mod:`deamtools.qc`, and any other consumer that needs to read a
UCSC BED file. See https://en.wikipedia.org/wiki/BED_(file_format).
"""

from __future__ import annotations

import io

import pandas as pd

# Per the UCSC BED spec a record has 3 required and up to 9 optional columns.
# https://en.wikipedia.org/wiki/BED_(file_format)
BED_COLUMNS: tuple[str, ...] = (
    "chrom", "start", "end", "name", "score", "strand",
    "thickStart", "thickEnd", "itemRgb",
    "blockCount", "blockSizes", "blockStarts",
)


def _load_regions(bed_path: str) -> pd.DataFrame:
    """Load a BED file and return non-overlapping merged intervals as a DataFrame.

    Recognises the UCSC BED format: the first three columns (``chrom``,
    ``start``, ``end``) are required and 0-based half-open; up to nine
    optional columns may follow. Lines starting with ``#``, ``track``, or
    ``browser`` and blank lines are skipped.

    Parameters
    ----------
    bed_path : str
        Path to a BED file. The file may contain any subset of the 12 columns
        defined in the UCSC BED specification; only the first three are
        consumed here.

    Returns
    -------
    pandas.DataFrame
        DataFrame with three columns — ``chrom`` (``str``), ``start``
        (``int``), and ``end`` (``int``) — sorted by chromosome then
        position. Overlapping or adjacent intervals on the same chromosome
        are merged so that the returned set is non-overlapping; this
        prevents reads from being processed twice by consumers that iterate
        the returned intervals. An empty DataFrame with the same columns is
        returned when the file contains no data rows.

    Raises
    ------
    ValueError
        If any record has a non-integer ``start`` or ``end`` value, or if
        ``start > end`` for any interval.

    Notes
    -----
    Columns beyond ``end`` (``name``, ``score``, ``strand``, ``thickStart``,
    ``thickEnd``, ``itemRgb``, ``blockCount``, ``blockSizes``,
    ``blockStarts``) are silently dropped. The full list is exposed as
    :data:`BED_COLUMNS` for callers that need to round-trip them.

    References
    ----------
    UCSC Genome Browser, BED format specification:
    https://en.wikipedia.org/wiki/BED_(file_format)

    Examples
    --------
    >>> _load_regions("peaks.bed")
      chrom  start    end
    0  chr1    100    500
    1  chr1   1200   1800
    2  chr2     50    300
    """
    with open(bed_path) as f:
        cleaned = "".join(
            line for line in f
            if line.strip() and not line.startswith(("#", "track", "browser"))
        )

    if not cleaned:
        return pd.DataFrame(columns=["chrom", "start", "end"])

    df = pd.read_csv(
        io.StringIO(cleaned),
        sep="\t",
        header=None,
        usecols=[0, 1, 2],
        names=["chrom", "start", "end"],
        dtype={"chrom": str, "start": "Int64", "end": "Int64"},
    )

    if df[["start", "end"]].isna().any().any():
        raise ValueError(f"BED file has non-integer start/end values: {bed_path}")
    df = df.astype({"start": int, "end": int})

    bad = df["start"] > df["end"]
    if bad.any():
        row = df.loc[bad].iloc[0]
        raise ValueError(
            f"BED file has interval with start > end: "
            f"{row['chrom']}:{row['start']}-{row['end']} ({bed_path})"
        )

    df = df.sort_values(["chrom", "start", "end"], kind="stable").reset_index(drop=True)

    merged: list[tuple[str, int, int]] = []
    for chrom, group in df.groupby("chrom", sort=False):
        chrom = str(chrom)
        cur_start, cur_end = int(group.iat[0, 1]), int(group.iat[0, 2])
        for s, e in zip(group["start"].iloc[1:], group["end"].iloc[1:]):
            s, e = int(s), int(e)
            if s <= cur_end:
                cur_end = max(cur_end, e)
            else:
                merged.append((chrom, cur_start, cur_end))
                cur_start, cur_end = s, e
        merged.append((chrom, cur_start, cur_end))

    return pd.DataFrame(merged, columns=["chrom", "start", "end"])
