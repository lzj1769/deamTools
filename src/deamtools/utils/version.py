from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version


def get_version() -> str:
    """Return the installed package version for deamtools.

    Returns
    -------
    str
        Installed package version string. If package metadata is not
        available, returns a fallback development version.
    """
    try:
        return version("deamtools")
    except PackageNotFoundError:
        return "0.0.0"