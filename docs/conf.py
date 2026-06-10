"""Sphinx configuration for the DeamTools documentation.

The docs are written in Markdown and parsed with MyST (myst-parser), and
rendered with the Read the Docs theme — the same Sphinx + sphinx_rtd_theme
setup used by RGT (reg-gen.readthedocs.io).
"""

from importlib.metadata import PackageNotFoundError, version as _pkg_version

# -- Project information ------------------------------------------------------

project = "DeamTools"
author = "Zhijian Li"
copyright = "2026, Zhijian Li"

try:
    release = _pkg_version("deamtools")
except PackageNotFoundError:
    release = "0.1.0"
version = release

# -- General configuration ----------------------------------------------------

extensions = [
    "myst_parser",
]

# Parse both Markdown (MyST) and reStructuredText sources.
source_suffix = {
    ".md": "markdown",
    ".rst": "restructuredtext",
}

# MyST features: ::: fenced directives (admonitions), definition lists, and
# auto-generated anchors so cross-page header links resolve.
myst_enable_extensions = [
    "colon_fence",
    "deflist",
]
myst_heading_anchors = 3

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- HTML output --------------------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_title = "DeamTools"
html_theme_options = {
    "navigation_depth": 2,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "prev_next_buttons_location": "bottom",
}
