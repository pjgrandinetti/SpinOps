# pylint: disable=missing-module-docstring,invalid-name,redefined-builtin

import datetime
import os
import sys

sys.path.insert(0, os.path.abspath("../.."))

now = datetime.datetime.now()
year = now.year

project = "spinOps"
copyright = "2025, Philip Grandinetti"
author = "Philip Grandinetti"
release = "0.1.1"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.doctest",
    "matplotlib.sphinxext.plot_directive",
    "sphinx.ext.mathjax",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx_copybutton",
    # "sphinxcontrib.bibtex",
    "breathe",
    "sphinxjp.themes.basicstrap",
    "sphinx.ext.intersphinx",
    "sphinx_tabs.tabs",
    "sphinx.ext.todo",
    "versionwarning.extension",
]

templates_path = ["_templates"]
exclude_patterns: list[str] = []

# Set the master document
master_doc = "index"

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
# html_static_path = ['_static']

autodoc_member_order = "bysource"
