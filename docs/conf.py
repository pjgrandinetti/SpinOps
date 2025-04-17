# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'spinOps'
copyright = '2025, Philip Grandinetti'
author = 'Philip Grandinetti'
release = '0.1.1'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',    # Automatically document modules
    'sphinx.ext.napoleon',  # Support for NumPy/Google-style docstrings
    'sphinx.ext.viewcode',  # Add links to source code
    # 'myst_parser',         # Uncomment if using Markdown files
]

templates_path = ['_templates']
exclude_patterns = []

# Set the master document
master_doc = 'index'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']  # Ensure this directory exists or comment it out

# Add the module path
import os
import sys
sys.path.insert(0, os.path.abspath('../spinOps'))

