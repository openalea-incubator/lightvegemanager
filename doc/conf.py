# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------
import os
import sys

sys.path.insert( 0, os.path.abspath("../src/lightvegemanager") )
sys.path.insert( 0, os.path.abspath("../src") )


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'LightVegeManager'
copyright = '2023, INRAE P3F'
author = 'INRAE P3F'
release = '0.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ["sphinx.ext.todo", "sphinx.ext.viewcode", "sphinx.ext.autodoc", "nbsphinx", "nbsphinx_link"]

templates_path = ['_templates']
exclude_patterns = ['doc/_build', 'Thumbs.db', '.DS_Store', 'README.rst']

autodoc_member_order = 'bysource'

locale_dirs = ['locale/']   # path is example but recommended.
gettext_compact = False     # optional.

language = 'en'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
