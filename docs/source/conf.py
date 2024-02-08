# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys

from recommonmark.transform import AutoStructify

sys.path.insert(0, os.path.abspath("../../"))
# from pytaser import __version__


# -- Project information -----------------------------------------------------

project = "PyTASER"
copyright = "2022, Savyasanchi Aggarwal"
author = "Savyasanchi Aggarwal"

# The full version, including alpha/beta/rc tags
release = "2.3.0"

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx_rtd_theme",
    "sphinx_book_theme",
    "sphinx.ext.autosummary",
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx_copybutton",
    "sphinx_toggleprompt",
    "myst_nb",  # for jupyter notebooks
]

source_suffix = {
    ".rst": "restructuredtext",
    ".ipynb": "myst-nb",
}

myst_enable_extensions = [
    "html_admonition",
]

# Add any paths that contain templates here, relative to this directory.
# templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_book_theme"  # "sphinx_rtd_theme"
# "sphinx_book_theme"  #

html_theme_options = {
    "repository_url": "https://github.com/WMD-group/PyTASER",
    # "repository_branch": "main",
    "path_to_docs": "docs/source",
    "use_repository_button": True,
    "use_issues_button": True,
    "use_edit_page_button": True,  # add button to suggest edits
    "home_page_in_toc": True,
}

html_logo = "_static/PyTASER.png"
html_title = "PyTASER"

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
html_use_smartypants = True

html_context = {
    "display_github": True,
    "github_user": "WMD-group",
    "github_repo": "PyTASER",
    "github_version": "master",
    "conf_py_path": "/docs_rst/",
}
# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {
    "python": ("https://docs.python.org/3.10", None),
    "numpy": ("http://docs.scipy.org/doc/numpy/", None),
    "pymatgen": ("http://pymatgen.org/", None),
    "matplotlib": ("http://matplotlib.org", None),
}

# -- Options for autodoc -----------------------------------------------------
autoclass_content = "both"

# -- Options for nb extension -----------------------------------------------
nb_execution_mode = "off"
myst_heading_anchors = 2
github_doc_root = "https://github.com/executablebooks/MyST-Parser/tree/master/docs/"


def setup(app):
    """Add configuration for MyST parser."""
    app.add_config_value(
        "myst_parser_config",
        {
            "url_resolver": lambda url: github_doc_root + url,
            "auto_toc_tree_section": "Contents",
        },
        True,
    )
    app.add_transform(AutoStructify)
