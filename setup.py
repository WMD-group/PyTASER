""""
PyTASER: transient absorption prediction tool.
"""


import pathlib

from setuptools import find_packages, setup

long_description = pathlib.Path("README.md").read_text()
setup(
    name="pytaser",
    version="2.3.0",
    description="TAS prediction tool",
    url="https://pytaser.readthedocs.io/en/latest/",
    author="Savyasanchi Aggarwal",
    author_email="savya10@gmail.com",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Information Technology",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering",
        "Topic :: Other/Nonlisted Topic",
        "Operating System :: OS Independent",
    ],
    keywords="dft tas semiconductor-physics solar-fuel materials-project metal-oxides",
    test_suite="nose.collector",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib>=3.7.1",
        "pymatgen>=2023.05.31",
        "mp-api!=0.34.0,!=0.34.1,!=0.34.2",
        # bug: boto3 added as an unnecessary requirement
        # (https://github.com/materialsproject/pymatgen/issues/3241,
        # https://github.com/materialsproject/api/pull/836)
    ],
    extras_require={
        "tests": [
            "pytest>=7.1.3",
            "pytest-mpl==0.15.1",
            "monty",
            "pathlib",
        ],
        "docs": [
            "sphinx",
            "sphinx-book-theme",
            "sphinx-rtd-theme",
            "sphinx_click",
            "sphinx_design",
            "nbsphinx",
            "nbsphinx_link",
            "sphinx_copybutton",
            "sphinx_toggleprompt",
            "sphinx-autobuild",
            "sphinx_minipres",
            "sphinx_tabs",
            "sphinx_togglebutton",
        ],
    },
    data_files=["LICENSE"],
)
