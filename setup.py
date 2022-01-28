""""
PyTASER: Post-DFT TAS spectrum classification tool.
"""

from setuptools import find_packages, setup

with open("README.md") as file:
    long_description = file.read()

setup(
    name="pytaser",
    version="0.1.0",
    description="Post-DFT TAS spectrum classification tool",
    url="https://github.com/WMD-group/PyTASER",
    author="Savyasanchi Aggarwal",
    author_email="sa13018@ic.ac.uk",
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
        "matplotlib",
        "pymatgen>=2017.12.30",
    ],
    data_files=["LICENSE"],
)
