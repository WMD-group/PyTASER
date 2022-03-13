[![Documentation Status](https://readthedocs.org/projects/pytaser/badge/?version=latest)](https://pytaser.readthedocs.io/en/latest/?badge=latest)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# PyTASER

***The library is undergoing changes before it is released*** 

[Documentation](https://pytaser.readthedocs.io/en/latest/)

PyTASER is a Python (3.9+) library and set of command-line tools for simulating transient absorption spectroscopy (TAS) features from electronic band structures. The goal of this library is to simulate TAS spectra for comparison with and interpretation of experimental spectra. The main features include:

* An interactive TAS spectrum for a pristine semiconducting crystal
* Isolated spectra for individual band transitions
* Spectra for different conditions: temperature and carrier concentrations 
* Consideration of non-magnetic and magnetic materials
* Capability to input calculated bandstructure and density of states inputs with support for the [Materials Project](https://materialsproject.org/)

# Installation

To install the module with pip (recommended): 

```
pip install --user pytaser
```
To install directly from the git repository:
```
pip install --user git+https://github.com/WMD-group/PyTASER
```
To do a manual build and installation:
```
python3 setup.py build
python3 setup.py install --user
```

PyTASER is currently compatible with Python 3.9+ and relies on a number of open-source python packages, specifically:

* [pymatgen](https://pymatgen.org/index.html) 
* [numpy](https://numpy.org/), [scipy](https://scipy.org/) for data structures and unit conversion
* [matplotlib](https://matplotlib.org/) for plotting the spectra

# Visualisation 

Once the library is installed, please setup a file as done in the examples provided. Then just run it as a python file:
```
python3 filename.py
```

# Contributing

We appreciate any contributions in the form of a pull request. 
Additional test cases/example spectra performed with PyTASER would be welcomed. 
Please feel free to reach out to us if there are any questions or suggestions. 

Future topics we plan to build on:

* Direct treatment of pump-probe time delay
* Incorporating spin-flip processes for spin-polarised systems
* Incorporating finite-temperature effects (particularly indirect bandgaps & phonons, and defects)
* Incorporating more complex optical processes (e.g. stimulated emission)
* Cleaning regions further away from the bandgap
* Implementing optical transition probabilities

# Acknowledgements

Developed by Savyasanchi Aggarwal, Alex Ganose and Liam Harnett-Caulfield. Aron Walsh helped to design the project. Thanks to group members at Imperial College and Yonsei University for all the interesting discussions and improvements!
