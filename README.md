# PyTASER

PyTASER is a Python (3.9+) library and set of command-line tools for classifying spectral features in bulk materials, post-DFT. The goal of this library is to provide qualitative comparisons for experimental TAS spectra - a complex and tedious process, especially for pristine materials. The main features include:

* An interactive TAS spectrum for a pristine semiconducting material
* Isolating spectra for individual band transitions from the overall TAS spectrum for the material.
* Spectra in different conditions: temperature, carrier concentrations (analogous to pump-probe time delay) 
* Identifying partial occupancies of valence and conduction bands, using the Fermi-Dirac distribution for different Quasi-Fermi levels. 
* Considers both non-magnetic and magnetic materials. 
* Taking DFT-calculated bandstructure and dos inputs, with primary support for the [Materials Project](https://materialsproject.org/).

# Installation

The recommended way to install PyTASER is in a conda environment. 

```
Installation method to be updated here
```
PyTASER is currently compatible with Python 3.9+ and relies on a number of open-source python packages, specifically:

* [pymatgen](https://pymatgen.org/index.html) 
* [numpy](https://numpy.org/), [scipy](https://scipy.org/) for data structures and unit conversion
* [matplotlib](https://matplotlib.org/), [plotly](https://plotly.com/) for plotting the spectra
* [AbiPy](https://abinit.github.io/abipy/) for the gaussian function implemented in JDOS_simple


# Visualisation 

Once the library is installed, please setup a file as done in the examples provided. Then just run it as a python file:
```
python3 filename.py
```

# Contributing

***The library is currently undergoing some final changes before it is finalised.*** However, once it is completed, we would greatly appreciate any contributions in the form of a pull request. 
Additionally, any test cases performed with PyTASER would be welcomed. 

Future topics we'd like to build on:

* Converting between carrier concentrations and pump-probe time delay (for a more quantitative analysis)
* Incorporating spin-change processes (e.g. moving from Spin.up to Spin.down and vice-versa) for spin-polarised systems
* Incorporating finite-temperature effects (particularly with indirect bandgaps and phonons, and defects)
* Incorporating more complex optical processes (e.g. Stimulated Emissions)
* Cleaning the regions further away from the bandgap
* Implementing the optical transition probabilities alongside the JDOS
* Creating a kinetics plot for TAS analysis.


# Acknowledgements

Developed by Savyasanchi Aggarwal, Alex Ganose and Liam Harnett-Caulfield. Aron Walsh designed and led the project. 

Thanks to the whole WMD group @ Imperial/Yonsei for all the interesting discussions and improvements they suggested!


