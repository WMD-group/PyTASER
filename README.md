# PyTASER

PyTASER is a Python (3.9+) library and set of command-line tools for classifying spectral features in bulk materials, post-DFT. The goal of this library is to provide simulated TAS spectra, for qualitative comparison with experimental spectra. This will simplify the complex and tedious processes of TAS prep and analysis, especially for pristine materials. The main features include:

* An interactive TAS spectrum for a pristine semiconducting material
* Isolated spectra for individual band transitions from the overall TAS spectrum for the material.
* Spectra in different conditions: temperature, carrier concentrations (analogous to pump-probe time delay) 
* Partial occupancies of valence and conduction bands, using the Fermi-Dirac distribution for different Quasi-Fermi levels. 
* Consideration of both non-magnetic and magnetic materials. 
* Capability to input DFT-calculated bandstructure and dos inputs, with primary support for the [Materials Project](https://materialsproject.org/).

# Installation

The recommended way to install PyTASER is in a conda environment. 

```
Installation method to be updated here
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

***The library is currently undergoing some final changes before it is finalised.*** However, once it is completed, we would greatly appreciate any contributions in the form of a pull request. 
Additionally, any test cases/example spectra performed with PyTASER would be welcomed. Please feel free to reach out to us if there are any questions or suggestions. 

Future topics we'd like to build on:

* Direct conversion between carrier concentrations and pump-probe time delay (for more quantitative simulated-experimental comparisons)
* Incorporating spin-change processes (e.g. moving from Spin.up to Spin.down and vice-versa) for spin-polarised systems
* Incorporating finite-temperature effects (particularly indirect bandgaps & phonons, and defects)
* Incorporating more complex optical processes (e.g. Stimulated Emission)
* Cleaning regions further away from the bandgap
* Implementing optical transition probabilities alongside the JDOS
* Creating a similarly-obtained kinetics plot for TAS analysis.
* Relating spectral features with associated optical processes


# Acknowledgements

Developed by Savyasanchi Aggarwal, Alex Ganose and Liam Harnett-Caulfield. Aron Walsh designed and led the project. 

Thanks to the WMD group @ Imperial/Yonsei for all the interesting discussions and improvements!


