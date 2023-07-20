[![Documentation Status](https://readthedocs.org/projects/pytaser/badge/?version=latest)](https://pytaser.readthedocs.io/en/latest/?badge=latest)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![tests](https://github.com/WMD-Group/PyTASER/actions/workflows/build_and_test.yml/badge.svg)](https://github.com/WMD-Group/PyTASER/actions)
[![pypi](https://img.shields.io/pypi/v/pytaser)](https://pypi.org/project/pytaser)
[![black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://img.shields.io/badge/code%20style-black-000000.svg)

# PyTASER

[![PyTASER Docs](docs/PyTASER.png)](https://pytaser.readthedocs.io/en/latest/)

[Official Documentation](https://pytaser.readthedocs.io/en/latest/)

PyTASER is a Python (3.9+) library for simulating differential absorption in general, and transient absorption spectroscopy (TAS) in particular, features from 
quantum chemical calculations of crystals. The goal of this library is to simulate spectra for comparison with and interpretation of 
experiments. The main features include:

* Spectra for insulating and semiconducting crystals
* Components from individual band-to-band transitions
* Spectra for different conditions: temperature and carrier concentrations
* Consideration of non-magnetic and magnetic materials
* Capability to input calculated bandstructure and density of states inputs with support for the [Materials Project](https://materialsproject.org)

# Installation

To install the module with pip (recommended): 

```
pip install pytaser
```
To install directly from the git repository:
```
pip install git+https://github.com/WMD-group/PyTASER
```
To do a manual build and installation:
```
python3 setup.py build
python3 setup.py install
```

PyTASER is compatible with Python 3.9+ and relies on a number of open-source packages, specifically:

* [pymatgen](https://pymatgen.org/index.html) 
* [numpy](https://numpy.org/), [scipy](https://scipy.org/) for data structures and unit conversion
* [matplotlib](https://matplotlib.org/) for plotting the spectra

# Developer’s installation (optional)

For development work, PyTASER can also be installed from a copy of the source directory:

Download PyTASER source code using the command:
```
git clone https://github.com/WMD-group/PyTASER
```
Navigate to root directory:
```
cd PyTASER
```
Install the code with the command:
```
pip install -e .
```
This command tries to obtain the required packages and their dependencies and install them automatically.

# Visualisation 

One method is to generate a Jupyter Notebook, as shown in the [examples](https://github.com/WMD-group/PyTASER/blob/main/examples) folder.
Alternatively, you can set up a file in Python to run it in the command line of the terminal:
```
python3 <filename.py>
```

If using your own, locally calculated data, please follow [this](https://github.com/WMD-group/PyTASER/blob/main/examples/PyTASER_DFT_Example.ipynb) workflow.
Otherwise, use the MP-integrated workflow [here](https://github.com/WMD-group/PyTASER/blob/main/examples/PyTASER_MP_Example.ipynb). 

# Contributing

We appreciate any contributions in the form of a pull request. Please see the [Contributing](https://pytaser.readthedocs.io/en/latest/contributing.html) documentation for more details.
Additional analysis/example spectra performed with `PyTASER` would be welcomed.

Please feel free to reach out to us via the [Issue Tracker](https://github.com/WMD-group/PyTASER/issues) if there are any questions or suggestions. 

# Testing

Unit tests are in the `tests` directory and can be run from the top directory using pytest. Automatic testing is run on the master and develop branches using GH Actions. 

Please run tests and add new tests for any new features whenever submitting pull requests.

# Future Work

Future topics we plan to build on:

* Incorporating finite-temperature effects (particularly for indirect band-to-band transitions)
* Direct treatment of pump-probe time delay and relaxation kinetics 
* Incorporating spin-flip processes for spin-polarised systems
* Description of defective crystals 

# Acknowledgements

The project has been developed by @savya10, @kavanase and @youngwonwoo under the supervision of @utf and @aronwalsh.
