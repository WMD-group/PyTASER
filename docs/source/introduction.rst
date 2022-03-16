============
Introduction
============

PyTASER is a Python (3.9+) library and set of command-line tools for simulating transient absorption spectroscopy (TAS) features from electronic band structures. The goal of this library is to simulate TAS spectra for comparison with and interpretation of experimental spectra. The main features include:

* An interactive TAS spectrum for a pristine semiconducting crystal
* Isolated spectra for individual band transitions
* Spectra for different conditions: temperature and carrier concentrations 
* Consideration of non-magnetic and magnetic materials
* Capability to input calculated bandstructure and density of states inputs with support for https://materialsproject.org.

Background
===========

TAS is a powerful pump-probe tool to characterise the excited states of materials. It can be used to understand microscopic processes in photochemical and electrochemical transformations, including phenomena such as electron trapping and carrier recombination. 

The drawback is that TAS spectra are difficult to interpret, especially for crystals where the specific valence and conduction band structure can give rise to complex features. Our goal here is to predict TAS features from first-principles starting from the most simple models of static excitations through to the kinetics of relaxation of the excited state back to the ground state.

Installation
============

To install the module with pip *(recommended)*: 

..code::

	```pip install --user pytaser```


To install directly from the git repository:

..code::

	`pip install --user git+https://github.com/WMD-group/PyTASER`


To do a manual build and installation:

..code::

	`python3 setup.py build
	python3 setup.py install --user`

Dependencies
============

PyTASER is currently compatible with Python 3.9+ and relies on a number of open-source python packages, specifically:

* `pymatgen <https://pymatgen.org/index.html>`__
* `numpy <https://numpy.org/>`__, 'scipy <https://scipy.org/>`__ for data structures and unit conversion
* `matplotlib <https://matplotlib.org/>`__ for plotting the spectra

Visualisation 
=============

Once the library is installed, please setup a file as done in the examples provided. Then just run it as a python file:
```
python3 filename.py
```

Contributing
============

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

Acknowledgements
================

Developed by Savyasanchi Aggarwal, Alex Ganose and Liam Harnett-Caulfield. Aron Walsh helped to design the project. Thanks to group members at Imperial College and Yonsei University for all the interesting discussions and improvements!
