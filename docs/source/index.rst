.. PyTASER documentation master file, created by
   sphinx-quickstart on Sun Mar  6 03:03:51 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to PyTASER's documentation!
===================================

.. image:: https://github.com/WMD-Group/PyTASER/actions/workflows/build_and_test.yml/badge.svg
 :target: https://github.com/WMD-Group/PyTASER/actions

.. image:: https://readthedocs.org/projects/pytaser/badge/?version=latest&style=flat
 :target: https://pytaser.readthedocs.io/en/latest/

.. image:: https://img.shields.io/pypi/v/pytaser
 :target: https://pypi.org/project/pytaser

.. image:: https://img.shields.io/badge/Made%20with-Python-1f425f.svg
 :target: https://www.python.org/

.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
 :target: https://opensource.org/licenses/MIT

.. image:: https://img.shields.io/badge/code%20style-black-000000.svg
 :target: https://img.shields.io/badge/code%20style-black-000000.svg

============
Introduction
============

``PyTASER`` is a ``Python`` (3.9+) library built for simulating transient absorption spectroscopy (TAS) features from
DFT calculations. The goal of this library is to simulate TAS spectra for comparison with and interpretation of
experimental spectra. The main features include:

* An interactive TAS spectrum for a pristine semiconducting crystal
* Isolated spectra for individual band transitions
* Spectra for different conditions: temperature and carrier concentrations
* Consideration of non-magnetic and magnetic materials
* Capability to input calculated bandstructures and density of states inputs with support for https://materialsproject.org.

Background
===========

TAS is a powerful pump-probe tool to characterise the excited states of materials. It can be used to understand microscopic processes in photochemical and electrochemical transformations, including phenomena such as electron trapping and carrier recombination.

The drawback is that TAS spectra are difficult to interpret, especially for crystals where the specific valence and conduction band structure can give rise to complex features. Our goal here is to predict TAS features from first-principles starting from the most simple models of static excitations through to the kinetics of relaxation of the excited state back to the ground state.

Installation
============

To install the module with ``pip`` *(recommended)*:

.. code::

	pip install pytaser

To install directly from the git repository:

.. code::

	pip install git+https://github.com/WMD-group/PyTASER

To do a manual build and installation:

.. code::

	python3 setup.py build
	python3 setup.py install

Developer’s installation (optional)
===================================

For development work, ``PyTASER`` can also be installed from a copy of the source directory:

Download ``PyTASER`` source code using the command:

.. code::

    git clone https://github.com/WMD-group/PyTASER

Navigate to root directory:

.. code::

    cd PyTASER

Install the code with the command:

.. code::

    pip install -e .

This command tries to obtain the required packages and their dependencies and install them automatically.

Dependencies
============

``PyTASER`` is currently compatible with ``Python 3.9+`` and relies on a number of open-source python packages, specifically:

* `pymatgen <https://pymatgen.org/index.html>`__
* `numpy <https://numpy.org/>`__, `scipy <https://scipy.org/>`__ for data structures and unit conversion
* `matplotlib <https://matplotlib.org/>`__ for plotting the spectra

Visualisation
=============

The preferred method is to generate a ``Jupyter`` Notebook, as shown in the `examples <https://github.com/WMD-group/PyTASER/blob/main/examples/>`__ folder.
Alternatively, you can setup a ``Python`` file to run it in the command line of the terminal:

.. code::

    python3 <filename.py>

Contributing
============

We appreciate any contributions in the form of a pull request.
Additional test cases/example spectra performed with ``PyTASER`` would be welcomed.
Please feel free to reach out to us if there are any questions or suggestions.

Future topics we plan to build on:

* Incorporating finite-temperature effects (particularly for indirect bandgaps)
* Description of more complex optical processes (e.g. stimulated emission)
* Direct treatment of pump-probe time delay
* Incorporating spin-flip processes for spin-polarised systems
* Description of defective crystals

Acknowledgements
================

Developed by Savyasanchi Aggarwal and Alex Ganose. Aron Walsh helped to design the project. Thanks to group members for their support, especially Seán Kavanagh, Youngwon Woo, Anahita Manchala and Liam Harnett-Caulfield.


.. toctree::
   :maxdepth: 1
   :caption: Usage Guide
   :hidden:

   dft_examples
   jdos_examples


.. toctree::
   :caption: Information
   :hidden:

   changelog_link
   contributors
   contributing
   license
   PyTASER on GitHub <https://github.com/WMD-group/PyTASER>

.. toctree::
   :hidden:
   :maxdepth: 4
   :caption: Function Reference

   API reference <pytaser>
   genindex


