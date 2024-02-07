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

``PyTASER`` is a ``Python`` library for simulating differential absorption spectra in compounds using
quantum chemical calculations, including transient (TAS) and differential (DAS) absorption spectroscopies.
The goal of this library is to simulate spectra for comparison with and interpretation of experiments. The main features include:

* Transient (TAS) and differential (DAS) absorption spectra for insulating and semiconducting crystals
  (using pre-calculated DFT, or database calculations)
* Spectra for different temperature and carrier density conditions
* Separation of spectral contributions from individual band-to-band transitions
* Integration with the Materials Project database, allowing support for non-locally calculated materials
* Ability to produce publication-ready figures, with flexibility in plotting.

PyTASER is designed for users with moderate experience in computational methods and optical techniques, enabled by the following features:

* Use of Python as the programming language (due to its low entry barrier, flexibility and popularity in the materials modelling field)
* Documentation and easy-to-follow workflows with complete unit-test coverage.
* Interfaced with the popular materials analysis package `pymatgen <https://pymatgen.org/index.html>`__.
* Currently compatible with VASP, while support for other electronic structure codes is planned – if you have a specific request, please open an issue!

==========
Background
==========

TAS is a powerful pump-probe tool to characterise the excited states of materials, while DAS can be used to represent how changes in a system affects its ground-state optical absorption spectra. These techniques can be used to understand microscopic processes in photochemical and electrochemical transformations, including phenomena such as electron trapping and carrier recombination.

The drawback is that TAS spectra are difficult to interpret, especially for crystals where the specific valence and conduction band structure can give rise to complex features. Our goal here is to predict TAS features from first-principles starting from the most simple models of static excitations through to the kinetics of relaxation of the excited state back to the ground state.

To achieve this, PyTASER identifies the allowed vertical optical transitions between electronic bands of the material to determine possible excitations that can occur from the ground 'dark' and excited 'light' electronic states.
This is done by calculating the effective absorption in each state - this is a product of the joint density of states (JDOS) and the transition probability for each band transition, both of which are based on post-processing ground state DFT calculations. Once calculated, PyTASER compares changes in electronic transitions between the dark and light states, as demonstrated in the figure below.

.. figure:: _static/jdos-tas-schematic.png
   :scale: 15 %
   :align: center

   Schematics of the ground and excited state electronic structures and optical profiles. The ground 'dark' state is at the top, showing full occupancy and unoccupancy (blue, orange) for the conduction and valence bands respectively. The excited 'light' state shows partial occupancy in a similar plot at the bottom. The overall DA plot is displayed to the right, the difference between the dark and light effective absorption plots.

============
Installation
============

To install the module with ``pip`` *(recommended)*:

.. code::

	pip install pytaser

To install directly from the git repository:

.. code::

	pip install git+https://github.com/WMD-group/PyTASER

Alternatively ``python setup.py install`` can also be used.

===================================
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

============
Dependencies
============

``PyTASER`` is currently compatible with ``Python 3.9+`` and relies on a number of open-source python packages, specifically:

* `pymatgen <https://pymatgen.org/index.html>`__
* `numpy <https://numpy.org/>`__, `scipy <https://scipy.org/>`__ for data structures and unit conversion
* `matplotlib <https://matplotlib.org/>`__ for plotting the spectra

=============
Visualisation
=============

The recommended approach is to use ``PyTASER`` within an interactive python environment (such as a Jupyter
Notebook), as shown in the `tutorials <https://pytaser.readthedocs.io/en/latest/dft_examples.html>`__ on
the docs and the `examples <https://github.com/WMD-group/PyTASER/blob/main/examples/>`__ folder.
Alternatively, you can use ``PyTASER`` through python scripts (e.g. with ``python <filename.py>``).

If using your own, locally calculated data, please follow `this workflow <https://github.com/WMD-group/PyTASER/blob/main/examples/PyTASER_DFT_Example.ipynb>`__. Otherwise, please use the MP-integrated workflow as shown `here <https://github.com/WMD-group/PyTASER/blob/main/examples/PyTASER_MP_Example.ipynb>`__.

============
Contributing
============

We appreciate any contributions in the form of a pull request. Please see the `Contributing <https://pytaser.readthedocs.io/en/latest/contributing.html>`__ Documentation for more details.
Additional analysis/example spectra performed with ``PyTASER`` would be welcomed.

Please feel free to reach out to us via the `Issue Tracker <https://github.com/WMD-group/PyTASER/issues>`__ if there are any questions or suggestions.

=======
Testing
=======

Unit tests are in the ``tests`` directory and can be run from the top directory using pytest. Automatic testing is run on the master and develop branches using Github Actions.

Please run tests and add new tests for any new features whenever submitting pull requests.

===========
Future Work
===========

Future topics we plan to build on:

* Incorporating finite-temperature effects (particularly for indirect bandgaps)
* Direct treatment of pump-probe time delay and relaxation kinetics
* Incorporating spin-flip processes for spin-polarised systems
* Description of defective crystals

================
Acknowledgements
================

Developed by Savyasanchi Aggarwal, Seán Kavanagh, Lucas Verga, Youngwon Woo and Alex Ganose. Aron Walsh
designed and led the project. Thanks to group members for their support, especially Anahita Manchala and
Liam Harnett-Caulfield.


.. toctree::
   :maxdepth: 1
   :caption: Usage Guide
   :hidden:

   dft_examples
   mp_example


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


