import generator
import plotter
## IMPORTS

#numpy
import numpy as np

#abipy
from abipy.tools.numtools import gaussian

#collections
from collections import OrderedDict

#pymatgen
from pymatgen.ext.matproj import MPRester
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.dos import f0, FermiDos
from pymatgen.core import Structure

#scipy
import scipy.constants as scpc

#converting from ev_to_lambda
def ev_to_lambda(ev):
    """Convert photon energies from eV to a wavelength in nm."""
    wavelength = ((scpc.h * scpc.c)/(ev*scpc.electron_volt))*10e8
    return wavelength


#converting from lambda to ev
def lambda_to_ev(lambda_float):
    """Convert photon energies from a wavelength in nm to eV."""
    electronvolts = (10e8*(scpc.h * scpc.c))/(lambda_float*scpc.electron_volt)
    return electronvolts



class tas:
    """
    A container class for all the generated data from TASgenerator.
    Inputs:
        total_TAS: overall TAS spectrum for a material under the specified conditions
        TAS_decomp: TAS spectrum across the energy mesh for a specific band transition i (initial) -> f (final) [dict]
        jdos_light_total: overall JDOS (pump-on) for a material under the specified conditions
        jdos_light_decomp: JDOS (pump-on) across the energy mesh for a specific band transition i (initial) -> f (final) [dict]
        jdos_dark_total: overall JDOS (pump-off) for a material under the specified conditions
        jdos_dark_decomp: JDOS (pump-off) across the energy mesh for a specific band transition i (initial) -> f (final) [dict]
        energy_mesh_ev: Energy mesh of spectra in eV, with an interval of 'step'.
    """
    def __init__(self, total_TAS, TAS_decomp, jdos_light_total, jdos_light_decomp, jdos_dark_total, jdos_dark_decomp, energy_mesh_ev):
        self.total_TAS = total_TAS
        # self.denom_total_TAS = denom_TAS_total
        self.TAS_decomp = TAS_decomp
        # self.denom_TAS_decomp = denom_TAS_decomp
        self.jdos_light_tot = jdos_light_total
        self.jdos_light_decomp = jdos_light_decomp
        self.jdos_dark_tot = jdos_dark_total
        self.jdos_dark_decomp = jdos_dark_decomp
        self.energy_mesh_ev = energy_mesh_ev
        self.wavelength_mesh = ev_to_lambda(self.energy_mesh_ev)
