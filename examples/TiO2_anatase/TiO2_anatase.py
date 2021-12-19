import generator
import TAS
import plotter
## IMPORTS

#numpy
import numpy as np

#abipy
import abipy.data as abidata
import abipy.abilab as abilab
import abipy.electrons as abielectrons
from abipy.abilab import abiopen
from abipy.tools.numtools import gaussian
abilab.enable_notebook()
from abipy.abilab import AbinitInput

#collections
from collections import OrderedDict

#pymatgen
from pymatgen.ext.matproj import MPRester
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.dos import f0, FermiDos
from pymatgen.core import Structure

#scipy
import scipy.constants as scpc

#######################

tio2a_temp = 298
tio2a_conc = 10e21
bandgap = 3.1
transitions_tio2a = [(8,11), (9,10), (4,14), (6,20), (6,8)]

# TiO2 is not spin-polarised, so we use the band_occupancies method instead of band_occupancies_sppol

tio2a_generator = generator.TASGenerator.from_mpid('cTHyvyRHqX34gOLMbM', 'mp-390', tio2a_temp, tio2a_conc, bandgap)
tio2a_tas = tio2a_generator.generate_tas(energy_min=0.0, energy_max=10.0)
tio2a_plotter = plotter.TASPlotter(tio2a_tas, bandgap, "TiO2 (anatase)", tio2a_temp, tio2a_conc)
plot_tio2a = tio2a_plotter.get_plot(transitions_tio2a, xaxis="wavelength", xmin = 350, xmax = 1400, yaxis = "TAS (deltaT)")
plot_tio2a.show()
