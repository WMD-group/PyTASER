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

gaas_temp = 298
gaas_conc = 10e21
bandgap = 1.46
transitions_gaas = [(8,11), (9,10), (4,14), (6,20), (6,8)]


# GaAs is not spin-polarised, so we use the band_occupancies method instead of band_occupancies_sppol

GaAs_generator = generator.TASGenerator.from_mpid('cTHyvyRHqX34gOLMbM', 'mp-2534', gaas_temp, gaas_conc, bandgap)
GaAs_tas = GaAs_generator.generate_tas(energy_min=0.0, energy_max=10.0)
GaAs_plotter = plotter.TASPlotter(GaAs_tas, bandgap, "GaAs", gaas_temp, gaas_conc)
plot_gaas = GaAs_plotter.get_plot(transitions_gaas, xaxis="wavelength", xmin = 350, xmax = 1400, yaxis = "TAS (deltaT)")
plot_gaas.show()
