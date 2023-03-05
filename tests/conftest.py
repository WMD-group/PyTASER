import pytest
from pytaser.generator import TASGenerator
from pytaser.plotter import TASPlotter
from monty.serialization import loadfn


@pytest.fixture(scope='package')
def conditions():  # return [temperature, carrier conc, bandgap] in that order
    return [298, 1e18, 1.5]


@pytest.fixture(scope='package')
def generated_class():
    gaas_bs = loadfn('data_gaas/new_gaas_2534_bs.json')
    gaas_dos = loadfn('data_gaas/new_gaas_2534_dos.json')
    kweights = loadfn('data_gaas/gaas_kweights.json')
    return TASGenerator(gaas_bs, kweights, gaas_dos)


@pytest.fixture(name='dark', scope='module')
def dark_occs_generated(generated_class, conditions):
    return generated_class.band_occupancies(conditions[0], conditions[1])


@pytest.fixture(name='light', scope='module')
def light_occs_generated(generated_class, conditions):
    return generated_class.band_occupancies(conditions[0], conditions[1], dark=False)


@pytest.fixture(scope='module')
def tas_object(generated_class, conditions, dark, light):
    return generated_class.generate_tas(conditions[0],
                                        conditions[1],
                                        energy_min=0,
                                        energy_max=4,
                                        step=0.15,
                                        light_occs=light,
                                        dark_occs=dark)


@pytest.fixture(scope='module')
def plotter_gaas(tas_object, conditions):
    return TASPlotter(tas_object, bandgap_ev=conditions[2], material_name='GaAs', temp=conditions[0],
                      conc=conditions[1])

