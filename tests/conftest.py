import pytest
from monty.serialization import loadfn

from pytaser.generator import TASGenerator
from pytaser.plotter import TASPlotter
from pathlib import Path


@pytest.fixture(scope="package")
def datapath_gaas():
    """Path to data_gaas folder"""
    return Path(__file__).parent / "data_gaas"


@pytest.fixture(scope="package")
def datapath_cdte():
    """Path to data_cdte folder"""
    return Path(__file__).parent / "data_cdte"


@pytest.fixture(scope="package")
def examplepath_cdte():
    """Path to CdTe examples folder"""
    return Path(__file__).parent.parent / "examples/CdTe"


@pytest.fixture(scope="package")
def conditions():  # return [temperature, carrier conc, bandgap] in that order
    return [298, 1e18, 1.5]


@pytest.fixture(scope="package")
def cdte_conditions():  # return [temperature, carrier conc, bandgap] in that order
    return [300, 1e23, 0.67]


@pytest.fixture(scope="package")
def generated_class(datapath_gaas):
    gaas_bs = loadfn(datapath_gaas / "new_gaas_2534_bs.json")
    gaas_dos = loadfn(datapath_gaas / "new_gaas_2534_dos.json")
    kweights = loadfn(datapath_gaas / "gaas_kweights.json")
    return TASGenerator(gaas_bs, kweights, gaas_dos)


@pytest.fixture(scope="package")
def cdte_generated_class(datapath_cdte):
    cdte_bs = loadfn(datapath_cdte / "cdte_444_bs.json")
    cdte_dos = loadfn(datapath_cdte / "cdte_444_dos.json")
    cdte_kweights = loadfn(datapath_cdte / "cdte_444_kweights.json")
    return TASGenerator(cdte_bs, cdte_kweights, cdte_dos)


@pytest.fixture(scope="package")
def cdte_vasp_generated_class(examplepath_cdte):
    return TASGenerator.from_vasp_outputs(vasprun_file=examplepath_cdte / "k666_Optics/vasprun.xml",
                                          waveder_file=examplepath_cdte / "k666_Optics/WAVEDER")


@pytest.fixture(scope="package")
def cdte_vasp_generated_class_vr_only(examplepath_cdte):
    return TASGenerator.from_vasp_outputs(vasprun_file=examplepath_cdte / "k666_Optics/vasprun.xml")


@pytest.fixture(name="dark", scope="module")
def dark_occs_generated(generated_class, conditions):
    return generated_class.band_occupancies(conditions[0], conditions[1])


@pytest.fixture(name="cdte_dark", scope="module")
def cdte_dark_occs_generated(cdte_generated_class, cdte_conditions):
    return cdte_generated_class.band_occupancies(
        cdte_conditions[0], cdte_conditions[1]
    )


@pytest.fixture(name="light", scope="module")
def light_occs_generated(generated_class, conditions):
    return generated_class.band_occupancies(
        conditions[0], conditions[1], dark=False
    )


@pytest.fixture(name="cdte_light", scope="module")
def cdte_light_occs_generated(cdte_generated_class, cdte_conditions):
    return cdte_generated_class.band_occupancies(
        cdte_conditions[0], cdte_conditions[1], dark=False
    )


@pytest.fixture(scope="module")
def tas_object(generated_class, conditions, dark, light):
    return generated_class.generate_tas(
        conditions[0],
        conditions[1],
        energy_min=0,
        energy_max=4,
        step=0.15,
        light_occs=light,
        dark_occs=dark,
    )


@pytest.fixture(scope="module")
def cdte_tas_object(
    cdte_generated_class, cdte_conditions, cdte_dark, cdte_light
):
    return cdte_generated_class.generate_tas(
        cdte_conditions[0],
        cdte_conditions[1],
        energy_min=0,
        energy_max=10,
        step=0.01,  # default
        light_occs=cdte_light,
        dark_occs=cdte_dark,
    )


@pytest.fixture(scope="module")
def cdte_vasp_tas_object(
    cdte_vasp_generated_class, cdte_conditions
):
    return cdte_vasp_generated_class.generate_tas(
        cdte_conditions[0],
        cdte_conditions[1],
        energy_min=0,
        energy_max=10,
        step=0.01,  # default
    )


@pytest.fixture(scope="module")
def cdte_vasp_tas_object_vr_only(
    cdte_vasp_generated_class_vr_only, cdte_conditions
):
    return cdte_vasp_generated_class_vr_only.generate_tas(
        cdte_conditions[0],
        cdte_conditions[1],
        energy_min=0,
        energy_max=10,
        step=0.01,  # default
    )


@pytest.fixture(scope="module")
def plotter_gaas(tas_object, conditions):
    return TASPlotter(
        tas_object,
        material_name="GaAs",
    )


@pytest.fixture(scope="module")
def plotter_cdte(cdte_tas_object, cdte_conditions):
    return TASPlotter(
        cdte_tas_object,
        material_name="CdTe",
    )

@pytest.fixture(scope="module")
def plotter_cdte_vasp(cdte_vasp_tas_object, cdte_conditions):
    return TASPlotter(
        cdte_vasp_tas_object,
        material_name="CdTe",
    )

@pytest.fixture(scope="module")
def plotter_cdte_vasp_vr_only(cdte_vasp_tas_object_vr_only, cdte_conditions):
    return TASPlotter(
        cdte_vasp_tas_object_vr_only,
        material_name="CdTe",
    )
