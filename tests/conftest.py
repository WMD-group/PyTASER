"""
Setup for pytest.
"""

from pathlib import Path

import pytest
from monty.serialization import loadfn

from pytaser.das_generator import DASGenerator
from pytaser.generator import TASGenerator
from pytaser.plotter import TASPlotter


@pytest.fixture(scope="package")
def datapath_gaas():
    """Path to data_gaas folder."""
    return Path(__file__).parent / "data_gaas"


@pytest.fixture(scope="package")
def datapath_cdte():
    """Path to data_cdte folder."""
    return Path(__file__).parent / "data_cdte"


@pytest.fixture(scope="package")
def examplepath_cdte():
    """Path to CdTe examples folder."""
    return Path(__file__).parent.parent / "examples/CdTe"


@pytest.fixture(scope="package")
def examplepath_gaas():
    """Path to GaAs examples folder."""
    return Path(__file__).parent.parent / "examples/GaAs"


@pytest.fixture(scope="package")
def examplepath_tio2_das():
    """Path to TiO2_DAS examples folder."""
    return Path(__file__).parent.parent / "examples/TiO2-DAS"


@pytest.fixture(scope="package")
def conditions():  # return [temperature, carrier conc, bandgap] in that order
    return [298, 1e18, 1.5]


@pytest.fixture(scope="package")
def cdte_conditions():  # return [temperature, carrier conc, bandgap] in that order
    return [300, 1e23, 0.67]


@pytest.fixture(scope="package")
def tio2_das_conditions():  # return [energy_min, energy_max, temperature, c shift] in that order
    return [0, 10, 300, 1e-3]


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
    return TASGenerator.from_vasp_outputs(
        vasprun_file=examplepath_cdte / "k666_Optics/vasprun.xml",
        waveder_file=examplepath_cdte / "k666_Optics/WAVEDER",
    )


@pytest.fixture(scope="package")
def cdte_vasp_bg_3_generated_class(examplepath_cdte):
    return TASGenerator.from_vasp_outputs(
        vasprun_file=examplepath_cdte / "k666_Optics/vasprun.xml",
        waveder_file=examplepath_cdte / "k666_Optics/WAVEDER",
        bg=3,
    )


@pytest.fixture(scope="package")
def cdte_vasp_generated_class_vr_only(examplepath_cdte):
    return TASGenerator.from_vasp_outputs(vasprun_file=examplepath_cdte / "k666_Optics/vasprun.xml")


@pytest.fixture(scope="package")
def tio2_das_vasp_generated_class_vr_only(examplepath_tio2_das):
    vasprun_new_sys = examplepath_tio2_das / "mp554278-monoclinic/vasprun.xml"
    vasprun_ref = examplepath_tio2_das / "mp2657-tetragonal/vasprun.xml"
    return DASGenerator.from_vasp_outputs(vasprun_new_sys, vasprun_ref)


@pytest.fixture(scope="package")
def tio2_das_vasp_generated_class(examplepath_tio2_das):
    vasprun_new_sys = examplepath_tio2_das / "mp554278-monoclinic/vasprun.xml"
    vasprun_ref = examplepath_tio2_das / "mp2657-tetragonal/vasprun.xml"
    waveder_new_sys = examplepath_tio2_das / "mp554278-monoclinic/WAVEDER"
    waveder_ref = examplepath_tio2_das / "mp2657-tetragonal/WAVEDER"

    return DASGenerator.from_vasp_outputs(vasprun_new_sys, vasprun_ref, waveder_new_sys, waveder_ref)


@pytest.fixture(scope="package")
def das_class_vr_only(tio2_das_vasp_generated_class_vr_only, tio2_das_conditions):
    return tio2_das_vasp_generated_class_vr_only.generate_das(
        energy_min=tio2_das_conditions[0],
        energy_max=tio2_das_conditions[1],
        temp=tio2_das_conditions[2],
        cshift=tio2_das_conditions[3],
        processes=1,
    )


@pytest.fixture(scope="package")
def das_class_with_waveder(tio2_das_vasp_generated_class, tio2_das_conditions):
    return tio2_das_vasp_generated_class.generate_das(
        energy_min=tio2_das_conditions[0],
        energy_max=tio2_das_conditions[1],
        temp=tio2_das_conditions[2],
        cshift=tio2_das_conditions[3],
        processes=1,
    )


@pytest.fixture(name="dark", scope="module")
def dark_occs_generated(generated_class, conditions):
    return generated_class.band_occupancies(conditions[0], conditions[1])


@pytest.fixture(name="cdte_dark", scope="module")
def cdte_dark_occs_generated(cdte_generated_class, cdte_conditions):
    return cdte_generated_class.band_occupancies(cdte_conditions[0], cdte_conditions[1])


@pytest.fixture(name="light", scope="module")
def light_occs_generated(generated_class, conditions):
    return generated_class.band_occupancies(conditions[0], conditions[1], dark=False)


@pytest.fixture(name="cdte_light", scope="module")
def cdte_light_occs_generated(cdte_generated_class, cdte_conditions):
    return cdte_generated_class.band_occupancies(cdte_conditions[0], cdte_conditions[1], dark=False)


@pytest.fixture(scope="package")
def tas_object(datapath_gaas):
    return loadfn(datapath_gaas / "tas_gaas.json")


@pytest.fixture(scope="module")
def cdte_tas_object(cdte_generated_class, cdte_conditions, cdte_dark, cdte_light):
    return cdte_generated_class.generate_tas(
        cdte_conditions[0],
        cdte_conditions[1],
        energy_min=0,
        energy_max=10,
        step=0.01,  # default
        light_occs=cdte_light,
        dark_occs=cdte_dark,
    )


@pytest.fixture(scope="package")
def das_object(examplepath_tio2_das):
    return loadfn(examplepath_tio2_das / "TiO2_das.json.gz")


@pytest.fixture(scope="module")
def cdte_vasp_tas_object(cdte_vasp_generated_class, cdte_conditions):
    return cdte_vasp_generated_class.generate_tas(
        cdte_conditions[0],
        cdte_conditions[1],
        energy_min=0,
        energy_max=10,
        step=0.01,  # default
    )


@pytest.fixture(scope="module")
def cdte_vasp_tas_object_vr_only(cdte_vasp_generated_class_vr_only, cdte_conditions):
    return cdte_vasp_generated_class_vr_only.generate_tas(
        cdte_conditions[0],
        cdte_conditions[1],
        energy_min=0,
        energy_max=10,
        step=0.01,  # default
    )


@pytest.fixture(scope="module")
def cdte_vasp_tas_bg_3_object(cdte_vasp_bg_3_generated_class, cdte_conditions):
    return cdte_vasp_bg_3_generated_class.generate_tas(
        cdte_conditions[0],
        cdte_conditions[1],
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
def plotter_tio2_das(das_object, tio2_das_conditions):
    return TASPlotter(
        das_object,
        material_name="TiO2",
        system_name="Monoclinic",
        reference_name="Tetragonal",
    )


@pytest.fixture(scope="module")
def plotter_cdte_bg_3_vasp(cdte_vasp_tas_bg_3_object, cdte_conditions):
    return TASPlotter(
        cdte_vasp_tas_bg_3_object,
        material_name="CdTe",
    )


@pytest.fixture(scope="module")
def plotter_cdte_vasp_vr_only(cdte_vasp_tas_object_vr_only, cdte_conditions):
    return TASPlotter(
        cdte_vasp_tas_object_vr_only,
        material_name="CdTe",
    )
