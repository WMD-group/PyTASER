import pytest
import numpy as np
from pytaser.kpoints import get_kpoint_weights
from pytaser.tas import Tas
from pytaser import generator
from generator import TASGenerator
from monty.serialization import loadfn

gaussian = loadfn('data_gaas/gaussian_123.json')
gaas_dos_bs_old = [loadfn('data_gaas/gaas_2534_dos.json'), loadfn('data_gaas/gaas_2534_bs.json')]
gaas_dos_bs_new = [loadfn('data_gaas/new_gaas_2534_dos.json'), loadfn('data_gaas/new_gaas_2534_bs.json')]
gaas_kweights = [loadfn('gaas_kweights.json')]


@pytest.mark.parametrize("gaussian_example", gaussian)
def test_gaussian(gaussian_example):
    x = np.array([1, 2, 3])
    width = 0.1
    assert generator.gaussian(x, width) == gaussian_example


@pytest.fixture
def conditions():  # return [temperature, carrier conc, bandgap] in that order
    return [298, 1e18, 1.5]


@pytest.mark.parametrize("gaas_dos", gaas_dos_bs_new[0])
@pytest.mark.parametrize("gaas_bs", gaas_dos_bs_new[1])
@pytest.mark.parametrize("kweights", gaas_kweights)
@pytest.fixture
def generated_class(gaas_dos, gaas_bs, kweights):
    return TASGenerator(gaas_bs, kweights, gaas_dos)


@pytest.fixture(name='dark')
def dark_occs_generated(generated_class, conditions):
    return generated_class.band_occupancies(conditions[0], conditions[1])


@pytest.fixture(name='light')
def light_occs_generated(generated_class, conditions):
    return generated_class.band_occupancies(conditions[0], conditions[1], dark=False)


@pytest.fixture(scope='package')
def tas_object(generated_class, conditions, dark, light):
    return generated_class.generate_tas(conditions[0],
                                        conditions[1],
                                        energy_min=0,
                                        energy_max=4,
                                        step=0.15,
                                        light_occs=light,
                                        dark_occs=dark)


@pytest.mark.parametrize("old", gaas_dos_bs_old)
@pytest.mark.parametrize("new", gaas_dos_bs_new)
def test_set_bandgap(old, new, conditions):
    bs_new, dos_new = generator.set_bandgap(gaas_dos_bs_old[0], gaas_dos_bs_old[1], conditions[2])
    assert bs_new.nb_bands == new[0].nb_bands
    assert bs_new.get_band_gap() == new[0].get_band_gap()
    assert dos_new.densities == new[1].densities


def test_jdos(generated_class, light, dark, tas_object):
    # considering only spin-up
    i = (generated_class.vb[Spin.up]) - 1
    f = (generated_class.cb[Spin.up]) + 1
    energy_mesh = np.arange(0, 4, 0.15)
    new_jdos_light = generator.jdos(generated_class.bs,
                                    f,
                                    i,
                                    light[Spin.up],
                                    energy_mesh,
                                    kweights=generated_class.kpoint_weights,
                                    gaussian_width=0.1,
                                    spin=Spin.up)

    new_jdos_dark = generator.jdos(generated_class.bs,
                                   f,
                                   i,
                                   dark[Spin.up],
                                   energy_mesh,
                                   kweights=generated_class.kpoint_weights,
                                   gaussian_width=0.1,
                                   spin=Spin.up)

    i_tas = i - (generated_class.vb[Spin.up])
    f_tas = f - (generated_class.cb[Spin.up]) + 1
    assert new_jdos_light == tas_object.jdos_light_decomp[(i_tas, f_tas)]
    assert new_jdos_dark == tas_object.jdos_dark_decomp[(i_tas, f_tas)]


def test_get_cbm_vbm_index(generated_class):
    bs = generated_class.bs
    cb = generator.get_cbm_vbm_index(bs)[1][Spin.up]
    vb = generator.get_cbm_vbm_index(bs)[0][Spin.up]

    assert cb == bs.get_cbm()['band_index'][Spin.up][0]
    assert vb == bs.get_vbm()['band_index'][Spin.up][0]
    assert int(cb) == int(vb) + 1


def test_band_occupancies(generated_class, light, conditions):
    l_occ_inp = light[Spin.up]
    cbm = generated_class.cb[Spin.up]
    vbm = generated_class.vb[Spin.up]
    dark_occs = generated_class.band_occupancies(conditions[0], conditions[1], dark=True)
    vb_dark_occs = dark_occs[Spin.up][:(vbm + 1)]
    cb_dark_occs = dark_occs[Spin.up][cbm:]

    for band in vb_dark_occs:
        occ_bool = (all(kpoint_occ == 1 for kpoint_occ in band))
        assert occ_bool == True
    for band in cb_dark_occs:
        occ_bool = (all(kpoint_occ == 0 for kpoint_occ in band))
        assert occ_bool == True

    light_occs = generated_class.band_occupancies(conditions[0], conditions[1], dark=False)
    assert light_occs == l_occ_inp


def test_generate_tas(generated_class, light, dark, tas_object, conditions):
    energy_min = 0
    energy_max = 4
    gaussian_width = 0.1
    step = 0.15

    tas_class = generated_class.generate_tas(conditions[0],
                                             conditions[1],
                                             energy_min,
                                             energy_max,
                                             gaussian_width,
                                             step,
                                             light,
                                             dark)

    assert tas_class.total_tas == tas_object.total_tas
    assert tas_class.tas_decomp == tas_object.tas_decomp
    assert tas_class.jdos_light_tot == tas_object.jdos_light_tot
    assert tas_class.jdos_dark_tot == tas_object.jdos_dark_tot
    assert tas_class.energy_mesh_ev == tas_object.energy_mesh_ev
    assert tas_class.bandgap_ev == tas_object.bandgap_ev

    assert tas_class.total_tas == (tas_class.jdos_light_tot - tas_class.jdos_dark_tot)

    number_combinations_if = int((tas_class.bs.nb_bands * (tas_class.bs.nb_bands - 1)) / 2)
    assert number_combinations_if == len(tas_class.tas_decomp.keys()) == len(
        tas_class.jdos_light_decomp.keys()) == len(tas_class.jdos_dark_decomp.keys())

    cbm = generated_class.cb[Spin.up]
    vbm = generated_class.vb[Spin.up]
    energy_mesh = np.arange(0, 4, 0.15)
    spin = Spin.up
    jdos_vbm_cbm = generator.jdos(object.bs,
                                  cbm,
                                  vbm,
                                  light,
                                  energy_mesh,
                                  generated_class.kpoint_weights,
                                  gaussian_width,
                                  spin)
    assert tas_class.jdos_light_decomp[0, 1] == jdos_vbm_cbm


def test_from_mpid(generated_class, conditions):
    gaas2534 = TASGenerator.from_mpid('mp-2534', conditions[2])
    assert gaas2534.bs == generated_class.bs
    assert gaas2534.kpoint_weights == generated_class.kpoint_weights
    assert gaas2534.dos == generated_class.dos
    assert gaas2534.bg_centre == generated_class.bg_centre
    assert gaas2534.vb == generated_class.vb
    assert gaas2534.cb == generated_class.cb
