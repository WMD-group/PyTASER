import pytest
import numpy as np
from pytaser.kpoints import get_kpoint_weights
from pytaser.tas import Tas
from pytaser import generator
from generator import TASGenerator
from monty.serialization import loadfn

sil_orig = [loadfn("data_sil/sil_bs.json"), loadfn("data_sil/sil_dos.json")]
sil_new = [loadfn("data_sil/sil_new_bs.json"), loadfn("data_sil/sil_new_dos.json")]
sil_obj = loadfn('data_sil/sil_object.yaml')
sil_light_occs = loadfn('data_sil/sil_light_occs.yaml')
sil_dark_occs = loadfn('data_sil/sil_dark_occs.yaml')
gaussian = loadfn('data_sil/gaussian_123.json')
jdos_sil_short = loadfn('data_sil/example_jdos_sil.json')
basic_sil_tas_class = loadfn('data_sil/example_tas.yaml')


@pytest.mark.parametrize("gaussian_example", gaussian)
def test_gaussian(gaussian_example):
    x = np.array([1, 2, 3])
    width = 0.1
    assert generator.gaussian(x, width) == gaussian_example


@pytest.mark.parametrize("old", sil_orig)
@pytest.mark.parametrize("new", sil_new)
def test_set_bandgap(old, new):
    bg = 1.11
    bs_new, dos_new = generator.set_bandgap(sil_orig[0], sil_orig[1], bg)
    assert bs_new.nb_bands == sil_new[0].nb_bands
    assert bs_new.get_band_gap() == sil_new[0].get_band_gap()
    assert dos_new.densities == sil_new[1].densities


@pytest.mark.parametrize("object", sil_obj)
@pytest.mark.parametrize("occs", sil_light_occs[Spin.up])
@pytest.mark.parametrize("jdos_inp", jdos_sil_short)
def test_jdos(object, occs, jdos_inp):
    init_band = (object.vb[Spin.up] - 1)
    final_band = (object.cb[Spin.up] + 1)
    assert final_band > init_band

    bs = object.bs
    energy_mesh = np.arange(0, 2, 0.25)
    kweights_sil = object.kpoint_weights
    gaussian_width = 0.1
    spin = Spin.up

    example_jdos = generator.jdos(bs, final_band, init_band, occs, energy_mesh, kweights_sil, gaussian_width, spin)
    assert example_jdos == jdos_inp


@pytest.mark.parametrize("bs", sil_new)
def test_get_cbm_vbm_index(bs):
    cb = generator.get_cbm_vbm_index(bs)[1][Spin.up]
    vb = generator.get_cbm_vbm_index(bs)[0][Spin.up]

    assert cb == bs.get_cbm()['band_index'][Spin.up][0]
    assert vb == bs.get_vbm()['band_index'][Spin.up][0]
    assert int(cb) == int(vb) + 1


@pytest.mark.parametrize("object", sil_obj)
@pytest.mark.parametrize("l_occ_inp", sil_light_occs[Spin.up])
def test_band_occupancies(object, l_occ_inp):
    temp, conc = 298, 10e21
    cbm = object.cb[Spin.up]
    vbm = object.vb[Spin.up]
    dark_occs = object.band_occupancies(temp, conc, dark=True)
    vb_dark_occs = dark_occs[Spin.up][:(vbm + 1)]
    cb_dark_occs = dark_occs[Spin.up][cbm:]

    for band in vb_dark_occs:
        occ_bool = (all(kpoint_occ == 1 for kpoint_occ in band))
        assert occ_bool == True
    for band in cb_dark_occs:
        occ_bool = (all(kpoint_occ == 0 for kpoint_occ in band))
        assert occ_bool == True

    light_occs = object.band_occupancies(temp, conc, dark=False)
    assert light_occs == l_occ_inp


@pytest.mark.parametrize("object", sil_obj)
@pytest.mark.parametrize("light_occs", sil_light_occs)
@pytest.mark.parametrize("dark_occs", sil_dark_occs)
@pytest.mark.parametrize("example_tas", basic_sil_tas_class)
def test_generate_tas(object, light_occs, dark_occs, example_tas):
    # must use something with fewer bands than silicon. Si has 96, must be using something with ~30 bands. if
    # implementing an input for the occupancies, they must be a dict form rather than array. it has to match the
    # output from the band_occupancies method.
    temp, conc = 298, 10e21
    energy_min = 0
    energy_max = 2
    gaussian_width = 0.1
    step = 0.25

    tas_class = object.generate_tas(temp, conc, energy_min, energy_max, gaussian_width, step, light_occs, dark_occs)

    assert tas_class.total_tas == example_tas.total_tas
    assert tas_class.tas_decomp == example_tas.tas_decomp
    assert tas_class.jdos_light_tot == example_tas.jdos_light_tot
    assert tas_class.jdos_dark_tot == example_tas.jdos_dark_tot
    assert tas_class.energy_mesh_ev == example_tas.energy_mesh_ev
    assert tas_class.bandgap_ev == example_tas.bandgap_ev

    assert tas_class.total_tas == (tas_class.jdos_light_tot - tas_class.jdos_dark_tot)

    number_combinations_if = int((object.bs.nb_bands * (object.bs.nb_bands - 1)) / 2)
    assert number_combinations_if == len(basic_tas_class.tas_decomp.keys()) == len(
        basic_tas_class.jdos_light_decomp.keys()) == len(basic_tas_class.jdos_dark_decomp.keys())

    cbm = object.cb[Spin.up]
    vbm = object.vb[Spin.up]
    energy_mesh = np.arange(0, 2, 0.25)
    spin = Spin.up
    jdos_vbm_cbm = generator.jdos(object.bs, cbm, vbm, light_occs, energy_mesh, object.kpoint_weights, gaussian_width,
                                  spin)
    assert tas_class.jdos_light_decomp[0, 1] == jdos_vbm_cbm


@pytest.mark.parametrize("object", sil_obj)
def test_from_mpid(object):
    bg = 1.11
    sil165 = TASGenerator.from_mpid('mp-165', bg)
    assert sil165.bs == object.bs
    assert sil165.kpoint_weights == object.kpoint_weights
    assert sil165.dos == object.dos
    assert sil165.bg_centre == object.bg_centre
    assert sil165.vb == object.vb
    assert sil165.cb == object.cb
