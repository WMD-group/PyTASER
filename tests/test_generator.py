import numpy as np
import pytest
from deepdiff import DeepDiff
from monty.serialization import loadfn
from pymatgen.electronic_structure.core import Spin
from unittest import mock

from pytaser import generator
from pytaser.generator import TASGenerator, set_bandgap

gaussian = loadfn("data_gaas/gaussian_123.json")


@pytest.mark.parametrize(
    "gaussian_example", loadfn("data_gaas/gaussian_123.json")
)
def test_gaussian(gaussian_example):
    x = np.array([1, 2, 3])
    width = 0.1
    assert generator.gaussian(x, width).all() == gaussian_example.all()


def test_set_bandgap(conditions):
    old_bs = loadfn("data_gaas/gaas_2534_bs.json")
    old_dos = loadfn("data_gaas/gaas_2534_dos.json")
    new_bs = loadfn("data_gaas/new_gaas_2534_bs.json")
    new_dos = loadfn("data_gaas/new_gaas_2534_dos.json")
    bs_new, dos_new = set_bandgap(old_bs, old_dos, conditions[2])
    assert bs_new.nb_bands == new_bs.nb_bands
    assert bs_new.get_band_gap() == new_bs.get_band_gap()
    assert dos_new.densities[Spin.up].all() == new_dos.densities[Spin.up].all()


def test_jdos(generated_class, light, dark, tas_object):
    # considering only spin-up
    i = (generated_class.vb[Spin.up]) - 1
    f = (generated_class.cb[Spin.up]) + 1
    energy_mesh = np.arange(0, 4, 0.15)
    new_jdos_light = generator.jdos(
        generated_class.bs,
        f,
        i,
        light[Spin.up],
        energy_mesh,
        kweights=generated_class.kpoint_weights,
        gaussian_width=0.1,
        spin=Spin.up,
    )

    new_jdos_dark = generator.jdos(
        generated_class.bs,
        f,
        i,
        dark[Spin.up],
        energy_mesh,
        kweights=generated_class.kpoint_weights,
        gaussian_width=0.1,
        spin=Spin.up,
    )

    i_tas = i - (generated_class.vb[Spin.up])
    f_tas = f - (generated_class.cb[Spin.up]) + 1
    assert (
            new_jdos_light.all()
            == tas_object.jdos_light_decomp[(i_tas, f_tas)].all()
    )
    assert (
            new_jdos_dark.all()
            == tas_object.jdos_dark_decomp[(i_tas, f_tas)].all()
    )


def test_get_cbm_vbm_index(generated_class):
    bs = generated_class.bs
    cb = generator.get_cbm_vbm_index(bs)[1][Spin.up]
    vb = generator.get_cbm_vbm_index(bs)[0][Spin.up]

    assert cb == bs.get_cbm()["band_index"][Spin.up][0]
    assert vb == bs.get_vbm()["band_index"][Spin.up][-1]
    assert int(cb) == int(vb) + 1


def test_band_occupancies(generated_class, light, conditions):
    l_occ_inp = light[Spin.up]
    cbm = generated_class.cb[Spin.up]
    vbm = generated_class.vb[Spin.up]
    dark_occs = generated_class.band_occupancies(
        conditions[0], conditions[1], dark=True
    )
    vb_dark_occs = dark_occs[Spin.up][: (vbm + 1)]
    cb_dark_occs = dark_occs[Spin.up][cbm:]

    for band in vb_dark_occs:
        occ_bool = all(kpoint_occ == 1 for kpoint_occ in band)
        assert occ_bool == True
    for band in cb_dark_occs:
        occ_bool = all(kpoint_occ == 0 for kpoint_occ in band)
        assert occ_bool == True

    light_occs = generated_class.band_occupancies(
        conditions[0], conditions[1], dark=False
    )
    assert light_occs[Spin.up].all() == l_occ_inp.all()


def test_generate_tas(generated_class, light, dark, tas_object, conditions):
    energy_min = 0
    energy_max = 4
    gaussian_width = 0.1
    step = 0.15

    tas_class = generated_class.generate_tas(
        conditions[0],
        conditions[1],
        energy_min,
        energy_max,
        gaussian_width,
        step,
        light,
        dark,
    )

    assert tas_class.total_tas.all() == tas_object.total_tas.all()
    assert DeepDiff(tas_class.tas_decomp, tas_object.tas_decomp) == {}
    assert (
            DeepDiff(tas_class.jdos_light_decomp, tas_object.jdos_light_decomp)
            == {}
    )
    assert tas_class.jdos_light_tot.all() == tas_object.jdos_light_tot.all()
    assert tas_class.jdos_dark_tot.all() == tas_object.jdos_dark_tot.all()
    assert tas_class.energy_mesh_ev.all() == tas_object.energy_mesh_ev.all()
    assert tas_class.bandgap == tas_object.bandgap

    assert (
            tas_class.total_tas.all()
            == (tas_class.jdos_light_tot - tas_class.jdos_dark_tot).all()
    )

    number_combinations_if = int(
        (generated_class.bs.nb_bands * (generated_class.bs.nb_bands - 1)) / 2
    )
    assert (
            number_combinations_if
            == len(tas_class.tas_decomp.keys())
            == len(tas_class.jdos_light_decomp.keys())
            == len(tas_class.jdos_dark_decomp.keys())
    )

    cbm = generated_class.cb[Spin.up]
    vbm = generated_class.vb[Spin.up]
    energy_mesh = np.arange(0, 4, 0.15)
    spin = Spin.up
    jdos_vbm_cbm = generator.jdos(
        generated_class.bs,
        cbm,
        vbm,
        light[Spin.up],
        energy_mesh,
        generated_class.kpoint_weights,
        gaussian_width,
        spin,
    )
    assert tas_class.jdos_light_decomp[0, 1].all() == jdos_vbm_cbm.all()


@mock.patch("pymatgen.ext.matproj.MPRester")
def test_from_mpid(mocker, generated_class, conditions):
    mock_mpr = mocker.Mock()
    mock_mpr.get_dos_by_material_id.return_value = loadfn("data_gaas/gaas_2534_dos.json")
    mock_mpr.get_bandstructure_by_material_id.return_value = loadfn("data_gaas/gaas_2534_bs.json")

    gaas2534 = TASGenerator.from_mpid("mp-2534", conditions[2], mpr=mock_mpr)
    mock_mpr.get_dos_by_material_id.assert_called_once_with('mp-2534')
    mock_mpr.get_bandstructure_by_material_id.assert_called_once_with('mp-2534', line_mode=False)

    assert (
            gaas2534.bs.projections[Spin.up].all()
            == generated_class.bs.projections[Spin.up].all()
    )
    assert (
            gaas2534.bs.is_spin_polarized == generated_class.bs.is_spin_polarized
    )
    assert gaas2534.bs.efermi == generated_class.bs.efermi
    assert type(gaas2534.bs) == type(generated_class.bs)

    assert gaas2534.dos.structure == generated_class.dos.structure
    assert type(gaas2534.dos) == type(generated_class.dos)

    assert (
            gaas2534.kpoint_weights.all() == generated_class.kpoint_weights.all()
    )
    assert gaas2534.bg_centre == generated_class.bg_centre
    assert gaas2534.vb == generated_class.vb
    assert gaas2534.cb == generated_class.cb
