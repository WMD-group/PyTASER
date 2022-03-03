import numpy as np
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.dos import FermiDos, f0
from pymatgen.ext.matproj import MPRester

from pytaser.kpoints import get_kpoint_weights
from pytaser.tas import Tas


def gaussian(x, width, center=0.0, height=None):
    """
    Returns the values of gaussian(x) where x is array-like.

    Args:
        x: Input array.
        width: Width of the gaussian.
        center: Center of the gaussian.
        height: height of the gaussian. If height is None, a normalized gaussian is returned.
    """
    x = np.asarray(x)
    if height is None:
        height = 1.0 / (width * np.sqrt(2 * np.pi))

    return height * np.exp(-(((x - center) / width) ** 2) / 2.0)


def set_bandgap(bandstructure, dos, bandgap):
    """
    Shifts all bands of a material to correct the DFT-underestimated bandgap according to
    the input experimental bandgap.

    Args:
        bandstructure: PMG bandstructure object
        dos: PMG dos object
        bandgap: experimental bandgap from literature (float)

    Returns:
        new bandstructure and dos objects for the same material, but with corrected
        bandgap.
    """
    from copy import deepcopy

    if abs(dos.efermi - bandstructure.efermi) > 0.001:
        raise ValueError("DOS and band structure are not from the same calculation")

    if bandstructure.is_metal():
        raise ValueError("System is a metal, cannot change bandgap")

    scissor = bandgap - bandstructure.get_band_gap()["energy"]
    midgap = (bandstructure.get_cbm()["energy"] + bandstructure.get_vbm()["energy"]) / 2

    new_bandstructure = deepcopy(bandstructure)
    for spin, spin_energies in bandstructure.bands.items():
        vb_mask = spin_energies <= midgap
        new_bandstructure.bands[spin][vb_mask] -= scissor / 2
        new_bandstructure.bands[spin][~vb_mask] += scissor / 2

    fermi_idx = np.argmin(np.abs(dos.energies - midgap))
    de = np.diff(dos.energies).mean()
    shift = int(scissor / 2 // de)
    new_dos = deepcopy(dos)
    for spin in dos.densities.keys():
        dens = np.zeros_like(dos.energies)
        if shift > 0:
            dens[: fermi_idx - shift] = dos.densities[spin][shift:fermi_idx]
            dens[fermi_idx + shift :] = dos.densities[spin][fermi_idx:-shift]
        else:
            dens[abs(shift) : fermi_idx] = dos.densities[spin][: fermi_idx + shift]
            dens[fermi_idx:+shift] = dos.densities[spin][fermi_idx - shift :]
        new_dos.densities[spin] = dens

    new_bandstructure.efermi = midgap
    new_dos.efermi = midgap

    return new_bandstructure, new_dos


def jdos(bs, f, i, occs, energies, kweights, gaussian_width, spin=Spin.up):
    """
    Args:
        bs: bandstructure object
        f: final band
        i: initial band
        occs: occupancies over all bands.
        energies: energy mesh (eV)
        kweights: k-point weights
        gaussian_width: width of gaussian plot.
        spin: Which spin channel to include.

    Returns:
        Cumulative JDOS value for a specific i->f transition, with consideration of
        partial occupancy and spin polarisation.
    """

    jdos = np.zeros(len(energies))
    for k in range(len(bs.bands[spin][i])):
        final_occ = occs[f][k]
        init_energy = bs.bands[spin][i][k]
        final_energy = bs.bands[spin][f][k]
        init_occ = occs[i][k]
        k_weight = kweights[k]
        factor = k_weight * (
            (init_occ * (1 - final_occ)) - (final_occ * (1 - init_occ))
        )
        jdos += factor * gaussian(
            energies, gaussian_width, center=final_energy - init_energy
        )
    return jdos


def get_cbm_vbm_index(bs):
    """
    Args:
        bs: bandstructure object

    Returns:
        Valence and Conduction band as an index number for different spins.
    """
    vbm_index = {}
    cbm_index = {}
    for spin, spin_bands in bs.bands.items():
        vbm_index[spin] = np.where(np.all(spin_bands <= bs.efermi, axis=1))[0].max()
        cbm_index[spin] = vbm_index[spin] + 1
    return vbm_index, cbm_index


class TASGenerator:
    """
    Class to generate a TAS spectrum (decomposed and cumulative) from a bandstructure and
    dos object.

    Args:
        bs: Pymatgen-based bandstructure object
        kpoint_weights: kpoint weights either found by the function or inputted.
        dos: Pymatgen-based dos object
        bandgap: Experimental or Inputted  bandgap of a material.
        temp: Temperature (K) of material we wish to investigate (affects the FD distribution)
        conc: Carrier concentration (cm^-3) of holes and electrons (both are equivalent).
            Inversely proportional to pump-probe time delay.

    Attributes:
        bs: Pymatgen-based bandstructure object
        kpoint_weights: kpoint weights either found by the function or inputted
        dos: Pymatgen-based dos object
        bg: Experimental or Inputted bandgap of a material.
        bg_centre: Energy (eV) of the bandgap centre.
        occs_light: Occupancies of the bands for a pump-on system, across all kpoints
        occs_dark: Occupancies of the bands for a pump-off system, across all kpoints
    """

    def __init__(self, bs, kpoint_weights, dos, bandgap, temp, conc):
        self.bs = bs
        self.kpoint_weights = kpoint_weights
        self.dos = dos
        self.bg_centre = (bs.get_cbm()["energy"] + bs.get_vbm()["energy"]) / 2

        if float(bandgap) < 0.01:
            raise ValueError("Bandgap is smaller than 0.01 eV; cannot compute TAS")

        self.occs_light = TASGenerator.band_occupancies(self, temp, conc, dark=False)
        self.occs_dark = TASGenerator.band_occupancies(self, temp, conc)
        self.vb = get_cbm_vbm_index(self.bs)[0]
        self.cb = get_cbm_vbm_index(self.bs)[1]

    def band_occupancies(self, temp, conc, dark=True):
        """
        Gives band occupancies.

        Args:
            temp: Temperature of material we wish to investigate (affects the FD
                distribution)
            conc: Carrier concentration of holes and electrons (both are the same).
                Inversely proportional to pump-probe time delay.
            dark: Bool; dark = True indicates pump is on.

        Returns:
            A dictionary of {Spin: occs} for all bands across all k-points.
        """
        # Calculate the quasi-Fermi levels
        fermidos = FermiDos(self.dos)
        q_fermi_e = fermidos.get_fermi(-conc, temp)  # quasi-electron fermi level
        q_fermi_h = fermidos.get_fermi(conc, temp)  # quasi-hole fermi level

        occs = {}
        if dark:
            for spin, spin_bands in self.bs.bands.items():
                hole_mask = spin_bands < self.bg_centre
                elec_mask = spin_bands > self.bg_centre

                # fully occupied hole mask, completely empty electron mask
                spin_occs = np.zeros_like(spin_bands)
                spin_occs[hole_mask] = 1
                spin_occs[elec_mask] = 0

                occs[spin] = spin_occs
        else:
            for spin, spin_bands in self.bs.bands.items():
                # Calculate the occupancies at the initial and final energy according to
                # the Fermi-Dirac distribution
                electron_occs = f0(spin_bands, q_fermi_e, temp)  # q-e-f
                hole_occs = f0(spin_bands, q_fermi_h, temp)  # q-h-f

                hole_mask = spin_bands < self.bg_centre
                elec_mask = spin_bands > self.bg_centre

                # set occupancies equal to computed by f0
                spin_occs = np.zeros_like(hole_occs)
                spin_occs[hole_mask] = hole_occs[hole_mask]
                spin_occs[elec_mask] = electron_occs[elec_mask]

                occs[spin] = spin_occs

        return occs

    def generate_tas(self, energy_min, energy_max, gaussian_width=0.2, step=0.1):

        """
        Generates TAS spectra based on inputted occupancies, and a specified energy mesh.

        Args:
            energy_min: Minimum band transition energy to consider for energy mesh (eV)
            energy_max: Maximum band transition energy to consider for energy mesh (eV)
            gaussian_width: Width of gaussian curve
            step: Interval between energy points in the energy mesh.

        Returns:
            TAS class containing the following inputs:
                tas_cumulative: overall deltaT TAS spectrum for a material under the
                    specified conditions
                tas_if: deltaT TAS spectrum across the energy mesh for a specific band
                    transition i (initial) -> f (final) [dict]
                jdos_light_cumulative: overall JDOS (pump-on) for a material under the
                    specified conditions
                jdos_light_if: JDOS (pump-on) across the energy mesh for a specific band
                    transition i (initial) -> f (final) [dict]
                jdos_dark_cumulative: overall JDOS (pump-off) for a material under the
                    specified conditions
                jdos_dark_if: JDOS (pump-off) across the energy mesh for a specific band
                    transition i (initial) -> f (final) [dict]
                energy_mesh_ev: Energy mesh of spectra in eV, with an interval of 'step'.
        """

        energy_mesh_ev = np.arange(energy_min, energy_max, step)
        jdos_light_if = {}
        jdos_dark_if = {}
        tas_if = {}
        tas_cumulative = np.zeros(len(energy_mesh_ev))
        jdos_dark_cumulative = np.zeros(len(energy_mesh_ev))
        jdos_light_cumulative = np.zeros(len(energy_mesh_ev))
        for spin, spin_bands in self.bs.bands.items():
            for i in range(len(spin_bands)):
                for f in range(len(spin_bands)):
                    if f > i:
                        jd_light = jdos(
                            self.bs,
                            f,
                            i,
                            self.occs_light[spin],
                            energy_mesh_ev,
                            self.kpoint_weights,
                            gaussian_width,
                            spin=spin,
                        )
                        jd_dark = jdos(
                            self.bs,
                            f,
                            i,
                            self.occs_dark[spin],
                            energy_mesh_ev,
                            self.kpoint_weights,
                            gaussian_width,
                            spin=spin,
                        )
                        tas = jd_light - jd_dark
                        jdos_dark_cumulative += jd_dark
                        jdos_light_cumulative += jd_light
                        tas_cumulative += tas

                        new_i = 0
                        new_f = 0
                        if i <= self.vb[spin]:
                            new_i = (i - self.vb[spin])
                        elif i > self.vb[spin]:
                            new_i = (i - self.cb[spin]) + 1
                        if f <= self.vb[spin]:
                            new_f = (f - self.vb[spin])
                        elif f > self.vb[spin]:
                            new_f = (f - self.cb[spin]) + 1

                        if self.bs.is_spin_polarized:
                            spin_str = "up" if spin == Spin.up else "down"
                            key = (new_i, new_f, spin_str)
                        else:
                            key = (new_i, new_f)

                        jdos_dark_if[key] = jd_dark
                        jdos_light_if[key] = jd_light
                        tas_if[key] = tas

        return Tas(
            tas_cumulative,
            tas_if,
            jdos_light_cumulative,
            jdos_light_if,
            jdos_dark_cumulative,
            jdos_dark_if,
            energy_mesh_ev,
        )

    @classmethod
    def from_mpid(cls, mpid, temperature, concentration, bg=None, api_key=None):
        """
        Import the desired bandstructure and dos objects from the Materials Project
        database.

        Args:
            mpid: The Materials Project ID of the desired material.
            temperature:Temperature (K) of material we wish to investigate (affects the
                FD distribution)
            concentration: Carrier concentration (cm^-3) of holes and electrons (both are
                equivalent). Inversely proportional to pump-probe time delay.
            bg: The experimental bandgap (eV) of the material to be implemented. If the
                user wants to use the DFT-calculated bandgap, omit.

        Returns:
            A TASGenerator class with a uniform mode bandstructure & dos object, k-weights
            and a corrected bandgap.
        """
        mpr = MPRester(api_key=api_key)
        mp_dos = mpr.get_dos_by_material_id(mpid)
        mp_bs = mpr.get_bandstructure_by_material_id(mpid, line_mode=False)
        if bg is None:
            bg = mp_bs.get_band_gap()["energy"]
        else:
            mp_bs, mp_dos = set_bandgap(mp_bs, mp_dos, bg)
        kweights = get_kpoint_weights(mp_bs)
        return TASGenerator(mp_bs, kweights, mp_dos, bg, temperature, concentration)
