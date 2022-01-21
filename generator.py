from collections import OrderedDict

import numpy as np
from abipy.tools.numtools import gaussian
from pymatgen.core import Structure
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.dos import FermiDos, f0
from pymatgen.ext.matproj import MPRester

from tas import Tas


def get_kpoint_weights(bandstructure, time_reversal=True, symprec=0.1):
    """
    Function to calculate the kpoint_weights for non-magnetic materials (non-metals).

    Args:
        bandstructure: PMG bandstructure object
        time_reversal:
        symprec: Symmetry precision in Angstrom.(Lower value is more precise, but
            computationally more expensive)
    Returns:
        k-point_weights
    """

    kpoints = get_kpoints_from_bandstructure(bandstructure)
    _, _, _, _, _, kp_mapping = expand_kpoints(
        bandstructure.structure,
        kpoints,
        symprec=symprec,
        time_reversal=time_reversal,
        return_mapping=True,
    )
    weights = np.unique(kp_mapping, return_counts=True)[1].astype(float)
    weights /= np.sum(weights)
    return weights


def get_kpoints_from_bandstructure(bandstructure, cartesian=False):
    if cartesian:
        kpoints = np.array([k.cart_coords for k in bandstructure.kpoints])
    else:
        kpoints = np.array([k.frac_coords for k in bandstructure.kpoints])

    return kpoints


def expand_kpoints(
    structure,
    kpoints,
    symprec=0.01,
    return_mapping=False,
    time_reversal=True,
):
    kpoints = np.array(kpoints).round(8)

    # due to limited input precision of the k-points, the mesh is returned as a float
    mesh, is_shifted = get_mesh_from_kpoint_diff(kpoints)

    if is_shifted:
        shift = np.array([1, 1, 1])
    else:
        shift = np.array([0, 0, 0])

    # to avoid issues to limited input precision, recalculate the input k-points
    # so that the mesh is integer and the k-points are not truncated
    # to a small precision
    addresses = np.rint((kpoints + shift / (mesh * 2)) * mesh)
    mesh = np.rint(mesh)
    kpoints = addresses / mesh - shift / (mesh * 2)

    rotations, translations, is_tr = get_reciprocal_point_group_operations(
        structure, symprec=symprec, time_reversal=time_reversal
    )
    n_ops = len(rotations)
    # rotate all-kpoints
    all_rotated_kpoints = []
    for r in rotations:
        all_rotated_kpoints.append(np.dot(r, kpoints.T).T)
    all_rotated_kpoints = np.concatenate(all_rotated_kpoints)

    # map to first BZ
    all_rotated_kpoints -= np.rint(all_rotated_kpoints)
    all_rotated_kpoints = all_rotated_kpoints.round(8)

    # zone boundary consistent with VASP not with spglib
    all_rotated_kpoints[all_rotated_kpoints == -0.5] = 0.5

    # Find unique points
    unique_rotated_kpoints, unique_idxs = np.unique(
        all_rotated_kpoints, return_index=True, axis=0
    )

    # find integer addresses
    unique_addresses = (unique_rotated_kpoints + shift / (mesh * 2)) * mesh
    unique_addresses -= np.rint(unique_addresses)
    in_uniform_mesh = (np.abs(unique_addresses) < 1e-5).all(axis=1)

    n_mapped = int(np.sum(in_uniform_mesh))
    n_expected = int(np.product(mesh))
    if n_mapped != n_expected:
        raise ValueError(f"Expected {n_expected} points but found {n_mapped}")

    full_kpoints = unique_rotated_kpoints[in_uniform_mesh]
    full_idxs = unique_idxs[in_uniform_mesh]

    if not return_mapping:
        return full_kpoints

    op_mapping = np.floor(full_idxs / len(kpoints)).astype(int)
    kp_mapping = (full_idxs % len(kpoints)).astype(int)

    return full_kpoints, rotations, translations, is_tr, op_mapping, kp_mapping


def get_mesh_from_kpoint_diff(kpoints, ktol=1e-5):
    kpoints = np.array(kpoints)

    # whether the k-point mesh is shifted or Gamma centered mesh
    is_shifted = np.min(np.linalg.norm(kpoints, axis=1)) > 1e-6

    unique_a = np.unique(kpoints[:, 0])
    unique_b = np.unique(kpoints[:, 1])
    unique_c = np.unique(kpoints[:, 2])

    if len(unique_a) == 1:
        na = 1
    else:
        # filter very small changes, with a tol of 5e-4 this means k-point meshes
        # denser than 2000x2000x2000 will be treated as numerical noise. Meshes
        # this dense are extremely unlikely
        diff = np.diff(unique_a)
        diff = diff[diff > ktol]
        na = 1 / np.min(diff[diff > ktol])

    if len(unique_b) == 1:
        nb = 1
    else:
        diff = np.diff(unique_b)
        nb = 1 / np.min(diff[diff > ktol])

    if len(unique_c) == 1:
        nc = 1
    else:
        diff = np.diff(unique_c)
        nc = 1 / np.min(diff[diff > ktol])

    # due to limited precision of the input k-points, the mesh is returned as a float
    return np.array([na, nb, nc]), is_shifted


def get_reciprocal_point_group_operations(
    structure: Structure,
    symprec: float = 0.01,
    time_reversal: bool = True,
):
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    sga = SpacegroupAnalyzer(structure, symprec=symprec)
    if sga.get_symmetry_dataset() is None:
        # sometimes default angle tolerance doesn't work as expected
        sga = SpacegroupAnalyzer(structure, symprec=symprec, angle_tolerance=-1)

    rotations = sga.get_symmetry_dataset()["rotations"].transpose((0, 2, 1))
    translations = sga.get_symmetry_dataset()["translations"]
    is_tr = np.full(len(rotations), False, dtype=bool)

    if time_reversal:
        rotations = np.concatenate([rotations, -rotations])
        translations = np.concatenate([translations, -translations])
        is_tr = np.concatenate([is_tr, ~is_tr])

        rotations, unique_ops = np.unique(rotations, axis=0, return_index=True)
        translations = translations[unique_ops]
        is_tr = is_tr[unique_ops]

    # put identity first and time-reversal last
    sort_idx = np.argsort(np.abs(rotations - np.eye(3)).sum(axis=(1, 2)) + is_tr * 10)

    return rotations[sort_idx], translations[sort_idx], is_tr[sort_idx]


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
    print(shift)
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


def jdos_simple(
    bs, f, i, occs, energies, k_index, kweights, gaussian_width=0.2, sppol=False
):
    """
    Similar to jdos_simple, but also accounts for Spin.down transitions (spin polarisation).

    Args:
        bs: bandstructure object
        f: final band
        i: initial band
        occs: occupancies over all bands.
        energies: energy mesh (eV)
        k_index: kpoint grid for each band.
        kweights: k-point weights
        gaussian_width: width of gaussian plot.
        sppol: bool identifying the material as spin polarised.

    Returns:
        Cumulative JDOS value for a specific i->f transition, with consideration of
        partial occupancy and spin polarisation.
    """

    jdos = np.zeros(len(energies))
    if sppol == False:
        #         spin_sign = +1
        for k in k_index:
            final_occ = occs[f][k]
            init_energy = bs.bands[Spin.up][i][k]
            final_energy = bs.bands[Spin.up][f][k]
            init_occ = occs[i][k]
            k_weight = kweights[k]
            # factor_old = k_weight * init_occ * (1 - final_occ)  THIS IS THE OLD EQUATION THAT WE USED.
            factor = k_weight * (
                (init_occ * (1 - final_occ)) - (final_occ * (1 - init_occ))
            )
            jdos += factor * gaussian(
                energies, gaussian_width, center=final_energy - init_energy
            )

    # Doesn't consider Spin.up -> Spin.down transitions (or vice-versa)
    else:
        #         spin_sign = -1
        for k in k_index:
            final_occ = occs[f][k]
            init_energy = bs.bands[Spin.down][i][k]
            final_energy = bs.bands[Spin.down][f][k]
            init_occ = occs[i][k]
            k_weight = kweights[k]
            factor = k_weight * init_occ * (1 - final_occ)
            jdos += factor * gaussian(
                energies, gaussian_width, center=final_energy - init_energy
            )

    return jdos


class TASGenerator:
    """
    Class to generate a TAS spectrum (decomposed and cumulative) from a bandstructure and
    dos object.

    Args:
        bs: Pymatgen-based bandstructure object
        kpoint_weights: kpoint weights either found by the function or inputted.
        dos:Pymatgen-based dos object
        bandgap: Experimental or Inputted  bandgap of a material.
        temp: Temperature (K) of material we wish to investigate (affects the FD distribution)
        conc: Carrier concentration (cm^-3) of holes and electrons (both are equivalent).
            Inversely proportional to pump-probe time delay.

    Attributes:
        bs: Pymatgen-based bandstructure object
        kpoint_weights: kpoint weights either found by the function or inputted
        dos: Pymatgen-based dos object
        bg: Experimental or Inputted bandgap of a material.
        c_min_energy: Energy (eV) of the conduction band minima.
        v_max_energy: Energy (eV) of the valence band maxima.
        bg_centre: Energy (eV) of the bandgap centre.
        number_bands: Number of bands of material
        c_min_index: Conduction minimum band number
        v_max_index: Valence maximum band number
        band_index: Ordered list index of all bands that are spin-up bands
        kpoint_index: Ordered list index of all kpoints in spin-up bands
        sppol: Boolean determining whether material is spin-polarised or not.
        band_index_down: Ordered list index of all bands that are spin-down. (Only
            relevant if material is spin-polarised)
        kpoint_index_down: Ordered list index of all kpoints in spin-down. (Only
            relevant if material is spin-polarised)
        occs_light: Occupancies of spin-up bands for a pump-on system, across all kpoints
        occs_dark: Occupancies of spin-up bands for a pump-off system, across all kpoints
        occs_light_down: Occupancies of spin-down bands for a pump-on system, across all
            kpoints
        occs_dark_down: Occupancies of spin-down bands for a pump-off system, across all
            kpoints
    """

    def __init__(self, bs, kpoint_weights, dos, bandgap, temp, conc):
        self.bs = bs
        self.kpoint_weights = kpoint_weights
        self.dos = dos
        self.bg = bandgap
        self.c_min_energy = bs.get_cbm()["energy"]
        self.v_max_energy = bs.get_vbm()["energy"]
        self.bg_centre = (self.c_min_energy + self.v_max_energy) / 2
        self.number_bands = bs.nb_bands
        self.c_min_index = bs.get_cbm()["band_index"]
        self.v_max_index = bs.get_vbm()["band_index"]

        band_index = []
        for band_indices, band_values in enumerate(self.bs.bands[Spin.up]):
            band_index += [band_indices]
        self.band_index = band_index
        kpoint_index = []
        for k_indices, k_values in enumerate(self.bs.bands[Spin.up][0]):
            kpoint_index += [k_indices]
        self.kpoint_index = kpoint_index

        if float(bandgap) < 0.01:
            raise ValueError("Bandgap is smaller than 0.01 eV; cannot compute TAS")

        if len(self.bs.bands) > 1:
            sppol = True
        else:
            sppol = False
        self.sppol = sppol

        band_index_down = []
        kpoint_index_down = []
        if self.sppol == True:
            for band_indices_down, band_values_down in enumerate(
                self.bs.bands[Spin.down]
            ):
                band_index_down += [band_indices_down]
            for k_indices_down, k_values_down in enumerate(self.bs.bands[Spin.down][0]):
                kpoint_index_down += [k_indices_down]
        self.band_index_down = band_index_down
        self.kpoint_index_down = kpoint_index_down

        light_occs = None
        dark_occs = None
        light_occs_down = None
        dark_occs_down = None
        if self.sppol == True:
            light_occs = TASGenerator.band_occupancies_sppol(
                self, temp, conc, dark=False
            )[0]
            light_occs_down = TASGenerator.band_occupancies_sppol(
                self, temp, conc, dark=False
            )[1]
            dark_occs = TASGenerator.band_occupancies_sppol(self, temp, conc)[0]
            dark_occs_down = TASGenerator.band_occupancies_sppol(self, temp, conc)[1]
        elif self.sppol == False:
            light_occs = TASGenerator.band_occupancies(self, temp, conc, dark=False)
            dark_occs = TASGenerator.band_occupancies(self, temp, conc)
        self.occs_light = light_occs
        self.occs_dark = dark_occs
        self.occs_light_down = light_occs_down
        self.occs_dark_down = dark_occs_down

    def band_occupancies(self, temp, conc, dark=True):
        """
        Gives band occupancies for a typical non spin-polarised system.

        Args:
            temp: Temperature (K) of material we wish to investigate (affects the FD
                distribution)
            conc: Carrier concentration (cm^-3) of holes and electrons (both are
                equivalent). Inversely proportional to pump-probe time delay.
            dark: Bool; dark = True indicates pump is on.

        Returns:
            occs: Band occupancies for all bands for a material, across all k-points.
        """
        if dark == True:
            hole_mask = self.bs.bands[Spin.up] < self.bg_centre
            elec_mask = self.bs.bands[Spin.up] > self.bg_centre

            occs = np.zeros_like(self.bs.bands[Spin.up])
            # fully occupied hole mask, completely empty electron mask
            occs[hole_mask] = 1
            occs[elec_mask] = 0

        else:
            fermidos = FermiDos(self.dos)
            # quasi-electron fermi level
            q_fermi_e = fermidos.get_fermi(-conc, temp)
            # quasi-hole fermi level
            q_fermi_h = fermidos.get_fermi(conc, temp)

            # Calculate the occupancies at the initial and final energy according to the
            # Fermi-Dirac distribution
            # q-e-f
            electron_occs = f0(self.bs.bands[Spin.up], q_fermi_e, temp)
            # q-h-f
            hole_occs = f0(self.bs.bands[Spin.up], q_fermi_h, temp)

            hole_mask = self.bs.bands[Spin.up] < self.bg_centre
            elec_mask = self.bs.bands[Spin.up] > self.bg_centre

            #   set occupancies equal to computed by f0
            occs = np.zeros_like(hole_occs)
            occs[hole_mask] = hole_occs[hole_mask]
            occs[elec_mask] = electron_occs[elec_mask]

        return occs

    def band_occupancies_sppol(self, temp, conc, dark=True):
        """
        Gives band occupancies for a spin-polarised system.

        Args:
            temp: Temperature of material we wish to investigate (affects the FD
                distribution)
            conc: Carrier concentration of holes and electrons (both are the same).
                Inversely proportional to pump-probe time delay.
            dark: Bool; dark = True indicates pump is on.

        Returns:
            occs: Band occupancies for all spin-up bands, across all k-points.
            occs_down: Band occupancies for the spin-down bands, across all k-points
        """

        # THIS FUNCTION IS ONLY TO BE USED FOR SPIN POLARISED SYSTEMS. IF NOT SPPOL, CALL
        # THE OTHER FUNCTION.

        #   Calculate the quasi-Fermi levels
        fermidos = FermiDos(self.dos)
        #   quasi-electron fermi level
        q_fermi_e = fermidos.get_fermi(-conc, temp)
        #   quasi-hole fermi level
        q_fermi_h = fermidos.get_fermi(conc, temp)

        if dark == True:
            # DARK
            hole_mask = self.bs.bands[Spin.up] < self.bg_centre
            elec_mask = self.bs.bands[Spin.up] > self.bg_centre

            occs = np.zeros_like(self.bs.bands[Spin.up])
            # fully occupied hole mask, completely empty electron mask
            occs[hole_mask] = 1
            occs[elec_mask] = 0

            # DOWN
            hole_mask_down = self.bs.bands[Spin.down] < self.bg_centre
            elec_mask_down = self.bs.bands[Spin.down] > self.bg_centre

            occs_down = np.zeros_like(self.bs.bands[Spin.down])
            # fully occupied hole mask, completely empty electron mask
            occs_down[hole_mask_down] = 1
            occs_down[elec_mask_down] = 0

        else:
            # LIGHT
            # Calculate the occupancies at the initial and final energy according to the
            # Fermi-Dirac distribution
            # q-e-f
            electron_occs = f0(self.bs.bands[Spin.up], q_fermi_e, temp)
            # q-h-f
            hole_occs = f0(self.bs.bands[Spin.up], q_fermi_h, temp)

            hole_mask = self.bs.bands[Spin.up] < self.bg_centre
            elec_mask = self.bs.bands[Spin.up] > self.bg_centre

            #   set occupancies equal to computed by f0
            occs = np.zeros_like(hole_occs)
            occs[hole_mask] = hole_occs[hole_mask]
            occs[elec_mask] = electron_occs[elec_mask]

            # DOWN

            # q-e-f
            electron_occs_down = f0(self.bs.bands[Spin.down], q_fermi_e, temp)
            # q-h-f
            hole_occs_down = f0(self.bs.bands[Spin.down], q_fermi_h, temp)

            hole_mask_down = self.bs.bands[Spin.down] < self.bg_centre
            elec_mask_down = self.bs.bands[Spin.down] > self.bg_centre

            #   set occupancies equal to computed by f0
            occs_down = np.zeros_like(hole_occs_down)
            occs_down[hole_mask_down] = hole_occs_down[hole_mask_down]
            occs_down[elec_mask_down] = electron_occs_down[elec_mask_down]

        return occs, occs_down

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
        jdos_light_if = OrderedDict()
        jdos_dark_if = OrderedDict()
        tas_if = OrderedDict()
        # denom_TAS_if = OrderedDict()
        tas_cumulative = np.zeros(len(energy_mesh_ev))
        jdos_dark_cumulative = np.zeros(len(energy_mesh_ev))
        jdos_light_cumulative = np.zeros(len(energy_mesh_ev))
        # denom_TAS_cumulative = np.zeros(len(energy_mesh_ev))
        if self.sppol == False:
            for i in self.band_index:
                for f in self.band_index:
                    if f > i:
                        jd_light = jdos_simple(
                            self.bs,
                            f,
                            i,
                            self.occs_light,
                            energy_mesh_ev,
                            self.kpoint_index,
                            self.kpoint_weights,
                        )
                        jd_dark = jdos_simple(
                            self.bs,
                            f,
                            i,
                            self.occs_dark,
                            energy_mesh_ev,
                            self.kpoint_index,
                            self.kpoint_weights,
                        )
                        tas = jd_light - jd_dark
                        jdos_dark_if[(i, f)] = jd_dark
                        jdos_light_if[(i, f)] = jd_light
                        jdos_dark_cumulative += jd_dark
                        jdos_light_cumulative += jd_light
                        tas_if[(i, f)] = tas
                        tas_cumulative += tas
                        # denom_TAS_if[(i,f, 'spin = up')] = tas/jd_light
                        # denom_TAS_cumulative = tas/jd_light

        elif self.sppol == True:
            for i in self.band_index:
                for f in self.band_index:
                    if f > i:
                        jd_light = jdos_simple(
                            self.bs,
                            f,
                            i,
                            self.occs_light,
                            energy_mesh_ev,
                            self.kpoint_index,
                            self.kpoint_weights,
                        )
                        jd_dark = jdos_simple(
                            self.bs,
                            f,
                            i,
                            self.occs_dark,
                            energy_mesh_ev,
                            self.kpoint_index,
                            self.kpoint_weights,
                        )
                        tas = jd_light - jd_dark
                        jdos_dark_if[(i, f, "spin = up")] = jd_dark
                        jdos_light_if[(i, f, "spin = up")] = jd_light
                        jdos_dark_cumulative += jd_dark
                        jdos_light_cumulative += jd_light
                        tas_if[(i, f, "spin = up")] = tas
                        tas_cumulative += tas
                        # denom_TAS_if[(i,f, 'spin = up')] = tas/jd_light
                        # denom_TAS_cumulative = tas/jd_light

            for i in self.band_index_down:
                for f in self.band_index_down:
                    if f > i:
                        jd_light_down = jdos_simple(
                            self.bs,
                            f,
                            i,
                            self.occs_light_down,
                            energy_mesh_ev,
                            self.kpoint_index_down,
                            self.kpoint_weights,
                            sppol=True,
                        )
                        jd_dark_down = jdos_simple(
                            self.bs,
                            f,
                            i,
                            self.occs_dark_down,
                            energy_mesh_ev,
                            self.kpoint_index_down,
                            self.kpoint_weights,
                            sppol=True,
                        )
                        tas_down = jd_light_down - jd_dark_down
                        jdos_dark_if[(i, f, "spin = down")] = jd_dark_down
                        jdos_light_if[(i, f, "spin = down")] = jd_light_down
                        jdos_dark_cumulative += jd_dark_down
                        jdos_light_cumulative += jd_light_down
                        tas_if[(i, f, "spin = down")] = tas_down
                        tas_cumulative += tas_down
                        # denom_TAS_if[(i,f, 'spin = down')] = tas_down/jd_light_down
                        # denom_TAS_cumulative = tas_down/jd_light_down

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
    def from_mpid(cls, api, mpid, temperature, concentration, bg=None):
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
        mpr = MPRester(api_key=api)
        mp_dos = mpr.get_dos_by_material_id(mpid)
        mp_bs = mpr.get_bandstructure_by_material_id(mpid, line_mode=False)
        if bg is None:
            bg = mp_bs.get_band_gap()["energy"]
        else:
            mp_bs, mp_dos = set_bandgap(mp_bs, mp_dos, bg)
        kweights = get_kpoint_weights(mp_bs)
        return TASGenerator(mp_bs, kweights, mp_dos, bg, temperature, concentration)
