import itertools
from tqdm import tqdm
import warnings
import numpy as np

from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.dos import FermiDos, f0
from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp import optics
from pymatgen.io.vasp.inputs import UnknownPotcarWarning
from pymatgen.io.vasp.outputs import Vasprun, Waveder

from pytaser.kpoints import get_kpoint_weights
from pytaser.tas import Tas


warnings.filterwarnings("ignore", category=RuntimeWarning)


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
        raise ValueError(
            "DOS and band structure are not from the same calculation"
        )

    if bandstructure.is_metal():
        raise ValueError("System is a metal, cannot change bandgap")

    scissor = bandgap - bandstructure.get_band_gap()["energy"]
    midgap = (
        bandstructure.get_cbm()["energy"] + bandstructure.get_vbm()["energy"]
    ) / 2

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
            dens[abs(shift) : fermi_idx] = dos.densities[spin][
                : fermi_idx + shift
            ]
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


def occ_dependent_alpha(dfc, occs, spin=Spin.up, sigma=None, cshift=None):
    """Calculate the expected optical absorption given the groundstate orbital derivatives and
    eigenvalues (via dfc) and specified band occupancies.
    Templated from pymatgen.io.vasp.optics.epsilon_imag().

    Args:
        dfc: DielectricFunctionCalculator object
        occs: Array of band occupancies with shape (nspin, nbands, nkpoints)
        spin: Which spin channel to include.
        sigma: Smearing width (in eV) for broadening of the dielectric function (see
            https://www.vasp.at/wiki/index.php/SIGMA). If not set, uses the value of SIGMA from the
            underlying VASP WAVEDER calculation.
        cshift: Complex shift in the Kramers-Kronig transformation of the dielectric function (see
            https://www.vasp.at/wiki/index.php/CSHIFT). If not set, uses the value of CSHIFT from
            the underlying VASP WAVEDER calculation.

    Returns:
        (alpha_dict, tdm_array) where alpha_dict is a dictionary of band-to-band absorption,
        stimulated emission and summed contributions to the total overall absorption coefficient
        under illumination in cm^-1, and tdm_array is an array of shape (nbands, nbands,
        nkpoints) with the transition dipole matrix elements for each band-to-band transition.
    """
    if sigma is None:
        sigma = dfc.sigma
    if cshift is None:
        cshift = dfc.cshift
    egrid = np.arange(0, dfc.nedos * dfc.deltae, dfc.deltae)
    dielectric_dict = {
        key: np.zeros_like(egrid, dtype=np.complex128)
        for key in ["absorption", "emission", "both"]
    }
    # array of shape equal to cder but without the last two dimensions (ispin, idir)
    tdm_array = np.zeros_like(dfc.cder[:, :, :, 0, 0])  # ib, jb, ik, ispin

    norm_kweights = np.array(dfc.kweights) / np.sum(dfc.kweights)
    eigs_shifted = dfc.eigs - dfc.efermi
    rspin = (
        3 - dfc.cder.shape[3]
    )  # 2 for ISPIN = 1, 1 for ISPIN = 2 (spin-polarised)
    min_band0, max_band0 = np.min(np.where(dfc.cder)[0]), np.max(
        np.where(dfc.cder)[0]
    )
    min_band1, max_band1 = np.min(np.where(dfc.cder)[1]), np.max(
        np.where(dfc.cder)[1]
    )

    _, _, nk, _ = dfc.cder.shape[:4]
    iter_idx = [
        range(min_band0, max_band0 + 1),
        range(min_band1, max_band1 + 1),
        range(nk),
    ]
    num_ = (max_band0 - min_band0) * (max_band1 - min_band1) * nk
    spin_string = "up" if spin == Spin.up else "down"
    light_dark_string = (
        "under illumination"
        if any(
            occs[b][k] not in [0, 1]
            for b in range(min_band0, max_band0 + 1)
            for k in range(nk)
        )
        else "dark"
    )
    for ib, jb, ik in tqdm(
        itertools.product(*iter_idx),
        total=num_,
        desc=f"Calculating oscillator strengths (spin {spin_string}, "
        f"{light_dark_string})",
    ):
        ispin = 0 if spin == Spin.up else 1
        init_energy = eigs_shifted[ib, ik, ispin]
        final_energy = eigs_shifted[jb, ik, ispin]
        if final_energy > init_energy:
            init_occ = occs[ib][ik]
            final_occ = occs[jb][ik]

            A = (
                sum(
                    dfc.cder[ib, jb, ik, ispin, idir]
                    * np.conjugate(dfc.cder[ib, jb, ik, ispin, idir])
                    for idir in range(3)
                )
                / 3
            )
            decel = dfc.eigs[jb, ik, ispin] - dfc.eigs[ib, ik, ispin]
            matrix_el_wout_occ_factor = np.abs(A) * norm_kweights[ik] * rspin
            tdm_array[ib, jb, ik] = np.abs(A) * rspin * decel  # kweight and occ factor
            # already accounted for with JDOS

            abs_occ_factor = init_occ * (1 - final_occ)
            em_occ_factor = (1 - init_occ) * final_occ
            both_occ_factor = init_occ - final_occ

            abs_matrix_el = abs_occ_factor * matrix_el_wout_occ_factor
            em_matrix_el = em_occ_factor * matrix_el_wout_occ_factor
            both_matrix_el = both_occ_factor * matrix_el_wout_occ_factor

            if dfc.ismear == 0:  # error in pymatgen, TODO: PR!
                ismear = -0.1
            else:
                ismear = dfc.ismear
            smeared_wout_matrix_el = optics.get_delta(
                x0=decel,
                sigma=sigma,
                nx=dfc.nedos,
                dx=dfc.deltae,
                ismear=ismear,
            )

            dielectric_dict["absorption"] += (
                smeared_wout_matrix_el * abs_matrix_el
            )
            dielectric_dict["emission"] += (
                smeared_wout_matrix_el * em_matrix_el
            )
            dielectric_dict["both"] += (
                smeared_wout_matrix_el * both_matrix_el
            )

    tdm_array = tdm_array.real  # real part of A is the TDM (imag part is zero after taking the
    # complex conjugate)
    alpha_dict = {}
    for key, dielectric in dielectric_dict.items():
        eps_in = dielectric * optics.edeps * np.pi / dfc.volume
        eps = optics.kramers_kronig(
            eps_in, nedos=dfc.nedos, deltae=dfc.deltae, cshift=cshift
        )
        eps += 1.0 + 0.0j
        dielectric_dict[key] = eps  # complex dielectric function

        # convert to alpha:
        n = np.sqrt(eps)  # complex refractive index
        alpha = (
            n.imag * egrid * 4 * np.pi / 1.23984212e-4
        )  # absorption coefficient in cm^-1
        alpha_dict[key] = alpha

    return alpha_dict, tdm_array


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
        vbm_index[spin] = np.where(np.all(spin_bands <= bs.efermi, axis=1))[
            0
        ].max()
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
        dfc: Pymatgen-based DielectricFunctionCalculator object (for computing oscillator strengths)

    Attributes:
        bs: Pymatgen bandstructure object
        kpoint_weights: k-point weights (degeneracies).
        dos: Pymatgen-based dos object
        dfc: Pymatgen-based DielectricFunctionCalculator object (for computing oscillator strengths)
        bg_centre: Energy (eV) of the bandgap centre.
        vb: Spin dict detailing the valence band maxima.
        cb: Spin dict detailing the conduction band minima
    """

    def __init__(self, bs, kpoint_weights, dos, dfc=None):
        self.bs = bs
        self.kpoint_weights = kpoint_weights
        self.dos = FermiDos(dos)
        self.bg_centre = (bs.get_cbm()["energy"] + bs.get_vbm()["energy"]) / 2
        self.dfc = dfc

        if self.bs.is_metal():
            raise ValueError("System is metallic, cannot compute TAS")

        self.vb = get_cbm_vbm_index(self.bs)[0]
        self.cb = get_cbm_vbm_index(self.bs)[1]

    @classmethod
    def from_vasp_objects(cls, vasprun_file, waveder_file=None):
        """Create a TASGenerator object from VASP output files."""
        warnings.filterwarnings("ignore", category=UnknownPotcarWarning)
        warnings.filterwarnings(
            "ignore", message="No POTCAR file with matching TITEL fields"
        )
        vr = Vasprun(vasprun_file)
        if waveder_file:
            waveder = Waveder.from_binary(waveder_file)
            dfc = optics.DielectricFunctionCalculator.from_vasp_objects(
                vr, waveder
            )
        else:
            dfc = None
        return cls(
            vr.get_band_structure(),
            vr.actual_kpoints_weights,
            vr.complete_dos,
            dfc,
        )

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
        q_fermi_e = self.dos.get_fermi(
            -conc, temp
        )  # quasi-electron fermi level
        q_fermi_h = self.dos.get_fermi(conc, temp)  # quasi-hole fermi level

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

    def generate_tas(
        self,
        temp,
        conc,
        energy_min=0,
        energy_max=5,
        gaussian_width=0.1,
        cshift=None,
        step=0.01,
        light_occs=None,
        dark_occs=None,
    ):
        """
        Generates TAS spectra based on inputted occupancies, and a specified energy mesh.

        Args:
            temp: Temperature (K) of material we wish to investigate (affects the FD distribution)
            conc: Carrier concentration (cm^-3) of holes and electrons (both are equivalent).
                Inversely proportional to pump-probe time delay.
            energy_min: Minimum band transition energy to consider for energy mesh (eV)
            energy_max: Maximum band transition energy to consider for energy mesh (eV)
            gaussian_width: Width of gaussian curve
            cshift: Complex shift in the Kramers-Kronig transformation of the dielectric function
                (see https://www.vasp.at/wiki/index.php/CSHIFT). If not set, uses the value of
                CSHIFT from the underlying VASP WAVEDER calculation. (only relevant if the
                TASGenerator has been generated from VASP outputs)
            step: Interval between energy points in the energy mesh.
            light_occs: Optional input parameter for occupancies of material under light [dict]
            dark_occs: Optional input parameter for occupancies of material in dark [dict]

        Returns:
            TAS class containing the following inputs;
                - tas_cumulative: overall deltaT TAS spectrum for a material under the
                    specified conditions.
                - tas_if: deltaT TAS spectrum across the energy mesh for a specific band
                    transition i (initial) -> f (final) [dict]
                - jdos_light_cumulative: overall JDOS (pump-on) for a material under the
                    specified conditions
                - jdos_light_if: JDOS (pump-on) across the energy mesh for a specific band
                    transition i (initial) -> f (final) [dict]
                - jdos_dark_cumulative: overall JDOS (pump-off) for a material under the
                    specified conditions
                - jdos_dark_if: JDOS (pump-off) across the energy mesh for a specific band
                    transition i (initial) -> f (final) [dict]
                - energy_mesh_ev: Energy mesh of spectra in eV, with an interval of 'step'.
                - bandgap_ev: Bandgap of the system, in eV, rounded to 2 decimal points
        """
        occs_light = light_occs
        occs_dark = dark_occs
        if light_occs is None:
            occs_light = self.band_occupancies(temp, conc, dark=False)
        if dark_occs is None:
            occs_dark = self.band_occupancies(temp, conc)

        bandgap_ev = round(self.bs.get_band_gap()["energy"], 2)
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
                            occs_light[spin],
                            energy_mesh_ev,
                            self.kpoint_weights,
                            gaussian_width,
                            spin=spin,
                        )
                        jd_dark = jdos(
                            self.bs,
                            f,
                            i,
                            occs_dark[spin],
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
                            new_i = i - self.vb[spin]
                        elif i > self.vb[spin]:
                            new_i = (i - self.cb[spin]) + 1
                        if f <= self.vb[spin]:
                            new_f = f - self.vb[spin]
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
            bandgap_ev,
        )

    @classmethod
    def from_mpid(cls, mpid, bg=None, api_key=None, mpr=None):
        """
        Import the desired bandstructure and dos objects from the legacy Materials Project
        database.

        Args:
            mpid: The Materials Project ID of the desired material.
            bg: The experimental bandgap (eV) of the material. If None, the band gap
                of the MP calculation will be used.
            api_key: The user's Materials Project API key.
            mpr: An MPRester object if already generated by user.

        Returns:
            A TASGenerator object.
        """
        if mpr is None:
            mpr = MPRester(api_key=api_key)
        mp_dos = mpr.get_dos_by_material_id(mpid)
        mp_bs = mpr.get_bandstructure_by_material_id(mpid, line_mode=False)
        if bg is not None:
            mp_bs, mp_dos = set_bandgap(mp_bs, mp_dos, bg)
        kweights = get_kpoint_weights(mp_bs)
        return TASGenerator(mp_bs, kweights, mp_dos)
