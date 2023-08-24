import warnings
from multiprocessing import Pool, cpu_count, Array

import numpy as np
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.dos import FermiDos, f0
from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp import optics
from pymatgen.io.vasp.inputs import UnknownPotcarWarning
from pymatgen.io.vasp.outputs import Vasprun, Waveder
from tqdm import tqdm

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
        factor = k_weight * (init_occ - final_occ)

        jdos += factor * gaussian(
            energies, gaussian_width, center=final_energy - init_energy
        )
    return jdos


def _calculate_oscillator_strength(args):
    """Calculates the oscillator strength of a single band-band transition"""
    if len(args) == 9:  # shared memory arrays
        ib, jb, ik, rspin, spin, sigma, nedos, deltae, ismear = args

        # Convert shared arrays back to numpy arrays
        cder = np.frombuffer(_cder.get_obj()).reshape(_cder_shape)
        occs = np.frombuffer(_occs.get_obj()).reshape(_occs_shape)
        eigs_shifted = np.frombuffer(_eigs_shifted.get_obj()).reshape(
            _eigs_shifted_shape
        )
        norm_kweights = np.frombuffer(_norm_kweights.get_obj()).reshape(
            _norm_kweights_shape
        )

    else:  # normal function call
        (
            ib,
            jb,
            ik,
            rspin,
            spin,
            sigma,
            nedos,
            deltae,
            ismear,
            cder,
            occs,
            eigs_shifted,
            norm_kweights,
        ) = args

    ispin = 0 if spin == Spin.up else 1

    A = np.sum(np.abs(cder[ib, jb, ik, ispin, :3] ** 2)) / 3
    decel = eigs_shifted[jb, ik, ispin] - eigs_shifted[ib, ik, ispin]
    matrix_el_wout_occ_factor = np.abs(A) * norm_kweights[ik] * rspin
    tdm = (
        np.abs(A) * rspin * decel
    )  # kweight and occ factor already accounted for with JDOS

    init_occ = occs[ib][ik]
    final_occ = occs[jb][ik]
    abs_occ_factor = init_occ * (1 - final_occ)
    em_occ_factor = (1 - init_occ) * final_occ
    both_occ_factor = init_occ - final_occ

    abs_matrix_el = abs_occ_factor * matrix_el_wout_occ_factor
    em_matrix_el = em_occ_factor * matrix_el_wout_occ_factor
    both_matrix_el = both_occ_factor * matrix_el_wout_occ_factor

    smeared_wout_matrix_el = optics.get_delta(
        x0=decel,
        sigma=sigma,
        nx=nedos,
        dx=deltae,
        ismear=ismear,
    )

    absorption = smeared_wout_matrix_el * abs_matrix_el
    emission = smeared_wout_matrix_el * em_matrix_el
    both = smeared_wout_matrix_el * both_matrix_el

    return absorption, emission, both, tdm, ib, jb, ik


def get_nonzero_band_transitions(
    occs,
    eigs_shifted,
    rspin,
    spin,
    sigma,
    nedos,
    deltae,
    ismear,
    min_band,
    max_band,
    nk,
):
    """Helper function to filter band transitions before (multi)processing"""
    ispin_idx = 0 if spin == Spin.up else 1

    ib_vals, jb_vals, ik_vals = np.meshgrid(
        np.arange(min_band, max_band + 1),
        np.arange(min_band, max_band + 1),
        np.arange(nk),
        indexing="ij",
    )

    init_energy_vals = eigs_shifted[ib_vals, ik_vals, ispin_idx]
    final_energy_vals = eigs_shifted[jb_vals, ik_vals, ispin_idx]
    init_occ_vals = occs[ib_vals, ik_vals]
    final_occ_vals = occs[jb_vals, ik_vals]

    condition = (
        final_energy_vals > init_energy_vals
    ) & (  # one-way transitions (emission is accounted for by occupancy factor)
        np.abs(init_occ_vals - final_occ_vals)
        > 0.01  # non-negligible occupancy factor
    )

    return list(
        zip(
            ib_vals[condition].ravel(),
            jb_vals[condition].ravel(),
            ik_vals[condition].ravel(),
            [rspin] * np.sum(condition),
            [spin] * np.sum(condition),
            [sigma] * np.sum(condition),
            [nedos] * np.sum(condition),
            [deltae] * np.sum(condition),
            [ismear] * np.sum(condition),
        )
    )


def init_shared_memory(cder, occs, eigs_shifted, norm_kweights):
    global _cder, _occs, _eigs_shifted, _norm_kweights
    global _cder_shape, _occs_shape, _eigs_shifted_shape, _norm_kweights_shape

    _cder_shape = cder.shape
    _cder = Array("d", cder.ravel())
    _occs_shape = occs.shape
    _occs = Array("d", occs.ravel())
    _eigs_shifted_shape = eigs_shifted.shape
    _eigs_shifted = Array("d", eigs_shifted.ravel())
    _norm_kweights_shape = norm_kweights.shape
    _norm_kweights = Array("d", norm_kweights.ravel())


def occ_dependent_alpha(
    dfc,
    occs,
    spin=Spin.up,
    sigma=None,
    cshift=None,
    processes=None,
    energy_max=6,
):
    """
    Calculate the expected optical absorption given the groundstate orbital derivatives and
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
        processes: Number of processes to use for multiprocessing. If not set, defaults to one
            less than the number of CPUs available.
        energy_max: Maximum band transition energy to consider (in eV). A minimum range of -6 eV to +6
            eV is considered regardless of `energy_max`, to ensure a reasonable dielectric function output.

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
    tdm_array = np.zeros_like(dfc.cder[:, :, :, 0, 0])  # ib, jb, ik

    norm_kweights = np.array(dfc.kweights) / np.sum(dfc.kweights)
    eigs_shifted = dfc.eigs - dfc.efermi
    rspin = (
        3 - dfc.cder.shape[3]
    )  # 2 for ISPIN = 1, 1 for ISPIN = 2 (spin-polarised)

    # set band range to consider, based on energy_max:
    # use -6 eV to +6 eV as minimum range, extended to -energy_max/+energy_max if energy_max > 6 eV
    max_band_energy = max(6, energy_max)
    min_band = np.min((eigs_shifted > -max_band_energy).nonzero()[0])
    max_band = np.max((eigs_shifted < max_band_energy).nonzero()[0])

    _, _, nk, _ = dfc.cder.shape[:4]
    spin_string = "up" if spin == Spin.up else "down"
    light_dark_string = (
        "under illumination"
        if any(
            occs[b][k] not in [0, 1]
            for b in range(min_band, max_band + 1)
            for k in range(nk)
        )
        else "dark"
    )

    nonzero_transition_args = get_nonzero_band_transitions(
        occs,
        eigs_shifted,
        rspin,
        spin,
        sigma,
        dfc.nedos,
        dfc.deltae,
        dfc.ismear,
        min_band,
        max_band,
        nk,
    )

    if processes is None:
        processes = cpu_count() - 1

    shared_memory_args = (dfc.cder, occs, eigs_shifted, norm_kweights)
    if (
        len(nonzero_transition_args) <= 2.5e6
    ) or processes < 2:  # don't use multiprocessing for small arrays
        # append shared memory arguments to the list of arguments for each process
        nonzero_transition_args = [
            (*arg, *shared_memory_args) for arg in nonzero_transition_args
        ]

    if (
        processes > 1 and len(nonzero_transition_args) > 2.5e6
    ):  # quicker without multiprocessing below this arg length
        with Pool(
            processes,
            initializer=init_shared_memory,
            initargs=shared_memory_args,
        ) as pool:
            results = pool.map(
                _calculate_oscillator_strength,
                tqdm(
                    nonzero_transition_args,
                    desc=f"Calculating oscillator strengths (spin {spin_string}, {light_dark_string})",
                ),
            )
    else:
        results = [
            _calculate_oscillator_strength(arg)
            for arg in tqdm(
                nonzero_transition_args,
                desc=f"Calculating oscillator strengths (spin {spin_string}, {light_dark_string})",
            )
        ]

    results_array = np.array(
        results, dtype=object
    )  # dtype=object ensures that data is preserved

    # Accumulate the results
    dielectric_dict["absorption"] += results_array[:, 0].sum()
    dielectric_dict["emission"] += results_array[:, 1].sum()
    dielectric_dict["both"] += results_array[:, 2].sum()

    indices = results_array[:, 4:7].astype(int).T  # ib, jb, ik
    tdm_values = results_array[:, 3]
    tdm_array[indices[0], indices[1], indices[2]] = tdm_values

    tdm_array = (
        tdm_array.real
    )  # real part of A is the TDM (imag part is zero after taking the complex conjugate)
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
    def from_vasp_outputs(cls, vasprun_file, waveder_file=None, bg=None):
        """
        Create a TASGenerator object from VASP output files.

        Args:
            vasprun_file: Path to vasprun.xml file (to generate bandstructure object).
            waveder_file: Path to WAVEDER file (to generate dielectric function calculator object,
            to compute orbital derivatives and transition dipole moments).
            bg: Experimental bandgap of the material, used to rigidly shift the DFT results to match
            this value. If None (default), the bandgap of the DFT calculation will be used.

        Returns:
            A TASGenerator object.
        """
        warnings.filterwarnings("ignore", category=UnknownPotcarWarning)
        warnings.filterwarnings(
            "ignore", message="No POTCAR file with matching TITEL fields"
        )
        vr = Vasprun(vasprun_file)
        if waveder_file:
            waveder = Waveder.from_binary(waveder_file)
            # check if LVEL was set to True in vasprun file:
            if not vr.incar.get("LVEL", False):
                lvel_error_message = (
                    "LVEL must be set to True in the INCAR for the VASP optics calculation to output the full "
                    "band-band orbital derivatives and thus allow PyTASer to parse the WAVEDER and compute oscillator "
                    "strengths. Please rerun the VASP calculation with LVEL=True (if you use the WAVECAR from the "
                    "previous calculation this should only require 1 or 2 electronic steps!"
                )
                if vr.incar.get("ISYM", 2) in [-1, 0]:
                    raise ValueError(lvel_error_message)
                raise ValueError(
                    f"ISYM must be set to 0 and {lvel_error_message}"
                )

            dfc = optics.DielectricFunctionCalculator.from_vasp_objects(
                vr, waveder
            )
        else:
            dfc = None

        bs = vr.get_band_structure()
        dos = vr.complete_dos

        if bg is not None:  # apply scissor shift
            if dfc is not None:
                eigs_shifted = dfc.eigs - dfc.efermi
                scissor = bg - bs.get_band_gap()["energy"]
                # shift dfc.eigs[eigs_shifted > 0] up by scissor:
                dfc.eigs[eigs_shifted > 0] += scissor

            bs, dos = set_bandgap(bs, dos, bg)

        return cls(
            bs,
            vr.actual_kpoints_weights,
            dos,
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
        processes=None,
    ):
        """
        Generates TAS spectra based on inputted occupancies, and a specified energy mesh. If the
        TASGenerator has not been generated from VASP outputs (and thus does not have a dfc
        attribute), then the output TAS is generated using the change in joint density of states
        (JDOS) under illumination, with no consideration of oscillator strengths.
        Otherwise, the output TAS is generated considering all contributions to the predicted TAS
        spectrum.

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
            light_occs: Optional input parameter for occupancies of material under light, otherwise
                automatically calculated based on input temperature (temp) and carrier concentration
                (conc) [dict]
            dark_occs: Optional input parameter for occupancies of material in dark, otherwise
                automatically calculated based on input temperature (temp) and carrier concentration
                (conc) [dict]
            processes: Number of processes to use for multiprocessing. If not set, defaults to one
                less than the number of CPUs available.

        Returns:
            TAS class containing the following inputs;
                - tas_total: overall deltaT TAS spectrum for a material under the
                    specified conditions.
                - jdos_diff_if: JDOS difference (from dark to light) across the energy mesh for a
                    specific band transition i (initial) -> f (final) [dict]
                - jdos_light_total: overall JDOS (pump-on) for a material under the
                    specified conditions
                - jdos_light_if: JDOS (pump-on) across the energy mesh for a specific band
                    transition i (initial) -> f (final) [dict]
                - jdos_dark_total: overall JDOS (pump-off) for a material under the
                    specified conditions
                - jdos_dark_if: JDOS (pump-off) across the energy mesh for a specific band
                    transition i (initial) -> f (final) [dict]
                - energy_mesh_ev: Energy mesh of spectra in eV, with an interval of 'step'.
                - bandgap: Bandgap of the system, in eV, rounded to 2 decimal points
                - temp: Temperature of the system, in K
                - conc: Carrier concentration of the system, in cm^-3
                - alpha_dark: Absorption coefficient of the material in the dark, in cm^-1 (only
                    calculated if the TASGenerator has been generated from VASP outputs)
                - alpha_light_dict: Dictionary of band-to-band absorption, stimulated emission
                    and summed contributions to the total overall absorption coefficient under
                    illumination, in cm^-1 (only calculated if the TASGenerator has been
                    generated from VASP outputs)
                - weighted_jdos_diff_if: JDOS difference (from dark to light) across the energy
                    mesh for a specific band transition i (initial) -> f (final), weighted by the
                    oscillator strength of the transition [dict]
                - weighted_jdos_light_if: JDOS (pump-on) across the energy mesh for a specific band
                    transition i (initial) -> f (final), weighted by the oscillator strength of
                    the transition [dict]
                - weighted_jdos_dark_if: JDOS (pump-off) across the energy mesh for a specific band
                    transition i (initial) -> f (final), weighted by the oscillator strength of
                    the transition [dict]
        """
        occs_light = light_occs
        occs_dark = dark_occs
        if occs_light is None:
            occs_light = self.band_occupancies(temp, conc, dark=False)
        if occs_dark is None:
            occs_dark = self.band_occupancies(temp, conc)

        bandgap = round(self.bs.get_band_gap()["energy"], 2)
        energy_mesh_ev = np.arange(energy_min, energy_max, step)
        jdos_light_if = {}
        jdos_dark_if = {}
        jdos_diff_if = {}
        weighted_jdos_light_if = {}
        weighted_jdos_dark_if = {}
        weighted_jdos_diff_if = {}
        tas_total = np.zeros(len(energy_mesh_ev))
        jdos_dark_total = np.zeros(len(energy_mesh_ev))
        jdos_light_total = np.zeros(len(energy_mesh_ev))

        if self.dfc is not None:
            egrid = np.arange(
                0, self.dfc.nedos * self.dfc.deltae, self.dfc.deltae
            )
            alpha_dark = np.zeros_like(egrid, dtype=np.complex128)
            alpha_light_dict = {
                key: np.zeros_like(egrid, dtype=np.complex128)
                for key in ["absorption", "emission", "both"]
            }

        for spin, spin_bands in self.bs.bands.items():
            if self.dfc is not None:
                # _nonzero_ tdm_array values should be the same in dark or light, but this calculation gets
                # skipped for band-band transitions where the occupancy factor is near-zero, which happens
                # for certain band transitions in the dark but not in light and vice versa, so get both
                # and take the nonzero values:
                alpha_dark_dict, tdm_array_dark = occ_dependent_alpha(
                    self.dfc,
                    occs_dark[spin],
                    sigma=gaussian_width,
                    cshift=cshift,
                    processes=processes,
                    energy_max=energy_max,
                )
                alpha_dark += alpha_dark_dict[
                    "both"
                ]  # stimulated emission should be zero in the dark
                (
                    calculated_alpha_light_dict,
                    tdm_array_light,
                ) = occ_dependent_alpha(
                    self.dfc,
                    occs_light[spin],
                    sigma=gaussian_width,
                    cshift=cshift,
                    processes=processes,
                    energy_max=energy_max,
                )
                tdm_array_light[
                    tdm_array_dark.nonzero()
                ] = 0  # zero out duplicate TDM array values
                tdm_array = tdm_array_dark + tdm_array_light

                for key, array in alpha_light_dict.items():
                    alpha_light_dict[key] += calculated_alpha_light_dict[key]

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
                        jdos_diff = jd_light - jd_dark
                        jdos_dark_total += jd_dark
                        jdos_light_total += jd_light
                        tas_total += jdos_diff

                        new_i = i - self.vb[spin]
                        new_f = f - self.vb[spin]

                        if self.bs.is_spin_polarized:
                            spin_str = "up" if spin == Spin.up else "down"
                            key = (new_i, new_f, spin_str)
                        else:
                            key = (new_i, new_f)

                        jdos_dark_if[key] = jd_dark
                        jdos_light_if[key] = jd_light
                        jdos_diff_if[key] = jdos_diff

                        if self.dfc is not None:
                            weighted_jd_light = jdos(
                                self.bs,
                                f,
                                i,
                                occs_light[spin],
                                energy_mesh_ev,
                                np.array(self.kpoint_weights)
                                * tdm_array[i, f, :],
                                gaussian_width,
                                spin=spin,
                            )
                            weighted_jd_dark = jdos(
                                self.bs,
                                f,
                                i,
                                occs_dark[spin],
                                energy_mesh_ev,
                                np.array(self.kpoint_weights)
                                * tdm_array[i, f, :],
                                gaussian_width,
                                spin=spin,
                            )
                            weighted_jdos_diff = (
                                weighted_jd_light - weighted_jd_dark
                            )

                            weighted_jdos_light_if[key] = weighted_jd_light
                            weighted_jdos_dark_if[key] = weighted_jd_dark
                            weighted_jdos_diff_if[key] = weighted_jdos_diff

        # need to interpolate alpha arrays onto JDOS energy mesh:
        if self.dfc is not None:
            alpha_dark = np.interp(energy_mesh_ev, egrid, alpha_dark)
            for key, array in alpha_light_dict.items():
                alpha_light_dict[key] = np.interp(energy_mesh_ev, egrid, array)

            tas_total = (
                alpha_light_dict["absorption"]
                - alpha_light_dict["emission"]
                - alpha_dark
            )

        return Tas(
            tas_total,
            jdos_diff_if,
            jdos_light_total,
            jdos_light_if,
            jdos_dark_total,
            jdos_dark_if,
            energy_mesh_ev,
            bandgap,
            temp,
            conc,
            alpha_dark if self.dfc is not None else None,
            alpha_light_dict if self.dfc is not None else None,
            weighted_jdos_light_if if self.dfc is not None else None,
            weighted_jdos_dark_if if self.dfc is not None else None,
            weighted_jdos_diff_if if self.dfc is not None else None,
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
            mpr = MPRester() if api_key is None else MPRester(api_key=api_key)
        mp_dos = mpr.get_dos_by_material_id(mpid)
        mp_bs = mpr.get_bandstructure_by_material_id(mpid, line_mode=False)
        if bg is not None:
            mp_bs, mp_dos = set_bandgap(mp_bs, mp_dos, bg)
        kweights = get_kpoint_weights(mp_bs)
        return TASGenerator(mp_bs, kweights, mp_dos)
