#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 16:42:36 2023

@author: lucasverga
"""

import itertools
from tqdm import tqdm
import warnings
import numpy as np
from multiprocessing import Pool, cpu_count

from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.dos import FermiDos, f0
from pymatgen.ext.matproj import MPRester
from pymatgen.io.vasp import optics
from pymatgen.io.vasp.inputs import UnknownPotcarWarning
from pymatgen.io.vasp.outputs import Vasprun, Waveder

from pytaser.kpoints import get_kpoint_weights
from pytaser.tas import Tas
import pytaser.generator as generator


warnings.filterwarnings("ignore", category=RuntimeWarning)


#duplicate all entries


class DASGenerator:
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

    def __init__(self, bs, kpoint_weights, dos, 
                 bs_ref, kpoint_weights_ref,dos_ref,
                 dfc=None, dfc_ref=None):
        
        self.bs = bs
        self.kpoint_weights = kpoint_weights
        self.dos = FermiDos(dos)
        self.bg_centre = (bs.get_cbm()["energy"] + bs.get_vbm()["energy"]) / 2
        self.dfc = dfc
        
        
        self.bs_ref = bs_ref
        self.kpoint_weights_ref = kpoint_weights_ref
        self.dos_ref = FermiDos(dos_ref)
        self.bg_centre_ref = (bs_ref.get_cbm()["energy"] + bs_ref.get_vbm()["energy"]) / 2
        self.dfc_ref = dfc_ref

        # Lucas
        # if self.bs.is_metal():
        #     raise ValueError("System is metallic, cannot compute TAS")

        self.vb = generator.get_cbm_vbm_index(self.bs)[0]
        self.cb = generator.get_cbm_vbm_index(self.bs)[1]        
        self.vb_ref = generator.get_cbm_vbm_index(self.bs_ref)[0]
        self.cb_ref = generator.get_cbm_vbm_index(self.bs_ref)[1]
        
        

    @classmethod
    def from_vasp_outputs(cls, vasprun_file,vasprun_file_ref, waveder_file=None,waveder_file_ref=None):
        """Create a TASGenerator object from VASP output files."""
        warnings.filterwarnings("ignore", category=UnknownPotcarWarning)
        warnings.filterwarnings("ignore", message="No POTCAR file with matching TITEL fields")
        vr = Vasprun(vasprun_file)
        vr_ref = Vasprun(vasprun_file_ref)
        
        
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
                if vr.incar.get("ISYM", 2) not in [-1, 0]:
                    isym_error_message = "ISYM must be set to 0 and "
                    raise ValueError(isym_error_message + lvel_error_message)
                else:
                    raise ValueError(lvel_error_message)
            dfc = optics.DielectricFunctionCalculator.from_vasp_objects(
                vr, waveder
            )
        else:
            dfc = None 
            
        #Repeat the same code to get the data for the reference.     
        if waveder_file_ref:
            waveder_ref = Waveder.from_binary(waveder_file_ref)
            # check if LVEL was set to True in vasprun file:
            if not vr_ref.incar.get("LVEL", False):
                lvel_error_message = (
                    "LVEL must be set to True in the INCAR for the VASP optics calculation to output the full "
                    "band-band orbital derivatives and thus allow PyTASer to parse the WAVEDER and compute oscillator "
                    "strengths. Please rerun the VASP calculation with LVEL=True (if you use the WAVECAR from the "
                    "previous calculation this should only require 1 or 2 electronic steps!"
                )
                if vr_ref.incar.get("ISYM", 2) not in [-1, 0]:
                    isym_error_message = "ISYM must be set to 0 and "
                    raise ValueError(isym_error_message + lvel_error_message)
                else:
                    raise ValueError(lvel_error_message)
            dfc_ref = optics.DielectricFunctionCalculator.from_vasp_objects(
                vr_ref, waveder_ref
            )
        else:
            dfc_ref = None
            
            
        return cls(
            vr.get_band_structure(),
            vr.actual_kpoints_weights,
            vr.complete_dos,
            vr_ref.get_band_structure(),
            vr_ref.actual_kpoints_weights,
            vr_ref.complete_dos,
            dfc,
            dfc_ref
        )
 

    @classmethod
    def from_mpid(cls, mpid, mpid_ref, bg=None, bg_ref=None, api_key=None, mpr=None,mpr_ref=None):
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
            if api_key is None:
                mpr = MPRester()
            else:
                mpr = MPRester(api_key=api_key)
        mp_dos = mpr.get_dos_by_material_id(mpid)
        mp_bs = mpr.get_bandstructure_by_material_id(mpid, line_mode=False)
        if bg is not None:
            mp_bs, mp_dos = generator.set_bandgap(mp_bs, mp_dos, bg)
        kweights = get_kpoint_weights(mp_bs)
        
        
        
        if mpr_ref is None:
            if api_key is None:
                mpr_ref = MPRester()
            else:
                mpr_ref = MPRester(api_key=api_key)
        mp_dos_ref = mpr_ref.get_dos_by_material_id(mpid_ref)
        mp_bs_ref = mpr_ref.get_bandstructure_by_material_id(mpid_ref, line_mode=False)
        if bg_ref is not None:
            mp_bs_ref, mp_dos_ref = generator.set_bandgap(mp_bs_ref, mp_dos_ref, bg_ref)
        kweights_ref = get_kpoint_weights(mp_bs_ref)
        
                
        
        return DASGenerator(mp_bs, kweights, mp_dos,mp_bs_ref, kweights_ref, mp_dos_ref)





    def band_occupancies(self):
        """
        Gives band occupancies.

        Returns:
            Two dictionary of {Spin: occ} for all bands across all k-points.
            One for the modified and one for the reference system.
        """
        occs = {}  
        for spin, spin_bands in self.bs.bands.items():
            hole_mask = spin_bands < self.bg_centre
            elec_mask = spin_bands > self.bg_centre

            # fully occupied hole mask, completely empty electron mask
            spin_occs = np.zeros_like(spin_bands)
            spin_occs[hole_mask] = 1
            spin_occs[elec_mask] = 0

            occs[spin] = spin_occs
            
        occs_ref={}    
        for spin, spin_bands in self.bs_ref.bands.items():
            hole_mask = spin_bands < self.bg_centre_ref
            elec_mask = spin_bands > self.bg_centre_ref

            # fully occupied hole mask, completely empty electron mask
            spin_occs = np.zeros_like(spin_bands)
            spin_occs[hole_mask] = 1
            spin_occs[elec_mask] = 0

            occs_ref[spin] = spin_occs    
        
        return occs,occs_ref

    def generate_das(
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
        if light_occs is None:
            occs_light = self.band_occupancies()
        if dark_occs is None:
            occs_dark = self.band_occupancies()

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

        for spin, spin_bands in self.bs.bands.items():
            if self.dfc is not None:
                alpha_dark_dict, tdm_array = generator.occ_dependent_alpha(
                    self.dfc,
                    occs_dark[spin],
                    sigma=gaussian_width,
                    cshift=cshift,
                    processes=processes,
                    energy_max=energy_max,
                )
                alpha_dark += alpha_dark_dict[
                    "both"
                ]  # stimulated emission should be
                # zero in the dark
                alpha_light_dict = alpha_dark_dict                

        

            for i in range(len(spin_bands)):
                for f in range(len(spin_bands)):
                    if f > i:
                        jd_dark = generator.jdos(
                            self.bs,
                            f,
                            i,
                            occs_dark[spin],
                            energy_mesh_ev,
                            self.kpoint_weights,
                            gaussian_width,
                            spin=spin,
                        )
                        jdos_dark_total += jd_dark
                        
                        new_i = i - self.vb[spin]
                        new_f = f - self.vb[spin]

                        if self.bs.is_spin_polarized:
                            spin_str = "up" if spin == Spin.up else "down"
                            key = (new_i, new_f, spin_str)
                        else:
                            key = (new_i, new_f)

                        jdos_dark_if[key] = jd_dark
                        
                        
                        if self.dfc is not None:
                            weighted_jd_dark = generator.jdos(
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
 
                            weighted_jdos_dark_if[key] = weighted_jd_dark
 
        # need to interpolate alpha arrays onto JDOS energy mesh:
        if self.dfc is not None:
            alpha_dark = np.interp(energy_mesh_ev, egrid, alpha_dark)

            tas_total = (
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

