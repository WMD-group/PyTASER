#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 14 17:52:28 2023

@author: lucasverga
"""

import ast
from monty.json import MontyDecoder

class Das:
    """
    A container class for the data from Dasgenerator.

    Args:
        das_total: overall DAS spectrum between newSystem and reference system.
        jdos_diff_if: JDOS difference (from reference to newsystem) across the energy mesh for a
            specific band transition i (initial) -> f (final) [dict]
        jdos_newSys_total: overall JDOS for the new system under the specified
            conditions
        jdos_newSys_if: JDOS for the new system across the energy mesh for a specific band
            transition i (initial) -> f (final) [dict]
        jdos_ref_total: overall JDOS for the reference system under the specified
            conditions
        jdos_ref_if: JDOS for the reference system across the energy mesh for a specific band
            transition i (initial) -> f (final) [dict]
        energy_mesh_ev: Energy mesh of spectra in eV, with an interval of 'step'.
        bandgap_newSys: Bandgap of the new system in electronvolts (eV).
        bandgap_ref: Bandgap of the reference system in electronvolts (eV).
        temp: Temperature (K) of material we wish to investigate (affects the FD distribution)
        alpha_newSys: Absorption coefficient of the new system in the dark, in cm^-1.
        alpha_ref: Absorption coefficient of the reference system in the dark, in cm^-1.
        weighted_jdos_newSys_if: JDOS weighted by the transition dipole matrix (TDM) for the new system 
            across the energy mesh for a specific band transition i (initial) -> f (final) [dict]
        weighted_jdos_ref_if: JDOS weighted by the transition dipole matrix (TDM) for the reference system
            across the energy mesh for a specific band transition i (initial) -> f (final) [dict]
        weighted_jdos_diff_if: Difference in JDOS weighted by the transition dipole matrix
            (TDM) from reference to new system across the energy mesh for a specific band transition
            i (initial) -> f (final) [dict]
    """


    def __init__(
        self,
        das_total,
        jdos_newSys_total,
        jdos_newSys_if,
        jdos_ref_total,
        jdos_ref_if,
        energy_mesh_ev,
        bandgap_newSys,
        bandgap_ref,
        temp,
        alpha_newSys=None,
        alpha_ref=None,
        weighted_jdos_newSys_if=None,
        weighted_jdos_ref_if=None,
    ):
        self.das_total=das_total
        self.jdos_newSys_total=jdos_newSys_total
        self.jdos_newSys_if=jdos_newSys_if
        self.jdos_ref_total=jdos_ref_total
        self.jdos_ref_if=jdos_ref_if
        self.energy_mesh_ev=energy_mesh_ev
        self.bandgap_newSys=bandgap_newSys
        self.bandgap_ref=bandgap_ref
        self.temp=temp
        self.alpha_newSys=alpha_newSys 
        self.alpha_ref=alpha_ref
        self.weighted_jdos_newSys_if=weighted_jdos_newSys_if 
        self.weighted_jdos_ref_if=weighted_jdos_ref_if

    def as_dict(self):
        """
        JSON-serializable dict representation of Das
        """
        json_dict = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "das_total": self.das_total,
            "jdos_newSys_total": self.jdos_newSys_total,
            "jdos_newSys_if": self.jdos_newSys_if,
            "jdos_ref_total": self.jdos_ref_total,
            "jdos_ref_if": self.jdos_ref_if,
            "energy_mesh_ev": self.energy_mesh_ev,
            "bandgap_newSys": self.bandgap_newSys,
            "bandgap_ref": self.bandgap_ref,
            "temp": self.temp,
            "alpha_newSys": self.alpha_newSys,
            "alpha_ref": self.alpha_ref,
            "weighted_jdos_newSys_if": self.weighted_jdos_newSys_if,
            "weighted_jdos_ref_if": self.weighted_jdos_ref_if,
        }
        for key, value in json_dict.items():
            if isinstance(value, dict):
                json_dict[key] = {
                    str(k): v for k, v in value.items()
                }  # decomp dicts, can't have tuples as keys
        return json_dict

    @classmethod
    def from_dict(cls, d):
        """
        Reconstructs Das object from a dict representation of Das created using
        Das.as_dict().

        Args:
            d (dict): dict representation of Das

        Returns:
            Das object
        """

        def convert_to_tuple(subdict):
            if isinstance(subdict, dict) and "@module" not in subdict:
                return {
                    ast.literal_eval(k) if "(" in k and ")" in k else k: v
                    for k, v in subdict.items()
                }
            return subdict

        def decode_dict(subdict):
            if isinstance(subdict, dict) and "@module" in subdict:
                return MontyDecoder().process_decoded(subdict)

            elif isinstance(subdict, dict) and "@module" not in subdict:
                for k, v in subdict.items():
                    if isinstance(v, dict) and "@module" in v:
                        subdict[k] = MontyDecoder().process_decoded(v)

            return subdict

        d_dec = {k: convert_to_tuple(v) for k, v in d.items()}
        d_decoded = {k: decode_dict(v) for k, v in d_dec.items()}
        return cls(**d_decoded)