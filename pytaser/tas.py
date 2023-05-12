import ast
from monty.json import MontyDecoder

class Tas:
    """
    A container class for the data from TASgenerator.

    Args:
        tas_total: overall TAS spectrum for a material under the specified conditions
        jdos_diff_if: JDOS difference (from dark to light) across the energy mesh for a
            specific band transition i (initial) -> f (final) [dict]
        jdos_light_total: overall JDOS (pump-on) for a material under the specified
            conditions
        jdos_light_if: JDOS (pump-on) across the energy mesh for a specific band
            transition i (initial) -> f (final) [dict]
        jdos_dark_total: overall JDOS (pump-off) for a material under the specified
            conditions
        jdos_dark_if: JDOS (pump-off) across the energy mesh for a specific band
            transition i (initial) -> f (final) [dict]
        energy_mesh_ev: Energy mesh of spectra in eV, with an interval of 'step'.
        bandgap: Bandgap of the system in electronvolts (eV).
        temp: Temperature (K) of material we wish to investigate (affects the FD distribution)
        conc: Carrier concentration (cm^-3) of holes and electrons (both are equivalent).
            Inversely proportional to pump-probe time delay.
        alpha_dark: Absorption coefficient of the material in the dark, in cm^-1.
        alpha_light_dict: Dictionary of band-to-band absorption, stimulated emission and summed
            contributions to the total overall absorption coefficient under illumination, in cm^-1.
        weighted_jdos_light_if: JDOS weighted by the transition dipole matrix (TDM) (pump-on)
            across the energy mesh for a specific band transition i (initial) -> f (final) [dict]
        weighted_jdos_dark_if: JDOS weighted by the transition dipole matrix (TDM) (pump-off)
            across the energy mesh for a specific band transition i (initial) -> f (final) [dict]
        weighted_jdos_diff_if: Difference in JDOS weighted by the transition dipole matrix
            (TDM) from dark to illumination across the energy mesh for a specific band transition
            i (initial) -> f (final) [dict]
    """

    def __init__(
        self,
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
        alpha_dark=None,
        alpha_light_dict=None,
        weighted_jdos_light_if=None,
        weighted_jdos_dark_if=None,
        weighted_jdos_diff_if=None,
    ):
        self.tas_total = tas_total
        self.jdos_diff_if = jdos_diff_if
        self.jdos_light_total = jdos_light_total
        self.jdos_light_if = jdos_light_if
        self.jdos_dark_total = jdos_dark_total
        self.jdos_dark_if = jdos_dark_if
        self.energy_mesh_ev = energy_mesh_ev
        self.bandgap = bandgap
        self.temp = temp
        self.conc = conc
        self.alpha_dark = alpha_dark
        self.alpha_light_dict = alpha_light_dict
        self.weighted_jdos_light_if = weighted_jdos_light_if
        self.weighted_jdos_dark_if = weighted_jdos_dark_if
        self.weighted_jdos_diff_if = weighted_jdos_diff_if

    def as_dict(self):
        """
        JSON-serializable dict representation of Tas
        """
        json_dict = {
            "@module": type(self).__module__,
            "@class": type(self).__name__,
            "tas_total": self.tas_total,
            "jdos_diff_if": self.jdos_diff_if,
            "jdos_light_total": self.jdos_light_total,
            "jdos_light_if": self.jdos_light_if,
            "jdos_dark_total": self.jdos_dark_total,
            "jdos_dark_if": self.jdos_dark_if,
            "energy_mesh_ev": self.energy_mesh_ev,
            "bandgap": self.bandgap,
            "temp": self.temp,
            "conc": self.conc,
            "alpha_dark": self.alpha_dark,
            "alpha_light_dict": self.alpha_light_dict,
            "weighted_jdos_light_if": self.weighted_jdos_light_if,
            "weighted_jdos_dark_if": self.weighted_jdos_dark_if,
            "weighted_jdos_diff_if": self.weighted_jdos_diff_if,
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
        Reconstructs Tas object from a dict representation of Tas created using
        Tas.as_dict().

        Args:
            d (dict): dict representation of Tas

        Returns:
            Tas object
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
