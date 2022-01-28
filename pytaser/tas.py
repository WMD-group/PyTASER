class Tas:
    """
    A container class for all the generated data from TASgenerator.

    Args:
        total_tas: overall TAS spectrum for a material under the specified conditions
        tas_decomp: TAS spectrum across the energy mesh for a specific band transition i
            (initial) -> f (final) [dict]
        jdos_light_total: overall JDOS (pump-on) for a material under the specified
            conditions
        jdos_light_decomp: JDOS (pump-on) across the energy mesh for a specific band
            transition i (initial) -> f (final) [dict]
        jdos_dark_total: overall JDOS (pump-off) for a material under the specified
            conditions
        jdos_dark_decomp: JDOS (pump-off) across the energy mesh for a specific band
            transition i (initial) -> f (final) [dict]
        energy_mesh_ev: Energy mesh of spectra in eV, with an interval of 'step'.
    """

    def __init__(
        self,
        total_tas,
        tas_decomp,
        jdos_light_total,
        jdos_light_decomp,
        jdos_dark_total,
        jdos_dark_decomp,
        energy_mesh_ev,
    ):
        self.total_tas = total_tas
        self.tas_decomp = tas_decomp
        self.jdos_light_tot = jdos_light_total
        self.jdos_light_decomp = jdos_light_decomp
        self.jdos_dark_tot = jdos_dark_total
        self.jdos_dark_decomp = jdos_dark_decomp
        self.energy_mesh_ev = energy_mesh_ev
