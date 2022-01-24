from pytaser import generator
from pytaser import plotter

tio2r_temp = 298
tio2r_conc = 10e21
bandgap = 3.1
transitions_tio2r = [(8, 11), (9, 10), (4, 14), (6, 20), (6, 28)]

tio2r_generator = generator.TASGenerator.from_mpid(
    "cTHyvyRHqX34gOLMbM", "mp-2657", tio2r_temp, tio2r_conc, bandgap
)
tio2r_tas = tio2r_generator.generate_tas(energy_min=0.0, energy_max=10.0)
tio2r_plotter = plotter.TASPlotter(
    tio2r_tas, bandgap, "TiO2 (rutile)", tio2r_temp, tio2r_conc
)
plot_tio2r = tio2r_plotter.get_plot(
    transitions_tio2r, xaxis="wavelength", xmin=350, xmax=1400, yaxis="TAS (deltaT)"
)
plot_tio2r.show()
