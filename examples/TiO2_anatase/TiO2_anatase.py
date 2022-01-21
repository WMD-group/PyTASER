import generator
import plotter

tio2a_temp = 298
tio2a_conc = 10e21
bandgap = 3.1
transitions_tio2a = [(8, 11), (9, 10), (4, 14), (6, 20), (6, 8)]

tio2a_generator = generator.TASGenerator.from_mpid(
    "cTHyvyRHqX34gOLMbM", "mp-390", tio2a_temp, tio2a_conc, bandgap
)
tio2a_tas = tio2a_generator.generate_tas(energy_min=0.0, energy_max=10.0)
tio2a_plotter = plotter.TASPlotter(
    tio2a_tas, bandgap, "TiO2 (anatase)", tio2a_temp, tio2a_conc
)
plot_tio2a = tio2a_plotter.get_plot(
    transitions_tio2a, xaxis="wavelength", xmin=350, xmax=1400, yaxis="TAS (deltaT)"
)
plot_tio2a.show()
