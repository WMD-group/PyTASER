from pytaser import generator
from pytaser import plotter

gaas_temp = 298
gaas_conc = 10e21
bandgap = 1.46
transitions_gaas = [(8, 11), (9, 10), (4, 14), (6, 20), (6, 8)]

GaAs_generator = generator.TASGenerator.from_mpid(
    "cTHyvyRHqX34gOLMbM", "mp-2534", gaas_temp, gaas_conc, bandgap
)
GaAs_tas = GaAs_generator.generate_tas(energy_min=0.0, energy_max=10.0)
GaAs_plotter = plotter.TASPlotter(GaAs_tas, bandgap, "GaAs", gaas_temp, gaas_conc)
plot_gaas = GaAs_plotter.get_plot(
    transitions_gaas, xaxis="wavelength", xmin=350, xmax=1400, yaxis="TAS (deltaT)"
)
plot_gaas.show()
