import generator
import plotter

fe2o3_temp = 298
fe2o3_conc = 10e21
bandgap = 2.2  # based on https://aip.scitation.org/doi/full/10.1063/1.2177426
transitions_fe2o3 = [
    (8, 11, "spin = down"),
    (9, 10, "spin = up"),
    (4, 14, "spin = up"),
    (6, 20, "spin = down"),
    (6, 28, "spin = up"),
]

# Fe2O3 is  spin-polarised, so we must define whether the relevant transitions we want to
# see are for spin-up or spin-down

fe2o3_generator = generator.TASGenerator.from_mpid(
    "cTHyvyRHqX34gOLMbM", "mp-565814", fe2o3_temp, fe2o3_conc, bandgap
)
fe2o3_tas = fe2o3_generator.generate_tas(energy_min=0.0, energy_max=10.0)
fe2o3_plotter = plotter.TASPlotter(fe2o3_tas, bandgap, "Fe2O3", fe2o3_temp, fe2o3_conc)
plot_fe2o3 = fe2o3_plotter.get_plot(
    transitions_fe2o3, xaxis="wavelength", xmin=350, xmax=1400, yaxis="TAS (deltaT)"
)
plot_fe2o3.show()
