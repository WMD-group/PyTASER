from pytaser import generator, plotter

api_key = None  # Use your API key here
gaas_temp = 298
gaas_conc = 10e21
bandgap = 1.46

GaAs_generator = generator.TASGenerator.from_mpid(
    "mp-2534", gaas_temp, gaas_conc, bandgap, api_key=api_key
)
GaAs_tas = GaAs_generator.generate_tas(energy_min=0.0, energy_max=10.0)
GaAs_plotter = plotter.TASPlotter(
    GaAs_tas, bandgap, material_name="GaAs", temp=gaas_temp, conc=gaas_conc
)
plot_gaas = GaAs_plotter.get_plot(
    xaxis="wavelength", transition_cutoff=0.03, xmin=350, xmax=1000, yaxis="TAS (deltaT)"
)
plot_gaas.show()
