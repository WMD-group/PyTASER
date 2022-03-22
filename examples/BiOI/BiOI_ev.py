from pytaser import generator, plotter

api_key = None  # Use your API key here
bioi_temp = 298
bioi_conc = 10e21
bandgap = 1.94

bioi_generator = generator.TASGenerator.from_mpid(
    "mp-22987", bioi_temp, bioi_conc, bandgap, api_key=api_key
)
bioi_tas = bioi_generator.generate_tas(energy_min=0.0, energy_max=10.0)
bioi_plotter = plotter.TASPlotter(
    bioi_tas, bandgap, material_name="BiOI", temp=bioi_temp, conc=bioi_conc
)
plot_bioi = bioi_plotter.get_plot(
    xaxis="electronvolts", transition_cutoff=0.03, xmin=1.7, xmax=2.8, yaxis="TAS (deltaT)"
)
plot_bioi.show()
