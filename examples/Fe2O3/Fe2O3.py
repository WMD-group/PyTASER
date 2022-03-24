from pytaser.generator import TASGenerator
from pytaser.plotter import TASPlotter

api_key = None  # Use your API key here
temp = 298
conc = 10e21
bandgap = 2.2  # based on https://aip.scitation.org/doi/full/10.1063/1.2177426

generator = TASGenerator.from_mpid("mp-565814", bandgap, api_key=api_key)
tas = generator.generate_tas(temp, conc)
plotter = TASPlotter(
    tas, bandgap_ev=bandgap, material_name="Fe2O3", temp=temp, conc=conc
)
plt = plotter.get_plot(xaxis="wavelength", xmin=350, xmax=1400, yaxis="TAS (deltaT)")
plt.show()
