from pytaser.generator import TASGenerator
from pytaser.plotter import TASPlotter

api_key = None  # Use your API key here
temp = 298
conc = 10e21
bandgap = 3.1

generator = TASGenerator.from_mpid("mp-390", bandgap, api_key=api_key)
tas = generator.generate_tas(temp, conc)
plotter = TASPlotter(
    tas, bandgap_ev=bandgap, material_name="TiO2 (anatase)", temp=temp, conc=conc
)
plt = plotter.get_plot(xaxis="wavelength", xmin=350, xmax=1400, yaxis="TAS (deltaT)")
plt.show()
