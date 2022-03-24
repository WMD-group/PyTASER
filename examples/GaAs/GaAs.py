from pytaser.generator import TASGenerator
from pytaser.plotter import TASPlotter

api_key = None  # Use your API key here
temp = 298
conc = 10e21
bandgap = 1.46

generator = TASGenerator.from_mpid("mp-2534", bandgap, api_key=api_key)
tas = generator.generate_tas(temp, conc)
plotter = TASPlotter(
    tas, bandgap_ev=bandgap, material_name="GaAs", temp=temp, conc=conc
)
plt = plotter.get_plot(xaxis="wavelength", xmin=350, xmax=1000, yaxis="TAS (deltaT)")
plt.show()
