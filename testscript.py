from src.waveguide import Waveguide
from src.mode import WaveguideModes


params = {'w':0.8,'h': 0.3, 'alpha': 60}
wg = Waveguide(params)
print(wg)

modes = WaveguideModes(wg,n_modes=1)
modes.solve()

modes.visualize()

wg.resolution = 0.09

modes2 = WaveguideModes(wg,n_modes=1)
modes2.solve()

# modes.visualize(1)
modes.coupling_coefficients(modes2)