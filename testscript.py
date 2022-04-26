from src.waveguide import Waveguide
from src.mode import WaveguideModes


params = {'w':0.8,'h': 0.3, 'alpha': 60}
wg = Waveguide(params)
print(wg)

modes = WaveguideModes(wg,n_modes=5)
modes.solve()

modes.visualize(1)
