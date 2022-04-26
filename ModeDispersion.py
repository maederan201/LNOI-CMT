from src.mode import WaveguideModes
from src.waveguide import Waveguide
import numpy as np

def run():
    print('Ran!')

sweep_parameter_name = 'w'
sweep_parameter_values = np.linspace(0.6,1.2,10)

fixed_params = {'h': 0.3, 'alpha': 65}

for sweep_parameter in sweep_parameter_values:

    wg = Waveguide(fixed_params | {sweep_parameter_name: sweep_parameter})
    print(wg)

# modes = WaveguideModes(wg)
# modes.solve()