from src.mode import WaveguideModes
from src.waveguide import Waveguide
from emopt.misc import NOT_PARALLEL
import numpy as np
import matplotlib.pyplot as plt
# plt.style.use('https://raw.githubusercontent.com/maederan201/MPLStyles/master/mystyle-std.mplstyle')

N = 200
n = 1
wavelength = 1.55  # um
fixed_params = {'h': 0.25, 'w': 0.8, 'h_tf': 0.4, 'alpha': 65}

sweep_parameter_name = 'resolution'
sweep_parameter_label = 'Resoltuion [μm]'
sweep_parameter_values = np.linspace(0.1,0.002,N)


result_array = np.zeros((N,n+1),dtype=np.complex128)

for i,sweep_parameter in enumerate(sweep_parameter_values):
    result_array[i,0] = sweep_parameter
    wg = Waveguide({**fixed_params, **{sweep_parameter_name: sweep_parameter}})
    modes = WaveguideModes(wg,n_modes=n,wavelength=wavelength)
    modes.solve()
    for j in range(n):
        result_array[i,j+1] = modes.get_neff(j)

if NOT_PARALLEL:
    fig, ax = plt.subplots()
    for j in range(n):
        ax.plot(result_array[:,0],result_array[:,j+1],marker='.')
    ax.set_xlabel(sweep_parameter_label)
    ax.set_ylabel('Effective Index')
else:
    np.savez_compressed('result_array.npz', a=result_array)
