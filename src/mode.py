import emopt
import numpy as np

class WaveguideModes():
    def __init__(self, waveguide, wavelength=1.55, n_modes=1):
        self.wg = waveguide
        self.n_modes = n_modes
        self.wl = wavelength
    @property
    def solver(self):
        return emopt.modes.ModeFullVector(self.wl, self.wg.eps, self.wg.mu, self.wg.domain, n0=np.sqrt(self.wg.n_wg[0]), neigs=self.n_modes)

    def solve(self):
        self.solver.build()
        self.solver.solve()

    def get_field(self, i, component):
        return self.solver.get_field(i, component)

    def get_field_interp(self, i, component):
        return self.solver.get_field_interp(i, component)

    def get_neff(self, i):
        return self.solver.neff[i]