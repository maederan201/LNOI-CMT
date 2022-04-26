import emopt
from emopt.misc import NOT_PARALLEL
import numpy as np
import copy
import matplotlib.pyplot as plt


class WaveguideModes():
    def __init__(self, waveguide, wavelength=1.55, n_modes=1):
        self.wg = waveguide
        self.n_modes = n_modes
        self.wl = wavelength

        self.solver = emopt.modes.ModeFullVector(self.wl, self.wg.eps, self.wg.mu, self.wg.domain, n0=self.wg.n_wg[0], neigs=self.n_modes)


    def solve(self):
        local_solver = self.solver
        local_solver.build()
        local_solver.solve()

        self.solver = local_solver

    def get_field(self, i, component):
        return self.solver.get_field(i, component)

    def get_field_interp(self, i, component):
        return self.solver.get_field_interp(i, component)

    def get_neff(self, i):
        return self.solver.neff[i]

    def visualize(self, index=0):
        Ex = self.get_field_interp(index, 'Ex')[0]
        Ey = self.get_field_interp(index, 'Ey')[0]
        Ez = self.get_field_interp(index, 'Ez')[0]
        Hx = self.get_field_interp(index, 'Hx')[0]
        Hy = self.get_field_interp(index, 'Hy')[0]
        Hz = self.get_field_interp(index, 'Hz')[0]

        normE = np.sqrt(Ex**2 + Ey**2 + Ez**2)
        normH = np.sqrt(Hx**2 + Hy**2 + Hz**2)

        x = self.wg.domain.x
        y = self.wg.domain.y

        Dx = self.wg.domain.xspan
        Dy = self.wg.domain.yspan

        nxx = np.sqrt(self.wg.eps[0].get_values_in(self.wg.domain))
        colormap = 'inferno'


        fig, ax = plt.subplots(ncols=3,nrows=3,figsize=(11,11))

        pnormE = ax[0,0].imshow(np.absolute(normE),extent=[-Dx/2,Dx/2,-Dy/2.,Dy/2.],origin='lower',cmap=colormap)
        fig.colorbar(pnormE, ax=ax[0,0])
        ax[0,0].set_title(r'$\vert E \vert$')


        pnormH = ax[0,1].imshow(np.absolute(normH),extent=[-Dx/2,Dx/2,-Dy/2.,Dy/2.],origin='lower',cmap=colormap)
        fig.colorbar(pnormH, ax=ax[0,1])
        ax[0,1].set_title(r'$\vert H \vert$')

        pnxx = ax[0,2].imshow(nxx.real,interpolation='none',extent=[-Dx/2,Dx/2,-Dy/2.,Dy/2.],origin='lower',cmap=colormap,vmax=2.5,vmin=1.0)
        fig.colorbar(pnxx, ax=ax[0,2])
        ax[0,2].set_title(r'$n_{xx}$')

        pEx = ax[1,0].imshow(np.absolute(Ex),extent=[-Dx/2,Dx/2,-Dy/2.,Dy/2.],origin='lower',cmap=colormap)
        fig.colorbar(pEx, ax=ax[1,0])
        ax[1,0].set_title(r'$\vert E_x \vert$')

        pEy = ax[1,1].imshow(np.absolute(Ey),extent=[-Dx/2,Dx/2,-Dy/2.,Dy/2.],origin='lower',cmap=colormap)
        fig.colorbar(pEy, ax=ax[1,1])
        ax[1,1].set_title(r'$\vert E_y \vert$')

        pEz = ax[1,2].imshow(np.absolute(Ez),extent=[-Dx/2,Dx/2,-Dy/2.,Dy/2.],origin='lower',cmap=colormap)
        fig.colorbar(pEz, ax=ax[1,2])
        ax[1,2].set_title(r'$\vert E_z \vert$')

        pHx = ax[2,0].imshow(np.absolute(Hx),extent=[-Dx/2,Dx/2,-Dy/2.,Dy/2.],origin='lower',cmap=colormap)
        fig.colorbar(pHx, ax=ax[2,0])
        ax[2,0].set_title(r'$\vert H_x \vert$')

        pHy = ax[2,1].imshow(np.absolute(Hy),extent=[-Dx/2,Dx/2,-Dy/2.,Dy/2.],origin='lower',cmap=colormap)
        fig.colorbar(pHy, ax=ax[2,1])
        ax[2,1].set_title(r'$\vert H_y \vert$')

        pHz = ax[2,2].imshow(np.absolute(Hz),extent=[-Dx/2,Dx/2,-Dy/2.,Dy/2.],origin='lower',cmap=colormap)
        fig.colorbar(pHz, ax=ax[2,2])
        ax[2,2].set_title(r'$\vert H_z \vert$')

        I_TE = np.sum(Ex*np.conj(Hy)/2.)
        I_TM = np.sum(-Ey*np.conj(Hx)/2.)
        qTE = I_TE/(I_TE + I_TM)
        qTM = I_TM/(I_TE + I_TM)

        neff = self.get_neff(index)

        fig.suptitle(f'Mode {index:g}, TE:TM Ratio: {qTE.real:.2f}:{qTM.real:.2f}, neff = {neff.real:.5f}',fontsize=14)

        fig.show()
        if NOT_PARALLEL:
            pass
        else:
            print('Warning: Mode not visualized due to being in parallel computing mode.')

    def overlap_integral(self, mode2, indices=(0,0), gap=0, weight=None):
        Ex1 = self.get_field_interp(indices[0], 'Ex')
        Ey1 = self.get_field_interp(indices[0], 'Ey')
        Ez1 = self.get_field_interp(indices[0], 'Ez')
        pass
