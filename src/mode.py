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

    def _integrate_field(self,field, x, y):
        pass

    def coupling_coefficients(self, mode2, indices=(0,0), gap=0.0):
        Ex1 = self.get_field_interp(indices[0], 'Ex')[0]
        Ey1 = self.get_field_interp(indices[0], 'Ey')[0]
        Ez1 = self.get_field_interp(indices[0], 'Ez')[0]
        Hx1 = self.get_field_interp(indices[0], 'Hx')[0]
        Hy1 = self.get_field_interp(indices[0], 'Hy')[0]
        Hz1 = self.get_field_interp(indices[0], 'Hz')[0]

        Ex2 = mode2.get_field_interp(indices[1], 'Ex')[0]
        Ey2 = mode2.get_field_interp(indices[1], 'Ey')[0]
        Ez2 = mode2.get_field_interp(indices[1], 'Ez')[0]
        Hx2 = mode2.get_field_interp(indices[1], 'Hx')[0]
        Hy2 = mode2.get_field_interp(indices[1], 'Hy')[0]
        Hz2 = mode2.get_field_interp(indices[1], 'Hz')[0]

        nxx1 = np.sqrt(mode2.wg.eps[0].get_values_in(mode2.wg.domain))
        nyy1 = np.sqrt(mode2.wg.eps[1].get_values_in(mode2.wg.domain))
        nzz1 = np.sqrt(mode2.wg.eps[2].get_values_in(mode2.wg.domain))
        nxx2 = np.sqrt(mode2.wg.eps[0].get_values_in(mode2.wg.domain))
        nyy2 = np.sqrt(mode2.wg.eps[1].get_values_in(mode2.wg.domain))
        nzz2 = np.sqrt(mode2.wg.eps[2].get_values_in(mode2.wg.domain))

        xspan1 = 1

        E1 = np.array([Ex1,Ey1,Ez1],dtype='complex128')
        E2 = np.array([Ex2,Ey2,Ez2],dtype='complex128')
        H1 = np.array([Hx1,Hy1,Hz1],dtype='complex128')
        H2 = np.array([Hx2,Hy2,Hz2],dtype='complex128')

        k_pq = np.zeros((2,2),dtype='complex128')
        c_pq = np.zeros((2,2),dtype='complex128')
        chi_p = np.zeros((2,),dtype='complex128')

        E_i = [E1, E2]
        H_i = [H1, H2]

        for p in range(2):
            for q in range(2):
                # E_p* . E_q
                scalar_product = (
                    np.multiply(np.conj(E_i[p][0,:,:]),E_i[q][0,:,:])
                    + np.multiply(np.conj(E_i[p][1,:,:]),E_i[q][1,:,:])
                    + np.multiply(np.conj(E_i[p][2,:,:]),E_i[q][2,:,:])
                )
                # u_z . (E_p* x H_p + E_p x H_p*)
                cross_product1 = (
                    np.multiply(np.conj(E_i[p][0,:,:]),H_i[p][1,:,:])
                    - np.multiply(np.conj(E_i[p][1,:,:]),H_i[p][0,:,:])
                    + np.multiply(E_i[p][0,:,:],np.conj(H_i[p][1,:,:]))
                    - np.multiply(E_i[p][1,:,:],np.conj(H_i[p][0,:,:]))
                )
                # u_z . (E_p* x H_q + E_q x H_p*)
                cross_product2 = (
                    np.multiply(np.conj(E_i[p][0,:,:]),H_i[q][1,:,:])
                    - np.multiply(np.conj(E_i[p][1,:,:]),H_i[q][0,:,:])
                    + np.multiply(E_i[q][0,:,:],np.conj(H_i[p][1,:,:]))
                    - np.multiply(E_i[q][1,:,:],np.conj(H_i[p][0,:,:]))
                )
