import emopt
import numpy as np

class Waveguide():
    def __init__(self, parameters=None, **kwargs):
        self.w = 1.0
        self.h = 0.3
        self.h_tf = 0.6
        self.alpha = 65
        self.h_ins = 2.0
        self.n_clad = 1.0 #1.444023632787416
        self.n_wg = (2.1375596497855565,2.211111008653574,2.211111008653574)
        self.n_bg = 1.0
        self.n_ins = 1.444023632787416
        self.h_clad = None

        self._sim_region = (6,6)  # W, H
        self.resolution = 0.1

        if not self.h_clad:
            self.h_clad = self._sim_region[0]/2.

        if parameters:
            ks = parameters.keys()
        else:
            ks = []

        for attr in self.__dict__.keys():
            if attr in ks:
                setattr(self,attr,parameters[attr])

    def __repr__(self):
        return f'Waveguide on {self.h_tf*1000:.0f} nm thin film thickness of {self.w*1000:.0f} nm width and {self.h*1000:.0f} nm height'

    @property
    def _vertices(self):
        W = self._sim_region[0]
        w_wg_bot = self.w + 2*self.h/np.tan(self.alpha*np.pi/180.)
        return np.array([
            [-W/2, self.h_tf/2 - self.h],
            [-w_wg_bot/2., self.h_tf/2 - self.h],
            [-self.w/2., self.h_tf/2.],
            [self.w/2., self.h_tf/2.],
            [w_wg_bot/2., self.h_tf/2. - self.h],
            [W/2., self.h_tf/2. - self.h],
            [W/2., - self.h_tf/2.],
            [-W/2., - self.h_tf/2.],
            [-W/2, self.h_tf/2 - self.h]])

    @property
    def eps(self):
        W = self._sim_region[0]
        H = self._sim_region[1]

        dx = self.resolution
        dy = dx

        lnoi = emopt.grid.Polygon(self._vertices[:,0],self._vertices[:,1])
        bottom_cladding = emopt.grid.Rectangle(0, -self.h_tf/2. - self.h_clad/2. , W , self.h_clad)
        eps_bg = emopt.grid.Rectangle(0, self.h_tf/2. + (H-self.h_clad-self.h_tf)/2., W, (H-self.h_clad-self.h_tf))

        lnoi.layer = 1; lnoi.material_value = self.n_wg[0]**2
        bottom_cladding.layer = 2; bottom_cladding.material_value = self.n_ins
        eps_bg.layer = 3; eps_bg.material_value = self.n_clad

        eps_x = emopt.grid.StructuredMaterial2D(W,H,dx,dy)
        eps_x.add_primitives([lnoi, bottom_cladding, eps_bg])

        lnoi = emopt.grid.Polygon(self._vertices[:,0],self._vertices[:,1])
        bottom_cladding = emopt.grid.Rectangle(0, -self.h_tf/2. - self.h_clad/2. , W , self.h_clad)
        eps_bg = emopt.grid.Rectangle(0, self.h_tf/2. + (H-self.h_clad-self.h_tf)/2., W, (H-self.h_clad-self.h_tf))

        lnoi.layer = 1; lnoi.material_value = self.n_wg[1]**2
        bottom_cladding.layer = 2; bottom_cladding.material_value = self.n_ins
        eps_bg.layer = 3; eps_bg.material_value = 1.0

        eps_y = emopt.grid.StructuredMaterial2D(W,H,dx,dy)
        eps_y.add_primitives([lnoi, bottom_cladding, eps_bg])

        lnoi = emopt.grid.Polygon(self._vertices[:,0],self._vertices[:,1])
        bottom_cladding = emopt.grid.Rectangle(0, -self.h_tf/2. - self.h_clad/2. , W , self.h_clad)
        eps_bg = emopt.grid.Rectangle(0, self.h_tf/2. + (H-self.h_clad-self.h_tf)/2., W, (H-self.h_clad-self.h_tf))

        lnoi.layer = 1; lnoi.material_value = self.n_wg[2]**2
        bottom_cladding.layer = 2; bottom_cladding.material_value = self.n_ins
        eps_bg.layer = 3; eps_bg.material_value = 1.0

        eps_z = emopt.grid.StructuredMaterial2D(W,H,dx,dy)
        eps_z.add_primitives([lnoi, bottom_cladding, eps_bg])

        return [eps_x,eps_y,eps_z]

    @property
    def mu(self):
        return emopt.grid.ConstantMaterial2D(1.0)

    @property
    def domain(self):
        W = self._sim_region[0]
        H = self._sim_region[1]
        return emopt.misc.DomainCoordinates(-W/2., W/2., -self.h_tf/2. - self.h_clad, self.h_tf/2. + (H - self.h_tf - self.h_clad), 0, 0, self.resolution, self.resolution, 1.0)
