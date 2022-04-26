from emopt import grid

class Waveguide():
    def __init__(self, parameters=None, **kwargs):
        self.w = 1.0
        self.h = 0.3
        self.h_tf = 0.6
        self.alpha = 65
        self.h_ins = 2.0
        self.cladding = False

        if parameters:
            ks = parameters.keys()
        else:
            ks = []
        
        for attr in self.__dict__.keys():
            if attr in ks:
                setattr(self,attr,parameters[attr])

    def __repr__(self):
        return f'Waveguide on {self.h_tf*1000:.0f} nm thin film thickness of {self.w*1000:.0f} nm width and {self.h*1000:.0f} nm height'



p = {'w':0.6,'h': 0.3, 'alpha': 60}
wg = Waveguide(p)
print(wg)