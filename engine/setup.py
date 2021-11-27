import numpy as np

class MFLHYSimulationSetup():
    def __init__(self, **kwargs):
        # Set up of the scattering lengths
        self.a11 = 1.
        self.a22 = kwargs.get("a22") / kwargs.get("a11")
        self.a12 = kwargs.get("a12") / kwargs.get("a11")
        self.m1  = 1. 
        self.m2  = kwargs.get("m2") / kwargs.get("m1")
        
class SimulationSetup():
    def __init__(self, **kwargs):        
        
        self.nxyz = np.array(kwargs.get("nxyz"), dtype=np.int)
        self.L_grid = np.array(kwargs.get("L_grid"), dtype=np.float32)
        L_grid_nobuffer_percentage = kwargs.get("L_grid_nobuffer_percentage")
        self.L_grid_nobuffer = self.L_grid * L_grid_nobuffer_percentage
        self.dtau = kwargs.get("dtau")
        self.dt = kwargs.get("dt")
        self.entol = kwargs.get("entol")
        self.printfreq = kwargs.get("printfreq")
        

        # Set up of the DFT parameters
        
