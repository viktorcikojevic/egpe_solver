class SimulationSetup():
    def __init__(self, **kwargs):
        
        # Set up of the scattering lengths
        self.a11 = kwargs.get("a11")
        self.a22 = kwargs.get("a22")
        self.a12 = kwargs.get("a12")
        
        
        # Set up of the DFT parameters
        if kwargs["theory"] == "MFLHY":
            x = 1.
            alpha = 0.5 * 4*np.pi*(1.+x*x + 2*x*sim.a12/np.sqrt(sim.a22)) / \
                                np.power(1 + x/np.sqrt(sim.a22), 2);
            beta = 0.5 * (512*np.sqrt(np.pi)/15) * \
                        np.power((1+x*np.sqrt(sim.a22)) / (1+x/np.sqrt(sim.a22)), 5./2);
            gamma = 3./2;
        elif kwargs["theory"] == "DFT":
            alpha = kwargs["alpha"]
            beta  = kwargs["beta"]
            gamma = kwargs["gamma"]
        else:
            raise KeyError("You need to put --theory either to MFLHY or DFT")
        
        
