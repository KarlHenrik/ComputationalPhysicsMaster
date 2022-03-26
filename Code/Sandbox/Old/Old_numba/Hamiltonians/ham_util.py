from Hamiltonians.harmonicoscillator import HarmonicOscillator
import numpy as np

def new(config):
    if config["name"] == "HO":
        return HarmonicOscillator(config["omega"], np.array(config["HO shape"], dtype = "float64"))
    else:
        raise NotImplementedError