from Hamiltonians.harmonicoscillator import HarmonicOscillator
import numpy as np

def getHam(config):
    ham_config = config["hamiltonian"]
    if ham_config["name"] == "HO":
        return HarmonicOscillator(ham_config["omega"], np.array(ham_config["HO shape"], dtype = "float64"))
    else:
        raise NotImplementedError