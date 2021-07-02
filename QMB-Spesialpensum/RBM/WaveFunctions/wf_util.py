from WaveFunctions.simplegaussian import SimpleGaussian
import numpy as np

def new(config):
    if config["name"] == "simpleGaussian":
        return SimpleGaussian(config["alpha"], np.array(config["HO shape"], dtype = "float64"))
    else:
        raise NotImplementedError