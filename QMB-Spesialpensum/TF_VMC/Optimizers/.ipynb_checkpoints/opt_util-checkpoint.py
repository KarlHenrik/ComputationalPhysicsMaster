from Optimizers.gridSearch import GridSearch
import numpy as np

def getOptimizer(config):
    if config["name"] == "HO":
        return HarmonicOscillator(config["omega"], np.array(config["HO shape"], dtype = "float64"))
    else:
        raise NotImplementedError

