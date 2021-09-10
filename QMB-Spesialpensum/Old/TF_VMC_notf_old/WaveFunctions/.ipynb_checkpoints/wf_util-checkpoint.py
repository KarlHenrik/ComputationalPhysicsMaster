from WaveFunctions.simplegaussian import SimpleGaussian
import numpy as np

def getWf(config):
    wf_config = config["wavefunc"]
    if wf_config["name"] == "simpleGaussian":
        return SimpleGaussian(wf_config["params"], np.array(wf_config["HO shape"], dtype = "float64"), config["particles"])
    else:
        raise NotImplementedError