from WaveFunctions import *


class Optimizer():
    def __init__(self, config):
        self.waveFunc = wf_util.getWf(config)
        self.params = config["wavefunc"]["params"]
    
    def optimize(self, gradient):
        return
        
    def getWaveFunc(self):
        return self.waveFunc