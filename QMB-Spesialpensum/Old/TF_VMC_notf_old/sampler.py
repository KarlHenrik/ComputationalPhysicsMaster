import numpy as np
from WaveFunctions import *
from Hamiltonians import *


class Sampler:
    def __init__(self, config):
        self.steps = config["steps"]["sampled"]
        
        self.energy = 0
        self.energySq = 0
        # The wfDer can have any shape, so that should be updated
        self.wfDer = 0
        self.wfDerEnergy = 0
        
    def sample(self, particles, waveFunc, hamiltonian):
        energy = hamiltonian.potential(particles) + waveFunc.kinetic(particles)
        wfDer = waveFunc.paramDer(particles.getPositions())
        self.energy += energy
        self.energySq += energy * energy
        self.wfDer += wfDer
        self.wfDerEnergy += wfDer * energy
        
    def getSamples(self):
        return self.energy / self.steps, self.energySq / self.steps, self.wfDer / self.steps, self.wfDerEnergy / self.steps
        
    def getGrad(self):
        return 2 * (self.wfDerEnergy  / self.steps - self.wfDer * self.energy  / self.steps  / self.steps)