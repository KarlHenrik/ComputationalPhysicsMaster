import numpy as np
from particles import Particles   
from Hamiltonians.hamiltonian import Hamiltonian

class HarmonicOscillator(Hamiltonian):
    def __init__(self, omega, HOshape):
        self._omega2 = omega * omega
        self._HOs = HOshape
    
    def potential(self, particles):
        positions = particles.getPositions()
        return np.sum(0.5 * self._omega2 * positions * positions * self._HOs)