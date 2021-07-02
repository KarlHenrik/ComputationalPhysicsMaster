import numpy as np
from particles import Particles   
from numba import int32, float64
from numba.experimental import jitclass
from Hamiltonians.hamiltonian import Hamiltonian

spec = [
    ("_omega2", float64),
    ("_HOs", float64[:])
]

@jitclass(spec)
class HarmonicOscillator(Hamiltonian):
    def __init__(self, omega, HOshape):
        self._omega2 = omega * omega
        self._HOs = HOshape
    
    def potential(self, particles):
        positions = particles.getPositions()
        return np.sum(0.5 * self._omega2 * positions * positions * self._HOs)