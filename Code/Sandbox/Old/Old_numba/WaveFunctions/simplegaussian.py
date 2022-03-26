import numpy as np
from particles import Particles   
from numba import int32, float64
from numba.experimental import jitclass
from WaveFunctions.wavefunction import WaveFunction

spec = [
    ("_alpha", float64),
    ("_HOs", float64[:])
]

@jitclass(spec)
class SimpleGaussian(WaveFunction):
    def __init__(self, alpha, HOshape):
        self._alpha = alpha
        self._HOs = HOshape
    
    def kinetic(self, particles):
        return self._alpha * np.sum(self._HOs - 2 * self._alpha * np.power(particles.getPositions() * self._HOs, 2))
    
    def QF(self, particles, p1):
        return -4 * self._alpha * particles.getPos(p1) * self._HOs
    
    def ratio(self, particles, p1, oldPos):
        """
        Old wavefunc value term: exp(-alpha * old_r2)
        New wavefunc value term: exp(-alpha * r2)
        All other terms are the same and cancel. Since they are both exponentials, we can subtract the exponents
        """
        r2 = np.sum(particles.getPos(p1) * particles.getPos(p1) * self._HOs)
        old_r2 = np.sum(oldPos * oldPos * self._HOs)
        return np.exp(-self._alpha * r2 + self._alpha * old_r2)
    
    def paramDer(self, particles):
        positions = particles.getPositions()
        num, _ = particles.getSize()
        return -np.sum(positions * positions * self._HOs) / num
    
    def updateParam(self, alpha):
        self._alpha = alpha