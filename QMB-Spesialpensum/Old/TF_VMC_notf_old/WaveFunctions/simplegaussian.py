import numpy as np
from particles import Particles   
from WaveFunctions.wavefunction import WaveFunction

class SimpleGaussian(WaveFunction):
    def __init__(self, alpha, HOshape, particles):
        self._alpha = alpha
        self._HOs = HOshape
        self._num = particles
    
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
    
    def paramDer(self, positions):
        return -np.sum(positions * positions * self._HOs) / self._num