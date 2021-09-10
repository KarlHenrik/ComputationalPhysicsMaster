import numpy as np
from Hamiltonians.hamiltonian import Hamiltonian
import tensorflow as tf

class HarmonicOscillator(Hamiltonian):
    def __init__(self, omega, HOshape):
        self._omega2 = tf.constant(omega * omega, dtype = tf.float64)
        self._HOs = tf.constant(HOshape, dtype = tf.float64)
    
    @tf.function
    def potential(self, positions):
        return tf.math.reduce_sum(0.5 * self._omega2 * positions * positions * self._HOs)