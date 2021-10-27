from WaveFunctions.wavefunction import WaveFunction
import tensorflow as tf

class SimpleGaussian(WaveFunction):
    def __init__(self, alpha, HOshape, num):
        self._alpha = tf.constant(alpha, dtype = tf.float64)
        self._HOs = tf.constant(HOshape, dtype = tf.float64)
        self._num = tf.constant(num, dtype = tf.float64)
    
    @tf.function
    def kinetic(self, positions):
        return self._alpha * tf.math.reduce_sum(self._HOs - 2 * self._alpha * tf.math.pow(positions * self._HOs, 2))
    
    @tf.function
    def QF(self, positions, p1):
        return -4 * self._alpha * positions[p1] * self._HOs
    
    @tf.function
    def ratio(self, positions, p1, old_positions):
        """
        Old wavefunc value term: exp(-alpha * old_r2)
        New wavefunc value term: exp(-alpha * r2)
        All other terms are the same and cancel. Since they are both exponentials, we can subtract the exponents
        """
        r2 = tf.math.reduce_sum(positions[p1]**2 * self._HOs)
        old_r2 = tf.math.reduce_sum(old_positions[p1]**2 * self._HOs)
        return tf.math.exp(-self._alpha * r2 + self._alpha * old_r2)
    
    @tf.function
    def paramDer(self, positions):
        return -tf.math.reduce_sum(positions * positions * self._HOs) / self._num