from WaveFunctions import *
from Hamiltonians import *
import tensorflow as tf

class Sampler:
    def __init__(self, config):
        self.steps = tf.constant(config["steps"]["sampled"], dtype = tf.float64)
        
        self.energy = tf.Variable(0, dtype = tf.float64)
        self.energySq = tf.Variable(0, dtype = tf.float64)
        # The wfDer can have any shape, so that should be updated
        self.wfDer = tf.Variable(0, dtype = tf.float64)
        self.wfDerEnergy = tf.Variable(0, dtype = tf.float64)
    
    @tf.function
    def sample(self, positions, waveFunc, hamiltonian):
        energy = hamiltonian.potential(positions) + waveFunc.kinetic(positions)
        wfDer = waveFunc.paramDer(positions)
        self.energy.assign_add(energy)
        self.energySq.assign_add(energy * energy)
        self.wfDer.assign_add(wfDer)
        self.wfDerEnergy.assign_add(wfDer * energy)
    
    @tf.function
    def getSamples(self):
        return self.energy / self.steps, self.energySq / self.steps, self.wfDer / self.steps, self.wfDerEnergy / self.steps
    
    @tf.function
    def getGrad(self):
        return 2 * (self.wfDerEnergy  / self.steps - self.wfDer * self.energy  / self.steps  / self.steps)