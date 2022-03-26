import numpy as np
from numba import int32, float64, deferred_type
from numba.experimental import jitclass

from WaveFunctions.simplegaussian import SimpleGaussian

def sampler_util(particles, wf, ham):
    
    particles_type = deferred_type()
    particles_type.define(Particles.class_type.instance_type)
    #dynamically using any wf or ham jitclass in this jitclass
    wf_type = deferred_type()
    wf_type.define(globals()[wf.__class__.__name__].class_type.instance_type)
    ham_type = deferred_type()
    ham_type.define(globals()[ham.__class__.__name__].class_type.instance_type)

    spec = [
        ('_particles', particles_type),
        ('_wf', wf_type),
        ('_ham', ham_type),
    ]

    @jitclass(spec)
    class Sampler:
        def __init__(self, wf):
            self._wf = wf

        def get(self):
            return self._wf