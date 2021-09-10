import numpy as np
from particles import Particles
#from numba import njit

def create_randomUniform(config):
    num = config["particles"]
    dims = config["dimensions"]
    radius = config["initial"]["radius"]
    
    #@njit
    def findAvaliable(positions, max_idx):
        invalid = True
        attempts = 0
        while invalid:
            attempts += 1
            if attempts > 100:
                raise ValueError("No more room for particles. Reduce the number of particles or the hard-shell radius")
            invalid = False
            candidate = (np.random.rand(dims) - 0.5) * 2
            for p1 in range(max_idx):
                if np.linalg.norm(positions[p1] - candidate) < radius:
                    invalid = True
        return candidate
    
    def initialize():
        positions = np.zeros((num, dims))
        
        positions[0] = (np.random.rand(dims) - 0.5) * 2
        for p1 in range(1, num):
            positions[p1] = findAvaliable(positions, p1)
            
        particles = Particles(num, dims)
        particles.setPositions(positions)
        return particles
    
    return initialize