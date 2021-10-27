import numpy as np
from Hamiltonians import *
from initialstate import initial_util
from WaveFunctions import *
from particles import Particles
from numba import njit
#from sampler import Sampler


def metro_util(config):
    # Creating the correct step function
    metro_name = config["metropolis"]["name"]
    if metro_name == "metro":
        step = create_metro_step(config)
    elif metro_name == "importance":
        step = create_importance_step(config)
    else:
        raise NotImplementedError
        
    return create_run_steps(config, step)

def create_run_steps(config, step):
    num         = config["particles"]
    dims        = config["dimensions"]
    equil_steps = config["metropolis"]["equilibrationsteps"]
    num_steps   = config["metropolis"]["steps"]
    initialize  = initial_util(config)
    
    @njit
    def run_steps(wf, ham):
        particles = Particles(num, dims)
        initialize(particles)
        #sampler = Sampler(particles, wf, ham)

        for i in range(equil_steps):
            step(particles, wf, ham)
        for i in range(num_steps):
            accepted = step(particles, wf, ham)
            #sampler.sample(accepted)
        return 2#sampler
    
    return run_steps
    
def create_metro_step(config):
    num = config["particles"]
    dims = config["dimensions"]
    step_length = config["metropolis"]["step_length"]
    @njit
    def metro_step(particles, wf, ham = None):
        p1 = np.random.randint(num)
        oldPos = particles.getPos(p1)

        adj_dir = np.random.randint(dims)
        adj_sign = (np.random.randint(2) - 0.5) * 2
        newPos = oldPos + 0
        newPos[adj_dir] += step_length * adj_sign

        particles.setPos(p1, newPos)
        ratio = wf.ratio(particles, p1, oldPos)
        if (np.random.random() < np.power(ratio, 2)):
            return True
        else:
            particles.setPos(p1, oldPos)
            return False
    return metro_step

def create_importance_step(config):
    num = config["particles"]
    dims = config["dimensions"]
    step_length = config["metropolis"]["step_length"]
    @njit
    def importance_step(particles, wf, ham = None):
        p1 = np.random.randint(num)
        oldPos = particles.getPos(p1)

        oldQF = wf.QF(particles, p1)
        move = 0.5 * oldQF * step_length + np.random.normal(dims) * np.sqrt(step_length)
        newPos = oldPos + move
        
        particles.setPos(p1, newPos)
        ratio = wf.ratio(particles, p1, oldPos)
        
        newQF = wf.QF(particles, p1)
        greensFuncRatio = np.sum( (oldQF + newQF) * (oldPos - newPos + step_length * (oldQF - newQF)) )
        greensFuncRatio = np.exp(0.5 * greensFuncRatio)
        if (np.random.random() < greensFuncRatio * np.power(ratio, 2)):
            oldQF = newQF
            return True
        else:
            particles.setPos(p1, oldPos)
            return False
    return importance_step