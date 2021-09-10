from InitialPositions import *
from WaveFunctions import *
from Hamiltonians import *
from Optimizers import *
from sampler import Sampler
from metropolis import metro_util

def run(config):
    initializePos = initPos_util.getInit(config)
    metro_step    = metro_util.getStep(config)
    hamiltonian   = ham_util.getHam(config)
    optimizer     = opt_util.getOptimizer(config)
    
    while optimizer.notFinished:
        sampler   = Sampler(config)
        particles = initializePositions()
        waveFunc  = optimizer.getWaveFunc()
        
        #tf function call?
        for i in range(config["steps"]["equil"]):
            metro_step(particles, waveFunc)
        #tf function call with gradienttape and two wavefunctions?
        for i in range(config["steps"]["sampled"]):
            metro_step(particles, waveFunc)
            sampler.sample(particles, waveFunc, hamiltonian)
            
        gradient = sampler.getGrad()
        optimizer.optimize(gradient)
        