from particles import Particles

class Hamiltonian(): 
    
    def potential(self, particles):
        """Calculates the potential energy of the system.
        The kinetic energy depends on the wavefunction,
        and is handled by the wavefunction class"""
        return