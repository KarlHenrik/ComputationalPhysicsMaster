import numpy as np

class Particles:
    def __init__(self, num, dims):
        self._num = num
        self._dims = dims
        
        self._positions = np.zeros((num, dims), dtype = "float64")
        self._norms = np.zeros(num, dtype = "float64")
        
        self._dists = np.zeros((num, num), dtype = "float64") # from, to
        self._units = np.zeros((num, num, dims), dtype = "float64") #from, to, dim
        
    # Getters
    def getPos(self, p1):
        return self._positions[p1]

    def getDist(self, p1, p2):
        return self._dists[p1, p2]

    def getUnit(self, p1, p2):
        return self._units[p1, p2]
    
    def getSize(self):
        return (self._num, self._dims)
    
    def getPositions(self):
        return self._positions
    
    # Setters
    def setPos(self, p1, pos):
        self._positions[p1] = pos
        self._update(p1)
        
    def adjustPos(self, p1, adjustment):
        self._positions[p1] += adjustment
        self._update(p1)
    
    def setPositions(self, positions):
        self._positions = positions
        for p1 in range(self._num):
            self._update(p1) # will be called rarely, and is therefore not optimized
        
    # Utility
    def _update(self, p1):
        self._norms[p1] = np.linalg.norm(self._positions[p1])
        for p2 in range(self._num):
            if p1 != p2:
                dist = np.linalg.norm( self._positions[p1] - self._positions[p2] )
                self._dists[p1, p2] = dist
                self._dists[p2, p1] = dist
                if dist != 0: # If the distance is zero, the unit vectors will be put to zero, which is fine!
                    dist = 1 / dist
                self._units[p1, p2] = ( self._positions[p2] - self._positions[p1] ) * dist
                self._units[p2, p1] = -self._units[p1, p2]