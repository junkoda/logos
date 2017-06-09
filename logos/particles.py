import math
import numpy as np

class Particles:
    """
    particles = Particles(n, boxsize)

    Attributes:
      n            : number of particles
      x (array[n]) : position
      u (array[n]) : velocity (RSD displacement)
    """
    def __init__(self, n, boxsize):
        self.boxsize = boxsize
        self.n = n
        self.x = np.random.rand(n)*boxsize
        self.u = np.zeros_like(self.x)

        
    def set_velocity(self, grid):
        """ Set particle velocity from a grid
        Arg:
           grid (array): a grid of velocity field
        """
        u = self.u
        n = self.n
        nc = grid.shape[0]
        dx_inv = nc/self.boxsize
        
        for i, x in enumerate(self.x):
            ix = math.floor(x*dx_inv)
            w0 = x*dx_inv - ix
            w1 = 1.0 - w0
            u[i] = w0*grid[ix % nc] + w1*grid[(ix + 1) % nc]
