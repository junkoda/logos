import math
import numpy as np

def periodic_wrapup(x, boxsize):
    x[x < 0.0] += boxsize
    x[x >= boxsize] -= boxsize

    return x

class Particles:
    """
    particles = Particles(n, boxsize)

    Attributes:
      n            : number of particles
      x (array[n]) : position
      u (array[n]) : velocity (RSD displacement)
    """
    def __init__(self, n, boxsize, *, pos='random'):
        self.boxsize = boxsize
        self.n = n

        if pos == 'random':
            self.x = np.random.rand(n)*boxsize
        elif pos == 'lattice':
            self.x = (0.5 + np.arange(n))/n*boxsize
        else:
            raise TypeError('Unknown pos parameter %s' % pos)
        
        self.u = np.zeros_like(self.x)

        
    def set_velocity_from_grid(self, grid):
        """ Set particle velocity `particles.u` from a grid
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

    def set_velocity_from_particles(self, x, v):
        """ Set particle velocity particles.u from of nearest particle

        Args:
           x (array): array of positions
           v (array): array of velocities
        """

        assert(x.shape == v.shape)
        periodic_wrapup(x, self.boxsize)
        
        i = np.argsort(x)

        # set sorted array with periodic wrapup
        n = x.shape[0]
        xp = np.zeros(n + 2)
        vp = np.zeros(n + 2)
        xp[0] = x[i[-1]] - self.boxsize
        vp[0] = v[i[-1]]
        xp[-1] = x[i[0]] + self.boxsize
        vp[-1] = v[i[0]]
        
        xp[1:-1] = x[i]
        vp[1:-1] = v[i]

        ix = 0

        self.x.sort()
        
        for i, x_i in enumerate(self.x):
            while xp[ix] < x_i:
                ix += 1

            assert(xp[ix - 1] <= x_i <= xp[ix])
            
            if x_i - xp[ix - 1] < xp[ix] - x_i:
                # left neighbbr is nearest
                self.u[i] = vp[ix - 1]
            else:
                # right neighbor is nearest
                self.u[i] = vp[ix]
