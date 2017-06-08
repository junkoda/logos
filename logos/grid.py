import math
import numpy as np
import logos.particles

class Grid:
    """
    grid = Grid(nc, boxsize)

    Attributes:
        x (array):  grid centre position
        fx (array): real-space grid
        fk (array): Fourier-space grid
        shifted (Grid): shifted grid for interlacing

        n (int): number of particles assigned to density
    """
    def __init__(self, nc, boxsize, *, interlacing=True):
        """Create a Grid from array a
        Arg:
          a (array)
        """
        self.boxsize = boxsize
        self.fx = np.zeros(nc)
        self.mode = 'real-space'
        self.x = (np.arange(nc) + 0.5)*(boxsize/nc)
        self.n = 0

        if interlacing:
            self.shifted = Grid(nc, boxsize)
        else:
            self.shifted = None

    @property
    def nc(self):
        return self.fx.shape[0]
    
    def fft_forward(self):
        """FFT fx to fk
        fk = \int f(x) exp(-ikx) dx = \sum f(x) exp(-ikx) dx
        
        Returns:
            fk (complex array)
        """
        assert(self.mode == 'real-space')
            
        dx = self.boxsize/self.nc
        self.fk = dx*np.fft.rfft(self.fx)
            
        return self.fk

    def fft_inverse(self):
        """FFT fk to fx
        f(x) = \int f(k) exp(ikx) dk/(2pi) = 1/L \sum f(k) exp(ikx)
        
        Returns:
            fx (array)
        """
        assert(self.mode == 'fourier-space')

        self.fx = np.fft.irfft(uk)/boxsize

        return self.fx

    def assign_density(self, a, *, shift=0.5):
        """Assign density of particles to the grid
        Args:
            particles [Particles]
            shift (float): shift grids in units of grid spacing

        Result:
            fx become density contract delta(x) of particle density

        Note:
            density is cleared to zero before assignment
        """

        self.n = a.shape[0]
        self.shotnoise = self.boxsize/self.n

        nc = self.nc
        dx_inv = nc / self.boxsize
        d = self.fx
        d.fill(0.0)
    
        # CIC mass assignment
        for x in a:
            xx = x*dx_inv - shift
            ix = math.floor(xx)
            w1 = xx - ix
            w0 = 1.0 - w1
            d[ix % nc] += w0
            d[(ix + 1) % nc] += w1
    
        nbar = len(a)/nc
        d = d/nbar - 1.0

        self.mode = 'real-space'



        return self

