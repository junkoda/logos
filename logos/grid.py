import numpy as np
import particles

class Grid:
    """
    Attributes:
        fx (array): real-space grid
        fk (array): Fourier-space grid
        shifted (Grid): shifted grid for interlacing
    """
    def __init__(nc, boxsize, *, interlacing=False):
        """Create a Grid from array a
        Arg:
          a (array)
        """
        self.boxsize = boxsize
        self.nc = nc
        self.fx = np.zeros(nc)
        self.mode = 'real-space'

        if interlacing:
            self.shifted = Grid(nc, boxsize)
        else:
            self.shifted = None

        def fft_forward(self):
            """FFT fx to fk
            fk = \int f(x) exp(-ikx) dx = \sum f(x) exp(-ikx) dx

            Returns:
                fk (complex array)
            """
            assert(self.mode == 'real-space')
            
            dx = self.boxsize/fx.shape[0]
            self.fk = dx*np.fft.rfft(fx)
            
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
            """

            nc = self.nc
            dx_inv = nc / self.boxsize
            d = grid.fx
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

