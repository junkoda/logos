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
        nc (int): number of grid points
    """
    def __init__(self, nc, boxsize, *, interlacing=True):
        """
        Arg:
          nc (int): number of grid points
          boxsize (float): length of the periodic segment
          interlacing (bool): use interlacing aliasing reduction (default: True)
        """
        
        self.n = 0
        self.boxsize = boxsize
        self.fx = np.zeros(nc)
        self.mode = 'real-space'
        self.x = (np.arange(nc) + 0.5)*(boxsize/nc)
        self.interlaced = False

        # wavenumber
        nck = nc // 2 + 1
        self.k = np.arange(nck)*(2.0*math.pi/self.boxsize)


        if interlacing:
            self.shifted = Grid(nc, boxsize, interlacing=False)
        else:
            self.shifted = None

    @property
    def nc(self):
        return self.fx.shape[0]
    
    def fft_forward(self):
        """FFT fx to fk
        fk = \int f(x) exp(-ikx) dx = \sum f(x) exp(-ikx) dx        
        """
        assert(self.mode == 'real-space')
            
        dx = self.boxsize/self.nc
        self.fk = dx*np.fft.rfft(self.fx)
        self.mode = 'fourier-space'
        
        if self.shifted is not None:
            self.shifted.fft_forward()
            
        return self

    def fft_inverse(self):
        """FFT fk to fx
        f(x) = \int f(k) exp(ikx) dk/(2pi) = 1/L \sum f(k) exp(ikx)
        """
        assert(self.mode == 'fourier-space')
        self.mode == 'real-space'

        self.fx = np.fft.irfft(uk)/boxsize
        self.mode == 'real-space'

        if self.shifted is not None:
            self.shifted.fft_inverse()

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

        self.mode = 'real-space'
        self.n = a.shape[0]
        self.shotnoise = self.boxsize/self.n
        self.interlaced = False

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
        
        self.fx = d/nbar - 1.0

        if self.shifted is not None:
            self.shifted.assign_density(a, shift=shift - 0.5)

        return self


    def interlace(self):
        if self.shifted is None:
            print('no interlacing')
            return
        
        if self.mode == 'real-space':
            self.fft_forward()

        if self.interlaced:
            raise RuntimeError('grid.interlace already done')

        nck = self.nc // 2 + 1
        ik = np.arange(nck)
        
        expkx = np.zeros(nck, dtype=complex)
        fac = math.pi/self.nc
        expkx.real = np.cos(fac*ik)
        expkx.imag = np.sin(fac*ik)

        assert(self.shifted.mode == 'fourier-space')

        self.fk = 0.5*(self.fk + expkx*self.shifted.fk)

        self.interlaced = True

        return self
