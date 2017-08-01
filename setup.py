from distutils.core import setup, Extension
import numpy as np

idirs = ['src', np.get_include()]

setup(name='logos',
      version='0.0',
      description='1D reshift-space distortion library',
      author='Jun Koda',
      author_email='junkoda@gmail.com',
      url='https://github.com/junkoda/logos',
      py_modules=['logos.particles', 'logos.gird', 'logos.power_spectrum'],
      ext_modules=[
          Extension('logos._logos',
                    ['src/corr.cpp', ],
                    include_dirs = idirs,
          )
      ],
      packages=['logos'],
     )
