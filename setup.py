#
# fail on 64bit win7 using: python setup.py build -c mingw32
#

from numpy.distutils.core import setup
from numpy.distutils.extension import Extension
import os.path
import distutils.sysconfig as sf
import numpy

use_fftpack = 1

if use_fftpack:

    pymutt_module = Extension(
        'pymutt',
        ['pymutt.c'],
        library_dirs = [os.path.join(numpy.__path__[0], 'fft')],
        define_macros = [('USE_FFTPACK', None)],
		extra_compile_args = ["--verbose",],
        )

else:

    pymutt_module = Extension(
        'pymutt',
        ['pymutt.c'],
        libraries = ['fftw3'],
        )

setup(name='pymutt',
      version='1.0',
      description = "numpy support for multi-taper fourier transforms",
      author = "Martin L. Smith",
      long_description =
'''
Implements a numpy interface to an implementation of Thomson's multi-taper
fourier transform algorithms.  The key C module (mtbase.c) is derived from 
code written and made freely available by J. M. Lees and Jeff Park.
''',
      ext_modules=[pymutt_module],
    )
#@-node:@file setup.py
#@-leo
