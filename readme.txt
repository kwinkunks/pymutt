About Pymutt
------------

Pymutt is an implementation of Thomson (1982) multi-taper fourier spectral
estimator plus a python interface.  The core code is due to Lees and
Park (1995) and uses the conventions of Percival and Walden (1993).


Installation
------------

Pymutt requires numpy.  Most of the examples also require matplotlib.

Pymutt builds and installs properly on
    - win7-64 using 32-bit python
    - linux-64 using 64-bit python
It almost surely installs properly on linux-32 and win-32.  It
fails on win7-64 using 64-bit python.

To install pymutt execute

    python setup.py build

and as root or administrator

    python setup.py install

On windows I used mingw and the mingw shell.


Usage
-----

Go to the examples subdirectory and execute

    python runme.py

and it will print out help for runme.py as well as the self docs for
pymutt.  To run a test case try

    python runme.py -N 3

and you should get a figure closely resembling Percival and Walden
figure 512.

Look at the code in examples/minimalexamples.py for some easy-to-copy examples.


Other Code
----------

For a versatile, extensive, and newer fortran system for mtft see

     http://uniandes.academia.edu/gprieto


References
----------

Lees, J.M., and J. Parks (1995), Multi-taper Spectral Analysis: A
    Stand-alone C Subroutine, Computers and Geology, 21(2), 199-236.
Percival, D.B., and A.T. Walden (1993) Spectral Analysis for Physical
    Applications, Cambridge University Press, Cambridge.
Thomson, D.J. (1982), Spectral Estimation and Harmonic Analysis, IEEE
    Proc., 70, 1055-1096.
