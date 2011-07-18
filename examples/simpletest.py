'''This is a simple, randomized test just to see if mtft runs at all.
The code generates a series comprising white noise and a sinusoid.  It
then applies an F-test for harmonic terms at the sinusoid frequency
and at a frequency four fourier steps away.  The first of these should
be significant and have an amplitude of 5.0 and zero phase (on
average); the second should be insignificant and have a vanishing
amplitude.

Below is a typical result (but note that individual runs can
vary appreciable):
    
    test series length padded from 716 to 4096
    averaging over 20 trials
    frequency    real     imag         abs          F       line var
    0.06983     4.991    0.104      4.9921     8511.9      2.083e-03
    0.07542    -0.013   -0.005      0.0141        0.0      1.438e-03
'''

import sys
from math import log, ceil, pi, pow
import numpy as np
from numpy.random import normal
import pymutt

def doit(verbose = 0):
    # pick a miscellaneous length and then find a much bigger power-of-2

    n = 716
    dt = 1.0
    df = 1.0 / (dt * n)
    pow2 = ceil(log(5 * n) / log(2.0))
    paddedlen = int(pow(2, pow2))

    print __doc__
    print "here are the actual results for this run:\n"
    print "    test series length padded from %d to %d" % (n, paddedlen)

    # line and noise components

    fline   = 50 * df
    aline   = 10.0
    noise_amplitude = 1.0

    t = np.arange(n) * dt

    lines = np.array([fline, fline + 4.0 * df])
    nlines = len(lines)

    suma = complex(0.0, 0.0) * np.zeros([nlines])
    sumF = None
    sumv = np.zeros([nlines])

    trials = 20

    print "    averaging over %d trials" % trials

    for i in range(trials):

        y = noise_amplitude * normal(size = n) \
            + aline * np.cos(2.0 * pi * fline * t)

        r = pymutt.mtft(series = y,
                        dt = dt,
                        npi = 2,
                        nwin = 8,
                        kind = 1,
                        paddedlen = paddedlen,
                        dodof = 1,
                        lines = lines,
                        )
        suma += r['linea']
        sumv += r['linevar']
        if sumF is None:
            sumF = r['F']
        else:
            sumF += r['F']

    print "     frequency    real     imag         abs          F      line var"
    for i, f in enumerate(lines):
        amp = suma[i] / trials
        var = sumv[i] / trials
        ifrq = int(round(f / r['df']))
        Ftest = sumF[ifrq] / trials
        print "    %8.5f  %8.3f %8.3f    %8.4f   %8.1f    %10.3e" \
              % (f, amp.real, amp.imag, abs(amp), Ftest, var)

    return 0
