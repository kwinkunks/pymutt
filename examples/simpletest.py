#@+leo-ver=4
#@+node:@file examples/simpletest.py
#@@language python
#!/usr/bin/env python

''' a simple test just to see if mtft runs at all.

    here's a typical result using fftw3:
    
        python example/simpletest.py
        test series length padded from 716 to 4096
        averaging over 20 trials
         frequency    real     imag         abs          F
         0.06983     4.996    0.095      4.9973     9957.5
         0.07263    -0.204   -0.060      0.2130        0.0
    '''

import sys
from math import log, ceil, pi, pow

import numpy as np
from numpy.random import normal
import pymutt

# pick a miscellaneous length and then find a much bigger power-of-2

n = 716
dt = 1.0
df = 1.0 / (dt * n)
pow2 = ceil(log(5 * n) / log(2.0))
paddedlen = int(pow(2, pow2))

print __doc__
print "test series length padded from %d to %d" % (n, paddedlen)

# line and noise components

fline   = 50 * df
aline   = 10.0
noise_amplitude = 1.0

t = np.arange(n) * dt

lines = np.array([fline, fline + 2.0 * df])
nlines = len(lines)

suma = complex(0.0, 0.0) * np.zeros([nlines])
sumF = None
sumv = np.zeros([nlines])

trials = 20

print "averaging over %d trials" % trials

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
    print >> sys.stderr, r['linevar']
    suma += r['linea']
    sumv += r['linevar']
    if sumF is None:
        sumF = r['F']
    else:
        sumF += r['F']

print " frequency    real     imag         abs          F      line var"
for i, f in enumerate(lines):
    amp = suma[i] / trials
    var = sumv[i] / trials
    ifrq = int(round(f / r['df']))
    Ftest = sumF[ifrq] / trials
    print "%8.5f  %8.3f %8.3f    %8.4f   %8.1f    %10.3e" \
        % (f, amp.real, amp.imag, abs(amp), Ftest, var)

#@-node:@file examples/simpletest.py
#@-leo
