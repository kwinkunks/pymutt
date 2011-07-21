'''This test analyzes a long-period seismic records (courtesy of Pinon
Flat Observatory (http://pfostrain.ucsd.edu/pfo/).  The data are about
a day in duration and follow a magnitude 8.5 near Chile.  In the
frequency range plotted (1.4-3.0 mHz) the data processing shows clear
evidence for a series of line spectra.

These correspond to the Earth\'s free oscillations.  It\'s interesting
to compare the spectral lines inferred here to the list of observed
and theoretical frequencies tabulated by Masters and Widmer in
http://igppweb.ucsd.edu/~guy/sio224/mode_lect1.pdf. Note that the
actual signals are decaying sinusoids that vary perceptibly in
amplitude over the course of the dataset.  This violates our
theoretical model of a pure sinusoid and may account in part for the
residual structure in the reshaped power spectrum.
'''

import sys
import pymutt
import matplotlib.pyplot as mpl
import numpy as np
from utilities import mtanalyze

def doit(verbose = 0):

    print __doc__
    fn = "chile.txt"

    sps = None
    with file(fn, "rb") as f:
        for l in f:
            if l.startswith(('!', '$')):
                continue
            c = l.find(":")
            if c < 0:
                break
            if l[:c] == "SPS":
                sps = float(l[(c+1):])

        a = [float(l.strip())]
        for l in f:
            a.append(float(l.strip()))
        a = np.array(a)
        f.close()

    if verbose:
        print >> sys.stderr, "read in %d values at %e sps" % (len(a), sps)

    r = mtanalyze(a,
                  dt = 1.0/sps,
                  kind = 1,
                  npi = 4,
                  nwin = 7,
                  padby = 32,
                  nlines = 9,
                  linedomain = [0.0015, 0.003],
                  doplot = 1,
                  verbose = verbose,
                  )

    return r
