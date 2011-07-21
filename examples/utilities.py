'''interface to pymutt.mtft that includes line detection and extraction
as well as simple plotting support if the pylab package is available.'''

import sys
from math import floor, ceil
import pymutt
try:
    import matplotlib.pyplot as mpl
except:
    mpl = None
import numpy as np

def mtanalyze(data,
              dt = 1.0,
              kind = 1,     # adaptive: 2 for vanilla hi-res
              npi = 3.0,
              nwin = 5,     # s.b. 2 * npi - 1
              padby = 8,    # pad to padby * power-of-2-holding-data
              nlines = 0,   # number of lines to find
              linedomain = None,  # could be a range of frequencies to search
              doplot = 0,
              title = None,
              title2 = None,
              verbose = 0,
              ):

    if padby and padby > 1:
        paddedlen = 2
        while paddedlen < len(data):
            paddedlen *= 2
        paddedlen *= padby
    else:
        paddedlen = 0

    if verbose:
        print >> sys.stderr, "original length: %d padded to: %d" \
              % (len(data), paddedlen)
        
    r = pymutt.mtft(data, dt = dt, npi = npi, nwin = nwin,
                    paddedlen = paddedlen, dodof = (nlines > 0))
    r['f'] = r['df'] * np.arange(len(r['power']))

    if nlines:
        # we operate on a copy of r['F'] so we can pass the original up
        # unaltered
        Fcpy = np.array(r['F'], copy = 1)
        flines = []
        Fmin = Fcpy.min() - 1.0
        # For each line, find the maximum current value of the F-test.
        # Then zero the F-test values in the region fline +- W/2 before
        # searching for the next largest value.
        if verbose:
            print >> sys.stderr, "\nspectral lines:"
        while len(flines) < nlines:
            idxfmax = np.argsort(Fcpy)[-1]
            fline = r['df'] * idxfmax
            if not linedomain or \
                   (fline >= linedomain[0] and fline <= linedomain[1]):
                flines.append(fline)
                if verbose:
                    print >> sys.stderr, "%5d  %13.6e" % (len(flines), fline)
            wl = max(0, int(idxfmax - round(0.5 * r['W'] / r['df'])))
            wr = min(r['n'] - 1, int(idxfmax + round(0.5 * r['W'] / r['df'])))
            Fcpy[wl:wr] = Fmin

        if len(flines) > 0:
            reducedr = pymutt.mtft(data,
                                   dt = dt,
                                   npi = npi,
                                   nwin = nwin,
                                   kind = kind,
                                   paddedlen = paddedlen,
                                   dodof = 1,
                                   lines = np.array(flines),
                                   )
            # added reshaped spectrum values to r
            r['reshaped'] = reducedr['reshaped']
            r['linea'] = reducedr['linea']
            r['linevar'] = reducedr['linevar']

    if doplot and mpl:
        if not linedomain:
            il = 0
            ir = r['n'] - 1
        else:
            il = int(floor(linedomain[0] / r['df']))
            ir = int(ceil(linedomain[1] / r['df']))
        if nlines:
            mpl.subplot(212)
            if title2:
                mpl.title(title2)
            mpl.plot(r['f'][il:ir], r['F'][il:ir])
            mpl.ylabel('F')
        mpl.xlabel('frequency')
        if nlines:
            mpl.subplot(211)
        if title:
            mpl.title(title)
        mpl.plot(r['f'][il:ir], 10.0 * np.log10(r['power'][il:ir]))
        mpl.ylabel('dB')
        if r.has_key('reshaped'):
            mpl.plot(r['f'][il:ir], 10.0 * np.log10(reducedr['reshaped'][il:ir]))
        mpl.show()


    return r


    
              
