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
import pylab as mpl
import numpy as np

def doit():

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

    # print >> sys.stderr, "read in %d values at %e sps" % (len(a), sps)

    # select a padded length that is a power of 2
    n2 = 2
    while n2 < len(a):
        n2 *= 2
    paddedlen = 32 * n2

    # print >> sys.stderr, "padded to %d" % paddedlen

    kind = 2
    npi = 4
    nwin = 7
    dt = 1.0 / sps
    nlines_to_find = 9

    mpl.subplot(211)
    r = pymutt.mtft(a,
                    dt = dt,
                    kind = kind,
                    npi = npi,
                    nwin = nwin,
                    paddedlen =  paddedlen,
                    dodof = 1,
                    )
    freq = 1000.0 * r['df'] * np.arange(len(r['power']))
    il = int(round(1.5 / (1000.0 * r['df'])))
    ir = 2 * il

    Fcpy = np.array(r['F'], copy = 1)
    flines = []
    # For each line, find the maximum current value of the F-test.
    # Then zero the F-test values in the region fline +- W/2 before
    # searching for the next largest value.
    Fmin = Fcpy.min()
    while len(flines) < nlines_to_find:
        idxfmax = np.argsort(Fcpy)[-1]
        if idxfmax >= il and idxfmax <= ir:
            fline = r['df'] * idxfmax
            flines.append(fline)
            print >> sys.stderr, "%3d  %13.6e" % (len(flines), fline)
        # print >> sys.stderr, idxfmax
        wl = max(0, int(idxfmax - round(0.5 * r['W'] / r['df'])))
        wr = min(r['n'] - 1, int(idxfmax + round(0.5 * r['W'] / r['df'])))
        Fcpy[wl:wr] = Fmin

    mpl.subplot(212)
    mpl.plot(freq[il:ir], r['F'][il:ir])
    mpl.xlabel('mHz')
    mpl.ylabel('F')

    reducedr = pymutt.mtft(series = a,
                           dt = dt,
                           npi = npi,
                           nwin = nwin,
                           kind = kind,
                           paddedlen = paddedlen,
                           dodof = 1,
                           lines = np.array(flines),
                           )

    mpl.subplot(211)
    mpl.plot(freq[il:ir], 10.0 * np.log10(r['power'][il:ir]))
    mpl.ylabel('dB')
    mpl.plot(freq[il:ir], 10.0 * np.log10(reducedr['reshaped'][il:ir]))
    for line in flines:
        mpl.plot((1000.0 * line, 1000.0 * line), (80.0, 120.0))
    mpl.show()
