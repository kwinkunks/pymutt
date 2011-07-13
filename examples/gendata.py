#!/usr/bin/env python
#@+leo-ver=4
#@+node:@file gendata.py
#@@first
#@@language python

'''create datasets for pymutt usage doc.
'''

import sys, traceback
from math import log, ceil, pi, pow

import getopt

import numpy as np
from numpy.random import normal
import numpy.fft as dft

import pymutt

#@<<code>>
#@+node:<<code>>
#@+others
#@+node:setheader
def setheader(setout, title = None, subtitle = None):

    global _set_id
    _set_id = 0

    if title is not None:
        print >> setout, '@title "%s"' % str(title)

    if subtitle is not None:
        print >> setout, '@subtitle  "%s"' % str(subtitle)

#@-node:setheader
#@+node:setdata
def setdata(setout, x, y, legend = None):

    global _set_id

    if legend is not None:    
        print >> setout, '@ s%d legend "%s"' % (_set_id, str(legend))
    _set_id += 1

    print >> setout, "@type xy"

    for i in range(len(x)):
        print >> setout, "%12.5e  %12.5e" % (x[i], y[i])
    print >> setout, "&"

    setout.flush()

#@-node:setdata
#@+node:bandlimited
def bandlimited(n, dt, fu, fs):

    BigN = 20 * n
    y = normal(size = BigN)
    fy = dft.rfft(y, BigN)
    nf = len(fy)
    df = 1.0 / (BigN * dt)
    ifu = int(round(fu / df))
    ifs = int(round(fs / df))

    for i in range(ifu, ifs):
        fy[i] = (i - float(ifs))/float(ifu - ifs) * fy[i]
    fy[ifs:] = complex(0, 0)

    y = dft.irfft(fy, BigN)

    subidx  = int(0.374 * BigN)

    return y[subidx:(subidx + n)]


#@-node:bandlimited
#@+node:filteredprocess
def filteredprocess():

    n = 937
    dt = 1e-3
    fu = 80.0
    fs = 90.0
    yr = bandlimited(n, dt, fu, fs)
    setout = file('ranproc.dat', 'w')
    setheader(setout,
              title = 'filtered normal series',
              subtitle = "fpass = %.1f  fstop = %.1f" % (fu, fs)
        )
    setdata(setout, dt * np.arange(n), yr, legend = 'y(t)')
    setout.close()


    npi = 4.0
    nwin = 10

    setout = file('ranprocmtft.dat', 'w')
    setheader(setout,
              title = 'filtered normal series',
              subtitle = "npi = %.0f  nwin = %d" % (npi, nwin)
        )

    r = pymutt.mtft(series = yr,
                    dt = dt,
                    npi = npi,
                    nwin = nwin,
                    kind = 1,
                    paddedlen = 4 * n,
                    dodof = 1,
        )
    f = r['df'] * np.arange(len(r['power']))
    for psdkey in ['power', 'dof']:
        if not r.has_key(psdkey):
            continue
        setdata(setout, f, r[psdkey], legend = 'hires ' + psdkey)

    r = pymutt.mtft(series = yr,
                    dt = dt,
                    npi = npi,
                    nwin = nwin,
                    kind = 2,
                    paddedlen = 4 * n,
                    dodof = 1,
        )
    for psdkey in ['power', 'dof']:
        if not r.has_key(psdkey):
            continue
        setdata(setout, f, r[psdkey], legend = 'adaptive ' + psdkey)

    yrp  = dt * abs(dft.rfft(yr, n = n)) ** 2
    df = 1.0 / (n * dt)
    f = df * np.arange(len(yrp))
    setdata(setout, f, yrp / n, legend = 'raw')
    tscale = 1.0 /np.hanning(n).mean()
    yrph = dt * abs(dft.rfft(tscale * np.hanning(n) * yr, n = n)) ** 2
    setdata(setout, f, yrph / n, legend = 'hanning taper')


    setout.close()



#@-node:filteredprocess
#@+node:linedetection
def linedetection():

    n = 716
    dt = 1.0e-3
    df = 1.0 / (dt * n)

    npi = 4.0
    nwin = 10

    # line and noise components

    fline   = 50.0
    aline   = 10.0
    noise_amplitude = 1.0

    t = np.arange(n) * dt

    lines = np.array([fline, fline + 2.0 * df])
    nlines = len(lines)

    y = noise_amplitude * normal(size = n) \
            + aline * np.cos(2.0 * pi * fline * t)

    r = pymutt.mtft(series = y,
            dt = dt,
            npi = npi,
            nwin = nwin,
            kind = 1,
            paddedlen = 10 * n,
            dodof = 1,
           )
    Fsortidx = r['F'].argsort()
    lidx = Fsortidx[-1]
    linefs = [r['df'] * lidx]
    # print >> sys.stderr, (linefs, r['df'], n * r['df'])

    r = pymutt.mtft(series = y,
            dt = dt,
            npi = npi,
            nwin = nwin,
            kind = 1,
            paddedlen = 10 * n,
            dodof = 1,
            lines = np.array(linefs),
           )

    f = r['df'] * np.arange(len(r['power']))
    setout = file('lineprocmtft.dat', 'w')
    setheader(setout,
              title = 'line + noise',
              subtitle = "npi = %.0f  nwin = %d" % (npi, nwin)
        )

    for psdkey in ['power', 'reshaped', 'F']:
        if not r.has_key(psdkey):
            continue
        setdata(setout, f, r[psdkey], legend = psdkey)


    setout.close()
#@-node:linedetection
#@-others
#@nonl
#@-node:<<code>>
#@nl

#@<<main>>
#@+node:<<main>>
#@+others
#@+node:main
def main():

    global verbose
    pscommands = []

    options = "vhL:P:"

    opts, args = getopt.getopt(sys.argv[1:], options, ["help"])
    for opt, v in opts:
        if opt == "-v":
            verbose += 1
        else:
            print >> sys.stderr, "unknown option: %s %s\n" \
                        % (str(opt), str(v))
            print >> sys.stderr, __doc__
            return 3

    filteredprocess()

    linedetection()

    return 0
#@-node:main
#@-others

if __name__ == "__main__":
    try:
        r = main()
        if r != 0:
            print >> sys.stderr, "Exit code %s" % r
    except Exception, e:
        r = "error: %s" \
            % traceback.format_exception(sys.exc_type,
                                         sys.exc_value,
                                         sys.exc_traceback)
        print >> sys.stderr, r

    sys.stdout.flush()
#@-node:<<main>>
#@nl
#@nonl
#@-node:@file gendata.py
#@-leo
