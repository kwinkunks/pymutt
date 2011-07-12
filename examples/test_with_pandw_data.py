#!/usr/bin/env python
'''
usage: python test_with_pandw_data.py [options] [datasetname]
    where
        datasetname: one of willamette, oceanwave, ar4, sinusoid [willamette]
        options:
            -v              be more verbose
            -h              print this
            -W f            window width [4.0]
            -K n            number of tapers [6]
            -k n            kind: 1 -> highres, 2 -> adaptive [2]
            -l n            lines to find [2]
            -z n            padding multiplier [8]
'''
import sys, math

import numpy as nx
import numpy.random as nxrn
import numpy.linalg as linalg
import pylab as mpl

import pymutt


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

class Error(Exception):
    def __init__(self, msg):
        self.msg = msg

def getseries(which):

    if which == 'willamette':
        import willamettedata as wd
        return [nx.array(wd.data), wd.dt, wd.dtunits]

    if which == 'oceanwave':
        import oceanwavedata as od
        return [nx.array(od.data), od.dt, od.dtunits]

    if which == 'ar4':
        import ar4data as od
        return [nx.array(od.data), od.dt, od.dtunits]

    if which == "sinusoid":
        n = 128
        dt = 0.1
        omega = 2.0 * math.pi * 1.0
        cc = 0.1
        sc = 0.8
        eta = 1.0
        v = cc * nx.cos(omega * dt * nx.arange(n)) \
            + sc * nx.sin(omega * dt * nx.arange(n)) \
            + eta * nxrn.normal(loc = 0.0, scale = 1.0, size = n)
        return [v, dt, "decisecond"]

    raise Error("no such dataset: %s" % which)


def plotseries(y, dx = 1.0, xunits = "", title = None):

    mpl.figure()
    mpl.plot(dx * nx.arange(len(y)), y)
    mpl.xlabel(xunits)
    if title is not None:
        mpl.title(title)


def mtaddlines(r):
    if r.has_key('linea'):
        mpl.plot(r['linef'], 20.0 * nx.log10(2.0 * abs(r['linea'])), 'D')

def mtpanel(r, title = None):

    plkeys = ('power', 'F', 'reshaped')
    npl = 0
    for k in plkeys:
        npl += r.has_key(k)

    if npl == 0:
        print >> sys.stderr, "skipping mt panel"
        return

    flen = len(r['power'])
    f = r['df'] * nx.arange(flen)

    mpl.figure()

    pl = 100 * npl + 10 + 1

    mpl.subplot(pl)
    mpl.plot(f, 10.0 * nx.log10(r['power']))
    mpl.ylabel('power (db)')
    mtaddlines(r)

    if title is not None:
        mpl.title(title)

    pl += 1

    if r.has_key('F'):
        mpl.subplot(pl)
        mpl.plot(f, r['F'])
        mpl.ylabel('F statistic')
        if r.has_key('dof'):
            alpha = 1.0 / r['nin']
            dof = r['dof']
            Fth = 0.5 * dof * (alpha ** (-2.0 / dof) - 1.0)
            mpl.plot(f, Fth)
        pl += 1

    if r.has_key('reshaped'):
        mpl.subplot(pl)
        mpl.plot(f, 10.0 * nx.log10(r['reshaped']))
        mpl.ylabel('power (db)')
        mtaddlines(r)

    mpl.xlabel("frequency")


def mtpowerpanel(r, title = None):

    flen = len(r['power'])
    f = r['df'] * nx.arange(flen)

    mpl.figure()
    mpl.plot(f, 10.0 * nx.log10(r['power']))
    mpl.ylabel('power (db)')
    mpl.xlabel("frequency")
    if title is not None:
        mpl.title(title)


def mtshortpanel(r, title = None):

    flen = len(r['power'])
    f = r['df'] * nx.arange(flen)

    mpl.figure()

    mpl.subplot(311)
    mpl.plot(f, 10.0 * nx.log10(r['power']))
    mpl.ylabel('power (db)')

    if title is not None:
        mpl.title(title)

    mpl.subplot(312)
    mpl.plot(f, r['dof'])
    mpl.ylabel('nu')

    mpl.subplot(313)
    mpl.plot(f[1:], r['F'][1:])
    mpl.ylabel('F')

    mpl.xlabel("frequency")


def mtlampspanel(r, title = None):

    if not r.has_key('linea'):
        print >> sys.stderr, "skipping lamps panel"
        return

    mpl.figure()

    mpl.subplot(211)
    mpl.plot(r['linef'], 20.0 * nx.log10(abs(r['linea'])), 'D')
    mpl.ylabel('line power (db)')

    if title is not None:
        mpl.title(title)

    flen = len(r['reshaped'])
    f = r['df'] * nx.arange(flen)

    mpl.subplot(212)
    mpl.plot(f, 10.0 * nx.log10(r['reshaped']))
    mpl.plot(f, 10.0 * nx.log10(r['power']))
    mpl.ylabel('power (db)')

    mpl.xlabel("frequency")


def mtweightspanel(r, title = None):


    a = r['weights']
    nw, nf = a.shape
    f = r['df'] * nx.arange(nf)

    mpl.figure()
    for k in range(nw):
        mpl.subplot(100 * nw + 10 + k + 1)
        mpl.plot(f, nw * a[k, :])
        mpl.ylabel('N * bk')
        if k == 0 and title is not None:
            mpl.title(title)
                
    mpl.xlabel("frequency")


def show_contents(r, header = ""):
    
    arykeys = ('power', 'dof', 'F', 'reshaped', 'weights')
    print >> sys.stderr, "keys for scalars: %s" % header
    for k in r.keys():
        if k not in arykeys:
            print >> sys.stderr, "  %10s  %s" % (k, r[k])
    print >> sys.stderr, "keys for sequences: %s" % header
    for k in arykeys:
        if r.has_key(k):
            print >> sys.stderr, "  %10s  %d in [%s, %s]" \
                  % (k, len(r[k].ravel()),
                     min(r[k].ravel()), max(r[k].ravel()))

def Fthresh(nu, pct):
    alpha = 1.0 - 0.01 * pct
    return nu * (1.0 - alpha**(2.0 / nu)) / (2.0 * alpha**(2.0 / nu))


def main(argv = None):
    import getopt
    import traceback

    if argv is None:
        argv = sys.argv

    try:

        npi = 4.0
        nwin = 6
        kind = 2
        nlines_to_find = 2
        padfactor = 8

        arglist = "vhW:K:k:l:z:"
        try:
            opts, args = getopt.getopt(argv[1:], arglist)
        except getopt.error, msg:
            raise Usage(msg)

        for opt, v in opts:
            if opt == "-v":
                pass
            elif opt == "-h":
                print >> sys.stdout, __doc__
                return 0
            elif opt == "-W":
                npi = float(v)
            elif opt == "-K":
                nwin = int(v)
            elif opt == "-k":
                kind = int(v)
            elif opt == "-z":
                padfactor = int(v)
            elif opt == "-l":
                nlines_to_find = int(v)
    
        lines = []
        dodof = 1
        doweights = 1
        seriesname = args[0] if len(args) > 0  else "willamette"

        ts, dt, dtu = getseries(seriesname)
        serieslength = len(ts)
        paddedlen = 2
        while paddedlen < serieslength:
            paddedlen *= 2
        paddedlen *= padfactor

        print >> sys.stdout, "dataset %s has %d (%d) samples at %s %s" \
              % (seriesname, len(ts), paddedlen, dt, dtu)
        print >> sys.stdout, "  npi: %.1f   K: %d   lines: %d   %s" \
              % (npi, nwin, nlines_to_find, ("?", "high-res", "adaptive")[kind])

        ts -= nx.mean(ts)

        plotseries(ts, dx = dt, xunits = dtu, title = seriesname)

        r = pymutt.mtft(series = ts,
                        dt = dt,
                        npi = npi,
                        nwin = nwin,
                        kind = kind,
                        paddedlen = paddedlen,
                        dodof = dodof,
                        doweights = doweights,
                        lines = nx.array(lines),
                        )

        show_contents(r, "original")
        basetitle = seriesname + (" W %.1f  K %d  a %d  n %d (%d)"
                         % (r['npi'], r['nwin'], r['kind'], len(ts), r['n']))
        mtpanel(r, title = basetitle + ": original")

        tlen = len(ts)
        t = dt * nx.arange(tlen)
        A = nx.zeros((tlen, 2 * nlines_to_find + 1))
        A[:,0] = 1
        flines = []

        Fcpy = nx.array(r['F'], copy = 1)

        for l in range(nlines_to_find):

            idxfmax = nx.argsort(Fcpy)[-1]
            fline = r['df'] * idxfmax
            flines.append(fline)

            ir = max(0, int(idxfmax + round(r['W'] / r['df'])))
            il = min(r['n'] - 1, int(idxfmax - round(r['W'] / r['df'])))
            Fcpy[il:ir] = 0.0
    

        reducedr = pymutt.mtft(series = ts,
                               dt = dt,
                               npi = npi,
                               nwin = nwin,
                               kind = kind,
                               paddedlen = paddedlen,
                               dodof = dodof,
                               lines = nx.array(flines))
        show_contents(reducedr, "reshaped")

        if nlines_to_find > 0:
            if 0:
                plotseries(ts, dx = dt, xunits = dtu,
                           title = seriesname + ": reshaped")
            mtpanel(reducedr, title = basetitle + ": reshaped")

        mpl.show()
        return 0

    except Usage, err:
        print >> sys.stderr, "Usage: %s" % err.msg
        print >> sys.stderr, "use -h for help"
        return 1

    except Error, e:
        print >> sys.stderr, "Error: %s" % e.msg
        return 2

    except Exception, e:
        traceback.print_exc()
        return 3

if __name__ == "__main__":
    sys.exit(main())
