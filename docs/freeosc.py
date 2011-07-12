'''
Suck up the compressed PFO data and plot the spectrum.
'''
import sys
import gzip
import pymutt
import pylab
import numpy as np

fn = "pfochile8.5.psn.txt.gz"

sps = None
with gzip.open(fn, "rb") as f:
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

print >> sys.stderr, "read in %d values at %.1f sps" % (len(a), sps)

# start with a target Nyquist frequency
fNY = 0.010
newsps = 2.0 * fNY
decim = int(round(sps / newsps))
newsps = sps / decim
newsamples = int(len(a)*newsps/sps)

if 0:
    # stepping sum as a lazy low-pass filter
    step = decim / 2
    m = []
    for i in range(0, len(a) - decim, step):
        m.append(a[i:(i + decim)].sum())
        a = np.array(m)
else:
    # scipy.signal fourier-domain resampler
    import scipy.signal as sig
    a = sig.resample(a, newsamples)

print >> sys.stderr, "resampled to %d values at %.3e sps" % (len(a), newsps)

# select a padded length that is a power of 2
n2 = 2
while n2 < len(a):
    n2 *= 2
paddedlen = 16 * n2

print >> sys.stderr, "padded to %d" % paddedlen

r = pymutt.mtft(a,
                dt = 1.0/newsps,
                kind = 1,
                paddedlen =  paddedlen,
                dodof = 1,
                )
freq = 1000.0 * r['df'] * np.arange(len(r['power']))
pylab.plot(freq, np.log10(r['power']))
pylab.show()
