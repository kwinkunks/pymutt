'''
Suck up the compressed PFO data and plot the spectrum.
'''
import sys
import gzip
import pymutt
import pylab
import numpy as np

fnin = "pfochile8.5.psn.txt.gz"
fnout = "chile.txt"

sps = None
with gzip.open(fnin, "rb") as f:
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

out = file(fnout, "wb")
print >> out, "SPS: %e" % newsps
for y in a:
	print >> out, "%e" % y
out.close()
