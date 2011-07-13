'''
Suck up the compressed PFO data and plot the spectrum.
'''
import sys
import pymutt
import pylab
import numpy as np

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

print >> sys.stderr, "read in %d values at %e sps" % (len(a), sps)

# select a padded length that is a power of 2
n2 = 2
while n2 < len(a):
    n2 *= 2
paddedlen = 16 * n2

print >> sys.stderr, "padded to %d" % paddedlen

r = pymutt.mtft(a,
                dt = 1.0/sps,
                kind = 1,
                paddedlen =  paddedlen,
                dodof = 1,
                )
freq = 1000.0 * r['df'] * np.arange(len(r['power']))
pylab.plot(freq, np.log10(r['power']))
pylab.show()
