The code below consists of a series of short segments that exercise the major features of the pymutt module.  The first chunks defines a dataset and remaining chunks process it in different ways.  Also included are bits of simple plotting code, mostly commented out.


```
import numpy as np
import pymutt
try:
    import matplotlib.pyplot as mpl
except:
    mpl = None

#
# standard dataset: sine wave at 11.5 Hz plus noise (but it could be anything)
#
dt = 0.01
n = 430
t = dt * np.arange(n)
data = np.random.randn(n) + 0.6 * np.sin(2.0 * 3.141592654 * 11.5 * t)
#mpl.plot(data)
#mpl.figure()


#
# simplest processing: assume dt = 1 so frequency domain is [0, 0.5]
#
r = pymutt.mtft(data)
# mpl.plot(r['power'])


#
# next simplest: use correct dt so frequency has true domain and use
#                correct frequency in the plot
#
r = pymutt.mtft(data,
                dt = dt
                )
f = np.arange(r['nspec']) * r['df']
#mpl.plot(f, r['power'])


#
# next: set the padded length, the resolution parameters
#
r = pymutt.mtft(data,
                dt = dt,
                npi = 4.0,
                nwin = 7,
                paddedlen = 8192,
                )
f = np.arange(r['nspec']) * r['df']
# mpl.plot(f, r['power'])


#
# next: look for spectral lines
#
r = pymutt.mtft(data,
                dt = dt,
                npi = 4.0,
                nwin = 7,
                paddedlen = 8192,
                dodof = 1,
                )
f = np.arange(r['nspec']) * r['df']
idxfmax = np.argsort(r['F'])[-1]
fline = r['df'] * idxfmax
print "first spectral line at f = %.3f" % fline

#mpl.subplot(211)
#mpl.plot(f, r['power'])
#mpl.subplot(212)
#mpl.plot(f, r['F'])
#mpl.figure()
#
# now compute the reshaped spectrum and plot before and after
#
r2 = pymutt.mtft(data,
                dt = dt,
                npi = 4.0,
                nwin = 7,
                paddedlen = 8192,
                lines = np.array([fline,]),
                )

#mpl.plot(f, r['power'])
#mpl.plot(f, r2['reshaped'])

if mpl:
    mpl.show()

```