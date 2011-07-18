'''
     WILLAMETTE RIVER TIME SERIES
     SOURCE: U. S. GEOLOGICAL SURVEY
     DELTA T: 1 MONTH
     SAMPLE SIZE: 395
'''
import numpy as np
import pylab as mpl
import pymutt

dt = 1.0 / 12.0
dtunits = "year"
data = [
8.95361,
9.49409,
10.1943,
10.9566,
11.0777,
10.9817,
10.4897,
10.2754,
10.1662,
9.23484,
8.56537,
8.51064,
10.1211,
11.1849,
10.9352,
11.1803,
10.9572,
10.5683,
10.1404,
9.99542,
9.14185,
8.61151,
8.38892,
8.38799,
9.97498,
10.268,
11.051,
10.5171,
10.9691,
10.4305,
10.4104,
10.0797,
9.58923,
9.40861,
8.58714,
8.53241,
8.40138,
8.29538,
9.69605,
11.4099,
11.2158,
10.417,
10.0469,
10.3744,
10.1345,
9.14213,
8.69693,
9.0466,
9.87368,
11.4195,
10.782,
11.2071,
10.2734,
10.2904,
9.49902,
9.64484,
9.20419,
8.79096,
9.01315,
9.22128,
9.58608,
9.80404,
10.5156,
10.024,
9.96579,
10.7536,
10.2445,
10.0957,
9.38059,
8.83684,
8.87995,
9.70587,
10.4744,
11.503,
11.5894,
10.6264,
10.6584,
10.5401,
10.2864,
10.0068,
9.29737,
8.78104,
8.83459,
8.91153,
10.1838,
10.641,
10.0331,
9.88032,
11.2777,
10.4286,
9.78695,
9.45786,
8.81086,
8.65054,
8.80213,
9.1636,
9.39508,
10.46,
11.1003,
11.1673,
10.4914,
10.0914,
10.0245,
9.52539,
9.02301,
8.73272,
8.82439,
9.07249,
10.437,
10.3871,
10.8521,
10.913,
9.93564,
10.2633,
10.0245,
9.45103,
8.81036,
8.61567,
8.89151,
9.70648,
9.59803,
9.7055,
9.63761,
10.6818,
10.6586,
10.7008,
10.4167,
10.1996,
8.88604,
8.6455,
8.79051,
9.01514,
10.007,
10.9098,
10.2848,
11.1658,
11.1648,
10.3303,
10.0887,
9.54927,
8.84145,
8.68172,
8.78352,
9.01477,
9.56618,
10.8817,
10.7,
9.95657,
9.98688,
10.7029,
10.2483,
9.71987,
8.95066,
8.7634,
8.78263,
9.68751,
9.87777,
10.8209,
9.56413,
10.4624,
9.88511,
10.7228,
10.6574,
9.38429,
9.05626,
8.7505,
8.753,
9.15518,
10.3805,
10.2445,
10.8118,
10.8408,
10.2186,
9.93273,
9.76269,
9.99321,
9.28532,
8.86512,
8.82709,
9.16321,
9.30811,
11.0251,
11.7443,
11.22,
9.68413,
9.43131,
9.62201,
9.16605,
8.67197,
8.6422,
8.70123,
9.23354,
9.51295,
9.61525,
11.2013,
9.96707,
10.4793,
10.1727,
9.55935,
9.08441,
8.5973,
8.58352,
8.84072,
9.09978,
10.1174,
10.8026,
10.3599,
10.8238,
9.86989,
9.72548,
9.59583,
9.30527,
8.83127,
8.60183,
8.80634,
9.15617,
10.0035,
10.1266,
10.4927,
10.4235,
10.4886,
9.65358,
9.11574,
9.3545,
8.75826,
8.79066,
9.26137,
9.74722,
10.7272,
11.101,
11.1263,
10.4894,
9.93148,
9.91216,
9.94062,
9.80192,
9.66781,
8.91826,
9.3194,
9.44833,
9.72328,
10.0364,
10.9488,
11.3657,
9.97863,
9.46072,
9.94878,
9.24589,
8.8927,
8.83407,
9.09087,
9.21523,
9.88821,
10.8848,
10.8903,
11.0465,
10.5911,
10.6297,
10.1305,
9.92123,
9.49897,
9.03218,
9.42769,
9.56236,
10.0218,
11.2228,
10.714,
11.1722,
11.2966,
10.5697,
10.13,
10.0642,
9.02958,
8.95316,
9.27128,
9.35669,
9.62552,
9.69233,
10.8822,
9.89604,
9.59268,
9.56941,
9.2226,
8.78973,
8.74666,
8.70755,
8.91549,
9.29437,
10.6784,
11.3277,
11.3213,
10.995,
10.9183,
10.7554,
9.97161,
10.0842,
9.41925,
8.90048,
9.05455,
9.51439,
9.50853,
10.228,
11.0821,
10.6595,
10.4136,
10.3892,
10.0589,
9.75097,
9.22638,
8.88871,
9.2838,
9.48822,
10.3356,
10.9635,
11.1282,
10.2808,
10.2653,
10.2091,
10.3746,
9.46273,
9.13572,
8.90775,
9.22453,
9.36962,
9.47592,
9.06025,
8.8114,
8.54531,
9.60221,
9.3758,
9.42163,
9.51944,
8.73708,
8.7092,
8.96477,
9.45995,
9.85825,
11.3312,
11.108,
10.4693,
9.51782,
9.27974,
9.7181,
9.23711,
8.81247,
8.71841,
9.28461,
9.50457,
9.33256,
10.4049,
9.96525,
10.3933,
10.578,
10.0307,
10.2466,
9.30782,
8.82286,
8.75935,
8.98145,
9.18741,
9.86568,
10.5189,
10.7585,
10.3574,
10.0056,
9.96115,
9.62088,
9.31194,
8.95841,
8.72567,
8.98905,
9.27822,
9.6282,
10.6132,
10.7569,
9.64776,
10.1544,
9.90799,
9.53015,
9.91901,
9.32266,
8.84974,
9.05606,
9.59628,
9.75846,
11.1339,
11.1351,
10.8419,
11.0652,
10.1313,
9.94447,
9.54658,
9.17291,
8.89762,
9.06933,
]

def doit():
    '''This code reproduces figure 512 in Percival and Walden.  It (1)
    analyzes the time series and applies the F-test for spectral
    lines, (2) extracts the two lines with largest F-test values, and
    (3) computes the reshaped spectrum.  It then displays the F-test
    and the original and reshaped spectra as in P&W figure 512.
    '''

    print doit.__doc__

    global dt, data

    ts = np.array(data)
    ts -= ts.mean()
    paddedlen = 1024
    r = pymutt.mtft(series = ts,
                    dt = dt,
                    npi = 4,
                    nwin = 5,
                    kind = 2,
                    paddedlen = paddedlen,
                    dodof = 1,
                    )

    Fcpy = np.array(r['F'], copy = 1)
    flines = []
    nlines_to_find = 2
    # For each line, find the maximum current value of the F-test.
    # Then zero the F-test values in the region fline +- W/2 before
    # searching for the next largest value.
    for l in range(nlines_to_find):
        idxfmax = np.argsort(Fcpy)[-1]
        fline = r['df'] * idxfmax
        flines.append(fline)
        ir = max(0, int(idxfmax + round(0.5 * r['W'] / r['df'])))
        il = min(r['n'] - 1, int(idxfmax - round(0.5 * r['W'] / r['df'])))
        Fcpy[il:ir] = 0.0

    reducedr = pymutt.mtft(series = ts,
                           dt = dt,
                           npi = 4,
                           nwin = 5,
                           kind = 2,
                           paddedlen = paddedlen,
                           dodof = 1,
                           lines = np.array(flines),
                           )

    f = r['df'] * np.arange(len(r['power']))
    mpl.subplot(211)
    mpl.title("Compare with Percival and Walden, figure 512")
    mpl.plot(f, r['F'])
    for level in (8.6, 18.5):
        y = level + 0.0 * r['power']
        mpl.plot(f, y)
    mpl.ylabel("F-test")
    mpl.subplot(212)
    mpl.plot(f, 10.0 * np.log10(r['power']))
    mpl.plot(f, 10.0 * np.log10(reducedr['reshaped']))
    mpl.xlabel("f (cycles/year)")
    mpl.ylabel("dB")
    mpl.show()
