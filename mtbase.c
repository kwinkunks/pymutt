#ifndef MT_STATIC
#define DO_MTAP_SPEC_TYPE static
#else
#define DO_MTAP_SPEC_TYPE
#endif

// #define USE_FFTPACK


static int mt_verbose = 0; // some debugging output

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h>
#ifndef USE_FFTPACK
#include <fftw3.h>
#endif

#include "mtbase.h"

#define DIAG1 0

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))

#include "jtinvit.inc"


    //
    // local r2c dft interface: only route to underlying dft library
    //    - handles zero-padding and some temporary allocation
    //    - input array is not altered
    //

#ifdef USE_FFTPACK

#include "static_fftpack.inc"

    //    numpy fftpack interface

//void rffti(int n, double* wsave);            // fftpack init function
//void rfftf(int n, double* r, double* wsave); // fftpack r fwd

static double* fftpack_work = 0;
static int     fftpack_len = -1;

static void r2cdft(const double *in, int in_len, int pad_len,
                   double complex *out)
{

        //   in_len is length of the input series (in)
        //   pad_len is the length of the zero-padded series passed
        //       to the dft.  if pad_len <= in_len, it is ignored.
        //   out holds the dft computed with the fftw3 convention:
        //       a simple sum over data * exp(- omega ...).
        //   the out array has (n/2) + 1 elements,
        //       where n = max(in_len, pad_len)

    int finallen = (pad_len > in_len) ? pad_len : in_len;

    if (fftpack_len != finallen) {
        if(fftpack_work != 0)
            free(fftpack_work);
        fftpack_len = finallen;
        fftpack_work = malloc(sizeof(double) * (2 * finallen + 15));
        rffti(finallen, fftpack_work);
    }

    double* dptr = (double *) out;
    int outlength = (finallen / 2) + 1;
    memset(dptr, 0, 2 * outlength * sizeof(in[0]));
    memcpy(dptr + 1, in, in_len * sizeof(in[0]));

    rfftf(finallen, dptr + 1, fftpack_work);

    dptr[0] = dptr[1];
    dptr[1] = 0.0;

}

#else

    //    fftw3 interface

static void r2cdft(const double *in, int in_len, int pad_len,
                   double complex *out)
{

        //   in_len is length of the input series (in)
        //   pad_len is the length of the zero-padded series passed
        //       to the dft.  if pad_len <= in_len, it is ignored.
        //   out holds the dft computed with the fftw3 convention:
        //       a simple sum over data * exp(- omega ...).
        //   the out array has (n/2) + 1 elements,
        //       where n = max(in_len, pad_len)

    int finallen = (pad_len > in_len) ? pad_len : in_len;
    int nc = (finallen / 2) + 1;
    double* rin = fftw_malloc(sizeof(double) * finallen);
    double complex* cout = fftw_malloc(sizeof(double complex) * nc);
    fftw_plan p = fftw_plan_dft_r2c_1d(finallen, rin, cout,
                                FFTW_ESTIMATE);

    memcpy(rin, in, in_len * sizeof(in[0]));
    if(finallen > in_len)
        memset(rin + in_len, 0, (finallen - in_len) * sizeof(rin[0]));

    fftw_execute(p);
    memcpy(out, cout, nc * sizeof(double complex));

    fftw_destroy_plan(p);
    fftw_free(rin);
    fftw_free(cout);
}

#endif

/* 
** currently returns tapers normalized so that h * h = 1
*/ 

static void mt_multitap(int num_points, int nwin, double *lam,
                        double npi, double *tapers, double *tapsum)
{
        /*
         * get the multitaper slepian functions:
         * 
         *    num_points = number of points in data stream
         *    nwin = number of windows
         *    lam= vector of eigenvalues
         *    npi = order of slepian functions
         *    tapsum = sum of each taper,
         *        saved for use in adaptive weighting
         *    tapers =  matrix of slepian tapers,
         *        packed in a 1D double array
         */

    int             i, k, kk;
    double          ww, cs, ai, an, eps, rlu = 0.0, rlb = 0.0;
    double          dfac, drat, gamma, bh, TWOPI, DPI;
    double         *diag, *offdiag, *offsq;

    double         *scratch1, *scratch2, *scratch3, *scratch4, *scratch6;

    double tapsq;
    double aa;

        /* need to initialize iwflag = 0 */

    double          anpi;
    double         *ell;
    int            *ip;
    double         *evecs;

    long            len;
    int             ierr;
    int             m11;
    DPI = (double) M_PI;
    TWOPI = (double) 2 *DPI;

    anpi = npi;
    an = (double) (num_points);
    ww = (double) (anpi) / an;	/* this corresponds to P&W's W value  */
    cs = cos(TWOPI * ww);


    ell = (double *) malloc((size_t) nwin * sizeof(double));

    diag = (double *) malloc((size_t) num_points * sizeof(double));

    offdiag = (double *) malloc((size_t) num_points * sizeof(double));
    offsq = (double *) malloc((size_t) num_points * sizeof(double));

    scratch1 = (double *) malloc((size_t) num_points * sizeof(double));
    scratch2 = (double *) malloc((size_t) num_points * sizeof(double));
    scratch3 = (double *) malloc((size_t) num_points * sizeof(double));
    scratch4 = (double *) malloc((size_t) num_points * sizeof(double));
    scratch6 = (double *) malloc((size_t) num_points * sizeof(double));

        /* make the diagonal elements of the tridiag matrix  */

    for (i = 0; i < num_points; i++) {
        ai = (double) (i);
        diag[i] = -cs * (((an - 1.) / 2. - ai)) * (((an - 1.) / 2. - ai));
        offdiag[i] = -ai * (an - ai) / 2.;
        offsq[i] = offdiag[i] * offdiag[i];
    }

    eps = 1.0e-13;
    m11 = 1;

    ip = (int *) malloc((size_t) nwin * sizeof(int));

        /* call the eispac routines to invert the tridiagonal system */

    jtridib_(&num_points, &eps, diag, offdiag, offsq, &rlb,
             &rlu, &m11, &nwin, lam,
             ip, &ierr, scratch1, scratch2);
#if DIAG1
    fprintf(stderr, "ierr=%d rlb=%.8f rlu=%.8f\n", ierr, rlb, rlu);

    fprintf(stderr, "eigenvalues for the eigentapers\n");

    for (k = 0; k < nwin; k++)
        fprintf(stderr, "%.20f ", lam[k]);
    fprintf(stderr, "\n");
#endif


    len = num_points * nwin;

    evecs = (double *) malloc((size_t) len * sizeof(double));



    jtinvit_(&num_points, &num_points, diag, offdiag, offsq,
             &nwin, lam, ip, evecs, &ierr,
             scratch1, scratch2, scratch3, scratch4, scratch6);

    free(scratch1);
    free(scratch2);
    free(scratch3);
    free(scratch4);
    free(scratch6);

        /*
         * we calculate the eigenvalues of the dirichlet-kernel problem i.e.
         * the bandwidth retention factors from slepian 1978 asymptotic
         * formula, gotten from thomson 1982 eq 2.5 supplemented by the
         * asymptotic formula for k near 2n from slepian 1978 eq 61 more
         * precise values of these parameters, perhaps useful in adaptive
         * spectral estimation, can be calculated explicitly using the
         * rayleigh-quotient formulas in thomson (1982) and park et al (1987)
         * 
         */

    dfac = (double) an *DPI * ww;
    drat = (double) 8. *dfac;


    dfac = (double) 4. *sqrt(DPI * dfac) * exp((double) (-2.0) * dfac);


    for (k = 0; k < nwin; k++) {
        lam[k] = (double) 1.0 - (double) dfac;
        dfac = dfac * drat / (double) (k + 1);

            /* fails as k -> 2n */

    }


    gamma = log((double) 8. * an * sin((double) 2. * DPI * ww))
        + (double) 0.5772156649;



    for (k = 0; k < nwin; k++) {
        bh = -2. * DPI * (an * ww - (double) (k) /
                          (double) 2. - (double) .25) / gamma;
        ell[k] = (double) 1. / ((double) 1. + exp(DPI * (double) bh));

    }

    for (i = 0; i < nwin; i++)
        lam[i] = max(ell[i], lam[i]);

        /*
        ** normalize the eigentapers so they have an inner product of unity
        ** tapsum is the average of the eigentaper, should be near zero for
        ** antisymmetric tapers
        */
    for (k = 0; k < nwin; k++) {
        kk = k * num_points;
        tapsum[k] = 0.;
        tapsq = 0.0;
        for (i = 0; i < num_points; i++) {
            double v = evecs[i + kk];
            tapsum[k] = tapsum[k] + v;
            tapsq = tapsq + v * v;
        }

        aa = 1.0 / sqrt(tapsq);

        for (i = 0; i < num_points; i++)
            tapers[i + kk] = aa * evecs[i + kk];

        tapsum[k] = aa * tapsum[k];

    }

        /* Free Memory */

    free(ell);
    free(diag);
    free(offdiag);
    free(offsq);
    free(ip);

    free(evecs);
}

// hires estimate: P&W 369a

static void hires(double *sqr_spec, double *el, int nwin,
                  int num_freq, double *ares)
{
    int i, j, w;
    double wsum;

    memset(ares, 0, num_freq * sizeof(ares[0]));

    if(nwin < 1)
        return;

    wsum = 0.0;

    for (w = 0; w < nwin; w++) {
        double* off_spec = &sqr_spec[0 + w * num_freq];
        for (j = 0; j < num_freq; j++)
            ares[j] = ares[j] + el[w] * off_spec[j];
        wsum += el[w];
    }

    for(i = 0; i < num_freq; i++)
        ares[i] /= wsum;

    return;
}

// thomson's algorithm for the adaptive spectrum estimate
// P&W 370a ...

static int adweight(double *sqr_spec, double *dcf,
                    double *el, int nwin, int num_freq,
                    double *ares, double *degf, double avar)
{

    double as, das, tol, a1, ax, fn, fx;
    double *spw, *bias;
    double test_tol, dif;
    int jitter, i,  k, kpoint, ifreq;
    double df, dfs;

        /*
         * set tolerance for iterative scheme exit
         */

    tol = 3.0e-4;
    jitter = 0;

    spw = (double *) malloc((size_t) nwin * sizeof(double));
    bias = (double *) malloc((size_t) nwin * sizeof(double));


    for (i = 0; i < nwin; i++)
        bias[i] = (1.00 - el[i]) * avar;

    for (ifreq = 0; ifreq < num_freq; ifreq++) {

        for (i = 0; i < nwin; i++)
            spw[i] = sqr_spec[ifreq + i * num_freq];

            /*
            ** first guess is the average of the two
            ** lowest-order eigenspectral estimates
            */

        as = (spw[0] + spw[1]) / 2.00;

        for (k = 0; k < 20; k++) {

            fn = fx = 0.00;

            for (i = 0; i < nwin; i++) {
                a1 = sqrt(el[i]) * as / (el[i] * as + bias[i]);
                a1 = a1 * a1;
                fn = fn + a1 * spw[i];
                fx = fx + a1;
            }
            ax = fn / fx;
            dif = ax - as;
            das = fabs(dif);
            as = ax;
            test_tol = das / as;
            if (test_tol < tol)
                break;
        }

        if (k >= 20)
            jitter++;

        ares[ifreq] = as;

        df = 0.0;
        dfs = 0.0;

        for (i = 0; i < nwin; i++) {
            kpoint = ifreq + i * num_freq;
            dcf[kpoint] = as / (el[i] * as + bias[i]);
            df  += dcf[kpoint] * dcf[kpoint] * el[i];
            dfs += dcf[kpoint] * dcf[kpoint] * dcf[kpoint] * dcf[kpoint];
        }
        degf[ifreq] = 2.0 * df * df / dfs;

        dfs = 0.0;
        for(i = 0; i < nwin; i++) {
            kpoint = ifreq + i * num_freq;
            dfs += el[i] * dcf[kpoint] * dcf[kpoint];
        }
        for(i = 0; i < nwin; i++) {
            kpoint = ifreq + i * num_freq;
            dcf[kpoint] = dcf[kpoint] * dcf[kpoint] * el[i] / dfs;
        }

    }

    if(mt_verbose && (jitter > 0))
        fprintf(stderr, "%d failed iterations\n", jitter);
    free(spw);
    free(bias);

    return jitter;
}

static void get_F_values(complex double *s, int nf, int nwin,
                         double *Fvalue, double *b)
{
        /*
         * b is fft of slepian eigentapers at zero freq, s are the
         * eigenspectra, amu contains line frequency estimates and f-test
         * parameter
         */
    double sum, sum2;
    complex double z;
    int i, j;
    complex double *amu;

    amu = (double complex *) malloc((size_t) nf * sizeof(double complex));

    sum = 0.0;
    for (i = 0; i < nwin; i++)
        sum = sum + b[i] * b[i];

    for (i = 0; i < nf; i++) {
        amu[i] = 0.;
        for (j = 0; j < nwin; j++)
            amu[i] += s[i + j * nf] * b[j];
        amu[i] /= sum;

        sum2 = 0.;
        for (j = 0; j < nwin; j++) {
            z = s[i + j * nf] - amu[i] * b[j];
            sum2 += z * conj(z);
        }
        Fvalue[i] = (double) (nwin - 1) * (amu[i] * conj(amu[i])) * sum / sum2;
    }
    free(amu);
    return;
}

DO_MTAP_SPEC_TYPE
void do_mtap_spec(double* data, int truelen, double dt, int kind,
                  int nwin, double npi, int paddedlen,
                  double* mtpwr, double* dof, double* Fvalues,
                  int nlines,
                  double* lines,
                  double complex* lamp,
                  double* flamp,
                  double* reshaped,
                  double* linevar,
                  double* tweights,
                  double complex* tspectra
                  )
{
        /*
        **    input
        **    
        ** data:       time series
        ** truelen:    number of points in data
        ** dt:         sample interval
        ** kind:       1 ==> high-res, 2 ==> adaptive
        ** nwin:       taper count
        ** npi:        slepian function order
        ** paddedlen:  padded series length to use
        ** nlines      number of lines; >= 0
        ** lines       array of line frequencies; may be 0 if nlines == 0
        **
        **    output   -- there are (paddedlen / 2) + 1 frequencies
        **    
        ** mtpwr:      power(f)
        ** dof:        dof(f); may be 0
        ** Fvalues:    F(f); may be 0
        ** lamp:       line spectral amps(f); may be 0
        ** flamp:      actual frequencies used for lines; may be 0
        ** reshaped    reshaped spectrum; may be 0
        ++ linevar     variance of line amps; may be 0
        ** tweights:   weights(f, win); may be 0
        ** tspectra:   spectra(f, win); may be 0
        */

    int i, j, k, l;
    int iwin, kk;
    int kf;

    double* b;
    double* amu;
    double* sqr_spec;

    double complex *amp = 0;
    double complex *Jk = 0;
    double complex *Hk = 0;

    double* Fv = 0;

    double *dcf, *degf;
    double avar, amean;
    double lsum;

    double dtsqrt = sqrt(dt);
    double df = 1.0 / (paddedlen * dt);

        /*
        ** lambda = vector of eigenvalues
        ** tapsum = sum of each taper, saved for use in adaptive weighting
        ** tapers =  matrix of slepian tapers, packed in a 1D double array
        */

    double* lambda = (double *) malloc((size_t) nwin * sizeof(double));
    double* tapsum = (double *) malloc((size_t) nwin * sizeof(double));

    int len_taps = truelen * nwin;
    double* tapers = (double *) malloc((size_t) len_taps * sizeof(double));

    int num_freqs = paddedlen / 2 + 1;
    int num_freq_tap = num_freqs * nwin;


    dcf = (double *) malloc((size_t) num_freq_tap * sizeof(double));
    degf = (double *) malloc((size_t) num_freqs * sizeof(double));

    if(reshaped || nlines)
        Hk = (double complex *) malloc((size_t) num_freq_tap
                                       * sizeof(double complex));


    if(paddedlen < truelen) {
        fprintf(stderr,
                "uft::do_mtap_spec  requested: %d  <  input: %d\n",
                truelen, paddedlen);
        exit(1);
    }

        /* fetch the slepian tapers */

    mt_multitap(truelen, nwin, lambda, npi, tapers, tapsum);

    if(mt_verbose)
        for(l = 0; l < nwin; l++)
            fprintf(stderr, "%-17s %4d  %9.2f\n",
                    (l == 0) ? "eigenvalues:" : "", l, lambda[l]);

    b = (double *) malloc((size_t) truelen * sizeof(double));

    amu = (double *) malloc(num_freqs * sizeof(double));
    sqr_spec = (double *) malloc(num_freq_tap * sizeof(double));

    amp = (double complex *) malloc(num_freqs * sizeof(double complex));
    Jk = (double complex *) malloc(num_freq_tap * sizeof(double complex));

        // create the family of complex spectra

    for (iwin = 0; iwin < nwin; iwin++) {
        kk = iwin * truelen;
        kf = iwin * num_freqs;

        for (j = 0; j < truelen; j++)
            b[j] = data[j] * tapers[j + kk];

        r2cdft(b, truelen, paddedlen, amp);

        for (i = 0; i < num_freqs; i++) {
            sqr_spec[i + kf] = amp[i] * conj(amp[i]);
            Jk[i + kf] = dtsqrt * amp[i];
        }

        if(Hk) {
            r2cdft(tapers + kk, truelen, paddedlen, Hk + kf);
            for (i = 0; i < num_freqs; i++)
                Hk[i + kf] = dt * Hk[i + kf];
        }
    }

    if(Fvalues)
        Fv = (double *) malloc((size_t) num_freqs * sizeof(double));

    switch (kind) {

            /* hires */

        case 1:

            hires(sqr_spec, lambda, nwin, num_freqs, amu);

            lsum = 0.0;
            for(iwin = 0; iwin < nwin; iwin++)
                lsum += lambda[iwin];

            for (i = 0; i < num_freqs; i++) {
                degf[i] = 2 * nwin;
                for(iwin = 0; iwin < nwin; iwin++)
                    dcf[i + num_freqs * iwin] = lambda[iwin] / lsum;
            }

            break;

                /* adaptive */

        case 2:

            avar = 0.0;
            amean = 0.0;
            for (i = 0; i < truelen; i++) {
                avar  += data[i] * data[i];
                amean += data[i];
            }
            amean /= truelen;
            avar = (avar / truelen) - amean * amean;

            adweight(sqr_spec, dcf, lambda, nwin, num_freqs,
                     amu, degf, avar);

                // fprintf(stderr, "adweight returned %d\n", i);

            break;

        default:
            fprintf(stderr, "bad kind value: %d\n", kind);
            exit(1);
    }

    if(Fvalues) {
        get_F_values(Jk, num_freqs, nwin, Fv, tapsum);
        memcpy(Fvalues, Fv, num_freqs * sizeof(Fv[0]));
    }

    for (i = 0; i < num_freqs; i++)
        mtpwr[i] = dt * amu[i];

    if(tweights)
        memcpy(tweights, dcf, num_freq_tap * sizeof(dcf[0]));

    if(dof)
        memcpy(dof, degf, num_freqs * sizeof(degf[0]));

    if(nlines > 0 && lines && lamp) {

            // estimate complex line amplitudes at specified frequencies
            // using P&W 499a

        double hksum = 0.0;
        for(j = 0; j < nwin; j += 2)
            hksum += Hk[0 + j * num_freqs] * Hk[0 + j * num_freqs];

        for(k = 0; k < nlines; k++) {

            double complex clsum = 0.0;

            i = (int)(lines[k] / df + 0.5);
            if(i >= num_freqs)
                i = num_freqs - 1;
            if(i < 0)
                i = 0;

            for(j = 0; j < nwin; j += 2)
                clsum += Hk[0 + j * num_freqs] * Jk[i + num_freqs * j];

            lamp[k] = dtsqrt * clsum / hksum;
            if(flamp)
                flamp[k] = i * df;

#ifdef NOT_DEFINED            
            clsum = 0.0;
            for(j = 0; j < nwin; j += 2)
                clsum += dcf[i + j * num_freqs] * Jk[i + j * num_freqs];
            lamp[k] = clsum;
#endif
        }
        
        if(reshaped || linevar) {

            // removed the specified lines and reshape the spectrum
            // P&W 500

            double* localreshaped = 0;
            if(reshaped == 0)
                localreshaped = (double *) malloc(num_freqs * sizeof(double));
            else
                localreshaped = reshaped;

            for(i = 0; i < num_freqs; i++) {
                double sigl = 0.0;
                double wsum = 0.0;
                for(iwin = 0; iwin < nwin; iwin++) {
                    int kl = i + num_freqs * iwin;
                    double complex sdf = Jk[kl];
                    double wco = dcf[kl];
                    for(k = 0; k < nlines; k++) {
                        int ic = (int)(lines[k] / df + 0.1);
                        int deltai = (i > ic) ? (i - ic) : (ic - i);
                        double complex hk = Hk[deltai + num_freqs * iwin];
                        if(i < ic)
                            hk = conj(hk);
                        sdf -= lamp[k] * hk / dtsqrt;
                    }
                    sigl += wco * (sdf * conj(sdf));
                    wsum += wco;
                }
                localreshaped[i] = sigl / wsum;
            }
            if(linevar) {
                for(k = 0; k < nlines; k++) {
                    int ic = (int)(lines[k] / df + 0.1);
                    linevar[k] = dtsqrt * localreshaped[ic] / hksum;
                }
            }
            if(reshaped == 0)
                free(localreshaped);
        }
    }

    free(dcf);
    free(degf);

    if(Fvalues)
        free(Fv);

    free(amu);
    free(sqr_spec);

    if(tspectra)
        memcpy(tspectra, Jk, num_freq_tap * sizeof(Jk[0]));
    free(Jk);

    free(amp);

    free(b);

    free(lambda);
    free(tapsum);
    free(tapers);
}
