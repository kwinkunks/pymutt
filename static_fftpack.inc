/*
fftpack.c : A set of FFT routines in C.
Algorithmically based on Fortran-77 FFTPACK by Paul N. Swarztrauber (Version 4, 1985).

*/

/* isign is +1 for backward and -1 for forward transforms */

#include <math.h>
#include <stdio.h>
#define DOUBLE

#ifdef DOUBLE
#define Treal double
#else
#define Treal float
#endif


#define ref(u,a) u[a]

#define MAXFAC 13    /* maximum number of factors in factorization of n */
#define NSPECIAL 4   /* number of factors for which we have special-case routines */

#ifdef __cplusplus
extern "C" {
#endif


  /* ----------------------------------------------------------------------
radf2,radb2, radf3,radb3, radf4,radb4, radf5,radb5, radfg,radbg.
Treal FFT passes fwd and bwd.
---------------------------------------------------------------------- */

static void radf2(int ido, int l1, const Treal cc[], Treal ch[], const Treal wa1[])
  {
    int i, k, ic;
    Treal ti2, tr2;
    for (k=0; k<l1; k++) {
      ch[2*k*ido] =
          ref(cc,k*ido) + ref(cc,(k + l1)*ido);
      ch[(2*k+1)*ido + ido-1] =
          ref(cc,k*ido) - ref(cc,(k + l1)*ido);
    }
    if (ido < 2) return;
    if (ido != 2) {
      for (k=0; k<l1; k++) {
        for (i=2; i<ido; i+=2) {
          ic = ido - i;
          tr2 = wa1[i - 2]*ref(cc, i-1 + (k + l1)*ido) + wa1[i - 1]*ref(cc, i + (k + l1)*ido);
          ti2 = wa1[i - 2]*ref(cc, i + (k + l1)*ido) - wa1[i - 1]*ref(cc, i-1 + (k + l1)*ido);
          ch[i + 2*k*ido] = ref(cc,i + k*ido) + ti2;
          ch[ic + (2*k+1)*ido] = ti2 - ref(cc,i + k*ido);
          ch[i - 1 + 2*k*ido] = ref(cc,i - 1 + k*ido) + tr2;
          ch[ic - 1 + (2*k+1)*ido] = ref(cc,i - 1 + k*ido) - tr2;
        }
      }
      if (ido % 2 == 1) return;
    }
    for (k=0; k<l1; k++) {
      ch[(2*k+1)*ido] = -ref(cc,ido-1 + (k + l1)*ido);
      ch[ido-1 + 2*k*ido] = ref(cc,ido-1 + k*ido);
    }
  } /* radf2 */


static void radf3(int ido, int l1, const Treal cc[], Treal ch[],
      const Treal wa1[], const Treal wa2[])
  {
    static const Treal taur = -0.5;
    static const Treal taui = 0.866025403784439;
    int i, k, ic;
    Treal ci2, di2, di3, cr2, dr2, dr3, ti2, ti3, tr2, tr3;
    for (k=0; k<l1; k++) {
      cr2 = ref(cc,(k + l1)*ido) + ref(cc,(k + 2*l1)*ido);
      ch[3*k*ido] = ref(cc,k*ido) + cr2;
      ch[(3*k+2)*ido] = taui*(ref(cc,(k + l1*2)*ido) - ref(cc,(k + l1)*ido));
      ch[ido-1 + (3*k + 1)*ido] = ref(cc,k*ido) + taur*cr2;
    }
    if (ido == 1) return;
    for (k=0; k<l1; k++) {
      for (i=2; i<ido; i+=2) {
        ic = ido - i;
        dr2 = wa1[i - 2]*ref(cc,i - 1 + (k + l1)*ido) +
            wa1[i - 1]*ref(cc,i + (k + l1)*ido);
        di2 = wa1[i - 2]*ref(cc,i + (k + l1)*ido) - wa1[i - 1]*ref(cc,i - 1 + (k + l1)*ido);
        dr3 = wa2[i - 2]*ref(cc,i - 1 + (k + l1*2)*ido) + wa2[i - 1]*ref(cc,i + (k + l1*2)*ido);
        di3 = wa2[i - 2]*ref(cc,i + (k + l1*2)*ido) - wa2[i - 1]*ref(cc,i - 1 + (k + l1*2)*ido);
        cr2 = dr2 + dr3;
        ci2 = di2 + di3;
        ch[i - 1 + 3*k*ido] = ref(cc,i - 1 + k*ido) + cr2;
        ch[i + 3*k*ido] = ref(cc,i + k*ido) + ci2;
        tr2 = ref(cc,i - 1 + k*ido) + taur*cr2;
        ti2 = ref(cc,i + k*ido) + taur*ci2;
        tr3 = taui*(di2 - di3);
        ti3 = taui*(dr3 - dr2);
        ch[i - 1 + (3*k + 2)*ido] = tr2 + tr3;
        ch[ic - 1 + (3*k + 1)*ido] = tr2 - tr3;
        ch[i + (3*k + 2)*ido] = ti2 + ti3;
        ch[ic + (3*k + 1)*ido] = ti3 - ti2;
      }
    }
  } /* radf3 */


static void radf4(int ido, int l1, const Treal cc[], Treal ch[],
      const Treal wa1[], const Treal wa2[], const Treal wa3[])
  {
    static const Treal hsqt2 = 0.7071067811865475;
    int i, k, ic;
    Treal ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
    for (k=0; k<l1; k++) {
      tr1 = ref(cc,(k + l1)*ido) + ref(cc,(k + 3*l1)*ido);
      tr2 = ref(cc,k*ido) + ref(cc,(k + 2*l1)*ido);
      ch[4*k*ido] = tr1 + tr2;
      ch[ido-1 + (4*k + 3)*ido] = tr2 - tr1;
      ch[ido-1 + (4*k + 1)*ido] = ref(cc,k*ido) - ref(cc,(k + 2*l1)*ido);
      ch[(4*k + 2)*ido] = ref(cc,(k + 3*l1)*ido) - ref(cc,(k + l1)*ido);
    }
    if (ido < 2) return;
    if (ido != 2) {
      for (k=0; k<l1; k++) {
        for (i=2; i<ido; i += 2) {
          ic = ido - i;
          cr2 = wa1[i - 2]*ref(cc,i - 1 + (k + l1)*ido) + wa1[i - 1]*ref(cc,i + (k + l1)*ido);
          ci2 = wa1[i - 2]*ref(cc,i + (k + l1)*ido) - wa1[i - 1]*ref(cc,i - 1 + (k + l1)*ido);
          cr3 = wa2[i - 2]*ref(cc,i - 1 + (k + 2*l1)*ido) + wa2[i - 1]*ref(cc,i + (k + 2*l1)*
              ido);
          ci3 = wa2[i - 2]*ref(cc,i + (k + 2*l1)*ido) - wa2[i - 1]*ref(cc,i - 1 + (k + 2*l1)*
              ido);
          cr4 = wa3[i - 2]*ref(cc,i - 1 + (k + 3*l1)*ido) + wa3[i - 1]*ref(cc,i + (k + 3*l1)*
              ido);
          ci4 = wa3[i - 2]*ref(cc,i + (k + 3*l1)*ido) - wa3[i - 1]*ref(cc,i - 1 + (k + 3*l1)*
              ido);
          tr1 = cr2 + cr4;
          tr4 = cr4 - cr2;
          ti1 = ci2 + ci4;
          ti4 = ci2 - ci4;
          ti2 = ref(cc,i + k*ido) + ci3;
          ti3 = ref(cc,i + k*ido) - ci3;
          tr2 = ref(cc,i - 1 + k*ido) + cr3;
          tr3 = ref(cc,i - 1 + k*ido) - cr3;
          ch[i - 1 + 4*k*ido] = tr1 + tr2;
          ch[ic - 1 + (4*k + 3)*ido] = tr2 - tr1;
          ch[i + 4*k*ido] = ti1 + ti2;
          ch[ic + (4*k + 3)*ido] = ti1 - ti2;
          ch[i - 1 + (4*k + 2)*ido] = ti4 + tr3;
          ch[ic - 1 + (4*k + 1)*ido] = tr3 - ti4;
          ch[i + (4*k + 2)*ido] = tr4 + ti3;
          ch[ic + (4*k + 1)*ido] = tr4 - ti3;
        }
      }
      if (ido % 2 == 1) return;
    }
    for (k=0; k<l1; k++) {
      ti1 = -hsqt2*(ref(cc,ido-1 + (k + l1)*ido) + ref(cc,ido-1 + (k + 3*l1)*ido));
      tr1 = hsqt2*(ref(cc,ido-1 + (k + l1)*ido) - ref(cc,ido-1 + (k + 3*l1)*ido));
      ch[ido-1 + 4*k*ido] = tr1 + ref(cc,ido-1 + k*ido);
      ch[ido-1 + (4*k + 2)*ido] = ref(cc,ido-1 + k*ido) - tr1;
      ch[(4*k + 1)*ido] = ti1 - ref(cc,ido-1 + (k + 2*l1)*ido);
      ch[(4*k + 3)*ido] = ti1 + ref(cc,ido-1 + (k + 2*l1)*ido);
    }
  } /* radf4 */


static void radf5(int ido, int l1, const Treal cc[], Treal ch[],
      const Treal wa1[], const Treal wa2[], const Treal wa3[], const Treal wa4[])
  {
    static const Treal tr11 = 0.309016994374947;
    static const Treal ti11 = 0.951056516295154;
    static const Treal tr12 = -0.809016994374947;
    static const Treal ti12 = 0.587785252292473;
    int i, k, ic;
    Treal ci2, di2, ci4, ci5, di3, di4, di5, ci3, cr2, cr3, dr2, dr3, dr4, dr5,
        cr5, cr4, ti2, ti3, ti5, ti4, tr2, tr3, tr4, tr5;
    for (k = 0; k < l1; k++) {
      cr2 = ref(cc,(k + 4*l1)*ido) + ref(cc,(k + l1)*ido);
      ci5 = ref(cc,(k + 4*l1)*ido) - ref(cc,(k + l1)*ido);
      cr3 = ref(cc,(k + 3*l1)*ido) + ref(cc,(k + 2*l1)*ido);
      ci4 = ref(cc,(k + 3*l1)*ido) - ref(cc,(k + 2*l1)*ido);
      ch[5*k*ido] = ref(cc,k*ido) + cr2 + cr3;
      ch[ido-1 + (5*k + 1)*ido] = ref(cc,k*ido) + tr11*cr2 + tr12*cr3;
      ch[(5*k + 2)*ido] = ti11*ci5 + ti12*ci4;
      ch[ido-1 + (5*k + 3)*ido] = ref(cc,k*ido) + tr12*cr2 + tr11*cr3;
      ch[(5*k + 4)*ido] = ti12*ci5 - ti11*ci4;
    }
    if (ido == 1) return;
    for (k = 0; k < l1; ++k) {
      for (i = 2; i < ido; i += 2) {
        ic = ido - i;
        dr2 = wa1[i - 2]*ref(cc,i - 1 + (k + l1)*ido) + wa1[i - 1]*ref(cc,i + (k + l1)*ido);
        di2 = wa1[i - 2]*ref(cc,i + (k + l1)*ido) - wa1[i - 1]*ref(cc,i - 1 + (k + l1)*ido);
        dr3 = wa2[i - 2]*ref(cc,i - 1 + (k + 2*l1)*ido) + wa2[i - 1]*ref(cc,i + (k + 2*l1)*ido);
        di3 = wa2[i - 2]*ref(cc,i + (k + 2*l1)*ido) - wa2[i - 1]*ref(cc,i - 1 + (k + 2*l1)*ido);
        dr4 = wa3[i - 2]*ref(cc,i - 1 + (k + 3*l1)*ido) + wa3[i - 1]*ref(cc,i + (k + 3*l1)*ido);
        di4 = wa3[i - 2]*ref(cc,i + (k + 3*l1)*ido) - wa3[i - 1]*ref(cc,i - 1 + (k + 3*l1)*ido);
        dr5 = wa4[i - 2]*ref(cc,i - 1 + (k + 4*l1)*ido) + wa4[i - 1]*ref(cc,i + (k + 4*l1)*ido);
        di5 = wa4[i - 2]*ref(cc,i + (k + 4*l1)*ido) - wa4[i - 1]*ref(cc,i - 1 + (k + 4*l1)*ido);
        cr2 = dr2 + dr5;
        ci5 = dr5 - dr2;
        cr5 = di2 - di5;
        ci2 = di2 + di5;
        cr3 = dr3 + dr4;
        ci4 = dr4 - dr3;
        cr4 = di3 - di4;
        ci3 = di3 + di4;
        ch[i - 1 + 5*k*ido] = ref(cc,i - 1 + k*ido) + cr2 + cr3;
        ch[i + 5*k*ido] = ref(cc,i + k*ido) + ci2 + ci3;
        tr2 = ref(cc,i - 1 + k*ido) + tr11*cr2 + tr12*cr3;
        ti2 = ref(cc,i + k*ido) + tr11*ci2 + tr12*ci3;
        tr3 = ref(cc,i - 1 + k*ido) + tr12*cr2 + tr11*cr3;
        ti3 = ref(cc,i + k*ido) + tr12*ci2 + tr11*ci3;
        tr5 = ti11*cr5 + ti12*cr4;
        ti5 = ti11*ci5 + ti12*ci4;
        tr4 = ti12*cr5 - ti11*cr4;
        ti4 = ti12*ci5 - ti11*ci4;
        ch[i - 1 + (5*k + 2)*ido] = tr2 + tr5;
        ch[ic - 1 + (5*k + 1)*ido] = tr2 - tr5;
        ch[i + (5*k + 2)*ido] = ti2 + ti5;
        ch[ic + (5*k + 1)*ido] = ti5 - ti2;
        ch[i - 1 + (5*k + 4)*ido] = tr3 + tr4;
        ch[ic - 1 + (5*k + 3)*ido] = tr3 - tr4;
        ch[i + (5*k + 4)*ido] = ti3 + ti4;
        ch[ic + (5*k + 3)*ido] = ti4 - ti3;
      }
    }
  } /* radf5 */


static void radfg(int ido, int ip, int l1, int idl1,
      Treal cc[], Treal ch[], const Treal wa[])
  {
    static const Treal twopi = 6.28318530717959;
    int idij, ipph, i, j, k, l, j2, ic, jc, lc, ik, is, nbd;
    Treal dc2, ai1, ai2, ar1, ar2, ds2, dcp, arg, dsp, ar1h, ar2h;
    arg = twopi / ip;
    dcp = cos(arg);
    dsp = sin(arg);
    ipph = (ip + 1) / 2;
    nbd = (ido - 1) / 2;
    if (ido != 1) {
      for (ik=0; ik<idl1; ik++) ch[ik] = cc[ik];
      for (j=1; j<ip; j++)
        for (k=0; k<l1; k++)
          ch[(k + j*l1)*ido] = cc[(k + j*l1)*ido];
      if (nbd <= l1) {
        is = -ido;
        for (j=1; j<ip; j++) {
          is += ido;
          idij = is-1;
          for (i=2; i<ido; i+=2) {
            idij += 2;
            for (k=0; k<l1; k++) {
              ch[i - 1 + (k + j*l1)*ido] =
                  wa[idij - 1]*cc[i - 1 + (k + j*l1)*ido] + wa[idij]*cc[i + (k + j*l1)*ido];
              ch[i + (k + j*l1)*ido] =
                  wa[idij - 1]*cc[i + (k + j*l1)*ido] - wa[idij]*cc[i - 1 + (k + j*l1)*ido];
            }
          }
        }
      } else {
        is = -ido;
        for (j=1; j<ip; j++) {
          is += ido;
          for (k=0; k<l1; k++) {
            idij = is-1;
            for (i=2; i<ido; i+=2) {
              idij += 2;
              ch[i - 1 + (k + j*l1)*ido] =
                  wa[idij - 1]*cc[i - 1 + (k + j*l1)*ido] + wa[idij]*cc[i + (k + j*l1)*ido];
              ch[i + (k + j*l1)*ido] =
                  wa[idij - 1]*cc[i + (k + j*l1)*ido] - wa[idij]*cc[i - 1 + (k + j*l1)*ido];
            }
          }
        }
      }
      if (nbd >= l1) {
        for (j=1; j<ipph; j++) {
          jc = ip - j;
          for (k=0; k<l1; k++) {
            for (i=2; i<ido; i+=2) {
              cc[i - 1 + (k + j*l1)*ido] = ch[i - 1 + (k + j*l1)*ido] + ch[i - 1 + (k + jc*l1)*ido];
              cc[i - 1 + (k + jc*l1)*ido] = ch[i + (k + j*l1)*ido] - ch[i + (k + jc*l1)*ido];
              cc[i + (k + j*l1)*ido] = ch[i + (k + j*l1)*ido] + ch[i + (k + jc*l1)*ido];
              cc[i + (k + jc*l1)*ido] = ch[i - 1 + (k + jc*l1)*ido] - ch[i - 1 + (k + j*l1)*ido];
            }
          }
        }
      } else {
        for (j=1; j<ipph; j++) {
          jc = ip - j;
          for (i=2; i<ido; i+=2) {
            for (k=0; k<l1; k++) {
              cc[i - 1 + (k + j*l1)*ido] =
                  ch[i - 1 + (k + j*l1)*ido] + ch[i - 1 + (k + jc*l1)*ido];
              cc[i - 1 + (k + jc*l1)*ido] = ch[i + (k + j*l1)*ido] - ch[i + (k + jc*l1)*ido];
              cc[i + (k + j*l1)*ido] = ch[i + (k + j*l1)*ido] + ch[i + (k + jc*l1)*ido];
              cc[i + (k + jc*l1)*ido] = ch[i - 1 + (k + jc*l1)*ido] - ch[i - 1 + (k + j*l1)*ido];
            }
          }
        }
      }
    } else {  /* now ido == 1 */
      for (ik=0; ik<idl1; ik++) cc[ik] = ch[ik];
    }
    for (j=1; j<ipph; j++) {
      jc = ip - j;
      for (k=0; k<l1; k++) {
        cc[(k + j*l1)*ido] = ch[(k + j*l1)*ido] + ch[(k + jc*l1)*ido];
        cc[(k + jc*l1)*ido] = ch[(k + jc*l1)*ido] - ch[(k + j*l1)*ido];
      }
    }

    ar1 = 1;
    ai1 = 0;
    for (l=1; l<ipph; l++) {
      lc = ip - l;
      ar1h = dcp*ar1 - dsp*ai1;
      ai1 = dcp*ai1 + dsp*ar1;
      ar1 = ar1h;
      for (ik=0; ik<idl1; ik++) {
        ch[ik + l*idl1] = cc[ik] + ar1*cc[ik + idl1];
        ch[ik + lc*idl1] = ai1*cc[ik + (ip-1)*idl1];
      }
      dc2 = ar1;
      ds2 = ai1;
      ar2 = ar1;
      ai2 = ai1;
      for (j=2; j<ipph; j++) {
        jc = ip - j;
        ar2h = dc2*ar2 - ds2*ai2;
        ai2 = dc2*ai2 + ds2*ar2;
        ar2 = ar2h;
        for (ik=0; ik<idl1; ik++) {
          ch[ik + l*idl1] += ar2*cc[ik + j*idl1];
          ch[ik + lc*idl1] += ai2*cc[ik + jc*idl1];
        }
      }
    }
    for (j=1; j<ipph; j++)
      for (ik=0; ik<idl1; ik++)
        ch[ik] += cc[ik + j*idl1];

    if (ido >= l1) {
      for (k=0; k<l1; k++) {
        for (i=0; i<ido; i++) {
          ref(cc,i + k*ip*ido) = ch[i + k*ido];
        }
      }
    } else {
      for (i=0; i<ido; i++) {
        for (k=0; k<l1; k++) {
          ref(cc,i + k*ip*ido) = ch[i + k*ido];
        }
      }
    }
    for (j=1; j<ipph; j++) {
      jc = ip - j;
      j2 = 2*j;
      for (k=0; k<l1; k++) {
        ref(cc,ido-1 + (j2 - 1 + k*ip)*ido) =
            ch[(k + j*l1)*ido];
        ref(cc,(j2 + k*ip)*ido) =
            ch[(k + jc*l1)*ido];
      }
    }
    if (ido == 1) return;
    if (nbd >= l1) {
      for (j=1; j<ipph; j++) {
        jc = ip - j;
        j2 = 2*j;
        for (k=0; k<l1; k++) {
          for (i=2; i<ido; i+=2) {
            ic = ido - i;
            ref(cc,i - 1 + (j2 + k*ip)*ido) = ch[i - 1 + (k + j*l1)*ido] + ch[i - 1 + (k + jc*l1)*ido];
            ref(cc,ic - 1 + (j2 - 1 + k*ip)*ido) = ch[i - 1 + (k + j*l1)*ido] - ch[i - 1 + (k + jc*l1)*ido];
            ref(cc,i + (j2 + k*ip)*ido) = ch[i + (k + j*l1)*ido] + ch[i + (k + jc*l1)*ido];
            ref(cc,ic + (j2 - 1 + k*ip)*ido) = ch[i + (k + jc*l1)*ido] - ch[i + (k + j*l1)*ido];
          }
        }
      }
    } else {
      for (j=1; j<ipph; j++) {
        jc = ip - j;
        j2 = 2*j;
        for (i=2; i<ido; i+=2) {
          ic = ido - i;
          for (k=0; k<l1; k++) {
            ref(cc,i - 1 + (j2 + k*ip)*ido) = ch[i - 1 + (k + j*l1)*ido] + ch[i - 1 + (k + jc*l1)*ido];
            ref(cc,ic - 1 + (j2 - 1 + k*ip)*ido) = ch[i - 1 + (k + j*l1)*ido] - ch[i - 1 + (k + jc*l1)*ido];
            ref(cc,i + (j2 + k*ip)*ido) = ch[i + (k + j*l1)*ido] + ch[i + (k + jc*l1)*ido];
            ref(cc,ic + (j2 - 1 + k*ip)*ido) = ch[i + (k + jc*l1)*ido] - ch[i + (k + j*l1)*ido];
          }
        }
      }
    }
  } /* radfg */


static void factorize(int n, int ifac[MAXFAC+2], const int ntryh[NSPECIAL])
  /* Factorize n in factors in ntryh and rest. On exit,
ifac[0] contains n and ifac[1] contains number of factors,
the factors start from ifac[2]. */
  {
    int ntry=3, i, j=0, ib, nf=0, nl=n, nq, nr;
startloop:
    if (j < NSPECIAL)
      ntry = ntryh[j];
    else
      ntry+= 2;
    j++;
    do {
      nq = nl / ntry;
      nr = nl - ntry*nq;
      if (nr != 0) goto startloop;
      nf++;
      ifac[nf + 1] = ntry;
      nl = nq;
      if (ntry == 2 && nf != 1) {
        for (i=2; i<=nf; i++) {
          ib = nf - i + 2;
          ifac[ib + 1] = ifac[ib];
        }
        ifac[2] = 2;
      }
    } while (nl != 1);
    ifac[0] = n;
    ifac[1] = nf;
  }




  /* ----------------------------------------------------------------------
rfftf1, rfftb1, rfftf, rfftb, rffti1, rffti. Treal FFTs.
---------------------------------------------------------------------- */

static void rfftf1(int n, Treal c[], Treal ch[], const Treal wa[], const int ifac[MAXFAC+2])
  {
    int i;
    int k1, l1, l2, na, kh, nf, ip, iw, ix2, ix3, ix4, ido, idl1;
    Treal *cinput, *coutput;
    nf = ifac[1];
    na = 1;
    l2 = n;
    iw = n-1;
    for (k1 = 1; k1 <= nf; ++k1) {
      kh = nf - k1;
      ip = ifac[kh + 2];
      l1 = l2 / ip;
      ido = n / l2;
      idl1 = ido*l1;
      iw -= (ip - 1)*ido;
      na = !na;
      if (na) {
        cinput = ch;
        coutput = c;
      } else {
        cinput = c;
        coutput = ch;
      }
      switch (ip) {
      case 4:
        ix2 = iw + ido;
        ix3 = ix2 + ido;
        radf4(ido, l1, cinput, coutput, &wa[iw], &wa[ix2], &wa[ix3]);
        break;
      case 2:
        radf2(ido, l1, cinput, coutput, &wa[iw]);
        break;
      case 3:
        ix2 = iw + ido;
        radf3(ido, l1, cinput, coutput, &wa[iw], &wa[ix2]);
        break;
      case 5:
        ix2 = iw + ido;
        ix3 = ix2 + ido;
        ix4 = ix3 + ido;
        radf5(ido, l1, cinput, coutput, &wa[iw], &wa[ix2], &wa[ix3], &wa[ix4]);
        break;
      default:
        if (ido == 1)
          na = !na;
        if (na == 0) {
          radfg(ido, ip, l1, idl1, c, ch, &wa[iw]);
          na = 1;
        } else {
          radfg(ido, ip, l1, idl1, ch, c, &wa[iw]);
          na = 0;
        }
      }
      l2 = l1;
    }
    if (na == 1) return;
    for (i = 0; i < n; i++) c[i] = ch[i];
  } /* rfftf1 */


static void rfftf(int n, Treal r[], Treal wsave[])
  {
    if (n == 1) return;
    rfftf1(n, r, wsave, wsave+n, (int*)(wsave+2*n));
  } /* rfftf */


static void rffti1(int n, Treal wa[], int ifac[MAXFAC+2])
  {
    static const Treal twopi = 6.28318530717959;
    Treal arg, argh, argld, fi;
    int i, j;
    int k1, l1, l2;
    int ld, ii, nf, ip, is;
    int ido, ipm, nfm1;
    static const int ntryh[NSPECIAL] = {
      4,2,3,5    }; /* Do not change the order of these. */
    factorize(n,ifac,ntryh);
    nf = ifac[1];
    argh = twopi / n;
    is = 0;
    nfm1 = nf - 1;
    l1 = 1;
    if (nfm1 == 0) return;
    for (k1 = 1; k1 <= nfm1; k1++) {
      ip = ifac[k1 + 1];
      ld = 0;
      l2 = l1*ip;
      ido = n / l2;
      ipm = ip - 1;
      for (j = 1; j <= ipm; ++j) {
        ld += l1;
        i = is;
        argld = (Treal) ld*argh;
        fi = 0;
        for (ii = 3; ii <= ido; ii += 2) {
          i += 2;
          fi += 1;
          arg = fi*argld;
          wa[i - 2] = cos(arg);
          wa[i - 1] = sin(arg);
        }
        is += ido;
      }
      l1 = l2;
    }
  } /* rffti1 */


static void rffti(int n, Treal wsave[])
  {
    if (n == 1) return;
    rffti1(n, wsave+n, (int*)(wsave+2*n));
  } /* rffti */

#ifdef __cplusplus
}
#endif
