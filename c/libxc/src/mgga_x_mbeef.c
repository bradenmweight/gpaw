#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "util.h"

#define XC_MGGA_X_MBEEF          207 /* mBEEF Exchange*/


/*changes static with const*/
const XC(func_info_type) XC(func_info_mgga_x_mbeef) = {
  XC_MGGA_X_MBEEF,
  XC_EXCHANGE,
  "mBEEF",
  XC_FAMILY_MGGA,
  "mBEEF", 
  XC_PROVIDES_EXC | XC_PROVIDES_VXC
};


void XC(mgga_x_mbeef_init)(XC(mgga_type) *p)
{
  p->info = &XC(func_info_mgga_x_mbeef);

  p->lda_aux = (XC(lda_type) *) malloc(sizeof(XC(lda_type)));
  XC(lda_x_init)(p->lda_aux, XC_UNPOLARIZED, 3, XC_NON_RELATIVISTIC);

}


void XC(mgga_x_mbeef_end)(XC(mgga_type) *p)
{
  free(p->lda_aux);
}


static void
mbeef_exchange(XC(mgga_type) *pt, double *rho, double sigma, double tau_,
	    double *energy, double *dedd, double *vsigma, double *dedtau)
{
  double gdms, s2, ds2dd, ds2dsigma;
  double k, xi, xj, tmp, tmp1, tmp2, tmp3, tmp4, dxids2, dxjdalpha;
  double alpha, tau, tauw, tau_lsda, aux;
  double dtau_lsdadd, dalphadsigma, dalphadtau, dalphadd;
  double exunif, vxunif;
  double Fx, dFdxi, dFdxj;

  /* HEG energy and potential */
  XC(lda_vxc)(pt->lda_aux, rho, &exunif, &vxunif);

  /* calculate |nabla rho|^2 */
  gdms = max(MIN_GRAD*MIN_GRAD, sigma);
  
  /* reduced density gradient in transformation t1(s) */
  s2 = gdms/(4.0*POW(3.0*M_PI*M_PI, 2.0/3.0)*POW(rho[0], 8.0/3.0));
  ds2dd = -(8.0/3.0)*s2/rho[0];
  ds2dsigma = 1.0/(4.0*POW(3.0*M_PI*M_PI, 2.0/3.0)*POW(rho[0], 8.0/3.0));
  k = 3.6;
  tmp = k + s2;
  xi = 2.0 * s2 / tmp - 1.0;
  dxids2 = 2.0 * k / POW(tmp, 2.0);

  /* kinetic energy densities */
  tauw = max(gdms/(8.0*rho[0]), 1.0e-12);
  tau = max(tau_, tauw);
  aux = (3./10.) * POW((3.0*M_PI*M_PI), 2.0/3.0);
  tau_lsda = aux * POW(rho[0], 5.0/3.0); 
  dtau_lsdadd = aux * 5.0/3.0 * POW(rho[0], 2.0/3.0);

  /* alpha in transformation t2(a) */
  alpha = (tau - tauw)/tau_lsda;
  assert(alpha >= 0.0);
  tmp1 = POW(1.0 - POW(alpha, 2.0), 3.0);
  tmp2 = 1.0 + POW(alpha, 3.0) + POW(alpha, 6.0);
  xj = -1.0 * tmp1 / tmp2;
  tmp3 = -6.0*alpha +12.0*POW(alpha, 3.0) -6.0*POW(alpha, 5.0);
  tmp4 = 3.0*POW(alpha, 2.0) +6.0*POW(alpha, 5.0);
  dxjdalpha = -1.0 * (tmp3*tmp2 - tmp1*tmp4) / POW(tmp2, 2.0);

  if(ABS(tau - tauw) < 1.0e-20)
    {
    dalphadsigma = 0.0;
    dalphadtau = 0.0;
    dalphadd = 0.0; 
    }
  else
    {
    dalphadtau = 1.0/tau_lsda;
    dalphadsigma = -1.0/(tau_lsda*8.0*rho[0]);
    dalphadd = (tauw/rho[0]*tau_lsda - (tau-tauw)*dtau_lsdadd) / POW(tau_lsda, 2.0); 
    }

  /* product exchange enhancement factor and derivatives */
  double coefs[25] = { 1.21069701e+00,   2.01715812e-01,  -7.82094390e-03,   2.02237667e-03,  -6.14012920e-04,
                      -3.30416811e-02,   3.19945649e-02,   6.20318509e-03,  -4.63540317e-03,   4.24450456e-04,
                       1.16718139e-02,  -1.09824006e-02,   7.23718253e-03,  -7.88144670e-04,  -6.56934508e-05,
                      -1.56077377e-03,   4.65937251e-03,  -9.37254135e-04,  -2.03034625e-05,   8.12822447e-05,
                      -5.74437939e-04,   2.98249118e-04,   4.56158143e-05,   6.87660911e-05,  -1.14058408e-05 };

  int order = 5;
  double Li[order];
  double dLi[order];
  double Lj[order];
  double dLj[order];

  /* initializing */
  Li[0] = 1.0;
  Li[1] = xi;
  dLi[0] = 0.0;
  dLi[1] = 1.0;
  Lj[0] = 1.0;
  Lj[1] = xj;
  dLj[0] = 0.0;
  dLj[1] = 1.0;
  Fx = 0.0;
  dFdxi = 0.0;
  dFdxj = 0.0;

  /* recursively building polynomia and their derivatives */
  for(int i = 2; i < order; i++)
    {
    Li[i] = 2.0 * xi * Li[i-1] - Li[i-2] - (xi * Li[i-1] - Li[i-2])/i;
    Lj[i] = 2.0 * xj * Lj[i-1] - Lj[i-2] - (xj * Lj[i-1] - Lj[i-2])/i;
    dLi[i] = i * Li[i-1] + xi * dLi[i-1];
    dLj[i] = i * Lj[i-1] + xj * dLj[i-1];
    }

  /* building enhancement factor and derivatives */
  int m = 0;
  for(int j = 0; j < order; j++)
    {
    for(int i = 0; i < order; i++)
      {
      Fx += coefs[m] * Li[i] * Lj[j];
      dFdxi += coefs[m] * dLi[i] * Lj[j];
      dFdxj += coefs[m] * dLj[j] * Li[i];
      m += 1;
      }
    }

  /* exchange energy */
  *energy = exunif * Fx * rho[0];

  /* exunif is energy per particle already
     so we multiply by n the terms with exunif*/

  *dedd = vxunif * Fx + exunif * rho[0] * (dFdxi * dxids2 * ds2dd + dFdxj * dxjdalpha * dalphadd);
  *vsigma = exunif * rho[0] * (dFdxi * dxids2 * ds2dsigma + dFdxj * dxjdalpha * dalphadsigma);
  *dedtau = exunif * rho[0] * dFdxj * dxjdalpha * dalphadtau;
}


void
XC(mgga_x_mbeef)(XC(mgga_type) *p, double *rho, double *sigma, double *tau,
	    double *e, double *dedd, double *vsigma, double *dedtau)
{
  if(p->nspin == XC_UNPOLARIZED)
    {
    double en;
    mbeef_exchange(p, rho, sigma[0], tau[0], &en, dedd, vsigma, dedtau);
    *e = en/(rho[0] + rho[1]);
    }
  else
    {
    /* The spin polarized version is handled using the exact spin scaling
          Ex[n1, n2] = (Ex[2*n1] + Ex[2*n2])/2
    */

    double e2na, e2nb, rhoa[2], rhob[2];
    double vsigmapart[3];

    rhoa[0] = 2.0 * rho[0];
    rhoa[1] = 0.0;
    rhob[0] = 2.0 * rho[1];
    rhob[1] = 0.0;

    mbeef_exchange(p, rhoa, 4.0*sigma[0], 2.0*tau[0], &e2na, &(dedd[0]), &(vsigmapart[0]), &(dedtau[0]));
    mbeef_exchange(p, rhob, 4.0*sigma[2], 2.0*tau[1], &e2nb, &(dedd[1]), &(vsigmapart[2]), &(dedtau[1]));

    *e = (e2na + e2nb) / (2.0*(rho[0] + rho[1]));
    vsigma[0] = 2.0 * vsigmapart[0];
    vsigma[2] = 2.0 * vsigmapart[2];
    }
}
