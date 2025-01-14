#include "../gpu.h"
#include "../gpu-complex.h"
#include "numpy/arrayobject.h"
#include "assert.h"

#define BETA   0.066725
#define GAMMA  0.031091
#define MU     0.2195164512208958 // PBE mod in libxc
//#define MU     0.2195149727645171 from libxc
#define C2     0.26053088059892404
#define C0I    0.238732414637843
#define C1    -0.45816529328314287
#define CC1    1.9236610509315362
#define CC2    2.5648814012420482
#define IF2    0.58482236226346462
#define C3     0.10231023756535741
#define C0     4.1887902047863905
#define THIRD  0.33333333333333333
#define NMIN   1.0E-10

__global__ void calculate_residual_kernel_complex(int nG, int nn,
						  gpuDoubleComplex* residual_nG,
						  double* eps_n,
						  gpuDoubleComplex* wf_nG)
{
    int n = threadIdx.x + blockIdx.x * blockDim.x;
    int g = threadIdx.y + blockIdx.y * blockDim.y;
    if ((g < nG) && (n < nn))
    {
	residual_nG[n*nG + g] = gpuCsub(residual_nG[n*nG + g],
					gpuCmulD(wf_nG[n*nG + g], eps_n[n]));
    }
}

__global__ void calculate_residual_kernel_real(int nG, int nn,
					       double* residual_nG,
					       double* eps_n,
					       double* wf_nG)
{
    int n = threadIdx.x + blockIdx.x * blockDim.x;
    int g = threadIdx.y + blockIdx.y * blockDim.y;
    if ((g < nG) && (n < nn))
    {
	residual_nG[n*nG + g] -= eps_n[n] * wf_nG[n*nG + g];
    }
}

extern "C"
void calculate_residual_launch_kernel(int nG,
				      int nn,
				      double* residual_nG,
				      double* eps_n,
				      double* wf_nG,
				      int is_complex)
{
    if (is_complex)
    {
	gpuLaunchKernel(calculate_residual_kernel_complex,
			dim3((nn+15)/16, (nG+15)/16),
			dim3(16, 16),
			0, 0,
			nG, nn,
			(gpuDoubleComplex*) residual_nG,
			eps_n,
			(gpuDoubleComplex*) wf_nG);
    }
    else
    {
	gpuLaunchKernel(calculate_residual_kernel_real,
			dim3((nn+15)/16, (nG+15)/16),
			dim3(16, 16),
			0, 0,
			nG, nn,
			residual_nG,
			eps_n,
			wf_nG);
    }
}


template <bool gga> __device__ double pbe_exchange(double n, double rs, double a2,
						   double* dedrs, double* deda2)
{
    double e = C1 / rs;
    *dedrs = -e / rs;
    if (gga)
    {
	double kappa = 0.804;
	double c = C2 * rs / n;
	c *= c;
	double s2 = a2 * c;
	double x = 1.0 + MU * s2 / kappa;
	double Fx = 1.0 + kappa - kappa / x;
	double dFxds2 = MU / (x * x);
	double ds2drs = 8.0 * c * a2 / rs;
	*dedrs = *dedrs * Fx + e * dFxds2 * ds2drs;
	*deda2 = e * dFxds2 * c;
	e *= Fx;
    }
    return e;
}


__device__ double G(double rtrs, double A, double alpha1,
		    double beta1, double beta2, double beta3, double beta4,
		    double* dGdrs)
{
  double Q0 = -2.0 * A * (1.0 + alpha1 * rtrs * rtrs);
  double Q1 = 2.0 * A * rtrs * (beta1 +
				rtrs * (beta2 +
					rtrs * (beta3 +
						rtrs * beta4)));
  double G1 = Q0 * log(1.0 + 1.0 / Q1);
  double dQ1drs = A * (beta1 / rtrs + 2.0 * beta2 +
		       rtrs * (3.0 * beta3 + 4.0 * beta4 * rtrs));
  *dGdrs = -2.0 * A * alpha1 * G1 / Q0 - Q0 * dQ1drs / (Q1 * (Q1 + 1.0));
  return G1;
}


template <bool gga, int nspin> __device__ double pbe_correlation(double n, double rs, double zeta, double a2,
						       double* dedrs, double* dedzeta, double* deda2)
{
  bool spinpol = nspin == 2;
  double rtrs = sqrt(rs);
  double de0drs;
  double e0 = G(rtrs, GAMMA, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294,
		&de0drs);
  double e;
  double xp = 117.0;
  double xm = 117.0;
  if (spinpol)
    {
      double de1drs;
      double e1 = G(rtrs, 0.015545, 0.20548, 14.1189, 6.1977, 3.3662,
		    0.62517, &de1drs);
      double dalphadrs;
      double alpha = -G(rtrs, 0.016887, 0.11125, 10.357, 3.6231, 0.88026,
			0.49671, &dalphadrs);
      dalphadrs = -dalphadrs;
      double zp = 1.0 + zeta;
      double zm = 1.0 - zeta;
      xp = pow(zp, THIRD);
      xm = pow(zm, THIRD);
      double f = CC1 * (zp * xp + zm * xm - 2.0);
      double f1 = CC2 * (xp - xm);
      double zeta2 = zeta * zeta;
      double zeta3 = zeta2 * zeta;
      double zeta4 = zeta2 * zeta2;
      double x = 1.0 - zeta4;
      *dedrs = (de0drs * (1.0 - f * zeta4) +
	       de1drs * f * zeta4 +
	       dalphadrs * f * x * IF2);
      *dedzeta = (4.0 * zeta3 * f * (e1 - e0 - alpha * IF2) +
		 f1 * (zeta4 * e1 - zeta4 * e0 + x * alpha * IF2));
      e = e0 + alpha * IF2 * f * x + (e1 - e0) * f * zeta4;
    }
  else
    {
      *dedrs = de0drs;
      e = e0;
    }
  if (gga)
    {
      double n2 = n * n;
      double t2;
      double y;
      double phi = 117.0;
      double phi2 = 117.0;
      double phi3 = 117.0;
      if (spinpol)
	{
	  phi = 0.5 * (xp * xp + xm * xm);
	  phi2 = phi * phi;
	  phi3 = phi * phi2;
	  t2 = C3 * a2 * rs / (n2 * phi2);
	  y = -e / (GAMMA * phi3);
	}
      else
	{
	  t2 = C3 * a2 * rs / n2;
	  y = -e / GAMMA;
	}
      double x = exp(y);
      double A;
      if (x != 1.0)
	A = BETA / (GAMMA * (x - 1.0));
      else
	A = BETA / (GAMMA * y);
      double At2 = A * t2;
      double nom = 1.0 + At2;
      double denom = nom + At2 * At2;
      double H = GAMMA * log( 1.0 + BETA * t2 * nom / (denom * GAMMA));
      double tmp = (GAMMA * BETA /
		    (denom * (BETA * t2 * nom + GAMMA * denom)));
      double tmp2 = A * A * x / BETA;
      double dAdrs = tmp2 * *dedrs;
      if (spinpol)
	{
	  H *= phi3;
	  tmp *= phi3;
	  dAdrs /= phi3;
	}
      double dHdt2 = (1.0 + 2.0 * At2) * tmp;
      double dHdA = -At2 * t2 * t2 * (2.0 + At2) * tmp;
      *dedrs += dHdt2 * 7 * t2 / rs + dHdA * dAdrs;
      *deda2 = dHdt2 * C3 * rs / n2;
      if (spinpol)
	{
	  double dphidzeta = (1.0 / xp - 1.0 / xm) / 3.0;
	  double dAdzeta = tmp2 * (*dedzeta -
				   3.0 * e * dphidzeta / phi) / phi3;
	  *dedzeta += ((3.0 * H / phi - dHdt2 * 2.0 * t2 / phi ) * dphidzeta +
		      dHdA * dAdzeta);
	  *deda2 /= phi2;
	}
      e += H;
    }
  return e;
}


template <int nspin, bool gga> __global__ void evaluate_ldaorgga_kernel(int ng,
									double* n_sg,
									double* v_sg,
									double* e_g,
									double* sigma_xg,
									double* dedsigma_xg)
{
    int g = threadIdx.x + blockIdx.x * blockDim.x;
    if (g >= ng)
    {
	return;
    }
    if (nspin == 1)
    {
	double n = n_sg[g];
	if (n < NMIN)
	  n = NMIN;
	double rs = pow(C0I / n, THIRD);
	double dexdrs;
	double dexda2;
	double ex;
	double decdrs;
	double decda2;
	double ec;
	if (gga)
	  {
	    double a2 = sigma_xg[g];
	    ex = pbe_exchange<gga>(n, rs, a2, &dexdrs, &dexda2);
	    ec = pbe_correlation<gga, nspin>(n, rs, 0.0, a2, &decdrs, 0, &decda2);
	    dedsigma_xg[g] = n * (dexda2 + decda2);
	  }
	else
	  {
	    ex = pbe_exchange<gga>(n, rs, 0.0, &dexdrs, 0);
	    ec = pbe_correlation<gga, nspin>(n, rs, 0.0, 0.0, &decdrs, 0, 0);
	  }
	e_g[g] = n * (ex + ec);
	v_sg[g] += ex + ec - rs * (dexdrs + decdrs) / 3.0;
    }
    else
    {
	const double* na_g = n_sg;
	double* va_g = v_sg;
	const double* nb_g = na_g + ng;
	double* vb_g = va_g + ng;

	const double* sigma0_g = 0;
	const double* sigma1_g = 0;
	const double* sigma2_g = 0;
	double* dedsigma0_g = 0;
	double* dedsigma1_g = 0;
	double* dedsigma2_g = 0;

	if (gga)
	{
	    sigma0_g = sigma_xg;
	    sigma1_g = sigma0_g + ng;
	    sigma2_g = sigma1_g + ng;
	    dedsigma0_g = dedsigma_xg;
	    dedsigma1_g = dedsigma0_g + ng;
	    dedsigma2_g = dedsigma1_g + ng;
	}

	double na = 2.0 * na_g[g];
	if (na < NMIN)
	  na = NMIN;
	double rsa = pow(C0I / na, THIRD);
	double nb = 2.0 * nb_g[g];
	if (nb < NMIN)
	  nb = NMIN;
	double rsb = pow(C0I / nb, THIRD);
	double n = 0.5 * (na + nb);
	double rs = pow(C0I / n, THIRD);
	double zeta = 0.5 * (na - nb) / n;
	double dexadrs;
	double dexada2;
	double exa;
	double dexbdrs;
	double dexbda2;
	double exb;
	double decdrs;
	double decdzeta;
	double decda2;
	double ec;
	if (gga)
	{
	    exa = pbe_exchange<gga>(na, rsa, 4.0 * sigma0_g[g],
			       &dexadrs, &dexada2);
	    exb = pbe_exchange<gga>(nb, rsb, 4.0 * sigma2_g[g],
				   &dexbdrs, &dexbda2);
	    double a2 = sigma0_g[g] + 2 * sigma1_g[g] + sigma2_g[g];
	    ec = pbe_correlation<gga, nspin>(n, rs, zeta, a2,
					     &decdrs, &decdzeta, &decda2);
	    dedsigma0_g[g] = 2 * na * dexada2 + n * decda2;
	    dedsigma1_g[g] = 2 * n * decda2;
	    dedsigma2_g[g] = 2 * nb * dexbda2 + n * decda2;
	}
	else
	{
	   exa = pbe_exchange<gga>(na, rsa, 0.0, &dexadrs, 0);
	   exb = pbe_exchange<gga>(nb, rsb, 0.0, &dexbdrs, 0);
	   ec = pbe_correlation<gga, nspin>(n, rs, zeta, 0.0,
					    &decdrs, &decdzeta, 0);
	}

	e_g[g] = 0.5 * (na * exa + nb * exb) + n * ec;
	va_g[g] += (exa + ec -
		    (rsa * dexadrs + rs * decdrs) / 3.0 -
		    (zeta - 1.0) * decdzeta);
	vb_g[g] += (exb + ec -
		    (rsb * dexbdrs + rs * decdrs) / 3.0 -
		    (zeta + 1.0) * decdzeta);
    }
}

// The define wrappers do not allow special characters for the first argument
// Hence, here defining an expression in such way, that the first argument can be
// a well defined identifier, and the preprocessor macro parses it correctly.
constexpr void(*LDA_SPINPAIRED)(int, double*, double*, double*, double*, double*) = &evaluate_ldaorgga_kernel<1, false>;
constexpr void(*LDA_SPINPOLARIZED)(int, double*, double*, double*, double*, double*) = &evaluate_ldaorgga_kernel<2, false>;
constexpr void(*PBE_SPINPAIRED)(int, double*, double*, double*, double*, double*) = &evaluate_ldaorgga_kernel<1, true>;
constexpr void(*PBE_SPINPOLARIZED)(int, double*, double*, double*, double*, double*) = &evaluate_ldaorgga_kernel<2, true>;

extern "C"
void evaluate_pbe_launch_kernel(int nspin, int ng,
				double* n,
				double* v,
				double* e,
				double* sigma,
				double* dedsigma)
{
    if (nspin == 1)
    {
	gpuLaunchKernel(PBE_SPINPAIRED,
			dim3((ng+255)/256),
			dim3(256),
			0, 0,
			ng,
			n, v, e, sigma, dedsigma);
    }
    else if (nspin == 2)
    {
	gpuLaunchKernel(PBE_SPINPOLARIZED,
			dim3((ng+255)/256),
			dim3(256),
			0, 0,
			ng,
			n, v, e, sigma, dedsigma);
    }
}

extern "C"
void evaluate_lda_launch_kernel(int nspin, int ng,
				double* n,
				double* v,
				double* e)
{
    if (nspin == 1)
    {
	gpuLaunchKernel(LDA_SPINPAIRED,
			dim3((ng+255)/256),
			dim3(256),
			0, 0,
			ng,
			n, v, e, NULL, NULL);
    }
    else if (nspin == 2)
    {
	gpuLaunchKernel(LDA_SPINPOLARIZED,
			dim3((ng+255)/256),
			dim3(256),
			0, 0,
			ng,
			n, v, e, NULL, NULL);
    }
}

__global__ void pw_insert_many_16(int nb,
				  int nG,
				  int nQ,
				  gpuDoubleComplex* c_nG,
				  npy_int32* Q_G,
				  double scale,
				  gpuDoubleComplex* tmp_nQ)
{
    int G = threadIdx.x + blockIdx.x * blockDim.x;
    int b = threadIdx.y + blockIdx.y * blockDim.y;
    __shared__ npy_int32 locQ_G[16];
    if (threadIdx.y == 0)
	locQ_G[threadIdx.x] = Q_G[G];
    __syncthreads();

    if ((G < nG) && (b < nb))
    {
	npy_int32 Q = locQ_G[threadIdx.x];
	tmp_nQ[Q + b * nQ] = gpuCmulD(c_nG[G + b * nG], scale);
    }
}

__global__ void add_to_density_8(int nb,
				  int nR,
				  double* f_n,
				  double* psit_nR,
				  double* rho_R)
{
    //int b = threadIdx.x + blockIdx.x * blockDim.x;
    int R = threadIdx.x + blockIdx.x * blockDim.x;
    if (R < nR)
    {
	double rho = 0.0;
	for (int b=0; b< nb; b++)
	{
	    int idx = b * nR + R;
	    rho += f_n[b] * (psit_nR[idx] * psit_nR[idx]);
	}
	rho_R[R] += rho;
    }
}


__global__ void add_to_density_16(int nb,
				  int nR,
				  double* f_n,
				  gpuDoubleComplex* psit_nR,
				  double* rho_R)
{
    //int b = threadIdx.x + blockIdx.x * blockDim.x;
    int R = threadIdx.x + blockIdx.x * blockDim.x;
    if (R < nR)
    {
	double rho = 0.0;
	for (int b=0; b< nb; b++)
	{
	    int idx = b * nR + R;
	    rho += f_n[b] * (psit_nR[idx].x * psit_nR[idx].x + psit_nR[idx].y * psit_nR[idx].y);
	}
	rho_R[R] += rho;
    }
}


__global__ void pw_insert_16(int nG,
			     int nQ,
			     gpuDoubleComplex* c_G,
			     npy_int32* Q_G,
			     double scale,
			     gpuDoubleComplex* tmp_Q)
{
    int G = threadIdx.x + blockIdx.x * blockDim.x;
    if (G < nG)
	tmp_Q[Q_G[G]] = gpuCmulD(c_G[G], scale);
}

extern "C" void gpawDeviceSynchronize()
{
    gpuDeviceSynchronize();
}

extern "C"
void add_to_density_gpu_launch_kernel(int nb,
				      int nR,
				      double* f_n,
				      gpuDoubleComplex* psit_nR,
				      double* rho_R,
				      int wfs_is_complex)
{
    if (wfs_is_complex)
    gpuLaunchKernel(add_to_density_16,
		    dim3((nR+255)/256),
		    dim3(256),
		    0, 0,
		    nb, nR,
		    f_n,
		    psit_nR,
		    rho_R);
    else
    gpuLaunchKernel(add_to_density_8,
		    dim3((nR+255)/256),
		    dim3(256),
		    0, 0,
		    nb, nR,
		    f_n,
		    (double*) psit_nR,
		    rho_R);
}

extern "C"
void pw_insert_gpu_launch_kernel(
			     int nb,
			     int nG,
			     int nQ,
			     double* c_nG,
			     npy_int32* Q_G,
			     double scale,
			     double* tmp_nQ)
{
    if (nb == 1)
    {
       gpuLaunchKernel(pw_insert_16,
		       dim3((nG+15)/16, (nb+15)/16),
		       dim3(16, 16),
		       0, 0,
		       nG, nQ,
		       (gpuDoubleComplex*) c_nG, Q_G,
		       scale,
		       (gpuDoubleComplex*) tmp_nQ);
    }
    else
    {
       gpuLaunchKernel(pw_insert_many_16,
		       dim3((nG+15)/16, (nb+15)/16),
		       dim3(16, 16),
		       0, 0,
		       nb, nG, nQ,
		       (gpuDoubleComplex*) c_nG,
		       Q_G,
		       scale,
		       (gpuDoubleComplex*) tmp_nQ);
    }
}


__global__ void pwlfc_expand_kernel_8(double* f_Gs,
				       gpuDoubleComplex *emiGR_Ga,
				       double *Y_GL,
				       int* l_s,
				       int* a_J,
				       int* s_J,
				       int* I_J,
				       double* f_GI,
				       int nG,
				       int nJ,
				       int nL,
				       int nI,
				       int natoms,
				       int nsplines,
				       bool cc)
{
    int G = threadIdx.x + blockIdx.x * blockDim.x;
    int J = threadIdx.y + blockIdx.y * blockDim.y;
    gpuDoubleComplex imag_powers[4] = {make_gpuDoubleComplex(1.0,0),
				       make_gpuDoubleComplex(0.0,-1.0),
				       make_gpuDoubleComplex(-1.0,0),
				       make_gpuDoubleComplex(0, 1.0)};
    if ((G < nG) && (J < nJ))
    {
	f_Gs += G*nsplines;
	emiGR_Ga += G*natoms;
	Y_GL += G*nL;
	f_GI += G*nI*2 + I_J[J];

	int s = s_J[J];
	int l = l_s[s];
	gpuDoubleComplex f1 = gpuCmulD(gpuCmul(emiGR_Ga[a_J[J]],
					       imag_powers[l % 4]),
				       f_Gs[s]);
	for (int m = 0; m < 2 * l + 1; m++) {
	    gpuDoubleComplex f = gpuCmulD(f1, Y_GL[l * l + m]);
	    f_GI[0] = f.x;
	    f_GI[nI] = cc ? -f.y : f.y;
	    f_GI++;
	}
    }
}

__global__ void pwlfc_expand_kernel_16(double* f_Gs,
				       gpuDoubleComplex *emiGR_Ga,
				       double *Y_GL,
				       int* l_s,
				       int* a_J,
				       int* s_J,
				       int* I_J,
				       double* f_GI,
				       int nG,
				       int nJ,
				       int nL,
				       int nI,
				       int natoms,
				       int nsplines,
				       bool cc)

{
    int G = threadIdx.x + blockIdx.x * blockDim.x;
    int J = threadIdx.y + blockIdx.y * blockDim.y;
    gpuDoubleComplex imag_powers[4] = {make_gpuDoubleComplex(1.0,0),
				       make_gpuDoubleComplex(0.0,-1.0),
				       make_gpuDoubleComplex(-1.0,0),
				       make_gpuDoubleComplex(0, 1.0)};
    if ((G < nG) && (J < nJ))
    {
	f_Gs += G*nsplines;
	emiGR_Ga += G*natoms;
	Y_GL += G*nL;
	f_GI += (G*nI + I_J[J])*2;
	int s = s_J[J];
	int l = l_s[s];
	gpuDoubleComplex f1 = gpuCmulD(gpuCmul(emiGR_Ga[a_J[J]],
					       imag_powers[l % 4]),
				       f_Gs[s]);
	for (int m = 0; m < 2 * l + 1; m++) {
	    gpuDoubleComplex f = gpuCmulD(f1, Y_GL[l * l + m]);
	    *f_GI++ = f.x;
	    *f_GI++ = cc ? -f.y : f.y;
	}
    }
}

// outP_ani[a] = \sum_A H_aii[a] P_ani[a]
__global__ void dH_aii_times_P_ani_16(int nA, int nn, int nI,
				      npy_int32* ni_a, double* dH_aii_dev,
				      gpuDoubleComplex* P_ani_dev,
				      gpuDoubleComplex* outP_ani_dev)
{
    int n1 = threadIdx.x + blockIdx.x * blockDim.x;
    if (n1 < nn) {
	double* dH_ii = dH_aii_dev;
	int I = 0;
	for (int a=0; a< nA; a++)
	{
	    int ni = ni_a[a];
	    int Istart = I;
	    for (int i=0; i< ni; i++)
	    {
		gpuDoubleComplex* outP_ni = outP_ani_dev + n1 * nI + I;
		gpuDoubleComplex result = make_gpuDoubleComplex(0.0, 0.0);
		gpuDoubleComplex* P_ni = P_ani_dev + n1 * nI + Istart;
		for (int i2=0; i2 < ni; i2++)
		{
		   gpuDoubleComplex item = gpuCmulD(*P_ni, dH_ii[i2 * ni + i]);
		   result.x += item.x;
		   result.y += item.y;
		   P_ni++;
		}
		outP_ni->x = result.x;
		outP_ni->y = result.y;
		I++;
	    }
	    dH_ii += ni * ni;
	}
    }
}

__global__ void dH_aii_times_P_ani_8(int nA, int nn, int nI,
				      npy_int32* ni_a, double* dH_aii_dev,
				      double* P_ani_dev,
				      double* outP_ani_dev)
{
    int n1 = threadIdx.x + blockIdx.x * blockDim.x;
    if (n1 < nn) {
	double* dH_ii = dH_aii_dev;
	int I = 0;
	for (int a=0; a< nA; a++)
	{
	    int ni = ni_a[a];
	    int Istart = I;
	    for (int i=0; i< ni; i++)
	    {
		double* outP_ni = outP_ani_dev + n1 * nI + I;
		double result = 0;
		double* P_ni = P_ani_dev + n1 * nI + Istart;
		for (int i2=0; i2 < ni; i2++)
		{
		   double item = *P_ni * dH_ii[i2 * ni + i];
		   result += item;
		   P_ni++;
		}
		*outP_ni = result;
		I++;
	    }
	    dH_ii += ni * ni;
	}
    }
}



extern "C"
void dH_aii_times_P_ani_launch_kernel(int nA, int nn,
				      int nI, npy_int32* ni_a,
				      double* dH_aii_dev,
				      gpuDoubleComplex* P_ani_dev,
				      gpuDoubleComplex* outP_ani_dev,
				      int is_complex)
{
    if (is_complex)
    gpuLaunchKernel(dH_aii_times_P_ani_16,
		    dim3((nn+255)/256),
		    dim3(256),
		    0, 0,
		    nA, nn, nI, ni_a, dH_aii_dev,
		    P_ani_dev, outP_ani_dev);
    else
    gpuLaunchKernel(dH_aii_times_P_ani_8,
		    dim3((nn+255)/256),
		    dim3(256),
		    0, 0,
		    nA, nn, nI, ni_a, dH_aii_dev,
		    (double*) P_ani_dev, (double*) outP_ani_dev);

}



extern "C"
void pwlfc_expand_gpu_launch_kernel(int itemsize,
				    double* f_Gs,
				    gpuDoubleComplex *emiGR_Ga,
				    double *Y_GL,
				    int* l_s,
				    int* a_J,
				    int* s_J,
				    double* f_GI,
				    int* I_J,
				    int nG,
				    int nJ,
				    int nL,
				    int nI,
				    int natoms,
				    int nsplines,
				    bool cc)
{
    if (itemsize == 16)
    {
	gpuLaunchKernel(pwlfc_expand_kernel_16,
			dim3((nG+15)/16, (nJ+15)/16),
			dim3(16, 16),
			0, 0,
			f_Gs,
			emiGR_Ga,
			Y_GL,
			l_s,
			a_J,
			s_J,
			I_J,
			f_GI,
			nG,
			nJ,
			nL,
			nI,
			natoms,
			nsplines,
			cc);
    }
    else
    {
	gpuLaunchKernel(pwlfc_expand_kernel_8,
			dim3((nG+15)/16, (nJ+15)/16),
			dim3(16, 16),
			0, 0,
			f_Gs,
			emiGR_Ga,
			Y_GL,
			l_s,
			a_J,
			s_J,
			I_J,
			f_GI,
			nG,
			nJ,
			nL,
			nI,
			natoms,
			nsplines,
			cc);
    }
    //gpuDeviceSynchronize();
}
