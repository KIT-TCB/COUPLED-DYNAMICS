#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"qm_dftb.h"

/*
c===================================================================
c     calculate gamma or gamma^h-function and Gamma-function for 
c     one atom pairing (gamma^h_{ij} and Gamma_{ij}) 
c     details see Gaus JCTC 2011.
c     
c     INPUT:
c     real*8  r           distance between the two atoms i and j
c     real*8  ui,uj       Hubbard parameters for atom i and j
c     real*8  udi         Hubbard derivative dU/dq for atom i
c     logical xhgammahlp  .true. if h=exp(-((ui+uj)/2^zeta*r^2)) 
c     real*8  zeta        parameter for gamma^h
c    
c     OUTPUT:
c     real*8  gval        value for gamma_{ij} or gamma^h_{ij}
c     real*8  gder        value for Gamma_{ij} = dgamma_{ij}/dq_i
c===================================================================
*/

void gammaall(double r, double ui, double uj, double udi, int xhgammahlp, double zeta, double *gval, double *gder)
{
	const double zero=1.e-4;
	double a,b,a2,a3,a4,a6,b2,b3,b4,b6,
	       s,ar,g,fab,fba,dgda,dsdu,dfabda,dfbada,
	       h,dhdu;

	a = 3.2 * ui;
	b = 3.2 * uj;

	if (a+b < zero) {
		*gval = 0.;
		*gder = 0.;
	} else {
		if (r < zero) {
			*gval = 0.15625 * (a+b);
			*gder = 0.5;
		} else {
			if (fabs(a-b) < 1.0e-5) {
				ar   = (a+b)/2. * r;
				g    = (48. + 33.*ar + 9.*ar*ar + ar*ar*ar) / (r*48.);
				s    = exp(-ar) * g;
				dgda = (33. + 18.*ar + 3.*ar*ar) / 48.;
				dsdu = 3.2 * exp(-ar) * (dgda - r*g);
			} else {
				a2   = a*a;
				a3   = a2*a;
				a4   = a2*a2;
				a6   = a4*a2;
				b2   = b*b;
				b3   = b2*b;
				b4   = b2*b2;
				b6   = b4*b2;
				fab  = a*b4 / (2.*SQR(a2-b2)) - (b6-3.*a2*b4) / (CUB(a2-b2)*r);
				fba  = b*a4 / (2.*SQR(b2-a2)) - (a6-3.*b2*a4) / (CUB(b2-a2)*r);
				s    = exp(-a*r) * fab + exp(-b*r) * fba;
				dfabda = - (b6 + 3.*a2*b4) / (2. * CUB(a2-b2)) - 12.*a3*b4 / (r*QRT(a2-b2));
				dfbada = 2.*b3*a3 / (CUB(b2-a2)) + 12.*a3*b4 / (r*QRT(b2-a2));
				dsdu = 3.2 * (exp(-a*r) * (dfabda-r*fab) + exp(-b*r)*dfbada);
			}
			/* check for gamma^h */
			if (xhgammahlp) {
				h    = exp(-pow(((a+b)*0.15625), zeta) *r*r);
				dhdu = -h*zeta*r*r * pow(((a+b)*0.15625), zeta-1.) * 0.5;
				*gval= 1.0/r - s*h;
				*gder= -(dsdu*h + s*dhdu);
			} else {
				*gval= 1./r - s;
				*gder = -dsdu;
			}
		}
		*gder *= udi;
	}
	return;
}

/*
c===================================================================
c     calculate dgamma/dr or dgamma^h/dr and dGamma/dr for
c     one atom pairing (gamma^h_{ij} and Gamma_{ij})
c     details see Gaus JCTC 2011.
c     r=|R_j-R_i|
c
c     INPUT:
c     real*8  r           distance between the two atoms i and j
c     real*8  ui,uj       Hubbard parameters for atom i and j
c     real*8  udi         Hubbard derivative dU/dq for atom i
c     logical xhgammahlp  .true. if h=exp(-((ui+uj)/2^zeta*r^2))
c     real*8  zeta        parameter for gamma^h
c
c     OUTPUT:
c     real*8  dcdr        dgamma_{ij}/dr 
c     real*8  dcdr3       dGamma_{ij}/dr = d^2gamma_{ij}/drdq_i
c===================================================================
*/

void gammaall1(double r, double ui, double uj, double udi, int xhgammahlp, double zeta, double *dcdr, double *dcdr3)
{
	const double zero=1.e-4;
	double r2,a,b,a2,a3,a4,b2,b3,b4,z,z2,zr,g,dgdr,dsdr,fab,fba,dfabdr,dfbadr,
	       dcdudr,dsdudr,dgdadr,dgda,dfabdadr,dfbadadr,dfabda,dfbada,h,dhdu,dhdr,dhdudr,s,dsdu;

	r2 = r*r;
	a  = 3.2 * ui;
	b  = 3.2 * uj;

	if ((a+b < zero) || (r < zero)) {
	       *dcdr = 0.;
	       *dcdr3= 0.;
	       r    = 99999999.9; /* WHY ??? */
	} else { /* here: 1/r-s */
		if (fabs(a-b) < 1.e-5) {
			z    = 0.5 * (a+b);
			z2   = z*z;
			zr   = z*r;
			g    = (48. + 33.*zr + 9.*zr*zr + zr*zr*zr) / (48.*r);
			dgdr = -1./r2 + 3.*z2/16. + z2*zr/24.;
			dsdr = exp(-zr) * (dgdr - z*g);
			dgda   = (33. + 18.*zr + 3.0*zr*zr) / 48.;
			dgdadr = 0.375*z + 0.125*z2*r;
			dsdudr = 3.2 * exp(-zr) * (g*(zr-1.) - z*dgda + dgdadr - r*dgdr);
			if (xhgammahlp) {
				s    = exp(-zr)*g;
				dsdu = 3.2 * exp(-zr) * (dgda-r*g);
			}
		} else {
			a2  = a*a;
			a3  = a2*a;
			a4  = a2*a2;
			b2  = b*b;
			b3  = b2*b;
			b4  = b2*b2;
			fab = a*b4 / (2.*SQR(a2-b2)) - (b4*b2-3.*a2*b4) / (CUB(a2-b2)*r);
			fba = b*a4 / (2.*SQR(b2-a2)) - (a4*a2-3.*b2*a4) / (CUB(b2-a2)*r);
			dfabdr    = (b4*b2 - 3.*a2*b4) / (CUB(a2-b2)*r2);
			dfbadr    = (a4*a2 - 3.*b2*a4) / (CUB(b2-a2)*r2);
			dsdr      = exp(-a*r) * (dfabdr-a*fab) + exp(-b*r) * (dfbadr-b*fba);
			dfabda    = - (b2*b4 + 3.*a2*b4) / (2.*CUB(a2-b2)) -12.*a3*b4 / (r*QRT(a2-b2));
			dfbada    = 2.*b3*a3 / CUB(b2-a2) +12.*a3*b4 / (r*QRT(b2-a2));
			dfabdadr   = 12.*a3*b4 / (r2*QRT(a2-b2));
			dfbadadr  = -12.*a3*b4 / (r2*QRT(b2-a2));
			dsdudr    = 3.2 * (exp(-a*r) * (fab*(a*r-1.) - a*dfabda + dfabdadr - r*dfabdr) + exp(-b*r) * (dfbadadr-b*dfbada));
			if (xhgammahlp) {
				s   = exp(-a*r) * fab + exp(-b*r) * fba;
				dsdu= 3.2 * (exp(-a*r) * (dfabda-r*fab) + exp(-b*r) * dfbada);
			}
		}
		/* check for gamma^h */
		if (xhgammahlp) {
			h      = exp(- pow((a+b)*0.15625, zeta)*r*r);
			dhdu   = -h*zeta*r*r * pow((a+b)*0.15625, zeta-1.) * 0.5;
			dhdr   = -h*2.*r * pow((a+b)*0.15625, zeta);
			dhdudr = h*zeta*r * pow((a+b)*0.15625, zeta-1.) * (r*r * pow((a+b)*0.15625, zeta) - 1.);
			*dcdr   = - 1./r2 - (dsdr*h + s*dhdr);
			dcdudr = - (dsdudr*h + dsdu*dhdr + dsdr*dhdu + s*dhdudr);
			*dcdr3  = dcdudr * udi;
		} else {
			*dcdr   = -1./r2 - dsdr;
			dcdudr = -dsdudr;
			*dcdr3  = dcdudr * udi;
		}
	} /* end 1/r-s */
	return;
}

/*
c=======================================================================
c get the gamma and Gamma contribution to the gradient
c -F_{kx}= 0.5d0*\Delta q_k\sum_{a!=k}\Delta q_a(dgamma_{ak}/dR_{kx}+
c          dgamma_{ka}/dR_{kx})+1/3 \Delta q_k\sum_{a!=k}\Delta q_a (
c          \Delta q_a dGamma_{ak}/dR_{kx}+\Delta q_k dGamma_{ak}/dR_{kx}
c          )
c
c INPUT:
c integer  nn            number of atoms (in one cell)
c real*8   x(3,NNDIM)    coordinates
c real*8   izp(NNDIM)    map from atoms to atom types      
c real*8   uhubb(MAXTYP) Hubbard parameters
c real*8   uhder(MAXTYP) Hubbard derivatives
c logical  izpxh(MAXTYP) .true. for atom types which need extra term
c                         in gamma^h if switched on
c real*8   zeta          parameter for gamma^h (see Gaus JCTC 2011)
c real*8   qdiff(NNDIM)  atomic net charges (Mulliken)
c character*1 sccmode    last term of DFT taylor series which 
c                        is included (e.g. 2=2nd order, 3=3rdorder)
c
c OUTPUT:
c real*8   hgrad(3,NNDIM) gradient contribution 
c
c======================================================================
*/

void gammagrad(int nn, dvec *x, int *izp, double *uhubb, double *uhder,
		int *izpxh, double zeta, double *qdiff, int sccmode, double **hgrad)
{
	int ix,k,j;
	dvec tmp, tmp3, r;
	double bond,dgdrkj,dgdr3kj, /* dgdr[NNDIM][NNDIM],dgdr3[NNDIM][NNDIM], */
	       dgdrjk,dgdr3jk;

	for (k=0; k<nn; k++)
		clear_dvec(hgrad[k]);
	/*  get dgamma/dr and dGamma/dr   (r=|R_j-R_k|)
	 *  change with respect to original DFTB3:
	 *  get dr/dR_{kx} and build gammagradient contribution
	 *  right away in this cycle!
	 */
	for (k=0; k<nn; k++) {
		clear_dvec(tmp);
		clear_dvec(tmp3); 
		for (j=0; j<nn; j++) if (j != k) {
			dvec_sub(x[k], x[j], r);
			bond = dnorm(r);
			gammaall1(bond, uhubb[izp[k]], uhubb[izp[j]], uhder[izp[k]], izpxh[izp[k]] || izpxh[izp[j]], zeta, &dgdrkj, &dgdr3kj);
			gammaall1(bond, uhubb[izp[j]], uhubb[izp[k]], uhder[izp[j]], izpxh[izp[k]] || izpxh[izp[j]], zeta, &dgdrjk, &dgdr3jk);
			/*dgdr[k][j]  = dgdrkj;
			  dgdr3[k][j] = dgdr3kj; */
			for (ix=0; ix<3; ix++) {
				tmp[ix] += qdiff[j] * (dgdrjk + dgdrkj) / bond * r[ix];
				tmp3[ix]+= qdiff[j] * (qdiff[j]*dgdr3jk + qdiff[k]*dgdr3kj) / bond * r[ix];
			}
		}
		if (sccmode == 3)
			for (ix=0; ix<3; ix++)
				hgrad[k][ix] = qdiff[k] * (0.5*tmp[ix] + tmp3[ix]/3.);
		else
			for (ix=0; ix<3; ix++)
				hgrad[k][ix] = qdiff[k] * 0.5 * tmp[ix];
	}
	return;
}
/*
c======================================================================
c   Build symmetric matrix gammamat
c   Build non-symmetric matrix gammader
c
c   INPUT:
c   integer nn            number of atoms
c   real*8  x(3,NNDIM)    position of atoms
c   integer izp(NNDIM)    map from atoms to atom types
c   real*8  uhubb(MAXTYP) Hubbard parameters
c   real*8  uhder(MAXTYP) Hubbard derivatives dU/dq
c   real*8  zeta          parameter for gamma^h (see Gaus JCTC 2011)
c   logical izpxh(MAXTYP) .true. for atom types which need extra term
c                         in gamma^h if switched on
c
c   OUTPUT:
c   real*8 gammamat(*,*) matrix containing gamma/gamma^h for all atom-pairings
c   real*8 gammader(*,*) matrix containing Gamma=dgamma/dq for DFTB3 
c                        for all atom-pairings
c
c   Note that code is made efficient (but still easily readable) 
c   for DFTB3, but allows also running DFTB2, therefore gammader 
c   is calculated by default in this function but of course may 
c   be zeroed or controlled by a subroutine calling get_gammamat 
c   or using the OUTPUT.
c
c======================================================================
*/

void get_gammamat(int nn, dvec *x, int *izp, double *uhubb, double *uhder, double zeta,
		int *izpxh, double **gammamat, double **gammader)
{
	int i, j;
	dvec r;

	for (i=0; i<nn; i++)
		for (j=0; j<nn; j++) {
			dvec_sub(x[i], x[j], r);
			gammaall(dnorm(r), uhubb[izp[i]], uhubb[izp[j]], uhder[izp[i]], izpxh[izp[i]] || izpxh[izp[j]],
					zeta, &(gammamat[i][j]), &(gammader[i][j]));
		}
}

double gams(double r, double ua, double ub)
{
	double a,b,ar,a2,a4,a6,b2,b4,b6,fab,fba;

	a = 3.2 * ua;
	b = 3.2 * ub;
	if (fabs(a-b) < 1.e-5) {
		ar   = (a+b)/2. * r;
		return exp(-ar) * (48. + 33.*ar + 9.*ar*ar + ar*ar*ar) / (r*48.);
	} else {
		a2   = a*a;
		a4   = a2*a2;
		a6   = a4*a2;
		b2   = b*b;
		b4   = b2*b2;
		b6   = b4*b2;
		fab  = a*b4 / (2.*SQR(a2-b2)) - (b6 - 3.*a2*b4) / (CUB(a2-b2)*r);
		fba  = b*a4 / (2.*SQR(b2-a2)) - (a6 - 3.*b2*a4) / (CUB(b2-a2)*r);
		return exp(-a*r) * fab + exp(-b*r) * fba;
	}
}

double gams1(double r, double ua, double ub)
{
	double a,b,z,z2,zr,g,dgdr,a2,a4,b2,b4,fab,fba,dfabdr,dfbadr;

	a = 3.2 * ua;
	b = 3.2 * ub;
	if (fabs(a-b) < 1.e-5) {
		z     = (a+b)/2.;
		z2    = z*z;
		zr    = z*r;
		g     = (48. + 33.*zr + 9.*zr*zr + zr*zr*zr) / (48.*r);
		dgdr  = -1./(r*r) + 3.*z2/16. + z2*zr/24.;
		return exp(-zr) * (dgdr-z*g);
	} else {
		a2 = a*a;
		a4 = a2*a2;
		b2 = b*b;
		b4 = b2*b2;
		fab   = a*b4/(2.*SQR(a2-b2)) - (b4*b2 - 3.*a2*b4) / (CUB(a2-b2)*r);
		fba   = b*a4/(2.*SQR(b2-a2)) - (a4*a2 - 3.*b2*a4) / (CUB(b2-a2)*r);
		dfabdr= (b4*b2 - 3.*a2*b4) / (CUB(a2-b2)*r*r);
		dfbadr= (a4*a2 - 3.*b2*a4) / (CUB(b2-a2)*r*r);
		return exp(-a*r) * (dfabdr-a*fab) + exp(-b*r) * (dfbadr-b*fba);
	}
}
