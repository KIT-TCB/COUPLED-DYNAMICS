#include<stdio.h>
#include<stdlib.h>
#define Sqrt3 (1.732050808)
#include"charge_transfer.h"
#include"external-and-prototypes.h"

void skss(double x[6], double x2[6], int i, int j, double r2, int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
       int (*iovpar)(int, int, double, double [13], int [DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][3],
         int [DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][DFTB_MAXTYPES]),
       double em[LDIM][LDIM], tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3])
{
  int id;
  double parm[13];

  id = iovpar(i, j, r2, parm, lmax, skstab, skhtab, skself, dim, dr);
  em[0][0] = parm[9];
  return;
}

void sksp(double x[6], double x2[6], int i, int j, double r2, int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
       int (*iovpar)(int, int, double, double [13], int [DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][3],
         int [DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][DFTB_MAXTYPES]),
       double em[LDIM][LDIM], double emt[LDIM][LDIM], tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3])
{
  int l, id;
  double parm[13];

  id = iovpar(i, j, r2, parm, lmax, skstab, skhtab, skself, dim, dr);
  for (l=0; l<3; l++) {
    em[0][1+l] = x[l] * parm[8];
    emt[1+l][0] = -em[0][1+l];
    // em[1+l][0] = x[l] * parm[8];
    // emt[0][1+l] = -em[1+l][0];
  }
  return;
}

void sksd(double x[6], double x2[6], int i, int j, double r2, int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
       int (*iovpar)(int, int, double, double [13], int [DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][3],
         int [DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][DFTB_MAXTYPES]),
       double em[LDIM][LDIM], double emt[LDIM][LDIM], tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3])
{
  int l, id;
  double parm[13], es[5], d4, d5;

  d4 = x2[ZZ] - 0.5 * (x2[XX] + x2[YY]);
  d5 = x2[XX] - x2[YY];
  id = iovpar(i, j, r2, parm, lmax, skstab, skhtab, skself, dim, dr);
  for (l=0; l<3; l++)
    es[l] = Sqrt3 * x[l] * x[l+1];
  es[3] = 0.5 * Sqrt3 * d5;
  es[4] = d4;
  for (l=0; l<5; l++)
    emt[4+l][0] =
      em[0][4+l] = es[l] * parm[7];
  return;
}

void skpp(double x[6], double x2[6], int i, int j, double r2, int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
       int (*iovpar)(int, int, double, double [13], int [DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][3], 
         int [DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][DFTB_MAXTYPES]),
       double em[LDIM][LDIM], tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3])
{
  int k, l, ii, ir, is, id;
  double parm[13], epp[6], hp, dm[6];

  id = iovpar(i, j, r2, parm, lmax, skstab, skhtab, skself, dim, dr);
  for (l=0; l<3; l++) {
    epp[l] = x2[l];
    epp[l+3] = x[l] * x[l+1];
  }
  for (l=0; l<3; l++) {
    hp = epp[l];
    dm[l] = hp * parm[5] + (1.0-hp) * parm[6];
  }
  for (l=3; l<6; l++) {
    dm[l] = epp[l] * (parm[5] - parm[6]);
  }
  for (ir=0; ir<3; ir++)
    for (is=0; is<=ir; is++) {
      ii = ir - is;
      k = 3*ii - (ii*(ii-1))/2 + is;
      em[1+is][1+ir] = em[1+ir][1+is] = dm[k];
    }
  return;
}

void skpd(double x[6], double x2[6], int i, int j, double r2, int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
       int (*iovpar)(int, int, double, double [13], int [DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][3],
         int [DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][DFTB_MAXTYPES]),
       double em[LDIM][LDIM], double emt[LDIM][LDIM], tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3])
{
  int k, l, m, id, ir, is;
  double parm[13], epd[13][2], dm[15], d3, d4, d5, d6;

  d3 = x2[XX] + x2[YY];
  d4 = x2[ZZ] - 0.5 * d3;
  d5 = x2[XX] - x2[YY];
  d6 = x[XX] * x[YY] * x[ZZ];
  id = iovpar(i, j, r2, parm, lmax, skstab, skhtab, skself, dim, dr);

  for (l=0; l<3; l++) {
	  epd[l][0] = Sqrt3 * x2[l] * x[l+1];
	  epd[l][1] = x[l+1] * (1. - 2. * x2[l]);
	  epd[l+4][0] = Sqrt3 * x2[l] * x[l+2];
	  epd[l+4][1] = x[l+2] * (1. - 2. * x2[l]);
	  epd[l+7][0] = 0.5 * Sqrt3 * x[l] * d5;
	  epd[l+10][0] = x[l] * d4;
  }
  epd[3][0] = Sqrt3 * d6;
  epd[3][1] = -2. * d6;
  epd[7][1] = x[XX] * (1. - d5);
  epd[8][1] = -x[YY] * (1. + d5);
  epd[9][1] = -x[ZZ] * d5;
  epd[10][1] = - Sqrt3 * x[XX] * x2[ZZ];
  epd[11][1] = - Sqrt3 * x[YY] * x2[ZZ];
  epd[12][1] = Sqrt3 * x[ZZ] * d3;
  for (l=0; l<15; l++)
	  dm[l] = 0.;
  for (m=0; m<2; m++) {
	  dm[0] += epd[0][m] * parm[m+3];
	  dm[1] += epd[5][m] * parm[m+3];
	  dm[2] += epd[3][m] * parm[m+3];
	  dm[4] += epd[1][m] * parm[m+3];
	  dm[5] += epd[6][m] * parm[m+3];
	  dm[6] += epd[4][m] * parm[m+3];
	  dm[8] += epd[2][m] * parm[m+3];
	  for (l=7; l<13; l++)
		  dm[l+2] += epd[l][m] * parm[m+3];
  }
  dm[3] = dm[2];
  dm[7] = dm[2];
  for (ir=0; ir<5; ir++)
	  for (is=0; is<3; is++) {
		  k = 3 * ir + is;
		  emt[4+ir][1+is] = -dm[k];
		  em[1+is][4+ir] = dm[k];
	  }
  return;
}

void skdd(double x[6], double x2[6], int i, int j, double r2, int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
       int (*iovpar)(int, int, double, double [13], int [DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][3],
         int [DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][DFTB_MAXTYPES]),
       double em[LDIM][LDIM], tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3])
{
  int k, l, m, id, ii, ir, is;
  double parm[13], e[15][3], dm[15], dd[3], d3, d4, d5, d6;

  d3 = x2[XX] + x2[YY];
  d4 = x2[ZZ] - 0.5 * d3;
  d5 = x2[XX] - x2[YY];
  id = iovpar(i, j, r2, parm, lmax, skstab, skhtab, skself, dim, dr);

  for (l=0; l<3; l++) {
	  e[l][0] = x2[l] * x2[l+1];
	  e[l][1] = x2[l] + x2[l+1] - 4. * e[l][0];
	  e[l][2] = x2[l+2] + e[l][0];
	  e[l][0] *= 3.;
  }
  e[3][0] = SQR(d5);
  e[3][1] = d3 - e[3][0];
  e[3][2] = x2[ZZ] + 0.25 * e[3][0];
  e[3][0] *= 0.75;
  e[4][0] = SQR(d4);
  e[4][1] = 3. * x2[ZZ] * d3;
  e[4][2] = 0.75 * SQR(d3);
  dd[0] = x[XX] * x[ZZ];
  dd[1] = x[YY] * x[XX];
  dd[2] = x[ZZ] * x[YY];
  for (l=0; l<2; l++) {
	  e[l+5][0] = 3. * x2[l+1] * dd[l];
	  e[l+5][1] = dd[l] * (1. - 4. * x2[l+1]);
	  e[l+5][2] = dd[l] * (x2[l+1] - 1.);
  }
  e[7][0] = dd[0] * d5 * 1.5;
  e[7][1] = dd[0] * (1. - 2. * d5);
  e[7][2] = dd[0] * (0.5 * d5 - 1.);
  e[8][0] = 0.5 * Sqrt3 * d5 * d4;
  e[8][1] = -Sqrt3 * d5 * x2[ZZ];
  e[8][2] = 0.25 * Sqrt3 * d5 * (1. + x2[ZZ]);
  e[9][0] = 3. * x2[XX] * dd[2];
  e[9][1] = 4. * (0.25 - x2[XX]) * dd[2];
  e[9][2] = (x2[XX] - 1.) * dd[2];
  e[10][0] = 1.5 * dd[2] * d5;
  e[10][1] = -dd[2] * (1. + 2. * d5);
  e[10][2] = dd[2] * (1. + 0.5 * d5);
  e[12][2] = 0.5 * d5 * dd[1];
  e[12][1] = -2. * dd[1] * d5;
  e[12][0] = 3. * e[12][2];
  e[11][0] = Sqrt3 * d4 * dd[0];
  e[13][0] = Sqrt3 * d4 * dd[2];
  e[14][0] = Sqrt3 * d4 * dd[1];
  e[14][1] = -2. * Sqrt3 * dd[1] * x2[ZZ];
  e[14][2] = 0.5 * Sqrt3 * (1. + x2[ZZ]) * dd[1];
  e[13][1] = Sqrt3 * dd[2] * (d3 - x2[ZZ]);
  e[13][2] = -0.5 * Sqrt3 * dd[2] * d3;
  e[11][1] = Sqrt3 * dd[0] * (d3 - x2[ZZ]);
  e[11][2] = -0.5 * Sqrt3 * dd[0] * d3;
  for (l=0; l<15; l++) {
	  dm[l] = 0.;
	  for (m=0; m<3; m++)
		  dm[l] += e[l][m] * parm[m];
  }
  for (ir=0; ir<5; ir++)
	  for (is=0; is<=ir; is++) {
		  ii = ir-is;
		  k = 5*ii - ii*(ii-1)/2 + is;
		  em[4+ir][4+is] = em[4+is][4+ir] = dm[k];
	  }
  return;
}

void selfs(int i, int j, double r2, int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
       int (*iovpar)(int, int, double, double [13], int [DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][3],
         int [DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][DFTB_MAXTYPES]),
       double em[LDIM][LDIM], tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3])
{
  int id;
  double parm[13];

  id = iovpar(i, j, r2, parm, lmax, skstab, skhtab, skself, dim, dr);
  em[0][0] = parm[12];
  return;
}

void selfp(int i, int j, double r2, int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
       int (*iovpar)(int, int, double, double [13], int [DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][3],
         int [DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][DFTB_MAXTYPES]),
       double em[LDIM][LDIM], tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3])
{
  int l, m, id;
  double parm[13];

  id = iovpar(i, j, r2, parm, lmax, skstab, skhtab, skself, dim, dr);
  for (l=0; l<3; l++) {
    for (m=0; m<3; m++)
      em[1+m][1+l] = 0.0;
    em[1+l][1+l] = parm[11];
  }
  return;
}

void selfd(int i, int j, double r2, int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
       int (*iovpar)(int, int, double, double [13], int [DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][3],
         int [DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][DFTB_MAXTYPES]),
       double em[LDIM][LDIM], tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3])
{
  int l, m, id;
  double parm[13];

  id = iovpar(i, j, r2, parm, lmax, skstab, skhtab, skself, dim, dr);
  for (l=0; l<5; l++) {
    for (m=0; m<5; m++)
      em[4+m][4+l] = 0.0;
    em[4+l][4+l] = parm[10];
  }
  return;
}

