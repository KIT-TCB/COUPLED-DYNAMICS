#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"charge_transfer.h"
#include"external-and-prototypes.h"

int skspar(int i, int j, double r2, double dd[13],
           int lmax[DFTB_MAXTYPES], tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3],
           int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES])
{
  int maxmax, minmax, in, ind, inu, mxind;
  double r, grdr, x0, x1, x2, f0, f1, f2, xh, hl;

  if (lmax[i] < lmax[j]) {
    maxmax = lmax[j];
    minmax = lmax[i];
  } else {
    maxmax = lmax[i];
    minmax = lmax[j];
  }

  switch (maxmax) {
    case 1:
      inu = 9;
      break;
    case 2:
      if (minmax == 1)
        inu = 8;
      else
        inu = 5;
      break;
    case 3:
      switch (minmax) {
        case 1:
          inu = 7;
          break;
        case 2:
          inu = 3;
          break;
        case 3:
          inu = 0;
      }
  }

  // mxind = Maximaler Index bis zu dem der Spline weitergefuehrt wird
  mxind = dim[i][j] + 0.3/dr[i][j] - 2;
  //mxind = dim[i][j] + (0.3/dr[i][j]);
  r = sqrt(r2);
  //ind = r/dr[i][j] + 1; // REALLY PLUS ONE???
  ind = r/dr[i][j];

  if (r2 < 1e-8) {
    for (in=0; in<3; in++)
      dd[in+10] = 1.0;
  } else {
    if (ind > dim[i][j]-3) {
      // FREE CUBIC SPLINE
      if (ind == dim[i][j]-2) {
        x0 = (dim[i][j] - 3) * dr[i][j];
        x1 = x0 + dr[i][j];
        x2 = x1 + dr[i][j];
        xh = r - x1;
        hl = x2 - x1;
        for (in=inu; in<10; in++) {
          f0 = skstab[i][j][dim[i][j]-3][in];
          f1 = skstab[i][j][dim[i][j]-2][in];
          f2 = skstab[i][j][dim[i][j]-1][in];
          dd[in] = cubicspline(f0, f1, f2, x0, x1, xh, hl, dr[i][j]);
        }
      } else {
        if (ind < mxind-1) {
	// 5TH DEGREE SPLINE
          x0 = (dim[i][j] - 3) * dr[i][j];
          x1 = x0 + dr[i][j];
          x2 = x1 + dr[i][j];
          //xh = r - (mxind-1) * dr[i][j];
          xh = r - mxind * dr[i][j];
	  for (in=inu; in<10; in++) {
            f0 = skstab[i][j][dim[i][j]-3][in];
            f1 = skstab[i][j][dim[i][j]-2][in];
            f2 = skstab[i][j][dim[i][j]-1][in];
            dd[in] = spline5th(f0, f1, f2, x0, x1, x2, xh, dr[i][j], mxind);
	  }
	} else {
        // ZERO
	  for (in=inu; in<10; in++)
	    dd[in] = 0.0;
        }
      }
    } else {
      grdr = (r - ind*dr[i][j]) / dr[i][j];
      for (in=inu; in<10; in++) {
        f0 = skstab[i][j][ind][in];
        f1 = skstab[i][j][ind+1][in];
        f2 = skstab[i][j][ind+2][in];
        dd[in] = f0 + (f1-f0)*grdr + (f2+f0-2.0*f1)*grdr*(grdr-1.0) / 2.0;
      }
    }
  }
  return 0;
}

int skhpar(int i, int j, double r2, double dd[13],
           int lmax[DFTB_MAXTYPES], tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3],
           int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES])
{
  int maxmax, minmax, in, ind, inu, mxind;
  double r, grdr, x0, x1, x2, f0, f1, f2, xh, hl;

  if (lmax[i] < lmax[j]) {
    maxmax = lmax[j];
    minmax = lmax[i];
  } else {
    maxmax = lmax[i];
    minmax = lmax[j];
  }

  switch (maxmax) {
    case 1:
      inu = 9;
      break;
    case 2:
      if (minmax == 1)
        inu = 8;
      else
        inu = 5;
      break;
    case 3:
      switch (minmax) {
        case 1:
          inu = 7;
          break;
        case 2:
          inu = 3;
          break;
        case 3:
          inu = 0;
      }
  }

  // mxind = Maximaler Index bis zu dem der Spline weitergefuehrt wird
  mxind = dim[i][j] + 0.3/dr[i][j] - 2;
  //mxind = dim[i][j] + (0.3/dr[i][j]);
  r = sqrt(r2);
  //ind = r/dr[i][j] + 1; // REALLY PLUS ONE???
  ind = r/dr[i][j];

  if (r2 < 1e-8) {
    for (in=0; in<3; in++)
      dd[in+10] = skself[i][in];
  } else {
    if (ind > dim[i][j]-3) {
      // FREE CUBIC SPLINE
      if (ind == dim[i][j]-2) {
        x0 = (dim[i][j] - 3) * dr[i][j];
        x1 = x0 + dr[i][j];
        x2 = x1 + dr[i][j];
        xh = r - x1;
        hl = x2 - x1;
        for (in=inu; in<10; in++) {
          f0 = skhtab[i][j][dim[i][j]-3][in];
          f1 = skhtab[i][j][dim[i][j]-2][in];
          f2 = skhtab[i][j][dim[i][j]-1][in];
          dd[in] = cubicspline(f0, f1, f2, x0, x1, xh, hl, dr[i][j]);
        }
      } else {
        if (ind < mxind-1) {
	// 5TH DEGREE SPLINE
          x0 = (dim[i][j] - 3) * dr[i][j];
          x1 = x0 + dr[i][j];
          x2 = x1 + dr[i][j];
          //xh = r - (mxind-1) * dr[i][j];
          xh = r - mxind * dr[i][j];
	  for (in=inu; in<10; in++) {
            f0 = skhtab[i][j][dim[i][j]-3][in];
            f1 = skhtab[i][j][dim[i][j]-2][in];
            f2 = skhtab[i][j][dim[i][j]-1][in];
            dd[in] = spline5th(f0, f1, f2, x0, x1, x2, xh, dr[i][j], mxind);
	  }
	} else {
        // ZERO
	  for (in=inu; in<10; in++)
	    dd[in] = 0.0;
        }
      }
    } else {
      grdr = (r - ind*dr[i][j]) / dr[i][j];
      for (in=inu; in<10; in++) {
        f0 = skhtab[i][j][ind][in];
        f1 = skhtab[i][j][ind+1][in];
        f2 = skhtab[i][j][ind+2][in];
        dd[in] = f0 + (f1-f0)*grdr + (f2+f0-2.0*f1)*grdr*(grdr-1.0) / 2.0;
      }
    }
  }
  return 0;
}

double cubicspline(double f0, double f1, double f2, double x0, double x1,
         double xh, double hl, double dr)
{
  double f1abl, f2abl, a, b, c, d;

  f2abl= (f2 + f0 - 2.0*f1) / dr*dr;
  f1abl= (f1 - f0)/dr + 0.5*f2abl*(x1-x0);
  a = f1;
  b = f1abl;
  c = f2abl/2.0;
  d = (f2-a)/(hl*hl*hl) - b/(hl*hl) - c/hl;

  return a + b*xh + c*xh*xh + d*xh*xh*xh;
}

double spline5th(double f0, double f1, double f2, double x0, double x1, double x2,
         double xh, double dr, int mxind)
{
  double hl, f1abl, f2abl, a, b, c, d, hsp, isp, jsp;

  f2abl = (f2+f0-2.0*f1) / dr*dr;
  f1abl = (f1-f0)/dr + 0.5*f2abl*(x1-x0);
  a = f1;
  b = f1abl;
  c = f2abl/2.0;
  hl = x2-x1;
  d = (f2-a)/(hl*hl*hl) - b/(hl*hl) - c/hl;

  f1abl = b + 2.0*c*hl + 3.0*d*hl*hl;
  f2abl = 2.0*c + 6.0*d*hl;

  hl = x2 - mxind*dr;
  hsp = 10.0*f2/(hl*hl*hl) - 4.0*f1abl/(hl*hl) + f2abl/(2.0*hl);
  isp = -15.0*f2/(hl*hl*hl*hl) + 7.0*f1abl/(hl*hl*hl) - f2abl/(hl*hl);
  jsp = 6.0*f2/(hl*hl*hl*hl*hl) - 3.0*f1abl/(hl*hl*hl*hl) + f2abl/(2.0*hl*hl*hl);

  hl=xh*xh*xh;
  return (hsp + isp*xh + jsp*xh*xh) * hl;
}

