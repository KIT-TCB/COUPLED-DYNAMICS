#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"charge_transfer.h"
#include"external-and-prototypes.h"

//void slkmatrices(int i, int j, double (*xat)[3],
void slkmatrices(int i, int j, dvec *xat,
         double ham[LDIM][LDIM], double over[LDIM][LDIM],
         int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
         int *izp, tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3])
{
  double dif0[3];
  int k, l1, l2;

  for (k=0; k<3; k++)
    dif0[k] = xat[j][k] - xat[i][k];

  for (l2=0; l2<LDIM; l2++)
    for (l1=0; l1<LDIM; l1++) {
      ham[l1][l2] = 0.0;
      over[l1][l2] = 0.0;
    }

  slkode(dif0, izp[i], izp[j], ham, lmax, dim, dr, &skhpar, skstab, skhtab, skself);
  slkode(dif0, izp[i], izp[j], over, lmax, dim, dr, &skspar, skstab, skhtab, skself);
  return;
}

void slkode(double dum[3], int i, int j, double em[LDIM][LDIM], int lmax[DFTB_MAXTYPES], int dim[DFTB_MAXTYPES][DFTB_MAXTYPES], double dr[DFTB_MAXTYPES][DFTB_MAXTYPES],
               int (*iovpar)(int, int, double, double [13], int [DFTB_MAXTYPES],
                 tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *[DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][3],
                 int [DFTB_MAXTYPES][DFTB_MAXTYPES], double [DFTB_MAXTYPES][DFTB_MAXTYPES]),
               tendoubles *skstab[DFTB_MAXTYPES][DFTB_MAXTYPES], tendoubles *skhtab[DFTB_MAXTYPES][DFTB_MAXTYPES], double skself[DFTB_MAXTYPES][3])
{
  double x[6], x2[6], dummy[LDIM][LDIM], r2, r2i, ri;
  int k, l, minmax, maxmax;

  r2 = 0.0;
  for (l=0; l<3; l++) {
    x[l] = dum[l];
    x2[l] = x[l]*x[l];
    r2 += x2[l];
  }

  if (r2 > 1.0e-8) {

    r2i = 1.0/r2;
    ri = sqrt(r2i);
    for (l=0; l<3; l++) {
      x[l] *= ri;
      x[l+3] = x[l];
      x2[l] *= r2i;
      x2[l+3] = x2[l];
    }

    if (lmax[i] < lmax[j]) {
      maxmax = lmax[j];
      minmax = lmax[i];
    } else {
      maxmax = lmax[i];
      minmax = lmax[j];
    }

    skss(x, x2, i, j, r2, lmax, dim, dr, iovpar, em, skstab, skhtab, skself);

    if (maxmax == 1) return;

    if (minmax >= 2) {
      skpp(x, x2, i, j, r2, lmax, dim, dr, iovpar, em, skstab, skhtab, skself);
      sksp(x, x2, i, j, r2, lmax, dim, dr, iovpar, em, em, skstab, skhtab, skself);
      if (i != j)
        sksp(x, x2, j, i, r2, lmax, dim, dr, iovpar, dummy, em, skstab, skhtab, skself);
    } else {
      if (lmax[j] >= 2)
        sksp(x, x2, i, j, r2, lmax, dim, dr, iovpar, em, em, skstab, skhtab, skself);
      else
        sksp(x, x2, j, i, r2, lmax, dim, dr, iovpar, dummy, em, skstab, skhtab, skself);
    }

    if (maxmax == 2) return;

    if (minmax == 3) {
      skdd(x, x2, i, j, r2, lmax, dim, dr, iovpar, em, skstab, skhtab, skself);
      sksd(x, x2, i, j, r2, lmax, dim, dr, iovpar, em, em, skstab, skhtab, skself);
      skpd(x, x2, i, j, r2, lmax, dim, dr, iovpar, em, em, skstab, skhtab, skself);
      if (i != j) {
	sksd(x, x2, j, i, r2, lmax, dim, dr, iovpar, dummy, em, skstab, skhtab, skself);
	skpd(x, x2, j, i, r2, lmax, dim, dr, iovpar, dummy, em, skstab, skhtab, skself);
      }
    } else {
      if (lmax[i] == 1) {
        sksd(x, x2, i, j, r2, lmax, dim, dr, iovpar, em, em, skstab, skhtab, skself);
      }
      if (lmax[i] == 2) {
        sksd(x, x2, i, j, r2, lmax, dim, dr, iovpar, em, em, skstab, skhtab, skself);
        skpd(x, x2, i, j, r2, lmax, dim, dr, iovpar, em, em, skstab, skhtab, skself);
      }
      if (lmax[j] == 1) {
        sksd(x, x2, j, i, r2, lmax, dim, dr, iovpar, dummy, em, skstab, skhtab, skself);
      }
      if (lmax[j] == 2) {
        sksd(x, x2, j, i, r2, lmax, dim, dr, iovpar, dummy, em, skstab, skhtab, skself);
        skpd(x, x2, j, i, r2, lmax, dim, dr, iovpar, dummy, em, skstab, skhtab, skself);
      }
    }
  } else {
    //if (i != j) return;

    for (k=0; k<LDIM; k++)
      for (l=0; l<LDIM; l++)
        em[k][l] = 0.0;

    selfs(i, j, r2, lmax, dim, dr, iovpar, em, skstab, skhtab, skself);

    if (lmax[i] == 1) return;

    selfp(i, j, r2, lmax, dim, dr, iovpar, em, skstab, skhtab, skself);

    if (lmax[i] == 2) return;

    selfd(i, j, r2, lmax, dim, dr, iovpar, em, skstab, skhtab, skself);
  }

  return;
}

