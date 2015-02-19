#include<stdio.h>
#include"charge_transfer.h"
#include"external-and-prototypes.h"

void mulliken(int nn, double *qmat, double *qmulli, double *qtot, int ndim,
                double dacc, double *occ, double **a, double **overl, int *ind,
                int lmax[DFTB_MAXTYPES], int *izp)
{
  int i, j, lj, m, n, jofn, mj;
  double q, qhelp, temp;

/*
  for (i=0; i<ndim; i++)
    if (occ[i] < dacc)
      break;

  adim = i;

  // matrix multiplication
  for (i=0; i<ndim; i++)
    for (j=0; j<ndim; j++) {
      temp = 0.0;
      for (k=0; k<ndim; k++)
        temp += overl[i][k] * a[k][j];
      afoo[i][j] = temp;
    }

  for (n=0; n<ndim; n++) {
    qmulli[n] = 0.0;
    for (i=0; i<adim; i++)
      qmulli[n] += occ[i] * afoo[n][i] * a[n][i];
    // printf("Qmulli(%d) = %f\n", n+1, qmulli[n]);
  }
*/

  for (m=0; m<ndim; m++) {
    qmulli[m] = 0.0;
    for (n=0; n<ndim; n++) {
      temp = 0.0;
      for (i=0; i<ndim; i++)
        temp += occ[i] * a[m][i] * a[n][i];
      qmulli[m] += overl[m][n] * temp;
    }
  }

  // INDEXING OK HERE!
  for (j=0; j<nn; j++) {
    q = 0.0;
    // printf("j = %d\n", j);
    for (lj=0; lj<lmax[izp[j]]; lj++) {
      jofn = ind[j] + lj*lj;
      // printf("lj = %d, jofn = %d\n", lj, jofn);
      qhelp = 0.0;
      for (mj=0; mj<=2*lj; mj++) {
        // printf("mj = %d\n", mj);
        qhelp += qmulli[jofn+mj];
      }
      q += qhelp;
    }
    qmat[j] = q;
  }

  *qtot = 0.0;
  for (j=0; j<nn; j++)
    *qtot += qmat[j];
  // printf("Mulliken - total charge = %f\n", *qtot);

  return;
}
