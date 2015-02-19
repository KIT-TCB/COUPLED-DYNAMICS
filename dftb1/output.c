#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"charge_transfer.h"
#include"external-and-prototypes.h"

void outeigenvectors(double **a, double *ev, int *ind, int nn, dftb_phase1_t dftb1)
{
  FILE *f1, *f2;
  int i, j, k, l, neig;
  char orbital[4][4];

  strcpy(orbital[0], "S  ");
  strcpy(orbital[1], "Px ");
  strcpy(orbital[2], "Py ");
  strcpy(orbital[3], "Pz ");
  neig = ind[nn];

  f2 = fopen("ao2fo", "w");
  for (i=0; i<neig; i++)
    for (j=0; j<neig; j++)
      fprintf(f2, "%16.12f\n", a[j][i]);
  fclose(f2);

  f1 = fopen("EVC.DAT", "w");
  fprintf(f1, "THE ONE-ELECTRON EIGENVALUES AND EIGENVECTORS\n\n");
  for (i=0; i<neig; i++) {
    fprintf(f2, "%3dth eigenvalue = %10.6f H = %10.6f eV\n", i+1, ev[i], ev[i]*27.2116);
    fprintf(f2, "Atom No.   Atom Type\n");
    for (j=0; j<nn; j++) {
      fprintf(f2, "%5d  %5d\n", j+1, dftb1.izp[j]);
      for (l=ind[j], k=0; l<ind[j+1]; l++, k++)
        fprintf(f2, "%s%8.5f\n", orbital[k], a[l][i]);
    }
    fprintf(f2, "\n");
  }
  fclose(f1);
  
  return;
}

void outspec(int nn, int ndim, int *ind, double *ev, double *occ,
               double efermi, double *qmat, double *qmulli, dftb_t *dftb, dftb_phase1_t dftb1)
{
  FILE *f;
  int i, ind1, ind2, j;

  f = fopen("SPE.DAT", "w");
  fprintf(f, "%20.12f%20.12f\n", efermi, efermi*27.2114);
  for (i=0; i<ndim; i++)
    fprintf(f, "%20.12f%20.12f%20.12f\n", ev[i], ev[i]*27.2114, occ[i]);
  fclose(f);

  f = fopen("CHR.DAT", "w");
  for (i=0; i<nn; i++) {
    ind1 = ind[i];
    ind2 = ind1 + SQR(dftb->lmax[dftb1.izp[i]]);
    fprintf(f, "%4d%12.6f", i+1, qmat[i]);
    for (j=ind1; j<ind2; j++)
      fprintf(f, "%12.6f", qmulli[j]);
    fprintf(f, "\n");
  }
  fclose(f);

  return;
}
