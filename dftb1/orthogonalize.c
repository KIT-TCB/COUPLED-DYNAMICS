#include "charge_transfer.h"

long orthogonalize(double **THamil, double **OverlF, double **THamilOr, long n, dftb_orthogo_t arrays){
/* This function performs a standard DFTB calculation of a single fragment */

// PARAMETERS:
// THamil   = (in) non-othogonal FO Hamiltonian 
// OverlF   = (in) corresponding overlap matrix in FO basis
// THamilOr = (out) orthogonalized FO Hamiltonian
// n        = (in) dimension of the matirces (Total number of considered FOs. Same as number of sites if we take only HOMOs into consideration)
// arrays   = (auxiliary/out) data structure that helps in the orthogonalization
////////////////////////////////////////////////////////////////////////
  extern void dpotrf_(char *, long *, double *, long *, long *);
  extern void dpotri_(char *, long *, double *, long *, long *);
  extern void dsyevr_(char *, char *, char *, long *, double *, long *,
                double *, double *, long *, long *, double *, long *, double *,
		double *, long *, double *, double *, long *, double *, long *, long *);
  double abstol=1.0e-8;
  char itype='U';
  char jobz='V';
  char range='A';
  char uplo='U';
  long i, j, k, m, info, eval_no;
  double void_double;
  long void_long;

  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      arrays.sij[i+j*n] = OverlF[i][j];
      arrays.tij[i+j*n] = THamil[i][j];
    }

/*
  // print the arrays
  printf("Orthogonalize - array sij:\n");
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) printf("%12.6f", arrays.sij[i+j*n]);
    printf("\n");
  }
  printf("Orthogonalize - array tij:\n");
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) printf("%12.6f", arrays.tij[i+j*n]);
    printf("\n");
  }
*/

  // Cholesky-factorize the overlap matrix
  dpotrf_(&itype, &n, arrays.sij, &n, &info);
  //if (info) return info;

  // Invert the overlap matrix
  dpotri_(&itype, &n, arrays.sij, &n, &info);
  //if (info) return info;

  for (i=0; i<n; i++)
    for (j=0; j<i; j++)
      arrays.sij[i+j*n] = arrays.sij[j+i*n];

/*
  // print the arrays
  printf("Orthogonalize - inverse of sij:\n");
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) printf("%12.6f", arrays.sij[i+j*n]);
    printf("\n");
  }
*/

  // Diagonalize the inverse of overlap matrix
  // evec - eigenvectors, eval - eigenvalues
  dsyevr_(&jobz, &range, &uplo, &n, arrays.sij, &n,
          &void_double, &void_double, &void_long, &void_long, &abstol, &eval_no, arrays.eval,
	  arrays.evec, &n, arrays.issupz, arrays.work, &(arrays.lwork), arrays.iwork, &(arrays.liwork), &info);
  //if (info) return info;

  // sij = evec * sqrt(diag) * evec-transp
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      arrays.sij[i+j*n] = 0.0;
      for (k=0; k<n; k++)
        arrays.sij[i+j*n] += arrays.evec[i+k*n] * arrays.evec[j+k*n] * sqrt(arrays.eval[k]);
    }

  // THamilOr = sij * tij * sij
  for (i=0; i<n; i++)
    for (j=0; j<n; j++) {
      THamilOr[i][j] = 0.0;
      for (k=0; k<n; k++)
        for (m=0; m<n; m++)
	  THamilOr[i][j] += arrays.sij[i+k*n] * arrays.tij[k+m*n] * arrays.sij[m+j*n];
    }

  return 0;
}
